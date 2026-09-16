#include "vtkSMCinemaColorImagingExtractWriterProxy.h"

#include <CinemaWriter.h>
#include <PipelineProvenance.h>
#include <pqActiveObjects.h>
#include <pqView.h>
#include <string>
#include <vector>
#include <vtkAlgorithm.h>
#include <vtkCamera.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPVDataInformation.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkSMDomain.h>
#include <vtkSMExtractsController.h>
#include <vtkSMInputProperty.h>
#include <vtkSMOutputPort.h>
#include <vtkSMPropertyHelper.h>
#include <vtkSMRenderViewProxy.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMUncheckedPropertyHelper.h>
#include <vtkSmartPointer.h>

vtkStandardNewMacro(vtkSMCinemaColorImagingExtractWriterProxy);

namespace {
vtkSMSourceProxy* GetCameraSource(vtkSMProxy* proxy, unsigned int& port) {
  port = 0;
  if (!proxy || !proxy->GetProperty("Input")) { return nullptr; }

  vtkSMPropertyHelper input(proxy, "Input");
  port = input.GetOutputPort();
  return vtkSMSourceProxy::SafeDownCast(input.GetAsProxy());
}

vtkSMRenderViewProxy* GetActiveRenderView() {
  pqView* activeView = pqActiveObjects::instance().activeView();
  return activeView ? vtkSMRenderViewProxy::SafeDownCast(activeView->getViewProxy()) : nullptr;
}

vtkSMRenderViewProxy* GetSelectedRenderView(vtkSMProxy* proxy, const char* requestedName, bool& usedFallback) {
  usedFallback = false;
  vtkSMSessionProxyManager* pxm = proxy ? proxy->GetSessionProxyManager() : nullptr;
  if (!pxm) { return nullptr; }

  if (requestedName && *requestedName) {
    if (auto* view = vtkSMRenderViewProxy::SafeDownCast(pxm->GetProxy("views", requestedName))) { return view; }
    usedFallback = true;
  }

  return GetActiveRenderView();
}

vtkSmartPointer<vtkPointSet> FetchCameras(vtkSMSourceProxy* source, unsigned int port, double time) {
  if (!source) { return nullptr; }

  source->UpdatePipeline(time);
  vtkPVDataInformation* information = source->GetDataInformation(port);
  const int dataType = information ? information->GetDataSetType() : -1;
  vtkSMSessionProxyManager* pxm = source->GetSessionProxyManager();
  if (dataType < 0 || !pxm) { return nullptr; }

  vtkSmartPointer<vtkSMSourceProxy> mover;
  mover.TakeReference(vtkSMSourceProxy::SafeDownCast(pxm->NewProxy("filters", "ClientServerMoveData")));
  if (!mover) { return nullptr; }

  auto* input = vtkSMInputProperty::SafeDownCast(mover->GetProperty("Input"));
  if (!input) { return nullptr; }

  input->RemoveAllProxies();
  input->AddInputConnection(source, port);
  vtkSMPropertyHelper(mover, "OutputDataType").Set(dataType);
  mover->UpdateVTKObjects();
  mover->UpdatePipeline(time);

  auto* algorithm = vtkAlgorithm::SafeDownCast(mover->GetClientSideObject());
  auto* fetched = algorithm ? vtkPointSet::SafeDownCast(algorithm->GetOutputDataObject(0)) : nullptr;
  if (!fetched) { return nullptr; }

  vtkSmartPointer<vtkPointSet> copy;
  copy.TakeReference(vtkPointSet::SafeDownCast(fetched->NewInstance()));
  if (copy) { copy->ShallowCopy(fetched); }
  return copy;
}

void AddTuple(vtkFieldData* target, vtkDataArray* source, vtkIdType tuple, const char* name = nullptr) {
  if (!target || !source || tuple < 0 || tuple >= source->GetNumberOfTuples()) { return; }

  vtkSmartPointer<vtkDataArray> value;
  value.TakeReference(vtkDataArray::SafeDownCast(source->NewInstance()));
  if (!value) { return; }

  value->SetName(name ? name : source->GetName());
  value->SetNumberOfComponents(source->GetNumberOfComponents());
  value->SetNumberOfTuples(1);
  value->SetTuple(0, tuple, source);
  target->AddArray(value);
}

void AddCameraFieldData(vtkPointSet* cameras, vtkImageData* image, vtkIdType tuple) {
  vtkFieldData* fields = image->GetFieldData();
  vtkPointData* pointData = cameras->GetPointData();
  for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) { AddTuple(fields, pointData->GetArray(i), tuple); }
  AddTuple(fields, cameras->GetPoints()->GetData(), tuple, "CameraPos");
}

bool GetCameraVector(vtkPointSet* cameras, const char* name, vtkIdType index, double value[3]) {
  vtkDataArray* array = cameras->GetPointData()->GetArray(name);
  if (!array || array->GetNumberOfComponents() < 3 || index >= array->GetNumberOfTuples()) { return false; }

  array->GetTuple(index, value);
  return true;
}

class ViewStateGuard {
public:
  ViewStateGuard(vtkSMRenderViewProxy* view, double time) : View(view), SavedCamera(vtkSmartPointer<vtkCamera>::New()) {
    this->View->SynchronizeCameraProperties();
    this->SavedCamera->DeepCopy(this->View->GetActiveCamera());
    vtkSMPropertyHelper(this->View, "ViewSize").Get(this->ViewSize, 2);

    if (this->View->GetProperty("ViewTime")) {
      vtkSMPropertyHelper timeProperty(this->View, "ViewTime");
      this->SavedTime = timeProperty.GetAsDouble();
      this->HasViewTime = true;
      timeProperty.Set(time);
      this->View->UpdateVTKObjects();
    }
  }

  ~ViewStateGuard() {
    this->View->GetActiveCamera()->DeepCopy(this->SavedCamera);
    this->View->SynchronizeCameraProperties();
    vtkSMPropertyHelper(this->View, "ViewSize").Set(this->ViewSize, 2);
    if (this->HasViewTime) { vtkSMPropertyHelper(this->View, "ViewTime").Set(this->SavedTime); }
    this->View->UpdateVTKObjects();
    this->View->StillRender();
  }

private:
  vtkSMRenderViewProxy* View;
  vtkSmartPointer<vtkCamera> SavedCamera;
  int ViewSize[2] = {0, 0};
  double SavedTime = 0.0;
  bool HasViewTime = false;
};
} // namespace

void vtkSMCinemaColorImagingExtractWriterProxy::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

bool vtkSMCinemaColorImagingExtractWriterProxy::CanExtract(vtkSMProxy* proxy) {
  if (!proxy) { return false; }

  unsigned int port = 0;
  if (auto* output = vtkSMOutputPort::SafeDownCast(proxy)) {
    port = output->GetPortIndex();
    proxy = output->GetSourceProxy();
  }

  vtkSMProperty* input = this->GetProperty("Input");
  if (!input) { return false; }

  vtkSMUncheckedPropertyHelper helper(input);
  helper.Set(proxy, port);
  const bool supported = input->IsInDomains() == vtkSMDomain::IN_DOMAIN;
  helper.SetNumberOfElements(0);
  return supported;
}

bool vtkSMCinemaColorImagingExtractWriterProxy::IsExtracting(vtkSMProxy* proxy) {
  unsigned int port = VTK_UNSIGNED_INT_MAX;
  if (auto* output = vtkSMOutputPort::SafeDownCast(proxy)) {
    port = output->GetPortIndex();
    proxy = output->GetSourceProxy();
  }

  vtkSMPropertyHelper input(this, "Input");
  return input.GetAsProxy() == proxy && (port == VTK_UNSIGNED_INT_MAX || input.GetOutputPort() == port);
}

void vtkSMCinemaColorImagingExtractWriterProxy::SetInput(vtkSMProxy* proxy) {
  if (!proxy) {
    vtkErrorMacro("Input cannot be nullptr.");
    return;
  }

  unsigned int port = 0;
  if (auto* output = vtkSMOutputPort::SafeDownCast(proxy)) {
    port = output->GetPortIndex();
    proxy = output->GetSourceProxy();
  }
  vtkSMPropertyHelper(this, "Input").Set(proxy, port);
}

vtkSMProxy* vtkSMCinemaColorImagingExtractWriterProxy::GetInput() {
  vtkSMPropertyHelper input(this, "Input");
  auto* source = vtkSMSourceProxy::SafeDownCast(input.GetAsProxy());
  return source ? source->GetOutputPort(input.GetOutputPort()) : nullptr;
}

bool vtkSMCinemaColorImagingExtractWriterProxy::Write(vtkSMExtractsController* extractor) {
  if (!extractor) { return false; }

  unsigned int cameraPort = 0;
  vtkSMSourceProxy* cameraSource = GetCameraSource(this, cameraPort);
  if (!cameraSource) {
    vtkErrorMacro("Cinema color imaging requires a camera input.");
    return false;
  }

  const char* requestedViewName = vtkSMPropertyHelper(this, "View").GetAsString();
  bool usedViewFallback = false;
  vtkSMRenderViewProxy* view = GetSelectedRenderView(this, requestedViewName, usedViewFallback);
  if (!view) {
    vtkErrorMacro("Cinema color imaging requires a valid RenderView.");
    return false;
  }
  if (usedViewFallback) { vtkWarningMacro("Selected RenderView is unavailable; using the active RenderView."); }
  if (view->GetSession() != cameraSource->GetSession()) {
    vtkErrorMacro("The selected RenderView and camera source belong to different sessions.");
    return false;
  }

  vtkSmartPointer<vtkPointSet> cameras = FetchCameras(cameraSource, cameraPort, extractor->GetTime());
  if (!cameras || cameras->GetNumberOfPoints() == 0) {
    vtkErrorMacro("Unable to fetch camera data to the client.");
    return false;
  }

  vtkDataArray* cameraUp = cameras->GetPointData()->GetArray("CameraUp");
  vtkDataArray* cameraDir = cameras->GetPointData()->GetArray("CameraDir");
  if (!cameraUp || !cameraDir) {
    vtkErrorMacro("Camera input must provide CameraUp and CameraDir point-data arrays.");
    return false;
  }

  const bool orthographic = vtkSMPropertyHelper(this, "ProjectionMode").GetAsInt() == 1;
  vtkDataArray* cameraHeight = cameras->GetPointData()->GetArray("CameraHeight");
  if (orthographic && !cameraHeight) {
    vtkErrorMacro("Orthographic projection requires a CameraHeight point-data array.");
    return false;
  }

  int resolution[2];
  vtkSMPropertyHelper(this, "ImageResolution").Get(resolution, 2);
  const double cameraAngle = vtkSMPropertyHelper(this, "CameraAngle").GetAsDouble();

  ViewStateGuard restore(view, extractor->GetTime());
  vtkSMPropertyHelper(view, "ViewSize").Set(resolution, 2);
  view->UpdateVTKObjects();

  vtkCamera* camera = view->GetActiveCamera();
  vtkNew<vtkMultiBlockDataSet> images;
  images->SetNumberOfBlocks(static_cast<unsigned int>(cameras->GetNumberOfPoints()));

  for (vtkIdType index = 0; index < cameras->GetNumberOfPoints(); ++index) {
    double position[3];
    double up[3];
    double direction[3];
    cameras->GetPoint(index, position);
    if (!GetCameraVector(cameras, "CameraUp", index, up) || !GetCameraVector(cameras, "CameraDir", index, direction)) {
      vtkErrorMacro("Invalid CameraUp/CameraDir tuple at camera " << index << ".");
      return false;
    }

    const double focalPoint[3] = {position[0] + direction[0], position[1] + direction[1], position[2] + direction[2]};

    camera->SetPosition(position);
    camera->SetViewUp(up);
    camera->SetFocalPoint(focalPoint);
    if (orthographic) {
      camera->ParallelProjectionOn();
      camera->SetParallelScale(0.5 * cameraHeight->GetComponent(index, 0));
    } else {
      camera->ParallelProjectionOff();
      camera->SetViewAngle(cameraAngle);
    }

    view->SynchronizeCameraProperties();

    vtkSmartPointer<vtkImageData> image;
    image.TakeReference(view->CaptureWindow(1));
    if (!image) {
      vtkErrorMacro("Failed to capture camera " << index << ".");
      return false;
    }

    double finalNearFar[2];
    view->GetActiveCamera()->GetClippingRange(finalNearFar);

    if (vtkDataArray* colors = image->GetPointData()->GetScalars()) {
      const int components = colors->GetNumberOfComponents();
      if (components == 3) {
        colors->SetName("rgb");
      } else if (components == 4) {
        colors->SetName("rgba");
      }
    }

    AddCameraFieldData(cameras, image, index);
    image->GetFieldData()->RemoveArray("CameraNearFar");
    vtkNew<vtkDoubleArray> nearFar;
    nearFar->SetName("CameraNearFar");
    nearFar->SetNumberOfComponents(2);
    nearFar->SetNumberOfTuples(1);
    nearFar->SetTuple(0, finalNearFar);
    image->GetFieldData()->AddArray(nearFar);

    images->SetBlock(static_cast<unsigned int>(index), image);
  }

  const ProvenanceMap provenance = ComputeVisibleRepresentationProvenance(view, std::vector<vtkSMProxy*>{cameraSource});

  const char* fileName = vtkSMPropertyHelper(this, "FileName").GetAsString();
  if (!fileName || !*fileName) {
    vtkErrorMacro("Missing FileName property.");
    return false;
  }

  const std::string outputDirectory =
    this->GenerateExtractsFileName(fileName, extractor->GetRealExtractsOutputDirectory());

  vtkNew<CinemaWriter> writer;
  writer->SetInputDataObject(images);
  writer->SetOutputDirectory(outputDirectory);
  writer->SetFormat(vtkSMPropertyHelper(this, "Format").GetAsInt());
  writer->SetCompressionLevel(vtkSMPropertyHelper(this, "CompressionLevel").GetAsInt());
  for (const auto& item : provenance) { writer->AddProvenanceEntry(item.first.c_str(), item.second.c_str()); }
  writer->Update();

  extractor->AddSummaryEntry(this, outputDirectory);
  return true;
}
