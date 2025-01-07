#include "CinemaColorImaging.h"

#include <vtkDataObject.h>
#include <vtkObjectFactory.h>
#include <vtkInformation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkImageData.h>

#include <pqApplicationCore.h>
#include <pqView.h>
#include <pqServerManagerModel.h>
#include <vtkSMViewProxy.h>
#include <vtkSMRenderViewProxy.h>
#include <QList>
#include <vtkSMSaveScreenshotProxy.h>
#include <vtkCamera.h>
#include <vtkVector.h>       // needed for vtkVector2i.
#include <vtkRenderWindow.h>       // needed for vtkVector2i.
#include <vtkGenericOpenGLRenderWindow.h>       // needed for vtkVector2i.

vtkStandardNewMacro(CinemaColorImaging);

template<typename PT>
PT* getPointer(vtkAbstractArray* array){
  return static_cast<PT*>(array->GetVoidPointer(0));
}

//------------------------------------------------------------------------------
CinemaColorImaging::CinemaColorImaging(){
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
};

int CinemaColorImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int CinemaColorImaging::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------
int AddFieldDataArray_(vtkFieldData *fd, vtkDataArray *array, int tupleIdx, const std::string name="") {
  if(!array)
    return 0;

  size_t nComponents = array->GetNumberOfComponents();

  vtkNew<vtkFloatArray> newArray;
  newArray->SetName(name.empty() ? array->GetName() : name.data());
  newArray->SetNumberOfComponents(nComponents);
  newArray->SetNumberOfTuples(1);

  if(newArray->GetDataType() == array->GetDataType()) {
    newArray->SetTuple(0, tupleIdx, array);
  } else {
    for(size_t i = 0; i < nComponents; i++)
      newArray->SetValue(
        i, array->GetVariantValue(tupleIdx * nComponents + i).ToDouble());
  }

  fd->AddArray(newArray);

  return 1;
};

int AddAllFieldDataArrays_(vtkPointSet *cameras, vtkImageData *image, int tupleIdx) {
  auto imageFD = image->GetFieldData();

  auto camerasPD = cameras->GetPointData();
  for(int i = 0; i < camerasPD->GetNumberOfArrays(); i++) {
    AddFieldDataArray_(imageFD, camerasPD->GetArray(i), tupleIdx);
  }
  AddFieldDataArray_(imageFD, cameras->GetPoints()->GetData(), tupleIdx, "CameraPos");

  return 1;
};

//------------------------------------------------------------------------------
CinemaColorImaging::~CinemaColorImaging() = default;

int CinemaColorImaging::RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector){


  // initialize camera data
  auto cameras = vtkPointSet::GetData(inputVector[0]);
  const int nCameras = cameras->GetNumberOfPoints();
  const auto camPos = getPointer<float>(cameras->GetPoints()->GetData());
  const auto camUp = getPointer<float>(cameras->GetPointData()->GetArray("CameraUp"));
  const auto camDir = getPointer<float>(cameras->GetPointData()->GetArray("CameraDir"));
  const auto camHeight = getPointer<float>(cameras->GetPointData()->GetArray("CameraHeight"));
  const auto camNearFar = getPointer<float>(cameras->GetPointData()->GetArray("CameraNearFar"));

  // get render window
  const auto* pv = pqApplicationCore::instance();
  const auto* smm = pv->getServerManagerModel();
  QList<pqView*> views = smm->findItems<pqView*>();
  auto view = static_cast<vtkSMRenderViewProxy*>(views[0]->getViewProxy());
  auto camera = view->GetActiveCamera();
  auto window = view->GetRenderWindow();
  const int prev_resX = window->GetSize()[0];
  const int prev_resY = window->GetSize()[1];

  // render images
  this->printMsg("# Color Imaging (#"+std::to_string(nCameras)+")");

  window->SetSize(this->Resolution[0],this->Resolution[1]);
  auto outputCollection = vtkMultiBlockDataSet::GetData(outputVector);
  for(int c=0; c<nCameras; c++){
    camera->SetPosition(camPos[c*3+0],camPos[c*3+1],camPos[c*3+2]);
    camera->SetViewUp(camUp[c*3+0],camUp[c*3+1],camUp[c*3+2]);
    camera->SetFocalPoint(camPos[c*3+0]+camDir[c*3+0],camPos[c*3+1]+camDir[c*3+1],camPos[c*3+2]+camDir[c*3+2]);
    vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::Take(view->CaptureWindow(1));
    outputCollection->SetBlock(c, image);
    AddAllFieldDataArrays_( cameras, image, c );
  }

  // reset window
  window->SetSize(prev_resX,prev_resY);
  window->Render();

  return 1;
};
