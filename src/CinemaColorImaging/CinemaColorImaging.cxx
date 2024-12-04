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

//----------------------------------------------------------------------------
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

bool CinemaColorImaging::GetNeedsUpdate(){
  if(this->Ready) return false;

  const auto* pv = pqApplicationCore::instance();
  if(!pv) return false;
  const auto* smm = pv->getServerManagerModel();
  if(!smm) return false;
  QList<pqView*> views = smm->findItems<pqView*>();
  if(views.size()<1) return false;

  auto view = views[0]->getViewProxy();
  if(!view) return false;
  auto window = static_cast<vtkGenericOpenGLRenderWindow*>(view->GetRenderWindow());
  if(!window) return false;

  const auto ready = window->GetReadyForRendering();
  // std::cout<<"READY: "<<ready<<std::endl;
  if(ready && !this->Ready){
    this->Ready = true;
    this->Modified();
    return true;
  } else return false;
};

//----------------------------------------------------------------------------
CinemaColorImaging::~CinemaColorImaging() = default;

int CinemaColorImaging::RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector){

  if(!this->Ready){
    return 1;
  }

  auto cameras = vtkPointSet::GetData(inputVector[0]);
  const int nCameras = cameras->GetNumberOfPoints();

  // Initialize Output
  const auto camPos = getPointer<float>(cameras->GetPoints()->GetData());
  const auto camUp = getPointer<float>(cameras->GetPointData()->GetArray("CameraUp"));
  const auto camDir = getPointer<float>(cameras->GetPointData()->GetArray("CameraDir"));
  const auto camHeight = getPointer<float>(cameras->GetPointData()->GetArray("CameraHeight"));
  const auto camNearFar = getPointer<float>(cameras->GetPointData()->GetArray("CameraNearFar"));

  auto outputCollection = vtkMultiBlockDataSet::GetData(outputVector);

  // for(int c=0; c<nCameras; c++){
  //   auto image = vtkSmartPointer<vtkImageData>::New();
  //   image->SetDimensions(this->Resolution[0], this->Resolution[1], 1);
  //   image->SetSpacing(1, 1, 1);
  //   image->SetOrigin(0, 0, 0);
  //   image->AllocateScalars(VTK_FLOAT, 1);
  //   outputCollection->SetBlock(c, image);
  // }

  // window->Initialize();
  // std::cout<<(window->GetInitialized()?'y':'n')<<'\n';

  const auto* pv = pqApplicationCore::instance();
  const auto* smm = pv->getServerManagerModel();
  QList<pqView*> views = smm->findItems<pqView*>();

  // views[0]->forceRender();
  auto view = static_cast<vtkSMRenderViewProxy*>(views[0]->getViewProxy());
  auto camera = view->GetActiveCamera();
  auto window = view->GetRenderWindow();
  const int resX = window->GetSize()[0];
  const int resY = window->GetSize()[1];

  window->SetSize(this->Resolution[0],this->Resolution[0]);

  for(int c=0; c<nCameras; c++){
    std::cout<<c<<std::endl;
    camera->SetPosition(camPos[c*3+0],camPos[c*3+1],camPos[c*3+2]);
    camera->SetViewUp(camUp[c*3+0],camUp[c*3+1],camUp[c*3+2]);
    camera->SetFocalPoint(camPos[c*3+0]+camDir[c*3+0],camPos[c*3+1]+camDir[c*3+1],camPos[c*3+2]+camDir[c*3+2]);
    outputCollection->SetBlock(c, view->CaptureWindow(1));
  }

  // reset size
  window->SetSize(resX,resY);

  // std::cout<<view.<<std::endl;

// #ifdef _OPENMP
//   this->printMsg("#Imaging ("+std::to_string(nCameras)+" images, "+std::to_string(omp_get_max_threads())+" threads)");
//   #pragma omp parallel for
// #endif
//   for(int c=0; c<nCameras; c++){
//     std::cout<<c<<std::endl;
//     // auto image = static_cast<vtkImageData*>(outputCollection->GetBlock(c));


//     // size_t const nPixels = this->Resolution[0] * this->Resolution[1];
//     // auto imagePD = image->GetPointData();

//     // auto depthBuffer = imagePD->GetArray(0);
//     // depthBuffer->SetName("Depth");
//     // auto depthBuffer_ = getPointer<float>(depthBuffer);

//     // auto primitiveId = vtkSmartPointer<vtkUnsignedIntArray>::New();
//     // primitiveId->SetName("PrimitiveId");
//     // primitiveId->SetNumberOfComponents(1);
//     // primitiveId->SetNumberOfTuples(nPixels);
//     // auto primitiveId_ = getPointer<unsigned int>(primitiveId);
//     // // imagePD->AddArray(primitiveId);

//     // auto barycentricCoordinates = vtkSmartPointer<vtkFloatArray>::New();
//     // barycentricCoordinates->SetName("BarycentricCoordinates");
//     // barycentricCoordinates->SetNumberOfComponents(2);
//     // barycentricCoordinates->SetNumberOfTuples(nPixels);
//     // auto barycentricCoordinates_ = getPointer<float>(barycentricCoordinates);
//     // // imagePD->AddArray(barycentricCoordinates);

//     // renderImage(
//     //   depthBuffer_,
//     //   primitiveId_,
//     //   barycentricCoordinates_,

//     //   scene,
//     //   this->Resolution,
//     //   &camPos[c*3],
//     //   &camDir[c*3],
//     //   &camUp[c*3],
//     //   camHeight[c],
//     //   &camNearFar[c*2]
//     // );

//     // MapPointAndCellData(
//     //   image,

//     //   inputObject_, primitiveId_, barycentricCoordinates_,
//     //   inputObjectConnectivityList
//     // );

//     // AddAllFieldDataArrays( cameras, image, c );
//   }

  return 1;
};
