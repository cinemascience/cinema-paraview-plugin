#include "pcCameraGrid.h"

#include <vtkMultiProcessController.h>

#include <vtkObjectFactory.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>

#include <vtkDataSet.h>
#include <vtkPolyData.h>

#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkAbstractTransform.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <cmath>

vtkStandardNewMacro(pcCameraGrid);

class AddBoundsListOperator : public vtkCommunicator::Operation {
  // Description:
  // Performs a "B.AddBounds(A)" operation.
  void Function(const void* A, void* B, vtkIdType length, int datatype) override
  {
    (void)datatype;
    assert((datatype == VTK_DOUBLE) && (length % 6 == 0));
    assert("pre: A vector is nullptr" && (A != nullptr));
    assert("pre: B vector is nullptr" && (B != nullptr));
    vtkBoundingBox box;
    const double* aPtr = reinterpret_cast<const double*>(A);
    double* bPtr = reinterpret_cast<double*>(B);
    for (vtkIdType idx = 0; idx < length; idx += 6)
    {
      box.SetBounds(&bPtr[idx]);
      box.AddBounds(&aPtr[idx]);
      box.GetBounds(&bPtr[idx]);
    }
  }

  // Description:
  // Sets Commutative to true for this operation
  int Commutative() override { return 1; }
};

//----------------------------------------------------------------------------
pcCameraGrid::pcCameraGrid(){
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

void pcCameraGrid::SetController(vtkMultiProcessController* c){
  this->Controller = c;
}

int pcCameraGrid::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int pcCameraGrid::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData" );
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
pcCameraGrid::~pcCameraGrid() = default;

int pcCameraGrid::RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector){

  auto input = vtkDataSet::GetData(inputVector[0]);
  double bounds[6];
  {
    input->GetBounds(bounds);
    if (this->Controller->GetNumberOfProcesses() > 1) {
      double reduced_bounds[6];
      int procid = this->Controller->GetLocalProcessId();
      AddBoundsListOperator operation;
      this->Controller->AllReduce(bounds, reduced_bounds, 6, &operation);
      memcpy(bounds, reduced_bounds, 6 * sizeof(double));
    }
  }

  const double d[3]{
    bounds[1]-bounds[0],
    bounds[3]-bounds[2],
    bounds[5]-bounds[4]
  };
  const double diameter = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
  const double center[3]{
    bounds[0] + d[0]/2.0,
    bounds[2] + d[1]/2.0,
    bounds[4] + d[2]/2.0
  };

  const auto radius = diameter/2.0*this->RadiusFactor;

  vtkNew<vtkSphereSource> sphereSource;
  sphereSource->GenerateNormalsOff();
  sphereSource->SetRadius(1.0);

  sphereSource->SetThetaResolution(this->GetThetaResolution());
  sphereSource->SetStartTheta(this->GetStartTheta());
  sphereSource->SetEndTheta(this->GetEndTheta());

  sphereSource->SetPhiResolution(this->GetPhiResolution());
  // these swaps are intentional for a more intuitive sphere parameterization
  sphereSource->SetStartPhi(-this->GetStartPhi()+90);
  sphereSource->SetEndPhi(-this->GetEndPhi()+90);

  vtkNew<vtkTransform> transform;
  transform->Translate(center);
  transform->Scale(radius,radius,radius);
  if(this->Axis==AXIS::X){
    transform->RotateY(90);
  } else if(this->Axis==AXIS::Y){
    transform->RotateX(-90);
  }

  vtkNew<vtkTransformFilter> transformFilter;
  transformFilter->SetTransform(transform);
  transformFilter->SetInputConnection(0,sphereSource->GetOutputPort());
  transformFilter->Update();

  auto cameraGrid = vtkPolyData::GetData(outputVector);
  cameraGrid->ShallowCopy(transformFilter->GetOutput());

  // compute phi / theta
  {

    auto prepArray = [=](std::string name, int nTuples, int nComponents){
      vtkNew<vtkFloatArray> array;
      array->SetName(name.data());
      array->SetNumberOfComponents(nComponents);
      array->SetNumberOfTuples(nTuples);
      cameraGrid->GetPointData()->AddArray(array);
      return static_cast<float*>(array->GetVoidPointer(0));
    };

    const size_t nPoints = cameraGrid->GetNumberOfPoints();

    auto phi = prepArray("Phi",nPoints,1);
    auto theta = prepArray("Theta",nPoints,1);
    auto camUp = prepArray("CameraUp",nPoints,3);
    float camUp_[3]{
      this->Axis==AXIS::X ? 1.0f : 0.0f,
      this->Axis==AXIS::Y ? 1.0f : 0.0f,
      this->Axis==AXIS::Z ? 1.0f : 0.0f
    };
    float camUp_singularity[3]{
      0.0f,
      this->Axis==AXIS::Z ? 1.0f : 0.0f,
      this->Axis==AXIS::Z ? 0.0f : 1.0f
    };
    auto camDir = prepArray("CameraDir",nPoints,3);
    auto camHeight = prepArray("CameraHeight",nPoints,1);
    auto camNearFar = prepArray("CameraNearFar",nPoints,2);
    const float camHeight_ = this->CamHeight<=0 ? diameter : this->CamHeight;
    const float camNearFar_[2] = {
      (float)this->NearFar[0],
      this->NearFar[1]==0.0 ? (float)diameter : (float)this->NearFar[1]
    };

    {
      const auto coords = static_cast<float *>(sphereSource->GetOutput()->GetPoints()->GetData()->GetVoidPointer(0));
      for(size_t i=0; i<nPoints; i++){
        const auto coords_ = &(coords[i*3]);

        phi[i] = coords_[2]<-0.999999 ? -90 : coords_[2]>0.999999 ? 90 : (180.0-acos(coords_[2])/M_PI*180.0)-90.0;
        theta[i] = std::fmod(atan2(coords_[0],coords_[1])*180.0/M_PI+180.0, 360.0);
        camHeight[i] = camHeight_;
        camNearFar[i*2+0] = camNearFar_[0];
        camNearFar[i*2+1] = camNearFar_[1];

        if(coords_[2]<-0.999999 || coords_[2]>0.999999){
          camUp[i*3+0] = camUp_singularity[0];
          camUp[i*3+1] = camUp_singularity[1];
          camUp[i*3+2] = camUp_singularity[2];
        } else {
          camUp[i*3+0] = camUp_[0];
          camUp[i*3+1] = camUp_[1];
          camUp[i*3+2] = camUp_[2];
        }
      }
    }
    {
      const auto coords = static_cast<float *>(cameraGrid->GetPoints()->GetData()->GetVoidPointer(0));
      for(size_t i=0; i<nPoints; i++){
        const auto coords_ = &(coords[i*3]);
        auto camDir_ = &(camDir[i*3]);
        camDir_[0] = center[0]-coords_[0];
        camDir_[1] = center[1]-coords_[1];
        camDir_[2] = center[2]-coords_[2];
      }
    }
  }

  return 1;
}
