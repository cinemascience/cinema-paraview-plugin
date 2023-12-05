#include <pcImageCompositing.h>

#include <vtkObjectFactory.h>
#include <vtkInformation.h>

#include <vtkMultiProcessController.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>

#include <cmath>

vtkStandardNewMacro(pcImageCompositing);

class CompositeDepthOp : public vtkCommunicator::Operation {
  void Function(const void* A, void* B, vtkIdType length, int datatype) override
  {
    const float* aPtr = reinterpret_cast<const float*>(A);
    float* bPtr = reinterpret_cast<float*>(B);
    for(vtkIdType i=0; i<length; i++){
      bPtr[i] = std::min(aPtr[i],bPtr[i]);
    }
  }

  int Commutative() override { return 0; }
};

class CompositeArrayOp : public vtkCommunicator::Operation {
  void Function(const void* A, void* B, vtkIdType length, int datatype) override
  {
    auto aPtr = reinterpret_cast<const float*>(A);
    auto bPtr = reinterpret_cast<float*>(B);
    for(vtkIdType i=0; i<length; i++){
      bPtr[i] = std::isnan(aPtr[i]) ? bPtr[i] : aPtr[i];
    }
  }

  int Commutative() override { return 0; }
};

pcImageCompositing::pcImageCompositing() {
  this->Controller = nullptr;
  this->SetController(vtkMultiProcessController::GetGlobalController());

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

void pcImageCompositing::SetController(vtkMultiProcessController* c){
  this->Controller = c;
}

pcImageCompositing::~pcImageCompositing() = default;

int pcImageCompositing::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else
    return 0;
  return 1;
}

int pcImageCompositing::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(pcAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

template<typename DT>
DT* getPointer(vtkAbstractArray* array){
  return static_cast<DT*>(array->GetVoidPointer(0));
}

int pcImageCompositing::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  // Get Input and Output
  auto inputMB = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  if (this->Controller->GetNumberOfProcesses() < 2) {
    outputMB->ShallowCopy(inputMB);
    return 1;
  }

  int procid = this->Controller->GetLocalProcessId();

  size_t nBlocks = inputMB->GetNumberOfBlocks();

  CompositeDepthOp compositeDepthOp;
  CompositeArrayOp compositeArrayOp;

  // std::array<char[40],12> arrayNames;

  for(size_t b=0; b<nBlocks; b++){
    auto iImage = vtkImageData::SafeDownCast(inputMB->GetBlock(b));
    if(!iImage) continue;

    auto oImage = vtkSmartPointer<vtkImageData>::New();
    oImage->DeepCopy(iImage);

    auto iImagePD = iImage->GetPointData();
    auto oImagePD = oImage->GetPointData();

    auto iDepth = iImagePD->GetArray("Depth");
    if(!iDepth) continue;

    auto oDepth = oImagePD->GetArray("Depth");

    // composite global depth
    auto iDepthData = getPointer<float>(iDepth);
    auto oDepthData = getPointer<float>(oDepth);
    this->Controller->AllReduce(iDepthData, oDepthData, iDepth->GetNumberOfValues(), &compositeDepthOp);

    // composite arrays
    for(size_t a=0; a<iImagePD->GetNumberOfArrays(); a++){
      auto iArray = iImagePD->GetArray(a);
      auto oArray = oImagePD->GetArray(a);
      if(!iArray || std::string(iArray->GetName()).compare("Depth")==0 || iArray->GetDataType()!=VTK_FLOAT)
        continue;

      auto iArrayData = getPointer<float>(iArray);
      auto oArrayData = getPointer<float>(oArray);

      const size_t nValues = iArray->GetNumberOfValues();

      // mask iArray
      const auto nan_value = std::numeric_limits<float>::quiet_NaN();
      for(size_t i=0; i<nValues; i++)
        iArrayData[i] = oDepthData[i]==iDepthData[i] ? iArrayData[i] : nan_value;

      this->Controller->Reduce(iArrayData, oArrayData, nValues, &compositeArrayOp, 0);
    }

    if(procid==0)
      outputMB->SetBlock(b,oImage);
  }

  // outputMB->ShallowCopy(inputMB);

  return 1;
}
