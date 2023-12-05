#include <pcAlgorithm.h>

#include <vtkDataSet.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationIntegerKey.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCompositeDataPipeline.h>

// Pass input type information key
#include <vtkInformationKey.h>
vtkInformationKeyMacro(pcAlgorithm, SAME_DATA_TYPE_AS_INPUT_PORT, Integer);

// Constructor / Destructor
vtkStandardNewMacro(pcAlgorithm);
pcAlgorithm::pcAlgorithm() = default;
pcAlgorithm::~pcAlgorithm() = default;

template <class vtkDataType>
int prepOutput(vtkInformation *info, const std::string &className) {
  auto output = vtkDataObject::GetData(info);
  if(!output || !output->IsA(className.data())) {
    auto newOutput = vtkSmartPointer<vtkDataType>::New();
    info->Set(vtkDataObject::DATA_OBJECT(), newOutput);
  }
  return 1;
}

vtkDataSet *pcAlgorithm::GetOutput() {
  return this->GetOutput(0);
}

vtkDataSet *pcAlgorithm::GetOutput(int port) {
  return vtkDataSet::SafeDownCast(this->GetOutputDataObject(port));
}

void pcAlgorithm::SetInputData(vtkDataSet *input) {
  this->SetInputData(0, input);
}

void pcAlgorithm::SetInputData(int index, vtkDataSet *input) {
  this->SetInputDataInternal(index, input);
}

void pcAlgorithm::AddInputData(vtkDataSet *input) {
  this->AddInputData(0, input);
}

void pcAlgorithm::AddInputData(int index, vtkDataSet *input) {
  this->AddInputDataInternal(index, input);
}

int pcAlgorithm::RequestDataObject(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  // for each output
  for(int i = 0; i < this->GetNumberOfOutputPorts(); ++i) {
    auto outInfo = outputVector->GetInformationObject(i);
    if(!outInfo) {
      return 0;
    }

    auto outputPortInfo = this->GetOutputPortInformation(i);

    // always request output type again for dynamic filter outputs
    if(!this->FillOutputPortInformation(i, outputPortInfo)) {
      return 0;
    }

    if(outputPortInfo->Has(pcAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT())) {
      // Set output data type to input data type at specified port
      auto inPortIndex
        = outputPortInfo->Get(pcAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT());
      if(inPortIndex < 0 || inPortIndex >= this->GetNumberOfInputPorts()) {
        return 0;
      }
      auto inInfo = inputVector[inPortIndex]->GetInformationObject(0);
      if(!inInfo) {
        return 0;
      }

      auto input = vtkDataObject::GetData(inInfo);
      auto output = vtkDataObject::GetData(outInfo);

      if(!output || !output->IsA(input->GetClassName())) {
        auto newOutput
          = vtkSmartPointer<vtkDataObject>::Take(input->NewInstance());
        outputPortInfo->Set(
          vtkDataObject::DATA_TYPE_NAME(), input->GetClassName());
        outInfo->Set(vtkDataObject::DATA_OBJECT(), newOutput);
      }
    } else {
      // Explicitly create output by data type name
      if(!outputPortInfo->Has(vtkDataObject::DATA_TYPE_NAME())) {
        return 0;
      }
      std::string const outputType
        = outputPortInfo->Get(vtkDataObject::DATA_TYPE_NAME());

      if(outputType == "vtkUnstructuredGrid") {
        prepOutput<vtkUnstructuredGrid>(outInfo, outputType);
      } else if(outputType == "vtkPolyData") {
        prepOutput<vtkPolyData>(outInfo, outputType);
      } else if(outputType == "vtkMultiBlockDataSet") {
        prepOutput<vtkMultiBlockDataSet>(outInfo, outputType);
      } else if(outputType == "vtkTable") {
        prepOutput<vtkTable>(outInfo, outputType);
      } else if(outputType == "vtkImageData") {
        prepOutput<vtkImageData>(outInfo, outputType);
      } else {
        return 0;
      }
    }
  }

  return 1;
}

//==============================================================================
int pcAlgorithm::ProcessRequest(vtkInformation *request,
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector) {
  // 1. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_OBJECT())) {
    return this->RequestDataObject(request, inputVector, outputVector);
  }

  // 2. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_INFORMATION())) {
    return this->RequestInformation(request, inputVector, outputVector);
  }

  // 3. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_TIME())) {
    return this->RequestUpdateTime(request, inputVector, outputVector);
  }

  // 4. Pass
  if(request->Has(
       vtkCompositeDataPipeline::REQUEST_TIME_DEPENDENT_INFORMATION())) {
    return this->RequestUpdateTimeDependentInformation(
      request, inputVector, outputVector);
  }

  // 5. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_UPDATE_EXTENT())) {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }

  // 6. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA_NOT_GENERATED())) {
    return this->RequestDataNotGenerated(request, inputVector, outputVector);
  }

  // 7. Pass
  if(request->Has(vtkCompositeDataPipeline::REQUEST_DATA())) {
#ifdef TTK_ENABLE_MPI
    if(ttk::hasInitializedMPI() && inputVector != nullptr) {
      if(this->updateMPICommunicator(vtkDataSet::GetData(inputVector[0], 0))) {
        return 1;
      };
    }
#endif // TTK_ENABLE_MPI
    return this->RequestData(request, inputVector, outputVector);
  }

  request->Print(cout);

  return 0;
};
