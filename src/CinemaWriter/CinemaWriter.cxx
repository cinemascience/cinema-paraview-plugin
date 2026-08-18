#include "CinemaWriter.h"

#include "H5Export.h"
#include "PipelineAnnotation.h"

#include <algorithm>
#include <vtkDataArray.h>
#include <vtkDataObject.h>
#include <vtkDirectory.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
// #include <H5Cpp.h>
#include <fstream>
#include <sstream>
#include <vtk_hdf5.h>

vtkStandardNewMacro(CinemaWriter);

//----------------------------------------------------------------------------
CinemaWriter::CinemaWriter() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
};
CinemaWriter::~CinemaWriter() = default;

int CinemaWriter::FillInputPortInformation(int port, vtkInformation* info) {
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    return 1;
  }
  return 0;
};

int CinemaWriter::FillOutputPortInformation(int port, vtkInformation* info) {
  if (port == 0) {
    info->Set(CinemaAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
};

int CinemaWriter::RequestData(vtkInformation* request, vtkInformationVector** inputVector,
                              vtkInformationVector* outputVector) {
  auto annotation = ComputeInputTree(this);

  std::cout << "Pipeline annotation:\n" << annotation << std::endl;

  return 1;
}
