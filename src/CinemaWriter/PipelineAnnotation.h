// PipelineAnnotation.h

#pragma once

#include <string>

class vtkObjectBase;

std::string ComputeInputTree(vtkObjectBase* clientObject);
