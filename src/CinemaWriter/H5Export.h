#pragma once

#include <string>

class vtkImageData;

int WriteImageHDF5(vtkImageData* image, const std::string& path, int compressionLevel);
