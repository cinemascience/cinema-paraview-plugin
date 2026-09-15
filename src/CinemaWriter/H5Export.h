#pragma once

#include <string>

class vtkImageData;

bool WriteImageHDF5(vtkImageData* image, const std::string& path, int compressionLevel);
