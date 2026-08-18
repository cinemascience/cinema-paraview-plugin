#include "H5Export.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <vtkDataArray.h>
#include <vtkFieldData.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtk_hdf5.h>

namespace {

bool AddH5DataSet(hid_t group, const std::string& name, const std::vector<hsize_t>& dims, const void* data,
                  hid_t dataType, int compressionLevel) {
  if (group < 0 || name.empty() || dims.empty() || !data) { return false; }

  for (hsize_t dim : dims) {
    if (dim == 0) {
      std::cerr << "Cannot create HDF5 dataset '" << name << "' with a zero-sized dimension." << std::endl;

      return false;
    }
  }

  const int rank = static_cast<int>(dims.size());

  hid_t dataspace = H5Screate_simple(rank, dims.data(), nullptr);

  if (dataspace < 0) {
    std::cerr << "Failed to create HDF5 dataspace for '" << name << "'." << std::endl;

    return false;
  }

  hid_t creationProperty = H5P_DEFAULT;
  hid_t chunkProperty = -1;

  /*
   * Compression in HDF5 requires chunked storage.
   *
   * Only create a dataset creation property list when compression
   * is actually requested.
   */
  if (compressionLevel > 0) {
    chunkProperty = H5Pcreate(H5P_DATASET_CREATE);

    if (chunkProperty < 0) {
      H5Sclose(dataspace);
      return false;
    }

    std::vector<hsize_t> chunkDims(dims);

    /*
     * Preserve the intent of the old chunking strategy:
     *
     *   tuples:  up to 256
     *   image Y: up to 32
     *
     * but never create a chunk larger than the dataset.
     */
    if (rank >= 1) { chunkDims[0] = std::min<hsize_t>(dims[0], 256); }

    if (rank >= 2) { chunkDims[1] = std::min<hsize_t>(dims[1], 32); }

    if (rank >= 3) {
      /*
       * Usually the component count, e.g. RGB/RGBA.
       * Keep the whole component dimension in one chunk.
       */
      chunkDims[2] = dims[2];
    }

    if (H5Pset_chunk(chunkProperty, rank, chunkDims.data()) < 0) {
      std::cerr << "Failed to configure HDF5 chunking for '" << name << "'." << std::endl;

      H5Pclose(chunkProperty);
      H5Sclose(dataspace);
      return false;
    }

    const unsigned int level = static_cast<unsigned int>(std::clamp(compressionLevel, 0, 9));

    if (H5Pset_deflate(chunkProperty, level) < 0) {
      std::cerr << "Failed to enable HDF5 compression for '" << name << "'." << std::endl;

      H5Pclose(chunkProperty);
      H5Sclose(dataspace);
      return false;
    }

    creationProperty = chunkProperty;
  }

  hid_t dataset = H5Dcreate(group, name.c_str(), dataType, dataspace, H5P_DEFAULT, creationProperty, H5P_DEFAULT);

  if (dataset < 0) {
    std::cerr << "Failed to create HDF5 dataset '" << name << "'." << std::endl;

    if (chunkProperty >= 0) { H5Pclose(chunkProperty); }

    H5Sclose(dataspace);
    return false;
  }

  const herr_t writeStatus = H5Dwrite(dataset, dataType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset);

  if (chunkProperty >= 0) { H5Pclose(chunkProperty); }

  H5Sclose(dataspace);

  if (writeStatus < 0) {
    std::cerr << "Failed to write HDF5 dataset '" << name << "'." << std::endl;

    return false;
  }

  return true;
}

template <typename DT>
bool WriteArray(hid_t group, int compressionLevel, hid_t dataType, vtkDataArray* array, hsize_t resX = 0,
                hsize_t resY = 0) {
  if (!array) { return false; }

  const vtkIdType numberOfTuples = array->GetNumberOfTuples();

  const int numberOfComponents = array->GetNumberOfComponents();

  if (numberOfTuples < 1 || numberOfComponents < 1) { return false; }

  const char* arrayName = array->GetName();

  if (!arrayName || !*arrayName) {
    std::cerr << "Skipping unnamed VTK array." << std::endl;

    return false;
  }

  const auto* rawData = static_cast<const DT*>(array->GetVoidPointer(0));

  if (!rawData) { return false; }

  /*
   * Image channel:
   *
   * VTK point ordering is flipped vertically before storage.
   */
  if (resX > 0 && resY > 0) {
    const vtkIdType expectedTuples = static_cast<vtkIdType>(resX * resY);

    if (numberOfTuples != expectedTuples) {
      std::cerr << "Array '" << arrayName << "' has " << numberOfTuples << " tuples, but image resolution " << resX
                << "x" << resY << " requires " << expectedTuples << "." << std::endl;

      return false;
    }

    const std::vector<hsize_t> dims{resY, resX, static_cast<hsize_t>(numberOfComponents)};

    std::vector<DT> data(static_cast<std::size_t>(resX * resY * numberOfComponents));

    for (hsize_t y = 0; y < resY; ++y) {
      for (hsize_t x = 0; x < resX; ++x) {
        const hsize_t inputIndex = y * resX + x;

        const hsize_t outputIndex = (resY - 1 - y) * resX + x;

        for (int component = 0; component < numberOfComponents; ++component) {
          data[outputIndex * numberOfComponents + component] = rawData[inputIndex * numberOfComponents + component];
        }
      }
    }

    return AddH5DataSet(group, arrayName, dims, data.data(), dataType, compressionLevel);
  }

  /*
   * Metadata / generic array.
   *
   * One component:
   *
   *   [tuples]
   *
   * Multiple components:
   *
   *   [tuples, components]
   */
  std::vector<hsize_t> dims;

  if (numberOfComponents == 1) {
    dims = {static_cast<hsize_t>(numberOfTuples)};
  } else {
    dims = {static_cast<hsize_t>(numberOfTuples), static_cast<hsize_t>(numberOfComponents)};
  }

  return AddH5DataSet(group, arrayName, dims, rawData, dataType, compressionLevel);
}

bool WriteSupportedArray(hid_t group, int compressionLevel, vtkDataArray* array, hsize_t resX = 0, hsize_t resY = 0) {
  if (!array) { return false; }

  switch (array->GetDataType()) {
  case VTK_FLOAT: return WriteArray<float>(group, compressionLevel, H5T_NATIVE_FLOAT, array, resX, resY);

  case VTK_DOUBLE: return WriteArray<double>(group, compressionLevel, H5T_NATIVE_DOUBLE, array, resX, resY);

  case VTK_UNSIGNED_CHAR: return WriteArray<std::uint8_t>(group, compressionLevel, H5T_NATIVE_UCHAR, array, resX, resY);

  case VTK_INT: return WriteArray<int>(group, compressionLevel, H5T_NATIVE_INT, array, resX, resY);

  case VTK_UNSIGNED_INT: return WriteArray<unsigned int>(group, compressionLevel, H5T_NATIVE_UINT, array, resX, resY);

  default:
    std::cerr << "Unsupported VTK array type for '" << (array->GetName() ? array->GetName() : "(unnamed)")
              << "': " << array->GetDataTypeAsString() << std::endl;

    return false;
  }
}

} // namespace

int WriteImageHDF5(vtkImageData* image, const std::string& path, int compressionLevel) {
  if (!image || path.empty()) { return 0; }

  int dims[3] = {0, 0, 0};
  image->GetDimensions(dims);

  if (dims[0] <= 0 || dims[1] <= 0) {
    std::cerr << "Cannot write image with invalid dimensions: " << dims[0] << " x " << dims[1] << " x " << dims[2]
              << std::endl;

    return 0;
  }

  compressionLevel = std::clamp(compressionLevel, 0, 9);

  hid_t root = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (root < 0) {
    std::cerr << "Failed to create HDF5 file: " << path << std::endl;

    return 0;
  }

  hid_t meta = H5Gcreate(root, "meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (meta < 0) {
    H5Fclose(root);
    return 0;
  }

  hid_t channels = H5Gcreate(root, "channels", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (channels < 0) {
    H5Gclose(meta);
    H5Fclose(root);
    return 0;
  }

  bool success = true;

  // ---------------------------------------------------------------------------
  // Resolution
  // ---------------------------------------------------------------------------

  {
    /*
     * Preserve your existing Cinema database convention of storing resolution
     * as float.
     */
    const float resolution[2] = {static_cast<float>(dims[0]), static_cast<float>(dims[1])};

    success = AddH5DataSet(meta, "resolution", {2}, resolution, H5T_NATIVE_FLOAT, 0) && success;
  }

  // ---------------------------------------------------------------------------
  // Field data -> /meta
  // ---------------------------------------------------------------------------

  if (vtkFieldData* fieldData = image->GetFieldData()) {
    for (int i = 0; i < fieldData->GetNumberOfArrays(); ++i) {
      vtkDataArray* array = fieldData->GetArray(i);

      if (!array) { continue; }

      success = WriteSupportedArray(meta, compressionLevel, array) && success;
    }
  }

  // ---------------------------------------------------------------------------
  // Point data -> /channels
  // ---------------------------------------------------------------------------

  if (vtkPointData* pointData = image->GetPointData()) {
    for (int i = 0; i < pointData->GetNumberOfArrays(); ++i) {
      vtkDataArray* array = pointData->GetArray(i);

      if (!array) { continue; }

      success = WriteSupportedArray(channels, compressionLevel, array, static_cast<hsize_t>(dims[0]),
                                    static_cast<hsize_t>(dims[1])) &&
                success;
    }
  }

  H5Gclose(channels);
  H5Gclose(meta);
  H5Fclose(root);

  return success ? 1 : 0;
}
