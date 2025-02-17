#include "CinemaDepthImageProjection.h"

#include <vtkObjectFactory.h>
#include <vtkDataObject.h>
#include <vtkInformation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellData.h>

#include <vtkIntArray.h>
#include <vtkDoubleArray.h>

#ifdef _OPENMP
#include <omp.h>
#endif

vtkStandardNewMacro(CinemaDepthImageProjection);

//----------------------------------------------------------------------------
CinemaDepthImageProjection::CinemaDepthImageProjection(){
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int CinemaDepthImageProjection::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int CinemaDepthImageProjection::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
  }
  return 0;
}

template<typename PT>
PT* getPointer(vtkAbstractArray* array){
  return static_cast<PT*>(array->GetVoidPointer(0));
}

void cross(float* out, const float* v, const float* u){
  out[0] = v[1] * u[2] - v[2] * u[1];
  out[1] = v[2] * u[0] - v[0] * u[2];
  out[2] = v[0] * u[1] - v[1] * u[0];
}

void normalize(float* v){
  const float temp = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= temp;
  v[1] /= temp;
  v[2] /= temp;
}

int computeProjection(
  // output
  float *pointCoordinates,
  double *triangleDistortions,
  int *connectivityList,
  int *offsetArray,

  // input
  const float *depthValues,
  const float *camPos,
  const float *camDir,
  const float *camUp,
  const float *camNearFar,
  const float *camHeight,
  const float *resolution
){

  const float camDirMag = std::sqrt(
    camDir[0] * camDir[0] + camDir[1] * camDir[1] + camDir[2] * camDir[2]);
  const float camDirN[3]{
    camDir[0] / camDirMag, camDir[1] / camDirMag, camDir[2] / camDirMag};

  const size_t resolutionST[2] = {(size_t)resolution[0], (size_t)resolution[1]};

  // Compute camera size
  const float camSize[2]
    = {resolution[0] / resolution[1] * camHeight[0], camHeight[0]};

  // Compute camRight = camDirN x CamUp
  float camRight[3];
  cross(camRight,camDirN,camUp);
  normalize(camRight);

  // Compute true camUp = camDirN x -camRight
  float camUpTrue[3];
  float camLeft[3]{-camRight[0],-camRight[1],-camRight[2]};
  cross(camUpTrue,camDirN,camLeft);
  normalize(camUpTrue);

  // -------------------------------------------------------------------------
  // Create Vertices
  // -------------------------------------------------------------------------
  {
    // Compute pixel size in world coordinates
    const float pixelWidthWorld = camSize[0] / resolution[0];
    const float pixelHeightWorld = camSize[1] / resolution[1];

    // Optimization: precompute half of the camera size to reduce the number of
    // operations in the for loop Include a half pixel offset (-0.5) to center
    // vertices at pixel centers
    const float camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
    const float camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

    // Compute depth delta
    const float delta = camNearFar[1] - camNearFar[0];

    // Optimization: reorient camera model to bottom left corner to reduce
    // operations in for loop
    const float camPosCorner[3] = {
      camPos[0] - camRight[0] * camWidthWorldHalf - camUpTrue[0] * camHeightWorldHalf,
      camPos[1] - camRight[1] * camWidthWorldHalf - camUpTrue[1] * camHeightWorldHalf,
      camPos[2] - camRight[2] * camWidthWorldHalf - camUpTrue[2] * camHeightWorldHalf
    };

    // Compute vertex coordinates while parallelizing over rows
    for(size_t y = 0; y < resolutionST[1]; y++) {
      const float v = ((float)y) * pixelHeightWorld;
      const float vTimesUp[3]
        = {v * camUpTrue[0], v * camUpTrue[1], v * camUpTrue[2]};

      const size_t yOffset = y * resolutionST[0];
      for(size_t x = 0; x < resolutionST[0]; x++) {
        const size_t pixelIndex = x + yOffset;

        // double d = (double)(depthValues[ pixelIndex ])*delta+camNearFar[0];
        const float depth = ((float)depthValues[pixelIndex]);
        // double d = depth > 0.98 ? 0 : depth * delta + camNearFar[0];
        const float d = depth * delta + camNearFar[0];
        const float u = ((float)x) * pixelWidthWorld;

        // compute vertex coordinate
        const size_t pointCoordinateOffset = pixelIndex * 3;
        pointCoordinates[pointCoordinateOffset]
          = camPosCorner[0] + u * camRight[0] + vTimesUp[0] + d * camDirN[0];
        pointCoordinates[pointCoordinateOffset + 1]
          = camPosCorner[1] + u * camRight[1] + vTimesUp[1] + d * camDirN[1];
        pointCoordinates[pointCoordinateOffset + 2]
          = camPosCorner[2] + u * camRight[2] + vTimesUp[2] + d * camDirN[2];
      }
    }
  }

  // -------------------------------------------------------------------------
  // Create Triangles
  // -------------------------------------------------------------------------
  {
    auto absDiff = [](const float &a, const float &b) {
      return a > b ? a - b : b - a;
    };
    auto isNaN = [](const float &a) { return std::isnan(a) || a >= 1.0; };

    const size_t nTriangles = 2 * (resolution[0] - 1) * (resolution[1] - 1);

    for(size_t t = 0; t < nTriangles; t++)
      offsetArray[t] = t * 3;
    offsetArray[nTriangles] = 3 * nTriangles;

    /* Index Structure:
    0 - 1
    | / |
    2 - 3
    */
    const size_t xl = resolutionST[0] - 1;
    const size_t yl = resolutionST[1] - 1;
    const size_t trianglesPerRow = xl * 2;

    const float myNan = std::numeric_limits<float>::quiet_NaN();

    for(size_t y = 0; y < yl; y++) {
      size_t const yOffset = y * resolutionST[0];
      size_t triangleIndexOffset = y * trianglesPerRow * 3;
      size_t triangleDistortionOffset = y * trianglesPerRow;

      for(size_t x = 0; x < xl; x++) {
        size_t const i0 = x + yOffset;
        size_t const i1 = i0 + 1;
        size_t const i2 = i0 + resolutionST[0];
        size_t const i3 = i2 + 1;

        connectivityList[triangleIndexOffset++] = i0;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i1;

        connectivityList[triangleIndexOffset++] = i1;
        connectivityList[triangleIndexOffset++] = i2;
        connectivityList[triangleIndexOffset++] = i3;

        const float i0Depth = (float)depthValues[i0];
        const float i1Depth = (float)depthValues[i1];
        const float i2Depth = (float)depthValues[i2];
        const float i3Depth = (float)depthValues[i3];

        // wow
        triangleDistortions[triangleDistortionOffset++]
          = isNaN(i0Depth) || isNaN(i2Depth) || isNaN(i1Depth)
              ? myNan
              : std::max(
                absDiff(i0Depth, i1Depth),
                std::max(absDiff(i1Depth, i2Depth), absDiff(i0Depth, i2Depth)));

        triangleDistortions[triangleDistortionOffset++]
          = isNaN(i1Depth) || isNaN(i2Depth) || isNaN(i3Depth)
              ? myNan
              : std::max(
                absDiff(i1Depth, i3Depth),
                std::max(absDiff(i3Depth, i2Depth), absDiff(i2Depth, i1Depth)));
      }
    }
  }

  return 1;
}

//----------------------------------------------------------------------------
CinemaDepthImageProjection::~CinemaDepthImageProjection() = default;

int CinemaDepthImageProjection::RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector){

  auto input = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto output = vtkMultiBlockDataSet::GetData(outputVector);

  size_t nBlocks = input->GetNumberOfBlocks();

  #ifdef _OPENMP
    this->printMsg("#Depth Image Projection ("+std::to_string(nBlocks)+" images, "+std::to_string(omp_get_max_threads())+" threads)");
  #else
    this->printMsg("#Depth Image Projection ("+std::to_string(nBlocks)+" images)");
  #endif

  // initialize output
  for(size_t b=0; b<nBlocks; b++)
    output->SetBlock(b, vtkSmartPointer<vtkPolyData>::New());

  bool failed = false;

  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for(long long b=0; b<nBlocks; b++){
    auto image = vtkImageData::SafeDownCast(input->GetBlock(b));
    if(!image) continue;

    int dims[3];
    image->GetDimensions(dims);
    size_t resX = dims[0];
    size_t resY = dims[1];
    float res[2]{(float)resX,(float)resY};

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(resX * resY);

    // Prepare output cell buffer that holds two triangles for every quad
    // consisting of 4 vertices
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    const size_t nTriangles = 2 * (resX - 1) * (resY - 1);

    auto triangleDistortions = vtkSmartPointer<vtkDoubleArray>::New();
    triangleDistortions->SetName("TriangleDistortion");
    triangleDistortions->SetNumberOfComponents(1);
    triangleDistortions->SetNumberOfTuples(nTriangles);

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    connectivityArray->SetNumberOfTuples(3 * nTriangles);

    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    offsetArray->SetNumberOfTuples(nTriangles + 1);

    // Compute Projection
    auto imageFD = image->GetFieldData();
    auto imagePD = image->GetPointData();
    int status = computeProjection(
      // Output
      getPointer<float>(points->GetData()),
      getPointer<double>(triangleDistortions),
      getPointer<int>(connectivityArray),
      getPointer<int>(offsetArray),

      // Input
      getPointer<float>(imagePD->GetArray("Depth")),
      getPointer<float>(imageFD->GetArray("CameraPos")),
      getPointer<float>(imageFD->GetArray("CameraDir")),
      getPointer<float>(imageFD->GetArray("CameraUp")),
      getPointer<float>(imageFD->GetArray("CameraNearFar")),
      getPointer<float>(imageFD->GetArray("CameraHeight")),
      res
    );
    if(!status)
      failed = true;

    // Create VTK Output
    auto mesh = vtkPolyData::SafeDownCast(output->GetBlock(b));
    {
      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);
      mesh->SetPoints(points);
      mesh->SetPolys(cellArray);
      mesh->GetCellData()->AddArray(triangleDistortions);
    }

    // copy image point data to output point data
    {
      auto meshPD = mesh->GetPointData();
      for(int p = 0; p < imagePD->GetNumberOfArrays(); p++)
        meshPD->AddArray(imagePD->GetAbstractArray(p));
    }
  }

  return failed ? 0 : 1;
}
