#include "CinemaImaging.h"

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

#include <embree3/rtcore.h>

#ifdef _OPENMP
#include <omp.h>
#endif

vtkStandardNewMacro(CinemaImaging);

const unsigned int INVALID_ID{std::numeric_limits<unsigned int>::max()};

template<typename PT>
PT* getPointer(vtkAbstractArray* array){
  return static_cast<PT*>(array->GetVoidPointer(0));
}

//----------------------------------------------------------------------------
CinemaImaging::CinemaImaging(){
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
};

int CinemaImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int CinemaImaging::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
CinemaImaging::~CinemaImaging() = default;

int CinemaImaging::PrepareCameraGrid(vtkPointSet* cameraGrid) const {
  return 1;
};

template <typename IT>
int initializeScene(
  RTCScene &scene,

  const RTCDevice &device,
  const size_t &nVertices,
  const float *vertexCoords,
  const size_t &nTriangles,
  const IT *connectivityList) {

  scene = rtcNewScene(device);

  RTCGeometry mesh = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

  // vertices
  {
    rtcSetSharedGeometryBuffer(mesh, RTC_BUFFER_TYPE_VERTEX, 0,
                               RTC_FORMAT_FLOAT3, (const void *)vertexCoords, 0,
                               3 * sizeof(float), nVertices);
  }

  // triangles
  {
    unsigned int *indices = (unsigned int *)rtcSetNewGeometryBuffer(
      mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
      3 * sizeof(unsigned int), nTriangles);

    for(size_t t = 0, tn = nTriangles * 3; t < tn; t++)
      indices[t] = (unsigned int)connectivityList[t];
  }

  rtcCommitGeometry(mesh);
  rtcAttachGeometry(scene, mesh);
  rtcReleaseGeometry(mesh);
  rtcCommitScene(scene);

  return 1;
};

int deallocateScene(RTCDevice &device, RTCScene &scene) {
  rtcReleaseScene(scene);
  rtcReleaseDevice(device);
  return 1;
};

int initializeDevice(RTCDevice &device) {
  // device = rtcNewDevice(NULL);
  device = rtcNewDevice("hugepages=1,threads=1");
  if(!device){
    printf("error %s\n", std::to_string(rtcGetDeviceError(nullptr)).c_str());
    return 0;
  }

  auto errorFunction
    = [](void *userPtr, enum RTCError error, const char *str) {
        printf("error %d: %s\n", error, str);
      };
  rtcSetDeviceErrorFunction(device, errorFunction, nullptr);
  return 1;
};

vtkCellArray* getCells(vtkPointSet* pointSet) {
  switch(pointSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
      return static_cast<vtkUnstructuredGrid *>(pointSet)->GetCells();
    case VTK_POLY_DATA:
      return static_cast<vtkPolyData *>(pointSet)->GetPolys();
  }
  return nullptr;
};

int renderImage(
  float *depthBuffer,
  unsigned int *primitiveIds,
  float *barycentricCoordinates,

  const RTCScene &scene,
  const int resolution[2],
  const float camPos[3],
  const float camDirRaw[3],
  const float camUp[3],
  const float &camHeight,
  const float camNearFar[2]
){
  const int resX = resolution[0];
  const int resY = resolution[1];

  // embree 3
  struct RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  // embree 4
  // struct RTCRayQueryContext context;
  // rtcInitRayQueryContext(&context);

  const auto normalize = [](float out[3], const float in[3]) {
    const float temp = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
    out[0] = in[0] / temp;
    out[1] = in[1] / temp;
    out[2] = in[2] / temp;
  };

  float camDir[3]{0, 0, 0};
  normalize(camDir, camDirRaw);

  // Compute camRight = camDir x CamUp
  float camRight[3]{camDir[1] * camUp[2] - camDir[2] * camUp[1],
                     camDir[2] * camUp[0] - camDir[0] * camUp[2],
                     camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  normalize(camRight, camRight);

  // Compute true up std::vector
  float camUpTrue[3]{camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
                      camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
                      camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  normalize(camUpTrue, camUpTrue);


  size_t pixelIndex = 0;
  size_t bcIndex = 0;
  const float nan = std::numeric_limits<float>::quiet_NaN();

  // Compute camera size
  const float resXF = static_cast<float>(resX);
  const float resYF = static_cast<float>(resY);
  const float aspect = resXF / resYF;
  const float camSize[2] = {aspect * camHeight, camHeight};

  // Compute pixel size in world coordinates
  const float pixelWidthWorld = camSize[0] / resXF;
  const float pixelHeightWorld = camSize[1] / resYF;

  // Optimization: precompute half of the camera size to reduce the number of
  // operations in the for loop. Include a half pixel offset (-0.5) to center
  // vertices at pixel centers
  const float camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
  const float camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

  // Optimization: reorient camera model to bottom left corner to reduce
  // operations in for loop
  const float camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                                    - camUpTrue[0] * camHeightWorldHalf,
                                  camPos[1] - camRight[1] * camWidthWorldHalf
                                    - camUpTrue[1] * camHeightWorldHalf,
                                  camPos[2] - camRight[2] * camWidthWorldHalf
                                    - camUpTrue[2] * camHeightWorldHalf};

  const float nearFarDelta = camNearFar[1]-camNearFar[0];

  for(int y = 0; y < resY; y++) {
    const float v = ((float)y) * pixelHeightWorld;

    for(int x = 0; x < resX; x++) {
      const float u = ((float)x) * pixelWidthWorld;

      struct RTCRayHit rayhit;

      // set origin
      rayhit.ray.org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
      rayhit.ray.org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
      rayhit.ray.org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];
      // std::cout<<"("<<rayhit.ray.org_x<<","<<rayhit.ray.org_y<<","<<rayhit.ray.org_z<<")"<<std::endl;

      // set dir
      rayhit.ray.dir_x = camDir[0];
      rayhit.ray.dir_y = camDir[1];
      rayhit.ray.dir_z = camDir[2];
      // std::cout<<"("<<rayhit.ray.dir_x<<","<<rayhit.ray.dir_y<<","<<rayhit.ray.dir_z<<")"<<std::endl;

      // compute hit
      rayhit.ray.tnear = 0;
      rayhit.ray.tfar = std::numeric_limits<float>::infinity();
      rayhit.ray.mask = -1;
      rayhit.ray.flags = 0;
      rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
      rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;

      // embree 3
      rtcIntersect1(scene, &context, &rayhit);

      // embree 4
      // rtcIntersect1(scene, &rayhit);

      // write depth
      const bool hitPrimitive = rayhit.hit.geomID != RTC_INVALID_GEOMETRY_ID;
      if(hitPrimitive) {
        depthBuffer[pixelIndex] = std::max(0.0f,std::min(1.0f,(rayhit.ray.tfar-camNearFar[0])/nearFarDelta));
        primitiveIds[pixelIndex] = rayhit.hit.primID;
        barycentricCoordinates[bcIndex++] = rayhit.hit.u;
        barycentricCoordinates[bcIndex++] = rayhit.hit.v;
      } else {
        depthBuffer[pixelIndex] = 1.0;
        primitiveIds[pixelIndex] = -1;
        barycentricCoordinates[bcIndex++] = nan;
        barycentricCoordinates[bcIndex++] = nan;
      }
      pixelIndex++;
    }
  }

  return 1;
};

template <typename DT, typename IT>
int interpolateArray(
  float *outputArray,
  const unsigned int *primitiveIds,
  const float *barycentricCoordinates,
  const IT *connectivityList,

  const DT *inputArray,
  const size_t &nTuples,
  const size_t &nComponents,
  const float missingValue = std::numeric_limits<float>::quiet_NaN()
){
  if(nComponents==1){
    for(size_t i = 0; i < nTuples; i++) {
      const unsigned int &cellId = primitiveIds[i];
      if(cellId == INVALID_ID) {
        outputArray[i] = missingValue;
        continue;
      }

      const size_t cellIndex = cellId * 3;
      const IT &v0 = connectivityList[cellIndex + 0];
      const IT &v1 = connectivityList[cellIndex + 1];
      const IT &v2 = connectivityList[cellIndex + 2];

      const size_t bcIndex = i * 2;
      const float &u = barycentricCoordinates[bcIndex + 0];
      const float &v = barycentricCoordinates[bcIndex + 1];
      const float w = 1 - u - v;

      outputArray[i]
        = static_cast<float>(w * inputArray[v0] + u * inputArray[v1] + v * inputArray[v2]);
    }
  } else {
    for(size_t i = 0; i < nTuples; i++) {
      const unsigned int &cellId = primitiveIds[i];
      if(cellId == INVALID_ID) {
        for(size_t c=0; c<nComponents; c++)
          outputArray[i*nComponents+c] = missingValue;
        continue;
      }

      const size_t cellIndex = cellId * 3;
      const IT &v0 = connectivityList[cellIndex + 0];
      const IT &v1 = connectivityList[cellIndex + 1];
      const IT &v2 = connectivityList[cellIndex + 2];

      const size_t bcIndex = i * 2;
      const float &u = barycentricCoordinates[bcIndex + 0];
      const float &v = barycentricCoordinates[bcIndex + 1];
      const float w = 1 - u - v;

      for(size_t c=0; c<nComponents; c++)
        outputArray[i*nComponents+c]
          = static_cast<float>(w * inputArray[v0] + u * inputArray[v1] + v * inputArray[v2]);
    }
  }

  return 1;
};

template <typename DT>
int lookupArray(
  float *outputArray,
  const unsigned int *primitiveIds,
  const DT *inputArray,
  const size_t &nTuples,
  const size_t &nComponents,
  const float missingValue = std::numeric_limits<float>::quiet_NaN()
){
  size_t offset = 0;
  if(nComponents==1){
    for(size_t i = 0; i < nTuples; i++) {
      const auto& cellId = primitiveIds[i];
      outputArray[offset++] = cellId==INVALID_ID ? missingValue : inputArray[offset];
    }
  } else {
    for(size_t i = 0; i < nTuples; i++) {
      const auto& cellId = primitiveIds[i];

      if(cellId==INVALID_ID) {
        for(size_t c=0; c<nComponents; c++)
          outputArray[offset++] = missingValue;
      } else {
        for(size_t c=0; c<nComponents; c++)
          outputArray[offset++] = inputArray[offset];
      }
    }
  }

  return 1;
};

int MapPointAndCellData(
  vtkImageData *image,

  vtkPointSet *inputObject,
  const unsigned int *primitiveIdArray,
  const float *barycentricCoordinates,
  const vtkIdType *inputObjectConnectivityList
){
  auto inputObjectPD = inputObject->GetPointData();
  auto inputObjectCD = inputObject->GetCellData();
  auto imagePD = image->GetPointData();
  int dim[3];
  image->GetDimensions(dim);
  size_t const nPixels = dim[0] * dim[1];

  const size_t nInputObjectPDArrays = inputObjectPD->GetNumberOfArrays();
  const size_t nInputObjectCDArrays = inputObjectCD->GetNumberOfArrays();

  int status = 0;

  // Map Point Data
  for(size_t j = 0; j < nInputObjectPDArrays; j++) {
    auto inputArray = inputObjectPD->GetArray(j);
    if(
      !inputArray
      || std::string(inputArray->GetName()).compare("vtkGhostType")==0
    ) continue;

    // auto outputArray = vtkSmartPointer<vtkAbstractArray>::Take(inputArray->NewInstance());
    auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    imagePD->AddArray(outputArray);

    switch(inputArray->GetDataType()) {
      vtkTemplateMacro(status = interpolateArray(
                        getPointer<float>(outputArray),

                        primitiveIdArray,
                        barycentricCoordinates,
                        inputObjectConnectivityList,

                        getPointer<VTK_TT>(inputArray),
                        nPixels, inputArray->GetNumberOfComponents()));
    }
  }

  // Map Cell Data
  for(size_t j = 0; j < nInputObjectCDArrays; j++) {
    auto inputArray = inputObjectPD->GetArray(j);
    if(
      !inputArray
      || std::string(inputArray->GetName()).compare("vtkGhostType")==0
    ) continue;

    auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    imagePD->AddArray(outputArray);

    switch(inputArray->GetDataType()) {
      vtkTemplateMacro(status = lookupArray(
                         getPointer<float>(outputArray),

                         primitiveIdArray,
                         getPointer<VTK_TT>(inputArray),
                         nPixels, inputArray->GetNumberOfComponents()));
    }
  }

  return 1;
};

int AddFieldDataArray(vtkFieldData *fd, vtkDataArray *array, int tupleIdx, const std::string name="") {
  if(!array)
    return 0;

  size_t nComponents = array->GetNumberOfComponents();

  vtkNew<vtkFloatArray> newArray;
  newArray->SetName(name.empty() ? array->GetName() : name.data());
  newArray->SetNumberOfComponents(nComponents);
  newArray->SetNumberOfTuples(1);

  if(newArray->GetDataType() == array->GetDataType()) {
    newArray->SetTuple(0, tupleIdx, array);
  } else {
    for(size_t i = 0; i < nComponents; i++)
      newArray->SetValue(
        i, array->GetVariantValue(tupleIdx * nComponents + i).ToDouble());
  }

  fd->AddArray(newArray);

  return 1;
};

int AddAllFieldDataArrays(vtkPointSet *cameras, vtkImageData *image, int tupleIdx) {
  auto imageFD = image->GetFieldData();

  auto camerasPD = cameras->GetPointData();
  for(int i = 0; i < camerasPD->GetNumberOfArrays(); i++) {
    AddFieldDataArray(imageFD, camerasPD->GetArray(i), tupleIdx);
  }
  AddFieldDataArray(imageFD, cameras->GetPoints()->GetData(), tupleIdx, "CameraPos");

  return 1;
};

int CinemaImaging::RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector){

  auto inputObject = vtkPointSet::GetData(inputVector[0]);
  auto cameras = vtkPointSet::GetData(inputVector[1]);
  if(!inputObject || !cameras) return 0;

  double bounds[6];
  inputObject->GetBounds(bounds);
  const double d[3]{
    bounds[1]-bounds[0],
    bounds[3]-bounds[2],
    bounds[5]-bounds[4]
  };
  const double diameter = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);

  const int nCameras = cameras->GetNumberOfPoints();
  auto inputObjectCells = getCells(inputObject);
  if(!inputObjectCells) return 0;
  auto inputObjectConnectivityList = getPointer<vtkIdType>(inputObjectCells->GetConnectivityArray());

  // init device and scene
  RTCDevice device;
  if(!initializeDevice(device)) return 0;

  RTCScene scene;
  if(!initializeScene<vtkIdType>(
      scene,
      device,
      inputObject->GetNumberOfPoints(),
      getPointer<float>(inputObject->GetPoints()->GetData()),
      inputObjectCells->GetNumberOfCells(), inputObjectConnectivityList
  )) return 0;

  // Initialize Output
  const auto camPos = getPointer<float>(cameras->GetPoints()->GetData());
  const auto camUp = getPointer<float>(cameras->GetPointData()->GetArray("CameraUp"));
  const auto camDir = getPointer<float>(cameras->GetPointData()->GetArray("CameraDir"));
  const auto camHeight = getPointer<float>(cameras->GetPointData()->GetArray("CameraHeight"));
  const auto camNearFar = getPointer<float>(cameras->GetPointData()->GetArray("CameraNearFar"));
  auto outputCollection = vtkMultiBlockDataSet::GetData(outputVector);
  for(int c=0; c<nCameras; c++){
    auto image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions(this->Resolution[0], this->Resolution[1], 1);
    image->SetSpacing(1, 1, 1);
    image->SetOrigin(0, 0, 0);
    image->AllocateScalars(VTK_FLOAT, 1);
    outputCollection->SetBlock(c, image);
  }

#ifdef _OPENMP
  this->printMsg("#Imaging ("+std::to_string(nCameras)+" images, "+std::to_string(omp_get_max_threads())+" threads)");
  #pragma omp parallel for
#endif
  for(int c=0; c<nCameras; c++){
    auto image = static_cast<vtkImageData*>(outputCollection->GetBlock(c));

    size_t const nPixels = this->Resolution[0] * this->Resolution[1];
    auto imagePD = image->GetPointData();

    auto depthBuffer = imagePD->GetArray(0);
    depthBuffer->SetName("Depth");
    auto depthBuffer_ = getPointer<float>(depthBuffer);

    auto primitiveId = vtkSmartPointer<vtkUnsignedIntArray>::New();
    primitiveId->SetName("PrimitiveId");
    primitiveId->SetNumberOfComponents(1);
    primitiveId->SetNumberOfTuples(nPixels);
    auto primitiveId_ = getPointer<unsigned int>(primitiveId);
    // imagePD->AddArray(primitiveId);

    auto barycentricCoordinates = vtkSmartPointer<vtkFloatArray>::New();
    barycentricCoordinates->SetName("BarycentricCoordinates");
    barycentricCoordinates->SetNumberOfComponents(2);
    barycentricCoordinates->SetNumberOfTuples(nPixels);
    auto barycentricCoordinates_ = getPointer<float>(barycentricCoordinates);
    // imagePD->AddArray(barycentricCoordinates);

    renderImage(
      depthBuffer_,
      primitiveId_,
      barycentricCoordinates_,

      scene,
      this->Resolution,
      &camPos[c*3],
      &camDir[c*3],
      &camUp[c*3],
      camHeight[c],
      &camNearFar[c*2]
    );

    MapPointAndCellData(
      image,

      inputObject, primitiveId_, barycentricCoordinates_,
      inputObjectConnectivityList
    );

    AddAllFieldDataArrays( cameras, image, c );
  }

  if(!deallocateScene(device,scene)) return 0;

  return 1;
};
