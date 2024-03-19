#pragma once

#include <CinemaAlgorithm.h>
#include <CinemaImagingModule.h>

class vtkPointSet;

class CINEMAIMAGING_EXPORT CinemaImaging : public CinemaAlgorithm
{
public:
  static CinemaImaging* New();
  vtkTypeMacro(CinemaImaging, CinemaAlgorithm);

  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

protected:
  CinemaImaging();
  ~CinemaImaging();

  int RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int PrepareCameraGrid(vtkPointSet* cameraGrid) const;

private:
  CinemaImaging(const CinemaImaging&) = delete;
  void operator=(const CinemaImaging&) = delete;

  int Resolution[2]{256, 256};
};
