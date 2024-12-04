#pragma once

#include <CinemaAlgorithm.h>
#include <CinemaColorImagingModule.h>

class vtkPointSet;

class CINEMACOLORIMAGING_EXPORT CinemaColorImaging : public CinemaAlgorithm
{
public:
  static CinemaColorImaging* New();
  vtkTypeMacro(CinemaColorImaging, CinemaAlgorithm);

  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

  bool GetNeedsUpdate();

protected:
  CinemaColorImaging();
  ~CinemaColorImaging();

  int RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  CinemaColorImaging(const CinemaColorImaging&) = delete;
  void operator=(const CinemaColorImaging&) = delete;

  int Resolution[2]{256, 256};
  bool NeedsUpdate{false};
  bool Ready{false};
};
