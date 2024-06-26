#pragma once

#include <CinemaAlgorithm.h>
#include <CinemaDepthImageProjectionModule.h>

class CINEMADEPTHIMAGEPROJECTION_EXPORT CinemaDepthImageProjection : public CinemaAlgorithm
{

public:
  static CinemaDepthImageProjection* New();
  vtkTypeMacro(CinemaDepthImageProjection, CinemaAlgorithm);

protected:
  CinemaDepthImageProjection();
  ~CinemaDepthImageProjection();

  int RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

private:
  CinemaDepthImageProjection(const CinemaDepthImageProjection&) = delete;
  void operator=(const CinemaDepthImageProjection&) = delete;

};
