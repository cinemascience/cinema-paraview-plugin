#pragma once

#include <CinemaAlgorithm.h>
#include <CinemaCameraGridModule.h>

class vtkMultiProcessController;

class CINEMACAMERAGRID_EXPORT CinemaCameraGrid : public CinemaAlgorithm
{

public:
  static CinemaCameraGrid* New();
  vtkTypeMacro(CinemaCameraGrid, CinemaAlgorithm);

  enum class AXIS {
    X = 0,
    Y = 1,
    Z = 2
  };

private:
  CinemaCameraGrid(const CinemaCameraGrid&) = delete;
  void operator=(const CinemaCameraGrid&) = delete;

  double RadiusFactor{1};
  double Center[3]{0,0,0};
  int ThetaResolution{8};
  int PhiResolution{3};
  double StartTheta{0};
  double EndTheta{360};
  double StartPhi{45};
  double EndPhi{135};
  double CamHeight{0};
  double NearFar[2]{0,0};
  AXIS Axis{AXIS::Z};

  vtkMultiProcessController* Controller{nullptr};

public:

  ///@{
  /**
   * Set and get the controller.
   */
  void SetController(vtkMultiProcessController*);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  ///@}

  pcSetEnumMacro(Axis, AXIS);
  vtkGetEnumMacro(Axis, AXIS);

  vtkSetMacro(CamHeight, double);
  vtkGetMacro(CamHeight, double);

  vtkSetVector2Macro(NearFar, double);
  vtkGetVector2Macro(NearFar, double);

  ///@{
  /**
   * Set the radius of sphere. Default is 0.5.
   */
  vtkSetClampMacro(RadiusFactor, double, 0.0, VTK_DOUBLE_MAX);
  vtkGetMacro(RadiusFactor, double);
  ///@}

  ///@{
  /**
   * Set the center of the sphere. Default is (0,0,0).
   */
  vtkSetVector3Macro(Center, double);
  vtkGetVectorMacro(Center, double, 3);
  ///@}

  ///@{
  /**
   * Set the number of points in the longitude direction (ranging from
   * StartTheta to EndTheta).
   */
  vtkSetClampMacro(ThetaResolution, int, 3, VTK_INT_MAX);
  vtkGetMacro(ThetaResolution, int);
  ///@}

  ///@{
  /**
   * Set the number of points in the latitude direction (ranging
   * from StartPhi to EndPhi).
   */
  vtkSetClampMacro(PhiResolution, int, 3, VTK_INT_MAX);
  vtkGetMacro(PhiResolution, int);
  ///@}

  ///@{
  /**
   * Set the starting longitude angle. By default StartTheta=0 degrees.
   */
  vtkSetClampMacro(StartTheta, double, 0.0, 360.0);
  vtkGetMacro(StartTheta, double);
  ///@}

  ///@{
  /**
   * Set the ending longitude angle. By default EndTheta=360 degrees.
   */
  vtkSetClampMacro(EndTheta, double, 0.0, 360.0);
  vtkGetMacro(EndTheta, double);
  ///@}

  ///@{
  /**
   * Set the starting latitude angle (0 is at north pole). By default
   * StartPhi=0 degrees.
   */
  vtkSetClampMacro(StartPhi, double, -90.0, 90.0);
  vtkGetMacro(StartPhi, double);
  ///@}

  ///@{
  /**
   * Set the ending latitude angle. By default EndPhi=180 degrees.
   */
  vtkSetClampMacro(EndPhi, double, -90.0, 90.0);
  vtkGetMacro(EndPhi, double);
  ///@}

protected:
  CinemaCameraGrid();
  ~CinemaCameraGrid();

  int RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
};
