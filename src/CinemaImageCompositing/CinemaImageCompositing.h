#pragma once

#include <CinemaAlgorithm.h>
#include <CinemaImageCompositingModule.h>

class vtkMultiProcessController;

class CINEMAIMAGECOMPOSITING_EXPORT CinemaImageCompositing : public CinemaAlgorithm {

private:
  vtkMultiProcessController* Controller{nullptr};

public:

  ///@{
  /**
   * Set and get the controller.
   */
  void SetController(vtkMultiProcessController*);
  vtkGetObjectMacro(Controller, vtkMultiProcessController);
  ///@}

  static CinemaImageCompositing *New();
  vtkTypeMacro(CinemaImageCompositing, CinemaAlgorithm);

protected:
  CinemaImageCompositing();
  ~CinemaImageCompositing() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
