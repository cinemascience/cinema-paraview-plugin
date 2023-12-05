#pragma once

#include <pcAlgorithm.h>
#include <pcImageCompositingModule.h>

class vtkMultiProcessController;

class PCIMAGECOMPOSITING_EXPORT pcImageCompositing : public pcAlgorithm {

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

  static pcImageCompositing *New();
  vtkTypeMacro(pcImageCompositing, pcAlgorithm);

protected:
  pcImageCompositing();
  ~pcImageCompositing() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
