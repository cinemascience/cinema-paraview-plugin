#pragma once

#include <pcAlgorithm.h>
#include <pcImagingModule.h>

class vtkPointSet;

class PCIMAGING_EXPORT pcImaging : public pcAlgorithm
{
public:
  static pcImaging* New();
  vtkTypeMacro(pcImaging, pcAlgorithm);

  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

protected:
  pcImaging();
  ~pcImaging();

  int RequestData(vtkInformation *request,
                 vtkInformationVector **inputVector,
                 vtkInformationVector *outputVector) override;
  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int PrepareCameraGrid(vtkPointSet* cameraGrid) const;

private:
  pcImaging(const pcImaging&) = delete;
  void operator=(const pcImaging&) = delete;

  int Resolution[2]{256, 256};
};
