#pragma once

#include "CinemaColorImagingModule.h"

#include <vtkSMExtractWriterProxy.h>

class vtkSMExtractsController;
class vtkSMProxy;

class CINEMACOLORIMAGING_EXPORT vtkSMCinemaColorImagingExtractWriterProxy : public vtkSMExtractWriterProxy {
public:
  static vtkSMCinemaColorImagingExtractWriterProxy* New();
  vtkTypeMacro(vtkSMCinemaColorImagingExtractWriterProxy, vtkSMExtractWriterProxy);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  bool Write(vtkSMExtractsController* extractor) override;
  bool CanExtract(vtkSMProxy* proxy) override;
  bool IsExtracting(vtkSMProxy* proxy) override;
  void SetInput(vtkSMProxy* proxy) override;
  vtkSMProxy* GetInput() override;

protected:
  vtkSMCinemaColorImagingExtractWriterProxy() = default;
  ~vtkSMCinemaColorImagingExtractWriterProxy() override = default;

private:
  vtkSMCinemaColorImagingExtractWriterProxy(const vtkSMCinemaColorImagingExtractWriterProxy&) = delete;
  void operator=(const vtkSMCinemaColorImagingExtractWriterProxy&) = delete;
};
