#pragma once

#include "CinemaWriterModule.h"

#include <vtkSMSourceProxy.h>

class CINEMAWRITER_EXPORT vtkSMCinemaWriterProxy : public vtkSMSourceProxy {
public:
  static vtkSMCinemaWriterProxy* New();
  vtkTypeMacro(vtkSMCinemaWriterProxy, vtkSMSourceProxy);

  void UpdateVTKObjects() override;

protected:
  vtkSMCinemaWriterProxy() = default;
  ~vtkSMCinemaWriterProxy() override = default;

private:
  vtkSMCinemaWriterProxy(const vtkSMCinemaWriterProxy&) = delete;
  void operator=(const vtkSMCinemaWriterProxy&) = delete;

  void CaptureProvenance();
  bool UpdatingVTKObjects = false;
};
