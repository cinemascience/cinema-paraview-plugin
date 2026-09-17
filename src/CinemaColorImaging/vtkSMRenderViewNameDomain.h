#pragma once

#include "CinemaColorImagingModule.h"

#include <vtkSMStringListDomain.h>

class vtkSMProperty;

class CINEMACOLORIMAGING_EXPORT vtkSMRenderViewNameDomain : public vtkSMStringListDomain {
public:
  static vtkSMRenderViewNameDomain* New();
  vtkTypeMacro(vtkSMRenderViewNameDomain, vtkSMStringListDomain);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void Update(vtkSMProperty* property) override;
  int SetDefaultValues(vtkSMProperty* property, bool useUncheckedValues) override;

protected:
  vtkSMRenderViewNameDomain() = default;
  ~vtkSMRenderViewNameDomain() override = default;

private:
  vtkSMRenderViewNameDomain(const vtkSMRenderViewNameDomain&) = delete;
  void operator=(const vtkSMRenderViewNameDomain&) = delete;
};
