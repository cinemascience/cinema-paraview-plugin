#include "vtkSMCinemaWriterProxy.h"

#include "PipelineAnnotation.h"

#include <vtkObjectFactory.h>
#include <vtkSMStringVectorProperty.h>

vtkStandardNewMacro(vtkSMCinemaWriterProxy);

void vtkSMCinemaWriterProxy::CaptureProvenance() {
  auto* property = vtkSMStringVectorProperty::SafeDownCast(this->GetProperty("Provenance"));
  if (!property) {
    vtkWarningMacro("CinemaWriter proxy has no StringVectorProperty named 'Provenance'.");
    return;
  }

  const AnnotationMap provenance = ComputePipelineAnnotations(this);
  property->SetNumberOfElements(static_cast<unsigned int>(2 * provenance.size()));

  unsigned int index = 0;
  for (const auto& item : provenance) {
    property->SetElement(index++, item.first.c_str());
    property->SetElement(index++, item.second.c_str());
  }
}

void vtkSMCinemaWriterProxy::UpdateVTKObjects() {
  // This is the reliable synchronization hook for both built-in and
  // client/server execution. Capture the current client-side pipeline state
  // before modified properties are pushed to the VTK object.
  if (!this->UpdatingVTKObjects) {
    this->UpdatingVTKObjects = true;
    this->CaptureProvenance();
    this->Superclass::UpdateVTKObjects();
    this->UpdatingVTKObjects = false;
    return;
  }

  this->Superclass::UpdateVTKObjects();
}
