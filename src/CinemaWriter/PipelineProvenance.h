#pragma once

#include <map>
#include <string>
#include <vector>

class vtkSMProxy;
class vtkSMViewProxy;

using ProvenanceMap = std::map<std::string, std::string>;

// Return the proxies connected to a proxy property such as "Input".
std::vector<vtkSMProxy*> GetPipelineInputs(vtkSMProxy* proxy, const char* propertyName = "Input");

// Return the data-producing proxies feeding visible data representations in a view.
std::vector<vtkSMProxy*> GetVisibleRepresentationInputs(vtkSMViewProxy* view);

// Collect persistent user-facing ServerManager state from one or more pipeline
// roots and their upstream inputs. All roots share one visited set and one set
// of type counters, so common ancestors and naming are unified deterministically.
ProvenanceMap ComputePipelineProvenance(vtkSMProxy* root);
ProvenanceMap ComputePipelineProvenance(const std::vector<vtkSMProxy*>& roots);

// Collect the unified provenance of all visible data representations in a view.
// Additional roots may be supplied for state that contributes to the extract
// but is not represented in the view, e.g. a camera-producing pipeline.
ProvenanceMap ComputeVisibleRepresentationProvenance(vtkSMViewProxy* view,
                                                     const std::vector<vtkSMProxy*>& additionalRoots = {});

// Merge source into destination. Existing destination entries win by default.
void MergeProvenance(ProvenanceMap& destination, const ProvenanceMap& source, bool overwriteExisting = false);
