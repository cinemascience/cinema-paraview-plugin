#pragma once

#include <map>
#include <string>

class vtkSMProxy;

using AnnotationMap = std::map<std::string, std::string>;

// Collects the current ServerManager state of all proxies upstream of the
// supplied writer proxy. This function is client-side only: it follows the
// writer's Input property and records persistent user-facing properties as
// flat columns such as "Contour0.Isosurfaces" and
// "Clip0.ClipType0.Origin".
AnnotationMap ComputePipelineAnnotations(vtkSMProxy* writerProxy);
