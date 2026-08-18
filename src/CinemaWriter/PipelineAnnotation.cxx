#include "PipelineAnnotation.h"

#include <sstream>
#include <string>
#include <vtkNew.h>
#include <vtkObjectBase.h>
#include <vtkSMDoubleVectorProperty.h>
#include <vtkSMInputProperty.h>
#include <vtkSMIntVectorProperty.h>
#include <vtkSMProperty.h>
#include <vtkSMPropertyIterator.h>
#include <vtkSMProxy.h>
#include <vtkSMProxyIterator.h>
#include <vtkSMProxyManager.h>
#include <vtkSMProxyProperty.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkSMSourceProxy.h>
#include <vtkSMStringVectorProperty.h>

namespace {

// Forward declaration because DumpProperty() may recurse into DumpProxy().
void DumpProxy(vtkSMProxy* proxy, int indent, bool followPipelineInputs, std::ostringstream& out);

// -----------------------------------------------------------------------------
// Dump one ServerManager property.
//
// Supported for now:
//   vtkSMIntVectorProperty
//   vtkSMDoubleVectorProperty
//   vtkSMStringVectorProperty
//   vtkSMProxyProperty
//
// vtkSMInputProperty is handled separately in DumpProxy(), since it represents
// pipeline topology rather than an ordinary configurable property.
// -----------------------------------------------------------------------------
void DumpProperty(const char* name, vtkSMProperty* property, int indent, std::ostringstream& out) {
  if (!name || !property) { return; }

  const std::string pad(indent, ' ');

  // ---------------------------------------------------------------------------
  // Integer vector property
  // ---------------------------------------------------------------------------
  if (auto* p = vtkSMIntVectorProperty::SafeDownCast(property)) {
    out << pad << name << " = ";

    for (unsigned int i = 0; i < p->GetNumberOfElements(); ++i) {
      if (i > 0) { out << ", "; }

      out << p->GetElement(i);
    }

    out << '\n';
    return;
  }

  // ---------------------------------------------------------------------------
  // Double vector property
  // ---------------------------------------------------------------------------
  if (auto* p = vtkSMDoubleVectorProperty::SafeDownCast(property)) {
    out << pad << name << " = ";

    for (unsigned int i = 0; i < p->GetNumberOfElements(); ++i) {
      if (i > 0) { out << ", "; }

      out << p->GetElement(i);
    }

    out << '\n';
    return;
  }

  // ---------------------------------------------------------------------------
  // String vector property
  // ---------------------------------------------------------------------------
  if (auto* p = vtkSMStringVectorProperty::SafeDownCast(property)) {
    out << pad << name << " = ";

    for (unsigned int i = 0; i < p->GetNumberOfElements(); ++i) {
      if (i > 0) { out << ", "; }

      const char* value = p->GetElement(i);

      out << "\"" << (value ? value : "") << "\"";
    }

    out << '\n';
    return;
  }

  // ---------------------------------------------------------------------------
  // Nested proxy property
  //
  // Examples:
  //
  //   ClipFunction -> Plane
  //   Locator      -> MergePoints
  //
  // These are NOT pipeline connections, so when recursively dumping the
  // contained proxy we explicitly disable traversal through Input properties.
  // ---------------------------------------------------------------------------
  if (auto* p = vtkSMProxyProperty::SafeDownCast(property)) {
    out << pad << name << " =" << '\n';

    for (unsigned int i = 0; i < p->GetNumberOfProxies(); ++i) {
      vtkSMProxy* subProxy = p->GetProxy(i);

      if (!subProxy) { continue; }

      DumpProxy(subProxy, indent + 2, false, out);
    }

    return;
  }

  // ---------------------------------------------------------------------------
  // Unsupported property type.
  //
  // Keep this visible during development so that we notice property classes
  // which should eventually be serialized as well.
  // ---------------------------------------------------------------------------
  out << pad << name << " = <unsupported " << property->GetClassName() << ">" << '\n';
}

// -----------------------------------------------------------------------------
// Dump a ServerManager proxy.
//
// followPipelineInputs:
//
//   true:
//       This proxy is part of the actual upstream pipeline. Its
//       vtkSMInputProperty entries are followed recursively.
//
//   false:
//       This proxy is a nested helper/property proxy such as Plane or Locator.
//       Its Input property is deliberately ignored to avoid accidentally
//       jumping back into the main pipeline.
// -----------------------------------------------------------------------------
void DumpProxy(vtkSMProxy* proxy, int indent, bool followPipelineInputs, std::ostringstream& out) {
  if (!proxy) { return; }

  const std::string pad(indent, ' ');

  const char* xmlGroup = proxy->GetXMLGroup();
  const char* xmlName = proxy->GetXMLName();

  out << pad << (xmlGroup ? xmlGroup : "(null)") << " / " << (xmlName ? xmlName : "(null)") << '\n';

  vtkSMPropertyIterator* piter = proxy->NewPropertyIterator();

  if (!piter) { return; }

  for (piter->Begin(); !piter->IsAtEnd(); piter->Next()) {
    const char* key = piter->GetKey();

    vtkSMProperty* property = piter->GetProperty();

    if (!property) { continue; }

    // -------------------------------------------------------------------------
    // Pipeline topology
    // -------------------------------------------------------------------------
    if (auto* input = vtkSMInputProperty::SafeDownCast(property)) {
      // Helper proxies sometimes expose Input properties as well. Those must
      // not be interpreted as pipeline traversal.
      if (!followPipelineInputs) { continue; }

      out << pad << "  " << (key ? key : "Input") << " =" << '\n';

      for (unsigned int i = 0; i < input->GetNumberOfProxies(); ++i) {
        vtkSMProxy* upstream = input->GetProxy(i);

        if (!upstream) { continue; }

        DumpProxy(upstream, indent + 4, true, out);
      }

      continue;
    }

    // -------------------------------------------------------------------------
    // Ordinary property or nested proxy property.
    // -------------------------------------------------------------------------
    DumpProperty(key ? key : "(unnamed)", property, indent + 2, out);
  }

  piter->Delete();
}

} // namespace

// -----------------------------------------------------------------------------
// Public API.
//
// clientObject is normally:
//
//     this
//
// from CinemaWriter::RequestData().
//
// The function:
//
//   1. gets the active ParaView ServerManager session,
//   2. finds the vtkSMSourceProxy whose client-side object equals clientObject,
//   3. gets that proxy's "Input" property,
//   4. recursively dumps all upstream pipeline proxies and their properties,
//   5. returns the result as one string.
// -----------------------------------------------------------------------------
std::string ComputeInputTree(vtkObjectBase* clientObject) {
  if (!clientObject) { return {}; }

  // ---------------------------------------------------------------------------
  // Get the active ServerManager proxy manager.
  // ---------------------------------------------------------------------------
  vtkSMProxyManager* pxm = vtkSMProxyManager::GetProxyManager();

  if (!pxm) { return {}; }

  vtkSMSessionProxyManager* spxm = pxm->GetActiveSessionProxyManager();

  if (!spxm) { return {}; }

  // ---------------------------------------------------------------------------
  // Find the ServerManager proxy corresponding to the supplied VTK object.
  // ---------------------------------------------------------------------------
  vtkSMSourceProxy* ownerProxy = nullptr;

  vtkNew<vtkSMProxyIterator> iter;
  iter->SetSessionProxyManager(spxm);

  for (iter->Begin(); !iter->IsAtEnd(); iter->Next()) {
    vtkSMSourceProxy* sourceProxy = vtkSMSourceProxy::SafeDownCast(iter->GetProxy());

    if (!sourceProxy) { continue; }

    vtkObjectBase* candidate = sourceProxy->GetClientSideObject();

    if (candidate == clientObject) {
      ownerProxy = sourceProxy;
      break;
    }
  }

  if (!ownerProxy) { return {}; }

  // ---------------------------------------------------------------------------
  // Get the writer/filter's pipeline Input property.
  // ---------------------------------------------------------------------------
  vtkSMInputProperty* inputProperty = vtkSMInputProperty::SafeDownCast(ownerProxy->GetProperty("Input"));

  if (!inputProperty) { return {}; }

  // ---------------------------------------------------------------------------
  // Serialize all direct inputs recursively.
  // ---------------------------------------------------------------------------
  std::ostringstream out;

  for (unsigned int i = 0; i < inputProperty->GetNumberOfProxies(); ++i) {
    vtkSMProxy* inputProxy = inputProperty->GetProxy(i);

    if (!inputProxy) { continue; }

    // Add a separator if there are multiple writer inputs.
    if (i > 0) { out << '\n'; }

    DumpProxy(inputProxy, 0, true, out);
  }

  return out.str();
}
