#include "PipelineAnnotation.h"

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vtkSMDoubleVectorProperty.h>
#include <vtkSMInputProperty.h>
#include <vtkSMIntVectorProperty.h>
#include <vtkSMProperty.h>
#include <vtkSMPropertyIterator.h>
#include <vtkSMProxy.h>
#include <vtkSMProxyProperty.h>
#include <vtkSMStringVectorProperty.h>

namespace {

std::string JsonQuote(const std::string& value) {
  std::ostringstream out;
  out << '"';
  for (char c : value) {
    if (c == '"' || c == '\\') { out << '\\'; }
    out << c;
  }
  out << '"';
  return out.str();
}

template <typename T> std::string VectorValue(T* property) {
  const unsigned int n = property->GetNumberOfElements();
  if (n == 1) {
    std::ostringstream out;
    out << property->GetElement(0);
    return out.str();
  }

  std::ostringstream out;
  out << '[';
  for (unsigned int i = 0; i < n; ++i) {
    if (i) { out << ','; }
    out << property->GetElement(i);
  }
  out << ']';
  return out.str();
}

// Excludes properties that do not represent persistent user-facing proxy state.
bool ShouldRecordProperty(vtkSMProperty* property) {
  return property && !property->GetInformationOnly() && !property->GetIsInternal();
}

std::string StringVectorValue(vtkSMStringVectorProperty* property) {
  const unsigned int n = property->GetNumberOfElements();
  if (n == 1) {
    const char* value = property->GetElement(0);
    return value ? value : "";
  }

  std::ostringstream out;
  out << '[';
  for (unsigned int i = 0; i < n; ++i) {
    if (i) { out << ','; }
    const char* value = property->GetElement(i);
    out << JsonQuote(value ? value : "");
  }
  out << ']';
  return out.str();
}

// Recursively flattens one proxy. Pipeline inputs are followed as peers;
// nested proxy properties are namespaced below their owning property.
void CollectProxy(vtkSMProxy* proxy, const std::string& prefix, bool followInputs, std::map<std::string, int>& counters,
                  std::set<vtkSMProxy*>& visited, AnnotationMap& result) {
  if (!proxy || !visited.insert(proxy).second) { return; }

  const std::string type = proxy->GetXMLName() ? proxy->GetXMLName() : proxy->GetClassName();
  const std::string name = prefix.empty() ? type + std::to_string(counters[type]++) : prefix;

  vtkSMPropertyIterator* iter = proxy->NewPropertyIterator();
  if (!iter) { return; }

  for (iter->Begin(); !iter->IsAtEnd(); iter->Next()) {
    vtkSMProperty* property = iter->GetProperty();
    const char* key = iter->GetKey();
    if (!property || !key) { continue; }

    if (auto* input = vtkSMInputProperty::SafeDownCast(property)) {
      if (followInputs) {
        for (unsigned int i = 0; i < input->GetNumberOfProxies(); ++i) {
          CollectProxy(input->GetProxy(i), "", true, counters, visited, result);
        }
      }
      continue;
    }

    if (!ShouldRecordProperty(property)) { continue; }

    const std::string column = name + "." + key;

    if (auto* p = vtkSMIntVectorProperty::SafeDownCast(property)) {
      result[column] = VectorValue(p);
    } else if (auto* p = vtkSMDoubleVectorProperty::SafeDownCast(property)) {
      result[column] = VectorValue(p);
    } else if (auto* p = vtkSMStringVectorProperty::SafeDownCast(property)) {
      result[column] = StringVectorValue(p);
    } else if (auto* p = vtkSMProxyProperty::SafeDownCast(property)) {
      for (unsigned int i = 0; i < p->GetNumberOfProxies(); ++i) {
        const std::string child = column + std::to_string(i);
        CollectProxy(p->GetProxy(i), child, false, counters, visited, result);
      }
    }
  }

  iter->Delete();
}

} // namespace

AnnotationMap ComputePipelineAnnotations(vtkSMProxy* writerProxy) {
  AnnotationMap result;
  if (!writerProxy) { return result; }

  auto* input = vtkSMInputProperty::SafeDownCast(writerProxy->GetProperty("Input"));
  if (!input) { return result; }

  std::map<std::string, int> counters;
  std::set<vtkSMProxy*> visited;
  for (unsigned int i = 0; i < input->GetNumberOfProxies(); ++i) {
    CollectProxy(input->GetProxy(i), "", true, counters, visited, result);
  }
  return result;
}
