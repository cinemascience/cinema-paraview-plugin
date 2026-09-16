#include "PipelineProvenance.h"

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
#include <vtkSMViewProxy.h>

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

bool IsVisibleRepresentation(vtkSMProxy* representation) {
  if (!representation) { return false; }
  auto* visibility = vtkSMIntVectorProperty::SafeDownCast(representation->GetProperty("Visibility"));
  return !visibility || visibility->GetNumberOfElements() == 0 || visibility->GetElement(0) != 0;
}

void CollectProxy(vtkSMProxy* proxy, const std::string& prefix, bool followInputs, std::map<std::string, int>& counters,
                  std::set<vtkSMProxy*>& visited, ProvenanceMap& result) {
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
        CollectProxy(p->GetProxy(i), column + std::to_string(i), false, counters, visited, result);
      }
    }
  }

  iter->Delete();
}

} // namespace

std::vector<vtkSMProxy*> GetPipelineInputs(vtkSMProxy* proxy, const char* propertyName) {
  std::vector<vtkSMProxy*> result;
  if (!proxy || !propertyName) { return result; }

  auto* input = vtkSMInputProperty::SafeDownCast(proxy->GetProperty(propertyName));
  if (!input) { return result; }

  result.reserve(input->GetNumberOfProxies());
  for (unsigned int i = 0; i < input->GetNumberOfProxies(); ++i) {
    if (auto* upstream = input->GetProxy(i)) { result.push_back(upstream); }
  }
  return result;
}

std::vector<vtkSMProxy*> GetVisibleRepresentationInputs(vtkSMViewProxy* view) {
  std::vector<vtkSMProxy*> result;
  if (!view) { return result; }

  auto* representations = vtkSMProxyProperty::SafeDownCast(view->GetProperty("Representations"));
  if (!representations) { return result; }

  std::set<vtkSMProxy*> unique;
  for (unsigned int i = 0; i < representations->GetNumberOfProxies(); ++i) {
    vtkSMProxy* representation = representations->GetProxy(i);
    if (!IsVisibleRepresentation(representation)) { continue; }

    auto* input = vtkSMInputProperty::SafeDownCast(representation->GetProperty("Input"));
    if (!input) { continue; }
    for (unsigned int j = 0; j < input->GetNumberOfProxies(); ++j) {
      vtkSMProxy* source = input->GetProxy(j);
      if (source && unique.insert(source).second) { result.push_back(source); }
    }
  }
  return result;
}

ProvenanceMap ComputePipelineProvenance(vtkSMProxy* root) {
  return ComputePipelineProvenance(std::vector<vtkSMProxy*>{root});
}

ProvenanceMap ComputePipelineProvenance(const std::vector<vtkSMProxy*>& roots) {
  ProvenanceMap result;
  std::map<std::string, int> counters;
  std::set<vtkSMProxy*> visited;
  for (vtkSMProxy* root : roots) { CollectProxy(root, "", true, counters, visited, result); }
  return result;
}

ProvenanceMap ComputeVisibleRepresentationProvenance(vtkSMViewProxy* view,
                                                     const std::vector<vtkSMProxy*>& additionalRoots) {
  auto roots = GetVisibleRepresentationInputs(view);
  roots.insert(roots.end(), additionalRoots.begin(), additionalRoots.end());
  return ComputePipelineProvenance(roots);
}

void MergeProvenance(ProvenanceMap& destination, const ProvenanceMap& source, bool overwriteExisting) {
  for (const auto& item : source) {
    if (overwriteExisting) {
      destination[item.first] = item.second;
    } else {
      destination.emplace(item.first, item.second);
    }
  }
}
