#include "vtkSMRenderViewNameDomain.h"

#include <pqActiveObjects.h>
#include <pqView.h>
#include <string>
#include <vector>
#include <vtkCollection.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkSMProperty.h>
#include <vtkSMProxy.h>
#include <vtkSMRenderViewProxy.h>
#include <vtkSMSessionProxyManager.h>
#include <vtkSMStringVectorProperty.h>

vtkStandardNewMacro(vtkSMRenderViewNameDomain);

void vtkSMRenderViewNameDomain::PrintSelf(ostream& os, vtkIndent indent) { this->Superclass::PrintSelf(os, indent); }

void vtkSMRenderViewNameDomain::Update(vtkSMProperty* property) {
  std::vector<std::string> names;

  vtkSMProxy* parent = property ? property->GetParent() : nullptr;
  vtkSMSessionProxyManager* pxm = parent ? parent->GetSessionProxyManager() : nullptr;
  if (pxm) {
    vtkNew<vtkCollection> views;
    pxm->GetProxies("views", views);
    for (vtkIdType i = 0; i < views->GetNumberOfItems(); ++i) {
      vtkSMProxy* proxy = vtkSMProxy::SafeDownCast(views->GetItemAsObject(i));
      if (!vtkSMRenderViewProxy::SafeDownCast(proxy)) { continue; }

      const char* name = pxm->GetProxyName("views", proxy);
      if (name && *name) { names.emplace_back(name); }
    }
  }

  this->SetStrings(names);
}

int vtkSMRenderViewNameDomain::SetDefaultValues(vtkSMProperty* property, bool useUncheckedValues) {
  this->Update(property);

  auto* svp = vtkSMStringVectorProperty::SafeDownCast(property);
  vtkSMProxy* parent = property ? property->GetParent() : nullptr;
  vtkSMSessionProxyManager* pxm = parent ? parent->GetSessionProxyManager() : nullptr;
  if (!svp || !pxm || this->GetNumberOfStrings() == 0) { return 0; }

  std::string defaultName;
  if (pqView* pqview = pqActiveObjects::instance().activeView()) {
    vtkSMProxy* activeProxy = pqview->getViewProxy();
    if (vtkSMRenderViewProxy::SafeDownCast(activeProxy) && activeProxy->GetSession() == parent->GetSession()) {
      if (const char* name = pxm->GetProxyName("views", activeProxy)) { defaultName = name; }
    }
  }

  if (defaultName.empty()) { defaultName = this->GetString(0); }

  if (useUncheckedValues) {
    svp->SetUncheckedElement(0, defaultName.c_str());
  } else {
    svp->SetElement(0, defaultName.c_str());
  }
  return 1;
}
