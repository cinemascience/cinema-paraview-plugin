set(classes pcGridLayout)

vtk_module_add_module(pcGridLayout
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS pcGridLayout.xml )
