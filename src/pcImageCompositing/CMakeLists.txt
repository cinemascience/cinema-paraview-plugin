set(classes pcImageCompositing)

vtk_module_add_module(pcImageCompositing
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS pcImageCompositing.xml )
