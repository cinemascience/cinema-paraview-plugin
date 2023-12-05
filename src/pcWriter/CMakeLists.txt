set(classes pcWriter)

vtk_module_add_module(pcWriter
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS pcWriter.xml )
