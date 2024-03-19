set(classes CinemaCameraGrid)

vtk_module_add_module(CinemaCameraGrid
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaCameraGrid.xml )
