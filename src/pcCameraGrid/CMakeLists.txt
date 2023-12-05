set(classes pcCameraGrid)

vtk_module_add_module(pcCameraGrid
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS pcCameraGrid.xml )
