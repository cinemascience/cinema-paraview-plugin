set(classes CinemaColorImaging)

vtk_module_add_module(CinemaColorImaging
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaColorImaging.xml )
