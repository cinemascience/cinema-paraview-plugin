set(classes CinemaWriter)

vtk_module_add_module(CinemaWriter
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaWriter.xml )
