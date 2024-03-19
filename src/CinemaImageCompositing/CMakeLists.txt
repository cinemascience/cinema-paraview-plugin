set(classes CinemaImageCompositing)

vtk_module_add_module(CinemaImageCompositing
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaImageCompositing.xml )
