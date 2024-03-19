set(classes CinemaGridLayout)

vtk_module_add_module(CinemaGridLayout
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaGridLayout.xml )
