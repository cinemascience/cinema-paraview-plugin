set(classes CinemaWriter)

vtk_module_add_module(CinemaWriter
  FORCE_STATIC
  CLASSES ${classes}
  SOURCES
    H5Export.cxx
    PipelineAnnotation.cxx
)

paraview_add_server_manager_xmls( XMLS CinemaWriter.xml )
