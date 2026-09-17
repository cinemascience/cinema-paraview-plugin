set(classes
  CinemaWriter
  vtkSMCinemaWriterProxy
)

vtk_module_add_module(CinemaWriter
  FORCE_STATIC
  CLASSES ${classes}
  SOURCES
    H5Export.cxx
    PipelineProvenance.cxx
  NOWRAP_HEADERS
    PipelineProvenance.h
)

paraview_add_server_manager_xmls(XMLS CinemaWriter.xml)
