set(classes CinemaDepthImageProjection)

vtk_module_add_module(CinemaDepthImageProjection
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaDepthImageProjection.xml )

set(openmp_cxx_lib)
if (TARGET OpenMP::OpenMP_CXX)
  set(openmp_cxx_lib OpenMP::OpenMP_CXX)
endif()
target_link_libraries(CinemaDepthImageProjection PUBLIC ${openmp_cxx_lib})
