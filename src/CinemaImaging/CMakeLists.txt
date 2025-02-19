set(classes CinemaImaging)

vtk_module_add_module(CinemaImaging
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS CinemaImaging.xml )

target_include_directories(CinemaImaging PUBLIC ${EMBREE_INCLUDE_DIRS})
set(openmp_cxx_lib)
if (TARGET OpenMP::OpenMP_CXX)
  set(openmp_cxx_lib OpenMP::OpenMP_CXX)
endif()
target_link_libraries(CinemaImaging PUBLIC ${EMBREE_LIBRARY} ${openmp_cxx_lib})

if (USE_EMBREE3)
  target_compile_definitions(CinemaImaging PUBLIC USE_EMBREE3)
endif()
