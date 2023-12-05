set(classes pcImaging)

vtk_module_add_module(pcImaging
  FORCE_STATIC
  CLASSES ${classes}
)

paraview_add_server_manager_xmls( XMLS pcImaging.xml )

target_include_directories(pcImaging PUBLIC ${EMBREE_INCLUDE_DIR})
target_link_libraries(pcImaging PUBLIC ${EMBREE_LIBRARY} OpenMP::OpenMP_CXX)
