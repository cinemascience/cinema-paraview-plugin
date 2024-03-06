set(table_of_contents
"<toc>
  <section title=\"CinemaExport Plugin\" ref=\"index.html\">
  </section>
</toc>")

paraview_add_plugin(CinemaExport
  REQUIRED_ON_CLIENT
  REQUIRED_ON_SERVER
  VERSION "1.0"
  DOCUMENTATION_DIR "${CMAKE_CURRENT_SOURCE_DIR}/documentation"
  DOCUMENTATION_ADD_PATTERNS
    "*.txt"
    "img/*.jpg"
    "js/*.js"
    "style/style.css"
  DOCUMENTATION_TOC "${table_of_contents}"
  MODULES pcAlgorithm pcImaging pcImageCompositing pcCameraGrid pcGridLayout pcWriter
  MODULE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/pcAlgorithm/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/pcImaging/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/pcImageCompositing/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/pcCameraGrid/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/pcGridLayout/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/pcWriter/vtk.module"
)

option(BUILD_TESTING "Build Plugin Testing" OFF)
if (BUILD_TESTING AND BUILD_SHARED_LIBS)
  add_subdirectory(Testing)
endif()
