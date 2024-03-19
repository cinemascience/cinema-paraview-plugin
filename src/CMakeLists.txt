set(table_of_contents
"<toc>
  <section title=\"CinemaExporter Plugin\" ref=\"index.html\">
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
  MODULES CinemaAlgorithm CinemaImaging CinemaImageCompositing CinemaCameraGrid CinemaGridLayout CinemaWriter
  MODULE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaAlgorithm/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaImaging/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaImageCompositing/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaCameraGrid/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaGridLayout/vtk.module"
    "${CMAKE_CURRENT_SOURCE_DIR}/CinemaWriter/vtk.module"
)

option(BUILD_TESTING "Build Plugin Testing" OFF)
if (BUILD_TESTING AND BUILD_SHARED_LIBS)
  add_subdirectory(Testing)
endif()
