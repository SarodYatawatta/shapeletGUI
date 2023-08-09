include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIRS}/wcslib)
include_directories(${WCSLIB_INCLUDE_DIRS})
include_directories(${Qt6Widgets_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

qt_add_executable(shapeletGUI
    optionsdialog.ui
    shapeletmain.ui
    textdialog.ui
    main.cpp
    mygraphicsview.cpp
    optionsdialog.cpp
    shapeletmain.cpp
    textdialog.cpp
    )

target_link_libraries(shapeletGUI PRIVATE 
  Qt6::Widgets 
  -lshapelet
  ${CFITSIO_LIB}
  ${OpenBLAS_LIB}
  ${LAPACK_LIBRARIES}
  ${WCSLIB_LIBRARIES}
  ${GLIB_PKG_LIBRARIES}
  ${FFTW_DOUBLE_LIB}
  ${FFTW_DOUBLE_THREADS_LIB}
  -lpthread
  -lm
  )
qt_add_resources(TARGET shapeletGUI SOURCES resources.qrc)
set_property(TARGET shapeletGUI PROPERTY AUTOMOC ON)
set_property(TARGET shapeletGUI PROPERTY AUTORCC ON)
set_property(TARGET shapeletGUI PROPERTY AUTOUIC ON)
install(TARGETS shapeletGUI DESTINATION bin)

