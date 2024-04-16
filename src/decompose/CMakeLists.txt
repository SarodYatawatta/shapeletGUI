include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIRS}/wcslib)
include_directories(${WCSLIB_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

link_directories(${LIBRARY_OUTPUT_PATH})
link_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

FILE(GLOB SRCFILES *.c)
add_executable(decompose ${SRCFILES})

target_link_libraries(decompose
 -lshapelet
 ${CFITSIO_LIB}
 ${OpenBLAS_LIB}
 ${LAPACK_LIBRARIES}
 ${WCSLIB_LIBRARIES}
 ${FFTW_DOUBLE_LIB}
 ${FFTW_DOUBLE_THREADS_LIB}
 -lpthread
 -lm
)

install(TARGETS decompose DESTINATION bin)
