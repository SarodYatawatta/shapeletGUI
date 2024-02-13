message("-------------------------------------")
message (STATUS "Configuring lib")

# Common objects
set (objects
  decom_fits.c
  fits_utils.c
  lapack.c
  myblas.c
  shapelet_lm.c
  fft.c
  hermite.c
  multi_fits.c
  shapelet_uv.c
  apc.c
)
message (STATUS "lib objects --> ${objects}")

find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIRS}/wcslib)
include_directories(${WCSLIB_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(shapelet SHARED ${objects} shapelet.h)
SET_TARGET_PROPERTIES(shapelet PROPERTIES SOVERSION 0.0)
target_link_libraries(shapelet PUBLIC m ${CMAKE_THREAD_LIBS_INIT}
       ${BLAS_LIBRARIES}
       ${LAPACK_LIBRARIES})
