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

find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIRS}/wcslib)
include_directories(${WCSLIB_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(HAVE_CUDA)
  message (STATUS "compiling with CUDA support")
  message (STATUS "lib objects --> ${objects}")
  set (extra_objects
    )
  set (extra_objects_cuda
    shapelet_cuda
    )

  set(objects ${objects} ${extra_objects} ${extra_objects_cuda})
  foreach (object ${objects})
        file(GLOB SRC_FILE ${object}.*)
        set(CUDA_SRC_FILES ${CUDA_SRC_FILES} ${SRC_FILE})
  endforeach (object)

  CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
  CUDA_ADD_LIBRARY(shapelet ${objects} ${CUDA_SRC_FILES} shapelet.h)

  message (STATUS "CUDA SRC --> ${CUDA_SRC_FILES}")
  SET_TARGET_PROPERTIES(shapelet PROPERTIES SOVERSION 0.0)
  SET_TARGET_PROPERTIES(shapelet PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.18")
      SET_TARGET_PROPERTIES(shapelet PROPERTIES CUDA_ARCHITECTURES native)
  endif()
  target_link_libraries(shapelet ${CMAKE_THREAD_LIBS_INIT}
       ${BLAS_LIBRARIES}
       ${LAPACK_LIBRARIES})

else()
  message (STATUS "lib objects --> ${objects}")

  add_library(shapelet SHARED ${objects} shapelet.h)
  SET_TARGET_PROPERTIES(shapelet PROPERTIES SOVERSION 0.0)
  target_link_libraries(shapelet PUBLIC m ${CMAKE_THREAD_LIBS_INIT}
       ${BLAS_LIBRARIES}
       ${LAPACK_LIBRARIES})
endif()
