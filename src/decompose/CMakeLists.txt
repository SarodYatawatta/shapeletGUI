include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIRS}/wcslib)
include_directories(${WCSLIB_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(HAVE_CUDA)
    CUDA_INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../lib)
    CUDA_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
endif()

link_directories(${LIBRARY_OUTPUT_PATH})
link_directories(${CMAKE_CURRENT_SOURCE_DIR})
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../lib)

FILE(GLOB SRCFILES *.c)
target_link_libraries(${target} shapelet)

if (HAVE_CUDA)
  cuda_add_executable(decompose ${SRCFILES})
  add_dependencies(decompose shapelet)
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
   ${CUDA_LIBRARIES}
   ${CUDA_CUBLAS_LIBRARIES}
   ${CUDA_CUFFT_LIBRARIES}
   ${CUDA_cusolver_LIBRARY}
   ${CUDA_cudadevrt_LIBRARY}
   ${NVML_LIB_PATH}
  )

  install(TARGETS decompose DESTINATION bin)

else()
  add_executable(decompose ${SRCFILES})
  add_dependencies(decompose shapelet)

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
endif()
