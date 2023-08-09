cmake_minimum_required(VERSION 3.10)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.20")
 cmake_policy(SET CMP0115 OLD)
endif()
#--------------------------------------- project
project (ShapeletGUI)
set(PROJECT_VERSION_MAJOR 0)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 1)
set(PROJECT_VERSION_REVISION 0)
set(PROJECT_VERSION
  "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")

message(STATUS "\n############################\n# ShapeletGUI v${PROJECT_VERSION}\n############################")


#--------------------------------------- languages
enable_language(C)
enable_language(Fortran)


#--------------------------------------- compiler
message(STATUS "\n############################\n# Building ShapeletGUI for ${CMAKE_HOST_SYSTEM_NAME} \n############################")

if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
  endif()
endif()

include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)
include(CheckFunctionExists)


#--------------------------------------- paths
SET(MAINFOLDER ${CMAKE_BINARY_DIR})

SET(EXECUTABLE_OUTPUT_PATH "${MAINFOLDER}/dist/bin")
SET(LIBRARY_OUTPUT_PATH "${MAINFOLDER}/dist/lib")

# Load CMAKE configuration from environment variables
set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH} )
set(CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} )

#--------------------------------------- cmake modules
list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/CMakeModules)


#--------------------------------------- user options
option (ENABLE_MPI "Enable MPI support" OFF)

#--------------------------------------- rpath setup for libraries

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")


#--------------------------------------- dependencies

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

message(STATUS "\n############################\n#   Searching ShapeletGUI dependencies\n############################")
# BLAS: use -DBLA_VENDOR=OpenBLAS for OpenBLAS
find_package(BLAS REQUIRED)
include_directories(${BLAS_INCLUDE_DIR})

if ("${BLA_VENDOR}" STREQUAL "OpenBLAS")
  message(STATUS "Using ${BLA_VENDOR}")
  # add preprocessor def
  add_definitions(-DHAVE_OPENBLAS)
else()
  message(STATUS "Using generic BLAS")
endif()


# LAPACK (will also find the generic BLAS)
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIR})

#wcslib
find_package(WcsLib REQUIRED)
include_directories(${WCSLIB_INCLUDE_DIR})
#cfitsio
find_package(CfitsIO REQUIRED)
include_directories(${CFITSIO_INCLUDE})

#glib
find_package(PkgConfig REQUIRED)
pkg_check_modules(GLIB_PKG glib-2.0)
if (GLIB_PKG_FOUND)
    include_directories(${GLIB_PKG_INCLUDE_DIRS})
endif()

find_package(FFTW REQUIRED)

find_package(Qt6 REQUIRED COMPONENTS Widgets Core)
#qt_standard_project_setup()

#--------------------------------------- check environment variables


#--------------------------------------- build parameters for all targets
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -pg -Wall -DNDEBUG -Wextra -Wcast-align -Wpointer-arith -Woverloaded-virtual -Wwrite-strings -Wcast-qual -Wctor-dtor-privacy -Wno-long-long-Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused")

set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#--------------------------------------- summary
message(STATUS "\n############################\n#   Configuration summary\n############################")
message (STATUS "CMAKE_SYSTEM .............. = ${CMAKE_SYSTEM}")
message (STATUS "CMAKE_INSTALL_PREFIX ...... = ${CMAKE_INSTALL_PREFIX}")
message (STATUS "CMAKE_BUILD_TYPE .......... = ${CMAKE_BUILD_TYPE}")
message (STATUS "CMAKE_CXX_COMPILER ........ = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS_RELEASE ... = ${CMAKE_CXX_FLAGS_RELEASE}")

message (STATUS "ENABLE_MPI ................ = ${ENABLE_MPI}")

message (STATUS "OpenBLAS_LIB .............. = ${OpenBLAS_LIB}")

message (STATUS "GLIB_PKG_INCLUDE_DIRS...... = ${GLIB_PKG_INCLUDE_DIRS}")
message (STATUS "GLIB_PKG_LIBRARIES......... = ${GLIB_PKG_LIBRARIES}")

message (STATUS "LAPACK_LIBRARIES........... = ${LAPACK_LIBRARIES}")

message (STATUS "CMAKE_Fortran_COMPILER..... = ${CMAKE_Fortran_COMPILER}")

message (STATUS "CFITSIO_ROOT_DIR........... = ${CFITSIO_ROOT_DIR}")
message (STATUS "CFITSIO_INCLUDE............ = ${CFITSIO_INCLUDE}")
message (STATUS "CFITSIO_LIB................ = ${CFITSIO_LIB}")

message (STATUS "WCSLIB_INCLUDE_DIRS........ = ${WCSLIB_INCLUDE_DIRS}")
message (STATUS "WCSLIB_LIBRARIES .......... = ${WCSLIB_LIBRARIES}")

#--------------------------------------- include directories
add_subdirectory(src)

#--------------------------------------- install

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist/ DESTINATION ${CMAKE_INSTALL_PREFIX})

install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/shapeletGUI/
        DESTINATION bin
        FILES_MATCHING PATTERN "*.py" )
