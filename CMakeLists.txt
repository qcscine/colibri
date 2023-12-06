cmake_minimum_required(VERSION 3.9)

# Set the project name. The chosen project name must be in CamelCase. Your CMake
# tree must then provide a properly namespaced target with the same name as
# your project.
project(Colibri
  VERSION 1.0.0
  DESCRIPTION "Library for the implementation of vibrational-structure calculations."
)

# Set the module path for universal cmake files inclusion
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake)

# Component setup
include(ComponentSetup)
scine_setup_component()

# Enable testing
if(SCINE_BUILD_TESTS)
  enable_testing()
endif()

# Different parallelization options
option(MPIPARALLEL "Specifies the parallelization type. If it is on, MPI will be used. Otherwise, OpenMP is used." OFF)
if(MPIPARALLEL)
  MESSAGE(STATUS "Compiling Colibri with MPI Parallelization")
  add_definitions(-DMPI_PARALLEL)
endif()

# Also build on-the-fly stuff
option(ENABLE_OTF "Will build On-the-Fly PES calculation with SCINE calculators based on external QC programs." OFF)
option(BUILD_SPARROW "Will download and build Sparrow (the Scine semi-empirical module)." OFF)
if(ENABLE_OTF)
  add_definitions(-DOTF)
  # Also build sparrow
  if(BUILD_SPARROW)
    include(ImportSparrow)
    import_sparrow()
  endif()
endif()

option(SCINE_USE_MKL "Use the optimized MKL library for linear algebra operations of Eigen" OFF)

# Subdirectories
add_subdirectory(src)

# Documentation
scine_component_documentation(UtilsOSDocumentation)

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            )
        if(NOT ${GIT_BRANCH} STREQUAL "")
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE SVN_REV
            ERROR_QUIET
            )

        if(NOT ${SVN_REV} STREQUAL "")
            string(STRIP ${SVN_REV} SVN_REV)
        endif()

        set(PROJECT_VERSION_BUILD "${GIT_REVISION} (${GIT_BRANCH}) ${SVN_REV}")
    endif()
endif()

if(PROJECT_VERSION_BUILD)
  set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}-${PROJECT_VERSION_BUILD}")
else()
  set(PROJECT_VERSION "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}")
endif()

set(VIB_PROJECT_VERSION_STRING "Colibri module version: ${PROJECT_VERSION}")
MESSAGE(STATUS "Colibri module version: ${PROJECT_VERSION}")

######################################################################
# Configure files
######################################################################

configure_file(src/Colibri/version.h.in ${PROJECT_BINARY_DIR}/src/Colibri/version.h)

# add a directory with header files to see the generated file
include_directories (${PROJECT_BINARY_DIR}/src/Colibri)

install(FILES ${PROJECT_BINARY_DIR}/src/Colibri/version.h
        DESTINATION src/Colibri COMPONENT headers)
