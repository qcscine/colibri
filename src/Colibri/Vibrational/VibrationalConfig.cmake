# Dependencies
include(CMakeFindDependencyMacro)
if(NOT TARGET Boost::log)
  set(Boost_USE_STATIC_LIBS OFF)
  unset(Boost_FOUND)
  find_dependency(Boost REQUIRED COMPONENTS log)
endif()
find_dependency(Eigen3 3.3.2 REQUIRED NO_MODULE)
find_dependency(yaml-cpp REQUIRED)
# We are interface-dependent on Core and UtilsOS, so those headers must also be available
find_dependency(Scine REQUIRED COMPONENTS Core UtilsOS)

include(${CMAKE_CURRENT_LIST_DIR}/UtilsOSTargets.cmake)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

check_required_components(UtilsOS)
