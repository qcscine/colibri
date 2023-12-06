project(VibrationalTests)
import_gtest()
find_package(Boost REQUIRED)


# Tests
add_executable(VSCFTest VSCFTest.cpp)
target_link_libraries(VSCFTest
                      PRIVATE
                      gtest_main
                      gmock
                      Scine::VibrationalUtils
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )
add_test(NAME VSCFTest COMMAND VSCFTest)

if(ENABLE_OTF)
  add_executable(OTFVSCFTest OTFVSCFTest.cpp)
  target_link_libraries(OTFVSCFTest
                        PRIVATE
                        gtest_main
                        gmock
                        Scine::VibrationalUtils
                        Scine::VSCF
                        Boost::boost
                        Eigen3::Eigen
                        )
  add_test(NAME OTFVSCFTest COMMAND OTFVSCFTest)
endif()

# Setup App tests if possible
if((SCINE_BUILD_TESTS) AND (TARGET Scine::Sparrow))
  include(TargetLibName)
  target_lib_directory(Scine::Sparrow _sparrow_module_path)
  message(STATUS "Sparrow module found at ${_sparrow_module_path}")
  target_lib_directory(Scine::UtilsOS _utils_module_path)
  message(STATUS "Utils module found at ${_utils_module_path}")
  set_tests_properties(OTFVSCFTest PROPERTIES
    ENVIRONMENT "PATH=$ENV{PATH}:${CMAKE_CURRENT_BINARY_DIR};SCINE_MODULE_PATH=${_sparrow_module_path}:${_utils_module_path}"
  )
  unset(_sparrow_module_path)
  unset(_utils_module_path)
endif()
