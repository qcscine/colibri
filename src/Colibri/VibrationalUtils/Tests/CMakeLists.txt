project(VibrationalUtilsTests)
import_gtest()

# Tests
add_executable(VibrationalUtilsTest VibrationalUtilsTest.cpp)
target_link_libraries(VibrationalUtilsTest
                      PRIVATE
                      gtest_main
                      gmock
                      Scine::VibrationalUtils
                      Eigen3::Eigen
                      )
                      
add_test(NAME VibrationalUtilsTest COMMAND VibrationalUtilsTest)
