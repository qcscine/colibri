project(VCITests)
import_gtest()

# Tests
add_executable(VCITest VCITest.cpp)
target_link_libraries(VCITest
                      PRIVATE
                      Scine::VCI
                      Scine::VibrationalUtils
                      Scine::VSCF
                      gtest_main
                      gmock
                      Eigen3::Eigen
                      )
                      
add_test(NAME VCITest COMMAND VCITest)
