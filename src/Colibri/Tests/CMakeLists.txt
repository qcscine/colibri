# Tests
project(ModuleTests)

import_gtest()
import_utils_os()
import_core()
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

if (MPIPARALLEL)
  find_package(MPI REQUIRED)
endif()

add_executable(ModuleTest ModuleTest.cpp)

if (MPIPARALLEL)
    target_link_libraries(ModuleTest
                        PRIVATE
                        gtest_main
                        gmock
                        Scine::Core
                        Scine::UtilsOS
                        MPI::MPI_CXX
                        )
else()
    target_link_libraries(ModuleTest
                        PRIVATE
                        gtest_main
                        gmock
                        Scine::Core
                        Scine::UtilsOS
                        )
endif()

add_test(NAME ModuleTest COMMAND ModuleTest)