project(VSCF)

find_package(Eigen3 REQUIRED)
find_package(Boost COMPONENTS timer REQUIRED)
if (MPIPARALLEL)
  find_package(MPI REQUIRED)
  find_package(OpenMP REQUIRED)
else()
  find_package(OpenMP REQUIRED)
endif()

set(headers
    VSCF/VSCF.h
    VSCF/VscfCalculator.h
    VSCF/VibrationalSCF.h
    )
set(sources
    VSCF/VSCF.cpp
    VSCF/VscfCalculator.cpp
    VSCF/VibrationalSCF.cpp
    )

add_library(VSCF STATIC ${headers} ${sources})
target_include_directories(VibrationalUtils PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include/Vibrational>
                           )

set_target_properties(VSCF PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(VSCF PUBLIC $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

target_include_directories(VSCF PUBLIC
                           ${Boost_INCLUDE_DIRS}
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include/Vibrational>
                           ${MPI_INCLUDE_PATH}
                           )

if (MPIPARALLEL)
  target_link_libraries(VSCF
                        PUBLIC
                        Scine::Core
                        Eigen3::Eigen
                        Boost::boost
                        Boost::timer
                        Scine::UtilsOS
                        OpenMP::OpenMP_CXX
                        MPI::MPI_CXX
                        PRIVATE
                        Scine::VibrationalUtils)
else()
  target_link_libraries(VSCF
                        PUBLIC
                        Scine::Core
                        Eigen3::Eigen
                        Boost::boost
                        Boost::timer
                        Scine::UtilsOS
                        OpenMP::OpenMP_CXX
                        PRIVATE
                        Scine::VibrationalUtils)
endif()

set_target_properties(VSCF PROPERTIES OUTPUT_NAME "vscf")

install(TARGETS VSCF EXPORT VSCF DESTINATION lib)
install(FILES ${headers} DESTINATION include/Vibrational/VSCF)

add_library(Scine::VSCF ALIAS VSCF)

if (SCINE_BUILD_TESTS)
  add_subdirectory(Tests)
endif()
