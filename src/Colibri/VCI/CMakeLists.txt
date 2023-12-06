project(VCI)

find_package(Eigen3 REQUIRED)
find_package(Boost  COMPONENTS timer REQUIRED)
find_package(OpenMP REQUIRED)

set(headers
    VCI/VCI.h
    VCI/VibrationalCI.h
    VCI/Davidson.h
    )
set(sources
    VCI/VCI.cpp
    VCI/VibrationalCI.cpp
    VCI/Davidson.cpp
    )

add_library(VCI STATIC ${headers} ${sources})
target_include_directories(VibrationalUtils PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include/Vibrational>
                           )

set_target_properties(VCI PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(VCI PUBLIC $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>)

target_include_directories(VCI PUBLIC
                           ${Boost_INCLUDE_DIRS}
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:include/Vibrational>
                           )

target_link_libraries(VCI
                      PUBLIC
                      Scine::Core
                      Eigen3::Eigen
                      Boost::boost
                      Boost::timer
                      Scine::UtilsOS
                      OpenMP::OpenMP_CXX
                      PRIVATE
                      Scine::VibrationalUtils)

set_target_properties(VCI PROPERTIES OUTPUT_NAME "vci")

install(TARGETS VCI EXPORT VCI DESTINATION lib)
install(FILES ${headers} DESTINATION include/Vibrational/VCI)

add_library(Scine::VCI ALIAS VCI)

if (SCINE_BUILD_TESTS)
  add_subdirectory(Tests)
endif()
