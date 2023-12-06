find_package(Eigen3 REQUIRED)

# Find/Import dependecies
include(ImportGTest)
include(ImportCore)
include(ImportUtilsOS)
import_core()
import_utils_os()

add_subdirectory(VibrationalUtils)
add_subdirectory(VSCF)
add_subdirectory(VCI)
add_subdirectory(Applications)

set(headers
    ColibriModule.h
    )
set(sources
    ColibriModule.cpp
    )

add_library(Colibri SHARED ${headers} ${sources})

set_target_properties(Colibri PROPERTIES
  OUTPUT_NAME colibri
)

target_link_libraries(Colibri
                      PUBLIC
                      Scine::VSCF
                      Scine::VCI
                      PRIVATE
                      Scine::VibrationalUtils
                      )

target_compile_options(Colibri
  PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)

add_library(Scine::Colibri ALIAS Colibri)
add_library(Scine::ColibriModule ALIAS Colibri)
add_library(Scine::VibrationalModule ALIAS Colibri)

## Create symlinks for modules
# Symlink on build
add_custom_command(TARGET Colibri POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink
                   "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}colibri${CMAKE_SHARED_LIBRARY_SUFFIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/colibri.module${CMAKE_SHARED_LIBRARY_SUFFIX}"
)
# Symlink on install
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}colibri${CMAKE_SHARED_LIBRARY_SUFFIX} \
        ${CMAKE_INSTALL_PREFIX}/lib/colibri.module${CMAKE_SHARED_LIBRARY_SUFFIX})"
)

if (SCINE_BUILD_TESTS)
  add_subdirectory(Tests)
endif()