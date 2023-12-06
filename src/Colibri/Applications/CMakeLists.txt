# Applications
project(Applications)

import_utils_os()
import_core()

find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

if (MPIPARALLEL)
  find_package(MPI REQUIRED)
endif()

add_executable(vscf vscf.cpp)
add_executable(vscfMC vscfMC.cpp)
add_executable(vscfMCSP vscfMCSP.cpp)
add_executable(vci vci.cpp)
add_executable(vciFromFcidump vciFromFcidump.cpp)
add_executable(vib_from_SPs vib_from_SPs.cpp)

if(ENABLE_OTF)
  add_executable(vib_otf vib_otf.cpp)
  add_executable(vib_hybrid vib_hybrid.cpp)
endif()


target_link_libraries(vscf
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )

target_link_libraries(vscfMC
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )

target_link_libraries(vscfMCSP
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )

target_link_libraries(vci
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VCI
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )

target_link_libraries(vciFromFcidump
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VCI
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )

if (MPIPARALLEL)
target_link_libraries(vib_from_SPs
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VCI
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      MPI::MPI_CXX
                      )
else()
target_link_libraries(vib_from_SPs
                      PRIVATE
                      Scine::VibrationalUtils
                      Scine::VCI
                      Scine::VSCF
                      Boost::boost
                      Eigen3::Eigen
                      )
endif()

if(ENABLE_OTF)
  if (MPIPARALLEL)
    target_link_libraries(vib_otf
                          PRIVATE
                          Scine::VibrationalUtils
                          Scine::VCI
                          Scine::VSCF
                          Boost::boost
                          Eigen3::Eigen
                          MPI::MPI_CXX
                          )
    target_link_libraries(vib_hybrid
                          PRIVATE
                          Scine::VibrationalUtils
                          Scine::VCI
                          Scine::VSCF
                          Boost::boost
                          Eigen3::Eigen
                          MPI::MPI_CXX
                          )
  else()
    target_link_libraries(vib_otf
                          PRIVATE
                          Scine::VibrationalUtils
                          Scine::VCI
                          Scine::VSCF
                          Boost::boost
                          Eigen3::Eigen
                          )
    target_link_libraries(vib_hybrid
                          PRIVATE
                          Scine::VibrationalUtils
                          Scine::VCI
                          Scine::VSCF
                          Boost::boost
                          Eigen3::Eigen
                          )
  endif()
endif()

