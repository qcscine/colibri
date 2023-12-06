project(VibrationalUtils)
import_core()
import_utils_os()
find_package(Eigen3 REQUIRED)
find_package(Boost  COMPONENTS program_options timer  REQUIRED)
find_package(OpenMP REQUIRED)
if (MPIPARALLEL)
  find_package(MPI REQUIRED)
endif()

set(headers
    VibrationalUtils/VibrationalMode.h
    VibrationalUtils/GaussHermite.h
    VibrationalUtils/VibrationalParameters.h
    VibrationalUtils/IntegralDumper.h
    VibrationalUtils/PESLibrary/PES.h
    VibrationalUtils/PESLibrary/PesQFF.h
    VibrationalUtils/PESLibrary/PesTaylorParam.h
    VibrationalUtils/PESLibrary/PesMCSP.h
    VibrationalUtils/PESLibrary/PesChristoffel.h
    VibrationalUtils/PESLibrary/PesCartesian.h
    VibrationalUtils/PESLibrary/NMAPesCalculator.h
    VibrationalUtils/PESLibrary/H2OHarmonicPesCalculator.h
    VibrationalUtils/PESLibrary/PesCoupledOscillator.h
    VibrationalUtils/PESLibrary/PesFromSPs.h
    VibrationalUtils/VibrationalBases/ModeGrid.h
    VibrationalUtils/VibrationalBases/ModalBasis.h
    VibrationalUtils/VibrationalBases/DistributedGaussians.h
    VibrationalUtils/VibrationalBases/DVR.h
    )

set(sources
    VibrationalUtils/GaussHermite.cpp
    VibrationalUtils/VibrationalParameters.cpp
    VibrationalUtils/IntegralDumper.cpp
    VibrationalUtils/PESLibrary/PesChristoffel.cpp
    VibrationalUtils/PESLibrary/PesCartesian.cpp
    VibrationalUtils/PESLibrary/PesTaylorParam.cpp
    VibrationalUtils/PESLibrary/PesMCSP.cpp
    VibrationalUtils/PESLibrary/PesQFF.cpp
    VibrationalUtils/PESLibrary/NMAPesCalculator.cpp
    VibrationalUtils/PESLibrary/H2OHarmonicPesCalculator.cpp
    VibrationalUtils/PESLibrary/PesCoupledOscillator.cpp
    VibrationalUtils/PESLibrary/PesFromSPs.cpp
    VibrationalUtils/VibrationalBases/ModeGrid.cpp
    VibrationalUtils/VibrationalBases/ModalBasis.cpp
    VibrationalUtils/VibrationalBases/DistributedGaussians.cpp
    VibrationalUtils/VibrationalBases/DVR.cpp
    )

if(ENABLE_OTF)
  set(headers ${headers} VibrationalUtils/PESLibrary/OnTheFlyPes.h
                         VibrationalUtils/PESLibrary/HybridOTFfromSPs.h)
  set(sources ${sources} VibrationalUtils/PESLibrary/OnTheFlyPes.cpp
                         VibrationalUtils/PESLibrary/HybridOTFfromSPs.cpp)
endif()

add_library(VibrationalUtils STATIC ${headers} ${sources})

target_include_directories(VibrationalUtils PUBLIC ${Boost_INCLUDE_DIRS})
set_target_properties(VibrationalUtils PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(VibrationalUtils PUBLIC
                       $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
                       $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>)

set(dirsToInclude
    $<INSTALL_INTERFACE:include/Vibrational>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

target_include_directories(VibrationalUtils PUBLIC ${dirsToInclude})

if (MPIPARALLEL)
  set(libsToLink
      Eigen3::Eigen
      Boost::boost
      Boost::program_options
      Boost::timer
      MPI::MPI_CXX
      Scine::Core
      Scine::UtilsOS
      $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
      )
else()
  set(libsToLink
      Eigen3::Eigen
      Boost::boost
      Boost::program_options
      Boost::timer
      Scine::Core
      Scine::UtilsOS
      $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
      )
endif()

target_link_libraries(VibrationalUtils PUBLIC ${libsToLink} -lgfortran)

set_target_properties(VibrationalUtils PROPERTIES OUTPUT_NAME "vibrationalutils")

install(TARGETS VibrationalUtils EXPORT VibrationalUtils DESTINATION lib)
install(FILES ${headers} DESTINATION include/Vibrational/VibrationalUtils)
add_library(Scine::VibrationalUtils ALIAS VibrationalUtils)

if (SCINE_BUILD_TESTS)
  add_subdirectory(Tests)
endif()
