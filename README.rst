SCINE - Colibri
===============

Introduction
------------

Colibri is a code for vibrational structure theory calculations.
It provides the following methods/tools:

- Normal mode analysis
- VSCF
- VCI
- On-the-fly PES construction
- FCIDUMP file generator for vDMRG calculations

Colibri can calculate harmonic, mean-field and correlated vibrational energies
for both the ground state (ZPVE) and vibrationally excited states. 

License and Copyright Information
---------------------------------

Colibri is distributed under the BSD 3-clause "New" or "Revised" License.
For more license and copyright information, see the file ``LICENSE.txt`` in the repository.

Installation and Usage
----------------------

Prerequisites
.............

For the installation, the following requirements are needed:

- a C++ compiler. Only `GCC <https://gcc.gnu.org>`_ has been tested.
- `CMake <https://cmake.org>`_
- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_


Installation
............

Colibri can be cloned and installed with:

.. code-block:: bash

   git clone <colibri-repo>
   cd colibri
   git submodule init && git submodule update
   mkdir build && cd build
   cmake -D<CMAKE_COMPILE_OPTIONS> ..
   make

where ``<CMAKE_COMPILE_OPTIONS>`` can be (optionally) one of the following:

- ``ENABLE_OTF``: enable on-the-fly PES construction. If this option is turned on, either SCINE Sparrow, Turbomole or Orca should be available to perform electronic structure calculations.
- ``BUILD_SPARROW``: automatically build SCINE Sparrow with Scine Colibri.
- ``MPIPARALLEL``: enable MPI parallelization. Requires OpenMPI to run in parallel.
- ``SCINE_BUILD_DOCS``: Compile the doxygen documentation. Requires doxygen.

How to Cite
-----------

When publishing results obtained with Colibri, please cite the corresponding
release as archived on Zenodo (please use the DOI of the respective release).

In addition, we kindly request you to cite the following article when using Colibri:

N. Glaser, A. Baiardi, M. Reiher "Flexible DMRG-based framework for anharmonic
vibrational calculations", *J. Chem. Theory Comput.*, **2023**, *https://doi.org/10.1021/acs.jctc.3c00902*.

Support and Contact
-------------------

In case you should encounter problems or bugs, please write a short message
to scine@phys.chem.ethz.ch.

Third-Party Libraries Used
--------------------------

SCINE Colibri makes use of the following third-party libraries:

- `Boost <https://www.boost.org/>`_
- `Eigen <http://eigen.tuxfamily.org>`_
- `Google Test <https://github.com/google/googletest>`_
- `pybind11 <https://github.com/pybind/pybind11>`_
- `yaml-cpp <https://github.com/jbeder/yaml-cpp>`_

