HEPfit - a Code for the Combination of Indirect and Direct Constraints on High Energy Physics Models.
===================================================================

HEPfit is a flexible tool which, given the Standard Model or any extension,
allows to:

  - fit the model parameters to a given set of experimental observables
  - obtain fit results for observables
  - obtain predictions for observables

The HEPfit library can be used:

  - to perform a Bayesian Markov Chain Monte Carlo analysis of the implemented model using BAT support
  - to perform a Bayesian Markov Chain Monte Carlo analysis of user defined models and/or observables using BAT support
  - to perform a user defined statistical analysis with the models and observables implemented in the library
  - to perform a user defined statistical analysis with user defined models and/or observables

The Markov Chain Monte Carlo is based on the [BAT (Bayesian Analysis Toolkit) library](https://www.mppmu.mpg.de/bat/).

Installation (for developers):
--------
- Public libraries
  - [gsl](https://www.gnu.org/software/gsl/)
  - [ROOT](https://root.cern/install/)
  - [boost](https://www.boost.io/doc/user-guide/getting-started.html)
  - [OpenMPI](https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html)
- Developer Libraries
  - [RGESolver](https://github.com/silvest/RGESolver)
    ```{sh}
    git clone git@github.com:silvest/RGESolver.git --recursive
    cd RGESolver
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=<PATH>
    make -j
    make install
    ```
  - [bat_MPI](https://github.com/silvest/bat_MPI)
    ```{sh}
    git clone git@github.com:silvest/bat_MPI.git
    cd bat_MPI
    ./configure CXX=mpic++ --prefix=<PATH>
    make -j
    make install
    ```
- HEPfit (using CMake)
  ```{sh}
  git clone git@github.com:silvest/HEPfit.git --recursive
  cd HEPfit
  mkdir build
  cd build
  cmake .. -DBAT_INSTALL=OFF
  ```
For install options for ```bat_MPI```, ```RGESolver```, and ```HEPfit``` look at the top of the CMakeLists.txt file in the root directory


Authors:
--------
See doc/CREDITS file for list of contributors to HEPfit.

Availability:
-------------
The HEPfit distribution package is available from the [HEPfit webpage](https://hepfit.roma1.infn.it/).
See doc/COPYING and doc/LICENSE for licensing terms.

Documentation:
--------------
The documentation of the master branch is available from the [HEPfit webpage](https://hepfit.roma1.infn.it/doc/master/index.html). It is meant for the developers and is updated periodically. Documentation can also be generated from the master using the pacakging script in the `Packaging` directory by running
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
$ cd Packaging
$ sh makepackage.sh --doxygen
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The documentation will be generated in `Packaging/HEPfit-master-<commit no.>/Doxygen/html`. An installation of [Doxygen](https://www.doxygen.nl/) is necessary.
