Installation   {#PageInstallation}
===================================================================

version 1.0
-----------

[TOC]

The installation of `HEPfit` requires the availability of `CMake` in the system. 
A description of `CMake` and the details of how to install it can be found in the
[`CMake` website](https://cmake.org/). Most package managers for Linux distributions 
should have a `CMake` package available for installation. For Mac users, it can be 
either installed from source or from a Unix port like [Darwin Ports](https://www.macports.org/)
or [Fink](http://www.finkproject.org/), or the installation package can be downloaded 
from the [CMake website](https://cmake.org/). We list below the dependencies that need 
to be satisfied to successfully install `HEPfit`:


Dependencies
------------
    
   * `GSL`:  The GNU Scientific Library (`GSL`) is a `C` library for numerical computations. 
   It can be found on the [GSL website](http://www.gnu.org/software/gsl/). Most
   Linux package managers will have a stable version as will any ports
   for Mac. `HEPfit` is compatible with `GSL` v1.16 or greater.

   * `ROOT v5` or greater: `ROOT` is an object oriented data analysis framework. You can obtain 
   it from the [`ROOT` website](http://root.cern.ch/). `BAT` requires `ROOT v5.34.19` or greater. 
   Both `HEPfit` and `BAT` are compatible with `ROOT v6`.  __NOTE__: If `GSL` is installed 
   before compiling `ROOT` from source, then `ROOT` builds by default the `MathMore` library, which depends
   on `GSL`. Hence it is recommended to install `GSL` before installing `ROOT`.

   * `BOOST`:  `BOOST` is a `C++` library which can be obtained from the 
   [`BOOST` website](http://www.boost.org) or from Linux package managers or Mac ports. 
   `HEPfit` only requires the `BOOST` headers, not the full libraries, so a header-only
   installation is sufficient. `HEPfit` has been tested to work with `BOOST` v1.53 and greater.

   * `MPI`: Optionally, `HEPfit` can be compiled with MPI for usage in parallelized clusters and 
   processors supporting multi-threading. In this case, the `HEPfit` installer will patch and
   compile `BAT` with MPI support as described below. To this purpose one needs 
   [OpenMPI](https://www.open-mpi.org/) which is also available through package managers in 
   Linux and ports on Mac.

   * `BAT v1.0` (not required for the Library mode): The [`BAT` website](https://www.mppmu.mpg.de/bat/) 
   offers the source code for `BAT` but it should __not__ be used with `HEPfit` since a patch is required 
   to integrate `BAT` with `HEPfit`. With the compilation option `-DBAT\_INSTALL=ON` explained below, 
   the `HEPfit` installation package will download, patch and install `BAT`. The parallelized version of `BAT` 
   compatible with the parallelized version of `HEPfit` can be installed with the additional option 
   `-DMPIBAT=ON` for which MPI must be installed (see "MPI Support" below).
    
Build System
------------    
`HEPfit` uses the [`CMake` build system](https://cmake.org/) to compile all the codes and package 
it into a library (`libHEPfit.a`) and a header file (`HEPfit.h`)

QUICK Installation Instructions
----------------------

In a nutshell, if all dependencies are satisfied, for a fully `MPI` compatible MCMC capable 
`HEPfit` installation from the tarball downloaded from the [`HEPfit` website](https://hepfit.roma1.infn.it/):

~~~~~~~~~~~~~~~~~~~~~~
  $ tar xvzf HEPfit-1.0.tar.gz
  $ mkdir HEPfit-1.0/build 
  $ cd HEPfit-1.0/build 
  $ cmake .. -DLOCAL_INSTALL_ALL=ON -DMPIBAT=ON
  $ make
  $ make install
~~~~~~~~~~~~~~~~~~~~~~

To run your first example:

~~~~~~~~~~~~~~~~~~~~~~
$ cd examples/MonteCarloMode/
$ make  
$ mpiexec -n 5 ./analysis ../config/StandardModel.conf MonteCarlo.conf
~~~~~~~~~~~~~~~~~~~~~~

This is all you will need for running a MCMC on 5 cores with the model, parameters and observables 
specified in the configuration files in `examples/config` directory with `HEPfit`. For variations 
please read what follows.

Detailed Installation Instructions
----------------------

Unpack the tarball containing the `HEPfit` source which you can obtain from the 
[`HEPfit` website](https://hepfit.roma1.infn.it/). A directory called 
`HEPfit-1.0` will be created containing the source code. To generate 
Makefiles, enter the source directory and run `CMake`:

~~~~~~~~~~~~~~~~~~~~~~
  $ cd HEPfit-1.0  
  $ cmake . <options>  
~~~~~~~~~~~~~~~~~~~~~~

**RECOMMENDED:** Alternatively, a directory separate from the source directory can be made for
building `HEPfit` (recommended as it allows for easy deletion of the build):

~~~~~~~~~~~~~~~~~~~~~~
  $ mkdir HEPfit-1.0/build  
  $ cd HEPfit-1.0/build  
  $ cmake .. <options>  
~~~~~~~~~~~~~~~~~~~~~~

where the available options are:

   * `-DLOCAL\_INSTALL\_ALL=ON`: to install `BAT` and `HEPfit` in the current directory (default: OFF). 
    This is equivalent to setting the combination of the options:
~~~~~~~~~~~~~~~~~~~~~~
-DCMAKE_INSTALL_PREFIX=./HEPfit -DBAT_INSTALL_DIR=./BAT -DBAT_INSTALL=ON
~~~~~~~~~~~~~~~~~~~~~~
   These variables cannot be modified individually when `-DLOCAL\_INSTALL\_ALL=ON` is set. 

   * `-DCMAKE\_INSTALL\_PREFIX=<HEPfit installation directory>`: the directory in which `HEPfit` 
   will be installed (default: `/usr/local`).  
  
   * `-DNOMCMC=ON`: to enable the mode without MCMC (default: OFF).

   * `-DDEBUG\_MODE=ON`: to enable the debug mode (default: OFF).

   * `-DBAT\_INSTALL\_DIR=<BAT installation directory>`: (default: `/usr/local`). This option 
   is overridden by `-DLOCAL\_INSTALL\_ALL=ON` .

   * `-DBAT\_INSTALL=ON` to download and install `BAT`. This is relevant only if
`-DNOMCMC=ON` is not set. Use `-DBAT\_INSTALL=OFF` only if you know your `BAT` installation
is already patched by `HEPfit` and is with or without MPI support as needed. (default: ON).

   * `-DMPIBAT=ON`: to enable support for MPI for both `BAT` and `HEPfit`
    (requires an implementation of MPI, default: OFF).

   * `-DMPI\_CXX\_COMPILER=<path to mpi>/mpicxx`: You can specify the MPI compiler with this option.

   * `-DBOOST\_INCLUDE\_DIR=<boost custom include path>/boost/`: if BOOST is not installed in the 
   search path then you can specify where it is with this option. The path must end with the `boost/` 
   directory which contains the headers.

   * `-DGSL\_CONFIG\_DIR=<path to gsl-config>`: `HEPfit` used `gsl-config` to get the GSL parameters. 
   If this is not in the search path, you can specify it with this option. 

   * `-DROOT\_CONFIG\_DIR=<path to root-config>`: `HEPfit` used `root-config` to get the `ROOT` parameters. 
   If this is not in the search path, you can specify it with this option. 

   * `-DINTEL\_FORTRAN=ON`: If you are compiling with INTEL compilers then this flag turns 
   on support for the compilers (default: OFF).

Setting the option `-DBAT_INSTALL=ON`, the `HEPfit` installer will download, 
compile and install the `BAT` libraries.

**NOTE:**  
If `BAT` libraries and headers are present in target directory for `BAT` they will be overwritten unless 
`-DBAT\_INSTALL=OFF` is set. This is done so that the correct patched version of `BAT` compatible with
`HEPfit` gets installed.

**No MCMC mode:**  
The generated Makefiles are used for building a `HEPfit` library. If
you do not perform a Bayesian statistical analysis with the MCMC, you can use the option 
`-DNOMCMC=ON`. For this case, `BAT` is not required.

`MPI` **Support:**  
If you want to perform an MCMC run with `MPI` support, you can specify
the option `-DMPIBAT=ON`. This option must not be accompanied with
`-DBAT_INSTALL=OFF` in order to enable the `HEPfit` installer to
download, patch and compile `BAT` and build `HEPfi` with MPI support:

~~~~~~~~~~~~~~~~~~~~~~
  $ cmake . -DMPIBAT=ON <other options>  
~~~~~~~~~~~~~~~~~~~~~~

Please make sure non-MPI `BAT` libraries are not already installed in the 
`BAT installation directory.

`ROOT`:  
`CMake` checks for `ROOT` availability in the system and fails if `ROOT` is
not installed. You can specify the path to `root-config` using the
option `-DROOT_CONFIG_DIR=<path to root-config>`.

`BOOST`:  
`CMake` also checks for `BOOST` availability in the system and fails if
`BOOST` is not installed. You can specify the path to the `BOOST` include
files with `-DBOOST_INCLUDE_DIR=<boost custom include path>/boost/`.

The **recommended installation** flags for a locally installed `HEPfit` with full MPI and 
MCMC support is:

~~~~~~~~~~~~~~~~~~~~~~
  $ cmake . -DLOCAL_INSTALL_ALL=ON -DMPIBAT=ON
~~~~~~~~~~~~~~~~~~~~~~

This will enable easy portability of all codes and easy upgrading to future version as nothing 
will be installed system wide.  Also, this is useful if you do not have root access and cannot 
install software in system folders. After successful CMake run, execute the build commands:

~~~~~~~~~~~~~~~~~~~~~~
$ make  
$ make install  
~~~~~~~~~~~~~~~~~~~~~~

to compile and install `HEPfit`, where using the command `make VERBOSE=1`
enables verbose output and `make -j` allows for parallel compilation (not recommended if you are
building `BAT` with `HEPfit`). Note that depending on the setting of installation prefix you might
need root privileges to be able to install `HEPfit` with `sudo make
install` instead of just `make install`.


Post Installation
------------
After the completion fo the installation with `make install` the following three files can be found 
in the installation location. The file `libHEPfit.h` is a combined header file corresponding to the 
library `libHEPfit.a`.

**Executable:** `<CMAKE_INSTALL_PREFIX>/bin/hepfit-config`  
**Library:** `<CMAKE_INSTALL_PREFIX>/lib/libHEPfit.a`  
**Combined Header:** `<CMAKE_INSTALL_PREFIX>/include/HEPfit/HEPfit.h`  

**Using hepfit-config:**
A `hepfit-config` script can be found in the `<CMAKE_INSTALL_PREFIX>/bin/`
directory, which can be invoked with the following options:

   * `--cflags`: to obtain the include path needed for compilation against the `HEPfit` library
   * `--libs`: to obtain the flags needed for linking against the `HEPfit` library

**Examples:**  
The example programs can be found in the `HEPfit` build directory:  

   * `examples/LibMode_config/`  
   * `examples/LibMode_header/` 
   * `examples/MonteCarloMode/`
   * `examples/EventGeneration/`
   * `examples/myModel/`

The first two demonstrate the usage of the `HEPfit` library, while 
the third one can be used for testing a Monte Carlo run with the `HEPfit` 
library. The fourth example can be used to generate values of observables with 
a sample of parameters drawn from the parameter space. The fifth one is an example 
implementation of a custom model and custom observables. To make an executable to 
run these examples:

~~~~~~~~~~~~~~~~~~~~~~
  $ cd examples/MonteCarloMode/
  $ make  
~~~~~~~~~~~~~~~~~~~~~~

This will produce an executable called `analysis` in the current directory that can be used 
to run `HEPfit`. The details are elaborated on in the [Usage Page](@ref PageUsage).



