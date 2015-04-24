Installation   {#PageInstallation}
===================================================================

version VERSIONNUMBER
-----------

You can find here a short description of what to do to be able to compile 
and use SusyFit on your computer.


Dependencies
------------

  * GSL:  

    The GNU Scientific Library (GSL) is a library for numerical computation. 

  * ROOT v5 or later:  

    ROOT is an object oriented data analysis framework. You can obtain 
    it from http://root.cern.ch/.

  * BOOST:  

    BOOST is a C++ library. You can obtain it from http://www.boost.org.

  * BAT v0.9.4 (not required for the Library mode):  

    Optionally, BAT (Bayesian Analysis Toolkit) can be obtained from 
    https://www.mppmu.mpg.de/bat/ for a Monte Carlo run. With the compilation 
    option `-DBAT_INSTALL=ON` explained below, the SusyFit installation package 
    will download and install BAT. Notice that at present SusyFit is not
    compatible with BAT compiled with the `--enable-parallelization` option.

  * MPI:  

    Optionally, SusyFit can be compiled with MPI for usage in parallel 
    clusters and processpors supporting multi-threading. In this case,
    the SusyFit installer will patch and compile BAT with MPI support.


Installation Procedure
----------------------
Unpack the tarball containing the SusyFit source. A directory called 
SusyFit-x.x will be created containing the source code. To generate 
Makefiles, enter the source directory and run CMake:

~~~~~~~~~~~~~~~~~~~~~~  
  $ cd SusyFit-x.x  
  $ cmake . <options>  
~~~~~~~~~~~~~~~~~~~~~~

Alternatively, a directory separate from the source directory can be made for
building SusyFit (recommended, as it allows for easy deletion of the build):

~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  $ mkdir SusyFit-x.x/build  
  $ cd SusyFit-x.x/build  
  $ cmake .. <options>  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

where the available options are:

  * `-DLOCAL_INSTALL_ALL=ON`  

    to install BAT and SusyFit in the current directory (default: OFF). 
    This is equivalent in setting the options `-DCMAKE_INSTALL_PREFIX=./SusyFit`
    `-BAT_INSTALL_DIR=./BAT` and `-DBAT_INSTALL=ON`, where these variables cannot 
    be modified individually when `-DLOCAL_INSTALL_ALL=ON` is set. 

  * `-DCMAKE_INSTALL_PREFIX=<SusyFit installation directory>`  

    (default: `/usr/local`)  
  
  * `-DNOMCMC=ON`  

    to enable the No MCMC mode (default: OFF)

  * `-DDEBUG_MODE=ON`  

    to enable the debug mode (default: OFF)

  * `-DBAT_INSTALL_DIR=<BAT installation directory>`  

    (default: `/usr/local`)  

  * `-DBAT_INSTALL=ON`  

    to download and install BAT (default: OFF)

  * `-DMPIBAT=ON`  

    to enable support for MPI
    (requires an implementation of MPI, default: OFF)

  * `-DMPI_CXX_COMPILER=<path to mpi>/mpicxx`

  * `-DBOOST_INCLUDE_DIR=<boost custom include path>/boost/`  

  * `-DGSL_CONFIG_DIR=<path to gsl-config>`  

  * `-DROOT_CONFIG_DIR=<path to root-config>`  

Setting the option `-DBAT_INSTALL=ON`, the SusyFit installer will download, 
compile and install the BAT libraries.

**NOTE:**
Please make sure that the BAT libraries are not already present in the
BAT installation directory. If it is present, the installer does not
build BAT, and uses the pre-installed one. 

**No MCMC mode:**
The generated Makefiles are used for building a SusyFit library. If
you do not perform a Bayesian statistical analysis with the Markov
Chain Monte Carlo (MCMC), you can use the option `-DNOMCMC=ON`. For
this case, BAT is not required. 

**MPI Support:**
If you want to perform a MCMC run with MPI support, you can specify
the option `-DMPIBAT=ON`. This option must be accompanied with
`-DBAT_INSTALL=ON` in order to enable the SusyFit installer to
download, patch and compile BAT with MPI support:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  $ cmake . -DBAT_INSTALL=ON -DMPIBAT=ON <other options>  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**ROOT:**
CMake checks for ROOT availability in the system and fails if ROOT is
not installed. You can specify the path to root-config using the
option `-DROOT_CONFIG_DIR=<path to root-config>`. 

**BOOST:**
CMake also checks for BOOST availability in the system and fails if
BOOST is not installed. You can specify the path to the BOOST include
files with `-DBOOST_INCLUDE_DIR=<boost custom include path>/boost/`. 

After successful CMake run, execute the build commands:

~~~~~~~~~~~~~~~~~  
  $ make  
  $ make install  
~~~~~~~~~~~~~~~~~

to compile and install SusyFit, where the command '`make VERBOSE=1`'
enables verbose output and '`make -j`' allows for parallel compilation.
Note that depending on the setting of installation prefix you might
need root priviledges to be able to install SusyFit with '`sudo make
install`' instead of just '`make install`'.


Post Install
------------

**Executable:** `<CMAKE_INSTALL_PREFIX>/bin/susyfit-config`  

**Library:** `<CMAKE_INSTALL_PREFIX>/lib/libSusyFit.a`  

**Combined Header:** `<CMAKE_INSTALL_PREFIX>/include/SusyFit/SusyFit.h`  

### Using susyfit-config

A susyfit-config script can be found in the `<CMAKE_INSTALL_PREFIX>/bin/`
directory, which can be invoked with the following options:

  * `--cflags`  

    to obtain the include path needed for compilation against the SusyFit library

  * `--libs`  

    to obtain the flags needed for linking against the SusyFit library

  * `--variable=parameters | sh`  

    to get a list of the mandatory model parameters sorted alphabetically and
    their default values as defined in the InputParameters class

### Examples

The example programs can be found in the SusyFit build directory:  

  * `examples/LibMode_config/`  
  * `examples/LibMode_header/` 
  * `examples/MonteCarloMode/`  

where the first two demonstrate the usage of the SusyFit library, while 
the last one can be used for testing a Monte Carlo run with the SusyFit 
executable.


Please read the documentation at http://susyfit.roma1.infn.it/ for further 
information.

