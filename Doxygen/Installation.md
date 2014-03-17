Installation   {#PageInstallation}
=============================================

This page describes how to compile and install SusyFit.

Downloads 
------------

The SusyFit source codes are available for download: 

* Latest test release (version 0.1, released on Mar. 17, 2014): 
[<b>SusyFit-0.1.tar.gz</b>] (http://susyfit.roma1.infn.it/SusyFit-0.1.tar.gz) (381 KB)



Dependencies
------------

* **GSL**:





* **ROOT**:

    ROOT is an object oriented data analysis framework. It can be obtained
    from http://root.cern.ch/.

* **BOOST**:

    BOOST is a c++ library. It can be obtained from from http://www.boost.org.

* **BAT v0.9.3** (not required for the Library mode):

    Optionally, BAT (Bayesian Analysis Toolkit) can be obtained from
    https://www.mppmu.mpg.de/bat/ for a Monte Carlo run. 
    With the compilation option -DBAT_INSTALL=ON explained below, the SusyFit
    installation package will download and install BAT. Notice that at present 

    __NOTE:__ SusyFit is not compatible with BAT compiled with the `--enable-parallelization`
    option.

* **MPI**:

    Optionally, SusyFit can be compiled with MPI for usage in parallel 
    clusters and processpors supporting multi-threading. In this case, the SusyFit installer will patch and compile
    BAT with MPI support.

Installation Procedure
----------------------
Unpack the tarball containing the SusyFit source. A directory called 
SusyFit-x.x will be created containing the source code. To generate 
Makefiles, enter the source directory and run CMake:

~~~~~~~~~~~~~~~
  $ cd SusyFit-x.x/src  
  $ cmake . <options>
~~~~~~~~~~~~~~~

Alternatively, a directory separate from the
source directory can be made for building SusyFit (recommended, as it allows for easy deletion of the build):

~~~~~~~~~~~~~~~
  $ mkdir SusyFit-x.x/build  
  $ cd SusyFit-x.x/build  
  $ cmake ../src <options>
~~~~~~~~~~~~~~~

where the available options are:

* <b>`-DSUSYFIT_INSTALL_DIR=<SusyFit installation directory>`</b>  
  (default: `/usr/local`)

* <b>`-DBAT_INSTALL_DIR=<BAT installation directory>`</b>  
  (default: `/usr/local`)

* <b>`-DBAT_INSTALL=ON`</b>  
  to download and install BAT (default: OFF)

* <b>`-DMPIBAT=ON`</b>  
  to enable support for MPI
  (requires an implementation of MPI, default: OFF)

* <b>`-DROOT_CONFIG_PATH=<path to root-config>`</b>  

* <b>`-DBOOST_INC=<boost custom include path>/boost/`</b>  

* <b>`-DBOOST_LIB=<boost custom library path>`</b>  

* <b>`-DONLY_LIB=ON`</b>  
  not to compile the SusyFit executable (default: OFF)

Setting the option `-DBAT_INSTALL=ON`, the SusyFit installer will download, 
compile and install the BAT libraries.

__NOTE:__ Please make sure that the BAT libraries are
 not already present in the BAT installation directory. If it is present,
the installer does not build BAT, and uses the pre-installed one. 

__MPI Support:__ If you want to use SusyFit with MPI support, you can specify the option
`-DMPIBAT=ON`. This option must be accompanied with `-DBAT_INSTALL=ON` in order 
to enable the SusyFit installer to download, patch and compile BAT with MPI 
support:

~~~~~~~~~~~~~~~
  $ cmake . -DBAT_INSTALL=ON -DMPIBAT=ON <other options>
~~~~~~~~~~~~~~~

__ROOT:__ CMake checks for ROOT availability in the system and fails if ROOT is
not installed. You can specify the path to root-config using the option 
`-DROOT_CONFIG_PATH=<path to root-config>`. 

__BOOST:__ CMake also checks for BOOST availability in the system and fails if BOOST
is not installed. You can specify the path to the BOOST include files with 
`-DBOOST_INC=<boost custom include path>/boost/` and the path to the BOOST 
library with `-DBOOST_LIB=<boost custom library path>`.

The generated Makefiles are used for building an executable and a library. 
If you need only the library, you can use the option `-DONLY_LIB=ON`. For 
building only the library, BAT is not required. 

After successful CMake run, execute the build commands:

~~~~~~~~~~~~~~~
  $ make  
  $ make install
~~~~~~~~~~~~~~~

to compile and install SusyFit, where the command '`make VERBOSE=1`' enables
verbose output and '`make -j`' allows for parallel compilation. Note that depending on the setting of installation prefix
you might need root priviledges to be able to install SusyFit with
'`sudo make install`' instead of just '`make install`'.

Post Install
------------


__Executable:__ `SUSYFIT_INSTALL_DIR/bin/analysis`

__Library:__ `SUSYFIT_INSTALL_DIR/lib/libSusyFit.a`

__Combined Header:__ `SUSYFIT_INSTALL_DIR/include/SusyFit/SusyFit.h`

### Using susyfit-config

A susyfit-config script can be found in the
`SUSYFIT_INSTALL_DIR/bin/` directory, which can be invoked with the 
following options:

* <b>`--cflags`</b>
  to obtain the include path needed for compilation against the SusyFit library

* <b>`--libs`</b>
  to obtain the flags needed for linking against the SusyFit library
  
* <b>`--variable=parameters | sh`</b>
  to get a list of the mandatory model parameters sorted alphabetically and their default values
  as defined in the InputParameters class

### Examples

The example program in the `examples/LibraryMode/` directory demonstrates 
the usage of the SusyFit library, while the files in `examples/MonteCarloMode/` 
can be used for testing a Monte Carlo run with the SusyFit executable.


