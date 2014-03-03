Installation   {#PageInstallation}
=============================================

You can find here a short description of what to do to be able to compile 
and use SusyFit on your computer.

Downloads 
------------

The SusyFit source codes are available for download: 

* Latest test release (version 0.1, released on Mar. 2, 2014): 
[<b>SusyFit-0.1.tar.gz</b>] (http://susyfit.roma1.infn.it/SusyFit-0.1.tar.gz) (381 KB)



Dependencies
------------

* **ROOT**:

    ROOT is an object oriented data analysis framework. You can obtain 
    it from http://root.cern.ch/.

* **BOOST**:

    BOOST is a c++ library. You can obtain it from http://www.boost.org.

* **BAT v0.9.3**:

    BAT (Bayesian Analysis Toolkit) can be obtained from https://www.mppmu.mpg.de/bat/. 
    With the compilation option -DBAT_INSTALL=ON explained below, the SusyFit
    installation package will download and install BAT. Notice that at present 
    SusyFit is not compatible with BAT compiled with the `--enable-parallelization` 
    option.

* **MPI**:

    Optionally, SusyFit can be compiled with MPI for usage in parallel 
    clusters. In this case, the SusyFit installer will patch and compile 
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

or, alternatively, you may separate a directory for building from the 
source directory: 

~~~~~~~~~~~~~~~
  $ mkdir SusyFit-x.x/build  
  $ cd SusyFit-x.x/build  
  $ cmake ../src <options>
~~~~~~~~~~~~~~~

where the available options are:

* <b>-DSUSYFIT_INSTALL_DIR=\<SusyFit installation directory\></b>  
  (default: `/usr/local`)

* <b>-DBAT_INSTALL_DIR=\<BAT installation directory\></b>  
  (default: `/usr/local`)

* <b>-DBAT_INSTALL=ON</b>  
  to download and install BAT (default: OFF)

* <b>-DMPIBAT=ON</b>  
  to enable support for MPI
  (requires an implementation of MPI, default: OFF)

* <b>-DROOT_CONFIG_PATH=\<path to root-config\></b>  

* <b>-DBOOST_INC=\<boost custom include path\>/boost/</b>  

* <b>-DBOOST_LIB=\<boost custom library path\></b>  

* <b>-DONLY_LIB=ON</b>  
  not to compile the SusyFit executable (default: OFF)

Setting the option -DBAT_INSTALL=ON, the SusyFit installer will download, 
compile and install the BAT library. Please make sure that the BAT library 
is not already present in the BAT installation directory. It it is present, 
the installer does not build BAT, and uses the pre-installed one. 

If you want to use SusyFit with MPI support, you can specify the option 
-DMPIBAT=ON. This option must be accompanied with -DBAT_INSTALL=ON in order 
to enable the SusyFit installer to download, patch and compile BAT with MPI 
support:

~~~~~~~~~~~~~~~
  $ cmake . -DBAT_INSTALL=ON -DMPIBAT=ON <other options>
~~~~~~~~~~~~~~~

CMake checks for ROOT availability in the system and fails if ROOT is 
not installed. You can specify the path to root-config using the option 
-DROOT_CONFIG_PATH=\<path to root-config\>. 

CMake also checks for BOOST availability in the system and fails if BOOST 
is not installed. You can specify the path to the BOOST include files with 
-DBOOST_INC=\<boost custom include path\>/boost/ and the path to the BOOST 
library with -DBOOST_LIB=\<boost custom library path\>.

The generated Makefiles are used for building an executable and a library. 
If you need only the library, you can use the option -DONLY_LIB=ON. For 
building only the library, BAT is not required. 

After successful CMake run, type the build commands:

~~~~~~~~~~~~~~~
  $ make  
  $ make install
~~~~~~~~~~~~~~~

to compile and install SusyFit, where the command '`make VERBOSE=1`' enables
verbose output. Note that depending on the setting of installation prefix
you might need root priviledges to be able to install SusyFit and run
'`sudo make install`' instead of plain '`make install`'. 

The SusyFit executable is in `SUSYFIT_INSTALL_DIR/bin/analysis`, while the
library is in `SUSYFIT_INSTALL_DIR/lib/libSusyFit.a`, which is used with
the header file `SUSYFIT_INSTALL_DIR/include/SusyFit/SusyFit.h`. 

After installation, a susyfit-config script can be found in the 
`SUSYFIT_INSTALL_DIR/bin/` directory, which can be invoked with the 
following options:

* <b>--cflags</b>
  to obtain the include path needed for compilation against the SusyFit library

* <b>--libs</b>
  to obtain the flags needed for linking against the SusyFit library

The example program in the `examples/LibraryMode/` directory demonstrates 
the usage of the SusyFit library, while the files in `examples/MonteCarloMode/` 
can be used for testing a Monte Carlo run with the SusyFit executable.


