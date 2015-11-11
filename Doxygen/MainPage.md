HEPfit - a Fitting Tool for the Standard Model and Beyond
===================================================================

version VERSIONNUMBER
-----------

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

Authors:
--------
See doc/CREDITS file for list of contributors to HEPfit.

Availability:
-------------
HEPfit is available from [the HEPfit webpage](http://hepfit.roma1.infn.it/).
See doc/COPYING and doc/LICENSE for licensing terms.

Contents:
---------
  * doc/          - directory containing information about HEPfit
  * examples-src/ - directory containing well commented example programs
  * %INSTALL.md    - information about how to install HEPfit on your system
  * %README.md     - basic information about HEPfit (this file)
  * %Usage.md      - information about how to use HEPfit

Files for building HEPfit:  

  * BAT_make_wrapper.sh.in
  * BAT_mpi_patch.txt
  * CMakeLists.txt
  * cmake_uninstall.cmake.in
  * HEPfit.h.in
  * HEPfit.pc.in
  * HEPfit_noMCMC.h.in
  * hepfit-config.in
  
Directories containing HEPfit source and header files:  

  * ComputeObservables/
  * EW/
  * Flavour/
  * gslpp/
  * HiggsExtensions/
  * InputParser/
  * LoopFunctions/
  * MonteCarlo/
  * NewPhysics/
  * Observables/
  * StandardModel/

Installation:
-------------
See the INSTALL.md file for installation instructions.

Support:
--------
For additional information and contacting the authors, please, consult
[the HEPfit webpage](http://hepfit.roma1.infn.it/).

