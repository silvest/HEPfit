Overview
===================================================================

version 1.0
-----------

[TOC]

`HEPfit` is a flexible tool which, given the Standard Model or any extension,
allows to:

  - fit the model parameters to a given set of experimental observables
  - obtain fit results for observables
  - obtain predictions for observables

The `HEPfit` library can be used:

  - to perform a Bayesian Markov Chain Monte Carlo analysis of the implemented model using `BAT` support
  - to perform a Bayesian Markov Chain Monte Carlo analysis of user defined models and/or observables using `BAT` support
  - to perform a user defined statistical analysis with the models and observables implemented in the library
  - to perform a user defined statistical analysis with user defined models and/or observables

The Markov Chain Monte Carlo is based on the 
[`BAT` (Bayesian Analysis Toolkit)](https://www.mppmu.mpg.de/bat/) library.

Authors:
--------
See the [developer section](https://hepfit.roma1.infn.it/developers.html) of 
the [`HEPfit` webpage](http://hepfit.roma1.infn.it/).

Availability:
-------------
`HEPfit` is available from the [`HEPfit` webpage](http://hepfit.roma1.infn.it/).
See doc/COPYING and doc/LICENSE for licensing terms.

Contents:
---------
  * doc/           - directory containing information on `HEPfit`
  * examples-src/  - directory containing well commented example programs
  * %INSTALL.md    - information about how to install HEPfit on your system
  * %README.md     - basic information about `HEPfit` (this package)
  * %Usage.md      - information about how to use HEPfit
  * VERSION        - version number for `HEPfit`

Files for building HEPfit:  

  * BAT-1.0.0_mpi_patch.txt
  * BAT-1.0.0_patch.txt
  * BAT_make_wrapper.sh.in
  * CMakeLists.txt
  * cmake_uninstall.cmake.in
  * HEPfit.h.in
  * HEPfit.pc.in
  * HEPfit_noMCMC.h.in
  * hepfit-config.in

  
Directories containing HEPfit source and header files:  

  * ComputeObservables/
  * EW/
  * EventGeneration/
  * Flavour/
  * FlavourWilsonCoefficient/
  * GeneralTHDM/
  * GeorgiMachacek/
  * gslpp
  * InputParser/
  * LeptonFlavour/
  * LoopFunctions/
  * MonteCarlo/
  * NewPhysics/
  * Observables/
  * SUSY/
  * StandardModel/
  * THDM/
  * THDMW/

Installation:
-------------
See the %INSTALL.md file for installation instructions or the [`HEPfit` documentation website](http://hepfit.roma1.infn.it/doc/v1.0/).

Support:
--------
For additional information and contacting the authors, please consult
the [`HEPfit` webpage](http://hepfit.roma1.infn.it/).

