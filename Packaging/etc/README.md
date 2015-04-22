SusyFit - a Fitting Tool for the Standard Model and Beyond
===================================================================

version VERSIONNUMBER
-----------

SusyFit is a flexible tool which, given the Standard Model or any extension,
allows to:

  - fit the model parameters to a given set of experimental observables
  - obtain fit results for observables
  - obtain predictions for observables

SusyFit can be used:

  - as an executable in Monte Carlo mode, to perform a Bayesian Markov Chain
    Monte Carlo analysis of the given model
  - as an executable in single event mode, to obtain predictions of observables
    for a given point in the parameter space of the model
  - in library mode, to obtain predictions of observables for a given point
    in the parameter space of the model, allowing our computational tool to
    be used in any statistical framework

The Monte Carlo mode is based on the BAT (Bayesian Analysis Toolkit) library.

Authors:
--------
See doc/CREDITS file for list of contributors to SusyFit.

Availability:
-------------
SusyFit is available from [the SusyFit webpage](http://susyfit.roma1.infn.it/).
See doc/COPYING and doc/LICENSE for licensing terms.

Contents:
---------
  * doc/          - directory containing information about SusyFit
  * examples-src/ - directory containing well commented example programs
  * INSTALL       - information about how to install SusyFit on your system
  * README.md     - basic information about SusyFit (this file)

Files for building SusyFit:  

  * BAT_make_wrapper.sh.in
  * BAT_mpi_patch.txt
  * CMakeLists.txt
  * cmake_uninstall.cmake.in
  * SusyFit.h.in
  * SusyFit.pc.in
  * SusyFit_noMCMC.h.in
  * susyfit-config.in
  
Directories containing SusyFit source and header files:  

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
See the INSTALL file for installation instructions.

Support:
--------
For additional information and contacting the authors, please, consult
[the SusyFit webpage](http://susyfit.roma1.infn.it/).

