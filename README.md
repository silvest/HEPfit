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
