Usage   {#PageUsage}
=============================================

The HEPfit installer generates the 
library `"libSufyFit.a"` along with header files including a combined
header file, HEPfit.h. To perform a Bayesian statistical analysis with the
Markov Chain Monte Carlo [__Library mode with MCMC__] one can use the given example implementation.
Alternatively, one can use the library to obtain predictions
of observables for a given point in the parameter space of the model, 
allowing our computational tool to be called from the user's own program
[__Library mode without MCMC__].


Monte Carlo mode 
-----------------

The Monte Carlo analysis is performed with the [BAT library](https://www.mppmu.mpg.de/bat/). First,
a text configuration file containing a list of model parameters,
model flags and observables to be analyzed has to be prepared. Another configuration
file for the Monte Carlo run has to be prepared, too.

### Step 1: %Model configuration file:

A configuration file for model parameters, model flags and
observables are written as follows: 

~~~~~~~~~~~~~~~~
StandardModel
# Model parameters:
ModelParameter  mtop        173.2       0.9         0.  
ModelParameter  mHl         125.6       0.3         0.  
...
CorrelatedGaussianParameters   V1_lattice 2
ModelParameter  a_0V    0.496   0.067   0.
ModelParameter  a_1V   -2.03    0.92    0.
1.00    0.86
0.86    1.00

<All the model parameters have to be listed here>

# Observables:
Observable  Mw        Mw        M_{W}      80.3290 80.4064 MCMC weight 80.385 0.015 0.  
Observable  GammaW    GammaW    #Gamma_{W} 2.08569 2.09249 MCMC weight 2.085  0.042 0.  
#
# Correlated observables:
CorrelatedGaussianObservables Zpole2 7  
Observable  Alepton   Alepton   A_{l}      0.143568 0.151850 MCMC weight 0.1513  0.0021  0.  
Observable  Rbottom   Rbottom   R_{b}      0.215602 0.215958 MCMC weight 0.21629 0.00066 0.  
Observable  Rcharm    Rcharm    R_{c}      0.172143 0.172334 MCMC weight 0.1721  0.0030  0.  
Observable  AFBbottom AFBbottom A_{FB}^{b} 0.100604 0.106484 MCMC weight 0.0992  0.0016  0.  
Observable  AFBcharm  AFBcharm  A_{FB}^{c} 0.071750 0.076305 MCMC weight 0.0707  0.0035  0.  
Observable  Abottom   Abottom   A_{b}      0.934320 0.935007 MCMC weight 0.923   0.020   0.  
Observable  Acharm    Acharm    A_{c}      0.666374 0.670015 MCMC weight 0.670   0.027   0.  
1.00   0.00   0.00   0.00   0.00   0.09   0.05
0.00   1.00  -0.18  -0.10   0.07  -0.08   0.04
0.00  -0.18   1.00   0.04  -0.06   0.04  -0.06
0.00  -0.10   0.04   1.00   0.15   0.06   0.01
0.00   0.07  -0.06   0.15   1.00  -0.02   0.04
0.09  -0.08   0.04   0.06  -0.02   1.00   0.11
0.05   0.04  -0.06   0.01   0.04   0.11   1.00  
#
# Output correlations: 
Observable2D  MwvsGammaW Mw M_{W} 80.3290 80.4064 noMCMC noweight GammaW #Gamma_{W} 2.08569 2.09249  
...
Observable2D Bd_Bsbar_mumu noMCMC noweight
Observable   BR_Bdmumu      BR(B_{d}#rightarrow#mu#mu)    1. -1.  1.05e-10   0.   0.
Observable   BRbar_Bsmumu   BR(B_{s}#rightarrow#mu#mu)    1. -1.  3.65e-9    0.   0.
...
Observable2D S5_P5 noMCMC noweight
BinnedObservable   S_5      S_5    1. -1.  0.   0.   0.   4.   6.
BinnedObservable   P_5      P_5    1. -1.  0.   0.   0.   4.   6.
#
# Including other configuration files
IncludeFile Flavour.conf
~~~~~~~~~~~~~~~~

where the lines beginning with '#' are commented out. 

Each line has to be written as follows: 

1. The first line must be __the name of the model__ to be analyzed,
   where the available models are listed in the page @ref PageModels.<br><br>
  
2. __A model parameter__ is given in the format:

  `%ModelParameter <name> <central value> <Gaussian error> <flat error>`

  where all the parameters in a given model (see @ref PageModels) have
  to be listed in the configuration file.

3. __A set of correlated model parameter__ is specified with "CorrelatedGaussianParameters name Npar"
   which initializes a set of Npars correlated parameters. It must be
   followed by exactly Npar %Parameter lines and then by Npar lines of
   Npar numbers for the correlation matrix. See the above example.

4. Optionally, one can set __a model flag__ in the format:

  `ModelFlag <name> <value>`

  where the available flags for a given model can be found in the page
  @ref PageModels.

5. __An observable__ to be computed is specified in one of the following formats:

  `%Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <central value> <Gaussian error> <flat error>`<br>
  `%Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC file <filename> <histoname>`<br>
  `%Observable <name> <obs label> <histolabel> <min> <max> noMCMC noweight`

  - `<name>` is user given name for different observables which must be unique for each observable
  - `<obs label>` is the theory label of the observable (see @ref PageModels).
  - `<histolabel>` is used for the label of the output ROOT histogram,
  while `<min>` and `<max>` represent the range of the histogram.
  - (no)MCMC is the flag specifying whether the observable should be included in the Monte Carlo prerun or not 
  - (no)weight specifies if the observable weight will be computed or not. If weigh is specified with noMCMC
    then the weight of the observable will be stored in the MCout*.root file
  - Experimental data is specified as `central value  Gaussian error  flat error`. If weight is specified
    both the errors cannot be 0.
  - A histogram in a ROOT file is specifed by the name of the root file (filename) and then the
    name of the histogram (histoname).
  <br><br>

6. __Correlations among the data of observables__ can be taken into
   account with the line "CorrelatedGaussianObservables name Nobs",
   which initializes a set of Nobs correlated observables. It must be
   followed by exactly Nobs %Observable lines and then by Nobs lines of
   Nobs numbers for the correlation matrix. See the above example.
   One can use the keyword noMCMC and noweight, instead of MCMC and weight.
   <br><br>

7. __A correlation between two observables__ can be obtained with: 

  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> noMCMC noweight <obs2 label> <histolabel2> <min2> <max2>`<br>
  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> MCMC file <filename> <histoname> <obs2 label> <histolabel2> <min2> <max2>`
    <br> or <br>
  `%Observable2D <name> (no)MCMC (no)weight`<br>
  `(Binned)%Observable <obs label 1> <histolabel 1> <min> <max> <central value> <Gaussian error> <flat error> (<bin_min> <bin_max>)`<br>
  `(Binned)%Observable <obs label 2> <histolabel 2> <min> <max> <central value> <Gaussian error> <flat error> (<bin_min> <bin_max>)`<br>
  `%Observable2D <name> MCMC file filename histoname`<br>
  `(Binned)%Observable <obs label 1> <histolabel 1> <min> <max> (<bin_min> <bin_max>)`<br>
  `(Binned)%Observable <obs label 2> <histolabel 2> <min> <max> (<bin_min> <bin_max>)`<br>

8. __Include configuration files__ with the IncludeFile directive. This is useful if one
   wants to separate the input configurations for better organization and flexibility.


### Step 2: %Monte Carlo configuration file:

The parameters and options of the Monte Carlo run are specified in
a configuration file, separate from the one for model parameters,
etc. Each line in the file has a pair of a label
and its value, separated by space(s) or tab(s). The available
parameters and options are: 

* __NChains__ : The number of chains in the Monte Carlo run. A minimum of 5 is suggested (default). 
If the theory space is complicated and/or the number of parameters are large then more
is necessary. The amount of statistics collected in the main run is proportional to the
numbe rof chains.
* __PrerunMaxIter__ : The maximum number of iterations that the prerun will go through (default: 1000000). 
The prerun ends automatically when the chains converge (R < 1.1) and all efficiencies 
are adjusted. While it is not necessary for the prerun to converge for a run to be 
successful, one should exercise caution in this case.
* __Iterations__ : The number of iterations in the main run. This run is for the purpose of
collecting statistics and is at the users discretion. (default: 100000)
* __Seed__ : The seed can be fixed for deterministic runs.
* __PrintAllMarginalized__ : All marginalized distributions will be printed in a pdf file
(MonteCarlo_plots_*.pdf).
* __PrintCorrelationMatrix__ : The parametric corellation will be printer in ParamCorrelations*.pdf
and ParamCorrelations*.tex.
* __PrintKnowledgeUpdatePlots__ : A comparison between prior and posterior knowledge will be 
printed in a plot stored in ParamUpdate*.pdf.
* __PrintParameterPlot__ : As summary of the paramters will be printed in ParamSummary*.pdf.
* __OrderParameters__ : This determines whether the parameters will be randomized one
at a time or as a whole set in the entire run. While for quick estimates it is not necessary
to order the prameters (false), it is suggested to do so for more accurate estimates.


* __FindModeWithMinuit__ : To find global mode with MINUIT starting
from the best fit parameters in the MCMC run.
* __MinimumEfficiency__  : This allows the setting of the minimum eifficiency of 
all the prameters. The default is 0.15.
* __WriteChain__         : The chains will be written in the root file MCout*.root. 
This can be used for analyzing the performance of the chains.
* __CalculateNormalization__ : The normalization of the theory space will be calculated
at the end of the Monte Carlo run. This is useful for model comparison.
* __WritePreRunData__ : The prerun data is written to a file for rerun.      
* __ReadPreRunData__  : The prerun data will be read from a previously stored prerun file. 

For example, a Monte Carlo configuration file is written as: 

~~~~~~~~~~~~~~~~
NChains                    10    (Default: 5)
PrerunMaxIter              50000 (Default: 1000000)
Iterations                 10000 (Default: 100000)
#Seed                       1
PrintAllMarginalized       true  (Default: false)
PrintCorrelationMatrix     true  (Default: false)    
PrintKnowledgeUpdatePlots  false (Default: false)
PrintParameterPlot         false (Default: false)
OrderParameters            true  (Default: true)
~~~~~~~~~~~~~~~~

where a '#' can be placed at the beginning of each line to comment it
out. These are the most commonly used settings. Some less commonly used settings:

~~~~~~~~~~~~~~~~
FindModeWithMinuit      (Default: false)
MinimumEfficiency       (Default: 0.15, set to 0. - 1.)
WriteChain              (Default: false)
CalculateNormalization  (Default: false)
WritePreRunData         (Mandatory: name of file)
ReadPreRunData          (Read existing prerun data file)
~~~~~~~~~~~~~~~~

     


### Step 3: Run:

Library mode with MCMC 
------------------------- 

An example can be found in examples/MonteCarloMode

~~~~~~~~~~~~~~~
  $ cd examples/MonteCarloMode
  $ make
~~~~~~~~~~~~~~~

After making the configuration files, run with the command:

~~~~~~~~~~~~~~~
  $ analysis <model conf> <Monte Carlo conf>
~~~~~~~~~~~~~~~

### Alternative: Run with MPI.

HEPfit allows for parallel processing of the MCMC run and the observable computations.
To allow for this HEPfit and BAT has to be compiled with MPI support as explained in the
@ref PageInstallation page. The command

~~~~~~~~~~~~~~~
  $ mpiexec -n N analysis <model conf> <Monte Carlo conf>
~~~~~~~~~~~~~~~

will launch analysis on `N` thread/cores/processors depending on the smallest
processing unit of the hardware used. Our MPI implementation allows for runs on multi-threaded single processors as
well as clusters with MPI support.

__NOTE:__ Our MPI implementation of HEPfit cannot be used with
BAT compiled with the `--enable-parallelization` option. It is
mandatory to use the MPI patched version of BAT as explained in the @ref PageInstallation page.


### Output Files:

* log.txt: This is the log file.
* MCout.root: This is the root file containing all the information of the run and
  the histogram objects.
* MonteCarlo_results.txt: This file contains the fits to the parameters.
* MonteCarlo_plots.pdf: This File contains the histograms for the parameters.
* Observables: This directory will contain the histograms of all the observables
  specified in the config file.
* Observables/HistoLog.txt: This file contains the information on over-run and under-run
  of the histogram filling.
* Observables/Statistics.txt: This file contains the compilation of the statistics
  extracted from the histograms.

Other files might be generated depending on the options specified in the Monte Carlo 
configuration file.


Single event mode
------------------

Using the model configuration file used in the Monte Carlo mode, one
can obtain predictions of observabes for the central values of the
model parameters with the command: 

~~~~~~~~~~~~~~~
  $ analysis <model conf> --noMC
~~~~~~~~~~~~~~~

where the results are printed on the standard output. 


Library mode without MCMC 
------------------------- 

The library mode allows for access to all the observables implemented in HEPfit
without a Monte Carlo run. The users can use one of our defined @ref PageModels and vary ModelParameters
according to their own algorithm and get the corresponding predictions for the observables. 

This is made possible through:

* a combined library: libHEPfit.a (installed in `HEPFIT_INSTALL_DIR/lib/`
* a combined header file: HEPfit.h (installed in `HEPFIT_INSTALL_DIR/include/HEPfit/`)


The HEPfit library allows for two different implementations of the access algorithm.

### Non-Minimal Mode:

In the non-minimal mode the user can use the SomeModel.conf file to pass the default value of
the model parameters. The following elements must be present in the user code to define
the parameters and access the observable. (For details of model parameters, observables etc. please lookup @ref PageModels.)

~~~~~~~~~~~~~~~{.c}
// Include the necessary header file. 
#include <HEPfit.h>

// Define the model configuration file. 
std::string ModelConf = "SomeModel.conf";

// Define a map for the observables. 
std::map<std::string, double> DObs;

// Define a map for the parameters to be varied.
std::map<std::string, double> DPars;

// Initialize the observables to be returned. 
DObs["Mw"] = 0;
DObs["GammaZ"] = 0.;
DObs["AFBbottom"] = 0.;

// Create and object of the class ComputeObservables. 
ComputeObservables CO(ModelConf, DObs);

// Vary the parameters that need to be varied in the analysis. 
DPars["Mz"] = 91.1875;
DPars["AlsMz"] = 0.1184;

// Get the map of observables with the parameter values defined above. 
DObs = CO.compute(DPars);

~~~~~~~~~~~~~~~


### Minimal Mode:

In the minimal mode the user can use the default values in InputParameters header file to define the
default values of the model parameters therefore not requiring any additional input files to be
parsed. (For details of model name, flags, parameters, observables etc. please lookup @ref PageModels.)

~~~~~~~~~~~~~~~{.c}
// Include the necessary header file. 
#include <HEPfit.h>

// Define a map for the observables. 
std::map<std::string, double> DObs;

// Define a map for the mandatory model parameters used during initializing a model. 
std::map<std::string, double> DPars_IN;

// Define a map for the parameters to be varied. 
std::map<std::string, double> DPars;

// Define a map for the model flags. 
std::map<std::string, std::string> DFlags;

// Define the name of the model to be used. 
std::string ModelName = "NPZbbbar";

// Create and object of the class InputParameters. 
InputParameters IP;

// Read a map for the mandatory model parameters. (Default values in InputParameters.h) 
DPars_IN = IP.getInputParameters(ModelName);

// Change the default values of the mandatory model parameters if necessary. 
// This can also be done with Dpars after creating an object of ComputeObservables 
DPars_IN["mcharm"] = 1.3;
DPars_IN["mub"] = 4.2;

// Initialize the observables to be returned. 
DObs["Mw"] = 0;
DObs["GammaZ"] = 0.;
DObs["AFBbottom"] = 0.;

// Initialize the model flags to be set. 
DFlags["NPZbbbarLR"] = "TRUE";

// Create and object of the class ComputeObservables. 
ComputeObservables CO(ModelName, DPars_IN, DObs);

// Set the flags for the model being used, if necessary. 
CO.setFlags(DFlags);

// Vary the parameters that need to be varied in the analysis. 
DPars["mtop"] = 170.0;
DPars["mHl"] = 126.0;

// Get the map of observables with the parameter values defined above. 
DObs = CO.compute(DPars);
~~~~~~~~~~~~~~~


### Use of hepfit-config:

If `'make install'` has been done, a hepfit-config script can be found in the
`HEPFIT_INSTALL_DIR/bin/` directory, which can be invoked with the 
following options:

~~~~~~~~~~~~~~~
Library and Library Path: hepfit-config --libs

Include Path: hepfit-config --cflags

Parameters List: hepfit-config --variable=parameters | sh
~~~~~~~~~~~~~~~

The last command lists all the mandatory parameters in all the models sorted alphabetically and their
default values as set in the class InputParameters.
