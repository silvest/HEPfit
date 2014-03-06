Usage   {#PageUsage}
=============================================

The SusyFit installer generates the executable `"analysis"` and the
library `"libSufyFit.a"` along with header files including a combined
header file, SusyFit.h. One can use
the executable to perform a Bayesian statistical analysis with the
Markov Chain Monte Carlo [__Monte Carlo mode__],
or to obtain predictions of observabes for a give point in the
parameter space of the model [__Single event mode__].
Alternatively, one can use the library to obtain predictions
of obsrvables for a given point in the parameter space of the model, 
allowing our computational tool to be called from the user's own program
[__Library mode__].


Monte Carlo mode 
-----------------

The Monte Carlo analysis is performed with the [BAT library](https://www.mppmu.mpg.de/bat/). First,
a text configuration file containing a list of model parameters,
model flags and observables to be analysed has to be prepared. Another configuration
file for the Monte Carlo run has to be prepared too.

### Step 1: %Model configuration file:

A configuration file for model parameters, model flags and
obaservables are written as follows: 

~~~~~~~~~~~~~~~~
StandardModel
# Model parameters:
ModelParameter  mtop        173.2       0.9         0.  
ModelParameter  mHl         125.6       0.3         0.  

<All the model parameters have to be listed here>

# Observables:
Observable  Mw        Mw        M_{W}      80.3290 80.4064 MCMC weight 80.385 0.015 0.  
Observable  GammaW    GammaW    #Gamma_{W} 2.08569 2.09249 MCMC weight 2.085  0.042 0.  
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
# Output correlations: 
Observable2D  MwvsGammaW Mw M_{W} 80.3290 80.4064 noMCMC noweight GammaW #Gamma_{W} 2.08569 2.09249  
ModelParaVsObs MtvsMw mtop m_{t} 168.7 177.7 Mw M_{W} 80.3290 80.4064
~~~~~~~~~~~~~~~~

where the lines beginning with '#' are commented out. 

Each line has to be written as follows: 

1. The first line must be __the name of the model__ to be analyzed,
   where the available models are listed in the page @ref PageModels.<br><br>
  
2. __A model parameter__ is given in the format:

  `%ModelParameter <name> <central value> <Gaussian error> <flat error>`

  where all the parameters in a given model (see @ref PageModels) have
  to be listed in the configuration file. 

3. Optionally, one can set __a model flag__ in the format:

  `ModelFlag <name> <value>`

  where the available flags for a given model can be found in the page
  @ref PageModels.

4. __An observable__ to be computed is specified in one of the following formats:

  `%Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC weight <central value> <Gaussian error> <flat error>`<br>
  `%Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC file <filename> <histoname>`<br>
  `%Observable <name> <obs label> <histolabel> <min> <max> noMCMC noweight`

  - `<name>` can be given arbitrarily and used only internally in the codes.
  - `<obs label>` is the label of the observable (see @ref PageModels).
  - `<histolabel>` is used for the label of the output ROOT histogram,
  while `<min>` and `<max>` represent the range of the histogram.
  - noMCMC is ... 
  - Experimenta data is ... 
  - A histogram in a ROOT file is ...
  <br><br>

5. __Correlations among the data of observables__ can be taken into
   account with the line "CorrelatedGaussianObservables name Nobs",
   which initializes a set of Nobs correlated observables. It must be
   followed by exactly Nobs %Observable lines and then by Nobs lines of
   Nobs numbers for the correlation matrix. See the above example.
   One can use the keyword noMCMC and noweight, instead of MCMC and weight.
   <br><br>

6. __A correlation between two observables__ can be obtained with: 

  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> noMCMC noweight <obs2 label> <histolabel2> <min2> <max2>`<br>
  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> MCMC file <filename> <histoname> <obs2 label> <histolabel2> <min2> <max2>`

7. __A correlation between a model parameter and an observable__ can
   be obtained with:

  `%ModelParaVsObs <name> <par name> <par histolabel> <par min> <par max> <obs label> <obs histolabel> <obs min> <obs max>`


### Step 2: %Monte Carlo configuration file:

The parameters and options of the Monte Carlo run are specified in
a configuration file, separate from the one for model parameters,
etc. Each line in the file has a pair of a label
and its value, separated by space(s) or tab(s). The available
parameters and options are: 

* __NChains__    --- the number of chains (default: 5)
* __PrerunMaxIter__ --- the maximum number of iterations in the pre-run (default: 1000000)
* __Iterations__ --- the number of iterations in the MCMC run (default: 100000)
* __Seed__ --- a particular seed number for the random number generator (default: use a random seed)
* __WriteChain__ --- write Markov Chain in `MCout.root` (default: false)
* __FindModeWithMinuit__  --- find global mode with MINUIT starging
  from the best fit parameters in the MCMC run (default: false) 
* __PrintAllMarginalized__ --- print all marginalized distributions of
  the input parameters into `MonteCarlo_plots.ps` (default: false) 
* __PrintCorrelationMatrix__ --- print the correlation matrix of the
  fitted parameters into `ParamCorrelations.eps` and
  `ParamCorrelations.tes` (default: false)
* __PrintKnowledgeUpdatePlots__ --- print comparisons of the prior
  knowledge to the posterior knowledge for the parameters into
  `ParamUpdate.ps` (default: false) 
* __PrintParameterPlot__ --- print an overview plot of the parameters
  into `ParamSummary.eps` (default: false) 

For example, a Monte Carlo configuration file is written as: 

~~~~~~~~~~~~~~~~
NChains                    10
PrerunMaxIter              500000
Iterations                 1000000
WriteChain                 false
#Seed                       1
FindModeWithMinuit         false
PrintAllMarginalized       true
PrintCorrelationMatrix     true
PrintKnowledgeUpdatePlots  false
PrintParameterPlot         false
~~~~~~~~~~~~~~~~

where a '#' can be placed at the beginnig of each line to comment it
out.


### Step 3: Run:

After making the configuration files, run with the command:

~~~~~~~~~~~~~~~
  $ analysis <model conf> <Monte Carlo conf> <options>
~~~~~~~~~~~~~~~

where the available options are:

* <b>`--rootfile=<filename>`</b> output root filename without extension (default: MCout)

* <b>`--job_tag=<arg>`</b> job tag, appended to output files (default: none)

* <b>`--thRange`</b> output the minimun and maximum of theory values of observables to the file `Observables/HistoLog.txt`

### Alternative: Run with MPI.

SusyFit allows for parallel processing of the MCMC run and the observable computations.
To allow for this SusyFit and BAT has to be compiled with MPI support as explained in the
@ref PageInstallation page. The command

~~~~~~~~~~~~~~~
  $ mpiexec -n N analysis <model conf> <Monte Carlo conf> <options>
~~~~~~~~~~~~~~~

will launch analysis on `N` thread/cores/processors depending on the smallest
processing unit of the hardware used. Our MPI implementation allows for runs on multi-threaded single processors as
well as clusters with MPI support.

__NOTE:__ Our MPI implementation of SusyFit cannot be used with
BAT compiled with the `--enable-parallelization` option. It is
mandatory to use the MPI patched version of BAT as explained in the @ref PageInstallation page.


### Output Files:


Single event mode
------------------

Using the model configuration file used in the Monte Carlo mode, one
can obtain predictions of observabes for the central values of the
model parameters with the command: 

~~~~~~~~~~~~~~~
  $ analysis <model conf> --noMC
~~~~~~~~~~~~~~~

where the results are printed on the standard output. 


Library mode 
------------- 

The library mode allows for access to all the observables implemented in SusyFit
without a Monte Carlo run. The users can use one of our defined @ref PageModels and vary ModelParameters
according to their own algorithm and get the corresponding predictions for the observables. 

This is made possible through:

* a combined library: libSusyFit.a (installed in `SUSYFIT_INSTALL_DIR/lib/`
* a combined header file: SusyFit.h (installed in `SUSYFIT_INSTALL_DIR/include/SusyFit/`)


The Susyfit library allows for two different implementations of the access algorithm.

### Non-Minimal Mode:

In the non-minimal mode the user can use the SomeModel.conf file to pass the default value of
the model parameters. The following elements must be present in the user code to define
the parameters and access the observable. (For details of model paramters, observables etc. please lookup @ref PageModels.)

~~~~~~~~~~~~~~~{.c}
// Include the necessary header file. 
#include <SusyFit.h>

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

// Vary the parameters that needs to be varied in the analysis. 
DPars["Mz"] = 91.1875;
DPars["AlsMz"] = 0.1184;

// Get the map of observables with the parameter values defined above. 
DObs = CO.compute(DPars);

~~~~~~~~~~~~~~~


### Minimal Mode:

In the minimal mode the user can use the default values in InputParameters header file to define the
default values of the model parameters therefore not requiring any additional input files to be
parsed. (For details of model name, flags, paramters, observables etc. please lookup @ref PageModels.)

~~~~~~~~~~~~~~~{.c}
// Include the necessary header file. 
#include <SusyFit.h>

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

// Vary the parameters that needs to be varied in the analysis. 
DPars["mtop"] = 170.0;
DPars["mHl"] = 126.0;

// Get the map of observables with the parameter values defined above. 
DObs = CO.compute(DPars);
~~~~~~~~~~~~~~~


### Use of susyfit-config:

If `'make install'` has been done, a susyfit-config script can be found in the
`SUSYFIT_INSTALL_DIR/bin/` directory, which can be invoked with the 
following options:

~~~~~~~~~~~~~~~
Library and Library Path: susyfit-config --libs

Include Path: susyfit-config --cflags

Parameters List: susyfit-config --variable=parameters | sh
~~~~~~~~~~~~~~~

The last command lists all the mandatory parameters in all the models sorted alphabetically and their
default values as set in the class InputParameters.
