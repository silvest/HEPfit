Usage   {#PageUsage}
=============================================

[TOC]

After the `HEPfit` installer generates the library `libHEPFit.a` along with header files included in a combined
header file, `HEPfit.h`, the given example implementation can be used to perform a MCMC based Bayesian statistical analysis.
Alternatively, the library can be used to obtain predictions of observables for a given point in the parameter space of a model, 
allowing `HEPfit` to be called from the user's own program. We explain both methods below. In addition `HEPfit` provides
the ability to the user to define custom models and observables as explained below. We give a brief description on how to
get started with custom models and observables.

Monte Carlo mode 
-----------------

The Monte Carlo analysis is performed with the [BAT library](https://www.mppmu.mpg.de/bat/). First,
a text configuration file containing a list of model parameters,
model flags and observables to be analyzed has to be prepared. Another configuration
file for the Monte Carlo run has to be prepared, too.

### Step 1: Model configuration file:

The configuration files are the primary way to control the behaviour of the code and to detail its input and output. 
While a lot of checks have been implemented in `HEPfit` to make sure the configuration files are of the right format, 
it is not possible to make it error-proof. Hence,  care should be taken in preparing these files. A configuration file 
for model parameters, model flags, and observables is written as follows: 

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
   
2. __Model flags__, if necessary, should be specified right after the model because some of them can control the way 
   the input parameters are read.   
~~~~~~~~~~~~~~~~
      ModelFlag <name> <value>
~~~~~~~~~~~~~~~~   

3. __A model parameter__ is given in the format:
~~~~~~~~~~~~~~~~
    %ModelParameter <name> <central value> <Gaussian error> <flat error>
~~~~~~~~~~~~~~~~
    where all the parameters in a given model (see @ref PageModels) have
    to be listed in the configuration file.  

4. __A set of correlated model parameter__ is specified with "
~~~~~~~~~~~~~~~~
    CorrelatedGaussianParameters name Npar
~~~~~~~~~~~~~~~~
   which initializes a set of `Npars` correlated parameters. It must be
   followed by exactly `Npar` %Parameter lines and then by `Npar` lines of
   `Npar` numbers for the correlation matrix. See the example above.

5. __An observable__ to be computed is specified in one of the following formats:
~~~~~~~~~~~~~~~~
    Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <central value> <Gaussian error> <flat error>
    Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC file <filename> <histoname>
    Observable <name> <obs label> <histolabel> <min> <max> noMCMC noweight
    Observable <name> <obs label> <histolabel> <min> <max> noMCMC writeChain
~~~~~~~~~~~~~~~~
    
   - `<name>` is a user given name for different observables which must be unique for each observable.
   - `<obs label>` is the theory label of the observable (see the online documentation).
   - `<histolabel>` is used for the label of the output `ROOT` histogram, while `<min>` and `<max>` represent the
       range of the histogram (if `<min>` <= `<max>` the range of the histogram is set automatically).
   - `(no)MCMC` is the flag specifying whether the observable should be included in likelihood used for the MCMC sampling.
   - `(no)weight` specifies if the observable weight will be computed or not. If `weight` is specified with `noMCMC`
     then a chain containing the weights for the observable will be stored in the `MCout*.root` file.  
   - `noMCMC noweight` is the combination to be used to get a prediction for an observable.
   - When the `weight` option is specified, at least one of the `<Gaussian error>` or the `<flat error>` must be
     nonvanishing, and the `<central value>` must of course be specified.
   - When using the `file` option, a histogram in a `ROOT` file must be specified by the name of the `ROOT` file 
   (`filename`) and then the name of the histogram (`histoname`) in the file (including, if needed, the directory).
   - The `writeChain` option allows one to write all the values of the observable generated during the main run of the 
        MCMC into the `ROOT` file.  

6. __A `BinnedObservable`__ is similar in construction to an `Observable` but with two extra arguments specifying the 
    upper and lower limit of the bin:
~~~~~~~~~~~~~~~~
    BinnedObservable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <central value> <Gaussian error> <flat error> <bin_min>       <bin_max>
    BinnedObservable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <filename> <histoname> <bin_min> <bin_max>
    BinnedObservable <name> <obs label> <histolabel> <min> <max> noMCMC writeChain 0. 0. 0. <bin_min> <bin_max>
~~~~~~~~~~~~~~~~
    Because of the order of parsing the `<central value> <Gaussian error> <flat error>` cannot be dropped out even in the
     `noMCMC noweight` case for a `BinnedObservable`.  
     
7. __A `FunctionObservable`__ is the same as a `BinnedObservable` but with only one extra argument that points to the 
    value at which the function is computed:
~~~~~~~~~~~~~~~~
    FunctionObservable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <central value> <Gaussian error> <flat error> <x_value>
    FunctionObservable <name> <obs label> <histolabel> <min> <max> (no)MCMC  (no)weight <filename> <histoname> <x_value>
    FunctionObservable <name> <obs label> <histolabel> <min> <max> noMCMC writeChain 0. 0. 0. <x_value>
~~~~~~~~~~~~~~~~     

8. __An asymmetric Gaussian__ constraint can be set using `AsyGausObservable`:
~~~~~~~~~~~~~~~~
    AsyGausObservable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight <central value> <left_error> <right_error>
~~~~~~~~~~~~~~~~  
   
9. __Correlations among observables__ can be taken into account with the line `CorrelatedGaussianObservables name Nobs`,
      which initializes a set of `Nobs` correlated observables. It must be followed by exactly `Nobs` `Observable` lines
       and then by `Nobs` rows of `Nobs` numbers for the correlation matrix (see the above example). One can use the 
       keywords `noMCMC` and `noweight`, instead of `MCMC` and `weight`.  
~~~~~~~~~~~~~~~~
    CorrelatedGaussianObservables <name> Nobs
    Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight    <central value> <Gaussian error> <flat error>
    Observable <name> <obs label> <histolabel> <min> <max> (no)MCMC (no)weight    <central value> <Gaussian error> <flat error>
    ...
    ...
    <Total of Nobs lines of Observables>
    <Nobs x Nobs correlation matrix>
~~~~~~~~~~~~~~~~  
   Any construction for `Observable` mentioned in item 5 of this list above can be used in a `CorrelatedGaussianObservables` 
   set. Also, `BinnedObservables` or `FunctionObservables` can be used instead of and alongside `Observable`. If 
   `noweight` is specified for any `Observable` then that particular `Observable` along with the corresponding row and 
   column of the correlation matrix is excluded from the set of `CorrelatedGaussianObservables`.
   
   In addition, the inverse covariance matrix of a set of `Nobs` `Observables` can be specified with the following:  
~~~~~~~~~~~~~~~~
    ObservablesWithCovarianceInverse <name> Nobs
    Observable <name> <obs label> <histolabel> <min> <max> MCMC weight    <central value> 0. 0.
    Observable <name> <obs label> <histolabel> <min> <max> MCMC weight    <central value> 0. 0.
    ...
    ...
    <Total of Nobs lines of Observables>
    <Nobs x Nobs inverse covariance matrix>
~~~~~~~~~~~~~~~~  

10. __A correlation between two observables__ can be obtained with: 
~~~~~~~~~~~~~~~~  
    %Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> noMCMC noweight <obs2 label> <histolabel2> <min2> <max2>
    %Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> MCMC file <filename> <histoname> <obs2 label> <histolabel2> <min2> <max2>
    or
    %Observable2D <name> (no)MCMC (no)weight
    (Binned)%Observable <obs label 1> <histolabel 1> <min> <max> <central value> <Gaussian error> <flat error> (<bin_min> <bin_max>)
    (Binned)%Observable <obs label 2> <histolabel 2> <min> <max> <central value> <Gaussian error> <flat error> (<bin_min> <bin_max>)
    %Observable2D <name> MCMC file filename histoname
    (Binned)%Observable <obs label 1> <histolabel 1> <min> <max> (<bin_min> <bin_max>)
    (Binned)%Observable <obs label 2> <histolabel 2> <min> <max> (<bin_min> <bin_max>)
~~~~~~~~~~~~~~~~  
    
9. __Include configuration files__ with the `IncludeFile` directive. This is useful if one
   wants to separate the input configurations for better organization and flexibility.


### Step 2: Monte Carlo configuration file:

The parameters and options of the Monte Carlo run are specified in a configuration file, separate from the one(s) for 
the model. Each line in the file has a pair of a label and its value, separated by space(s) or tab(s). The available
parameters and options are: 

* `NChains`: The number of chains in the Monte Carlo run. A minimum of 5 is suggested (default). If the theory space is 
    complicated and/or the number of parameters is large then more chains are necessary. The amount of statistics 
    collected in the main run is proportional to the number of chains.

* `PrerunMaxIter`: The maximum number of iterations that the prerun will go through (Default: 1000000). The prerun ends 
    automatically when the chains converge (by default R$<$1.1, see below) and all efficiencies are adjusted. While it 
    is not necessary for the prerun to converge for a run to be completed, one should exercise caution if convergence 
    is not attained.

* `NIterationsUpdateMax`: The maximum number of iterations after which the proposal functions are updated in the
    pre-run and convergence is checked. (Default: 1000)
    
* `Seed`: The seed can be fixed for deterministic runs. (Default: 0, corresponding to a random seed initialization)

* `Iterations`: The number of iterations in the main run. This run is for the purpose of collecting statistics and is 
    at the users discretion. (Default: 100000)
    
* `MinimumEfficiency`: This allows setting the minimum efficiency of all the parameters to be attained in the prerun. 
    (Default: 0.15)

* `MaximumEfficiency`: This allows setting the maximum efficiency of all the parameters to be attained in the prerun. 
    (Default: 0.35)

* `RValueForConvergence`: The $R$-value for which convergence is considered to be attained in the prerun can be set 
    with this flag. (Default: 1.1)

* `WriteParametersChains`: The chains will be written in the `ROOT` file MCout*.root. This can be used for analyzing 
    the performance of the chains and/or to use the sampled pdf for postprocessing. (Default: false)

* `FindModeWithMinuit`: To find the global mode with `MINUIT` starting from the best fit parameters in the MCMC run. 
    (Default: false)
    
* `RunMinuitOnly`: To run a `MINUIT` minimization only, without running the MCMC. (Default: false)

* `CalculateNormalization`: Whether the normalization of the posterior pdf will be calculated at the end of the Monte 
    Carlo run. This is useful for model comparison. (Default: false)

* `NIterationNormalizationMC`: The maximum number of iterations used to compute the normalization. (Default: 1000000)

* `PrintAllMarginalized`: Whether all marginalized distributions will be printed in a pdf file (MonteCarlo\_plots\_*.pdf).
    (Default: true)

* `PrintCorrelationMatrix`: Whether the parametric correlation will be printed in ParamCorrelations*.pdf and ParamCorrelations*.tex. 
    (Default: false)

* `PrintKnowledgeUpdatePlots`: Whether comparison between prior and posterior knowledge will be printed in a plot stored
    in ParamUpdate*.pdf. (Default: false)

* `PrintParameterPlot`: Whether a summary of the parameters will be printed in ParamSummary*.pdf. (Default: false)

* `PrintTrianglePlot`: Whether a triangle plot of the parameters will be printed. (Default: false)

* `WritePreRunData`: Whether the prerun data is written to a file. Useful to exploit a successful prerun for multiple runs. 
    (Default: false)

* `ReadPreRunData`: Whether the prerun data will be read from a previously stored prerun file. (Name of the file, default: empty)

* `MultivariateProposal`: Whether the proposal function will be multivariate or uncorrelated. (Default: true)

* `Histogram1DSmooth`: Sets the number of iterative smoothing of 1D histograms. (Default: 0)

* `Histogram2DType`: Sets the type of 2D histograms: 1001 -> Lego (default), 101 -> Filled, 1 -> Contour.

* `MCMCInitialPosition`: The initial distribution of chains over the parameter space. (Options: Center, RandomPrior, 
    RandomUniform (default))

* `PrintLogo`: Toggle the printing of the `HEPfit` logo on the histograms. (Default: true)

* `NoHistogramLegend`: Toggle the printing of the histogram legend. (Default: false)

* `PrintLoglikelihoodPlots`: Whether to print the 2D histograms for the parameters vs. loglikelihood. (Default: false)

* `WriteLogLikelihoodChain`: Whether to write the value of loglikelihood in a chain. (Default: false)

* `Histogram2DAlpha`: Control the transparency of the 2D histograms. This does not work with all 2D histogram types. 
    (Default: 1)

* `NBinsHistogram1D`: The number of bins in the 1D histograms. (Default: 100, 0 sets default)

* `NBinsHistogram2D`: The number of bins in the 2D histograms. (Default: 100, 0 sets default)

* `InitialPositionAttemptLimit`: The maximum number of attempts made at getting a valid logarithm of the likelihood for 
    all chains before the pre-run starts. (Default: 10000, 0 sets default)

* `SignificantDigits` The number of significant digits appearing in the Statistics file. (Default: computed based on 
    individual results, 0 sets default)

{\bf \texttt{HistogramBufferSize}}: The memory allocated to each
histogram. Also determines the number of events collected before
setting automatically the histogram range. (Default: 100000)\\

For example, a Monte Carlo configuration file is written as: 

~~~~~~~~~~~~~~~~
NChains                    10    
PrerunMaxIter              50000 
Iterations                 10000 
#Seed                       1
PrintAllMarginalized       true  
PrintCorrelationMatrix     true      
PrintKnowledgeUpdatePlots  false 
PrintParameterPlot         false 
MultivariateProposal       true  
~~~~~~~~~~~~~~~~

where a '#' can be placed at the beginning of each line to comment it
out. These are the most commonly used settings. Some less commonly used settings:

~~~~~~~~~~~~~~~~
FindModeWithMinuit      (Default: false)
MinimumEfficiency       (Default: 0.15, set to 0. - 1.)
WriteParametersChaina   (Default: false)
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

* `log.txt`: This is the log file.
* `MCout.root`: This is the root file containing all the information of the run and
  the histogram objects.
* `MonteCarlo_results.txt`: This file contains the fits to the parameters.
* `MonteCarlo_plots.pdf`: This File contains the histograms for the parameters.
* `Observables`: This directory will contain the histograms of all the observables
  specified in the config file.
* `Observables/HistoLog.txt`: This file contains the information on over-run and under-run
  of the histogram filling.
* `Observables/Statistics.txt`: This file contains the compilation of the statistics
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

Even Generation Mode
------------------

Using the model configuration file used in the Monte Carlo mode, one can obtain predictions of observables. 
An example can be found in `examples/EventGeneration` folder:

~~~~~~~~~~~~~~~
  $ cd examples/EventGeneration
  $ make
~~~~~~~~~~~~~~~

After making the configuration files, run with the command:
~~~~~~~~~~~~~~~
  $ ./analysis <model conf> <number of iterations> [output folder]
~~~~~~~~~~~~~~~

The `<number of iterations>` defines the number of random points in the parameter space that will be evaluated. Setting 
this to 0 gives the value of the observables at the central value of all the parameters. If the `[output folder]` is not
 specified everything is printed on the screen and no data is saved. Alternately, one can specify the output folder and 
 the run will be saved if `<number of iterations>` > 0. The output folder can be found in `./GeneratedEvents`. The 
 structure of the output folder is as follows.


###Output folder structure:
* `CGO: Contains any correlated Gaussian observables that might have been listed in the model configuration files.
* `Observables`: Contains any observables that might have been listed in the model configuration files.
* `Parameters`: Contains all the parameters that were varied in the model configuration files.
* `Summary.txt`: Contains a list of the model used, the parameters varied, the observables computed and the number of 
   events generated. This can be used, for example, to access all the files from a third party program.
   
The parameters and the observables are stored in the respective directories in files that are named after the same. 
For example, the parameter `lambda` will be saved in the file `lambda.txt` in the `Parameters` folder.


Library mode without MCMC 
------------------------- 

The library mode allows for access to all the observables implemented in `HEPfit`
without a Monte Carlo run. The users can use one of our defined @ref PageModels and vary ModelParameters
according to their own algorithm and get the corresponding predictions for the observables. 

This is made possible through:

* a combined library: libHEPfit.a (installed in `HEPFIT_INSTALL_DIR/lib/`
* a combined header file: HEPfit.h (installed in `HEPFIT_INSTALL_DIR/include/HEPfit/`)


The `HEPfit` library allows for two different implementations of the access algorithm.

### Non-Minimal Mode:

In the non-minimal mode the user can use the model configuration file to pass the default value of
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

A hepfit-config script can be found in the `HEPFIT_INSTALL_DIR/bin/` directory, which can be invoked with the 
following options:

~~~~~~~~~~~~~~~
Library and Library Path: hepfit-config --libs
Include Path: hepfit-config --cflags
~~~~~~~~~~~~~~~

Custom models and observables 
------------------------- 
A very useful feature of `HEPfit` is that it allows the user to create custom models and observables. We have already 
provided a template that can be found in the `examples/myModel` directory which can be used as a starting point. 
Below we describe how to implement both custom models and custom observables.

### Custom Models:
The idea of a custom model is to define an additional set of parameters over and above what is defined in any model in
`HEPfi`. Typically the starting point is the `StandardModel`, as in the template present in the `HEPfit` package. 
Going by this template in the `examples/myModel` directory, to create a model one has to define the following:

* In the `myModel.h` header file:
   1. Define the number of additional parameters:
~~~~~~~~~~~~~~~
static const int NmyModelvars = 4;
~~~~~~~~~~~~~~~
   2. Define the variables corresponding to the parameters:
~~~~~~~~~~~~~~~
double c1, c2, c3, c4;
~~~~~~~~~~~~~~~
   3. Define getters for all the parameters:
~~~~~~~~~~~~~~~
double getc1() const { return c1; }
double getc2() const { return c2; }
double getc3() const { return c3; }
double getc4() const { return c4; }
~~~~~~~~~~~~~~~        
    
* In the `myModel.cpp` file:
    1. Define the names of the parameters (they can be different from the variable names):
~~~~~~~~~~~~~~~
const std::string myModel::myModelvars[NmyModelvars] = {"c1", "c2", "c3", "c4"};
~~~~~~~~~~~~~~~
    2. Link the parameter name to the variable containing it for all the parameters:  
~~~~~~~~~~~~~~~
ModelParamMap.insert(std::make_pair("c1", std::cref(c1)));
ModelParamMap.insert(std::make_pair("c2", std::cref(c2)));
ModelParamMap.insert(std::make_pair("c3", std::cref(c3)));
ModelParamMap.insert(std::make_pair("c4", std::cref(c4)));
~~~~~~~~~~~~~~~
    3. Link the names of the parameters to the corresponding variables in the \texttt{setParameter} method:  
~~~~~~~~~~~~~~~
        if(name.compare("c1") == 0)
                c1 = value;
            else if(name.compare("c2") == 0)
                c2 = value;
            else if(name.compare("c3") == 0)
                c3 = value;
            else if(name.compare("c4") == 0)
                c4 = value;
            else
                StandardModel::setParameter(name,value);
~~~~~~~~~~~~~~~        

This completes the definition of the model. One can also define flags that will control certain aspects of the model, 
but since this is an advanced and not so commonly used feature we will not describe it here. There is an implementation 
in the template for the user to follow should it be needed. Finally, the custom model needs to be added with a name to 
the `ModelFactory` in the main function as is done in `examples/myModel/myModel_MCMC.cpp`.

~~~~~~~~~~~~~~~
ModelF.addModelToFactory("myModel", boost::factory<myModel*>() );
~~~~~~~~~~~~~~~

### Custom Observables

The definition of custom observables does not depend on having defined a custom model or not. A custom observable can be
 any observable that has not been defined in `HEPfit`. It can be a function of parameters already defined in a `HEPfit` 
 model or in a custom model or a combination of the two. However, a custom observable has to be explicitly added to the 
 `ThObsFactory` in the main function as is done in `examples/myModel/myModel_MCMC.cpp`.
~~~~~~~~~~~~~~~
ThObsF.addObsToFactory("BIN1", boost::bind(boost::factory<yield*>(), _1, 1) );
ThObsF.addObsToFactory("BIN2", boost::bind(boost::factory<yield*>(), _1, 2) );
ThObsF.addObsToFactory("BIN3", boost::bind(boost::factory<yield*>(), _1, 3) );
ThObsF.addObsToFactory("BIN4", boost::bind(boost::factory<yield*>(), _1, 4) );
ThObsF.addObsToFactory("BIN5", boost::bind(boost::factory<yield*>(), _1, 5) );
ThObsF.addObsToFactory("BIN6", boost::bind(boost::factory<yield*>(), _1, 6) );
ThObsF.addObsToFactory("C_3", boost::factory<C_3*>() );
ThObsF.addObsToFactory("C_4", boost::factory<C_4*>() );
~~~~~~~~~~~~~~~

The first 6 observables require an argument and hence needed `boost::bind`. The last two do not need an argument. The
implementation of these observables can be found in `examples/myModel/src/myObservables.cpp` and the corresponding
header file. In this template the `myObservables` class inherits from the `THObservable` class and the observables
called `yield`, `C_3` and `C_4` inherit from the former. Passing an object of the `StandardModel` class as a
reference is mandatory as is the overloading of the `computeThValue` method by the custom observables, which is used to 
compute the value of the observable at each iteration.
