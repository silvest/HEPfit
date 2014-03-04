Usage   {#PageUsage}
=============================================

The SusyFit installer generates an executable `"analysis"` and a
library `"libSufyFit.a"` accompanied with header files. You can use
the executable to perform a Bayesian statistical analysis with the
Markov Chain Monte Carlo [__Monte Carlo mode__],
or to obtain predictions of observabes for a give point in the
parameter space of the model [__Single event mdoe__].
Alternatively, you can use the library to obtain predictions
of obsrvables for a given point in the parameter space of the model, 
allowing our computational tool to be called from your own program
[__Library mode__].


Monte Carlo mode 
-----------------

The Monte Carlo analysis is performed with the BAT library. You first
have to prepare a text configuration file containing a list of model parameters,
model flags and observables to be analysed, and another confifuration
file for the Monte Carlo run. 

### %Model configuration file:

A configuration file for model parameters, model flags and
obaservables are written as follwos: 

~~~~~~~~~~~~~~~~
StandardModel
# Model parameters:
ModelParameter  mtop        173.2       0.9         0.  
ModelParameter  mHl         125.6       0.3         0.  

<You have to list all the model parameters here>

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

3. Optionally, you can set __a model flag__ in the format: 

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
   You can use the keyword noMCMC and noweight, instead of MCMC and weight. 
   <br><br>

6. __A correlation between two observables__ can be obtained with: 

  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> noMCMC noweight <obs2 label> <histolabel2> <min2> <max2>`<br>
  `%Observable2D <name> <obs1 label> <histolabel1> <min1> <max1> MCMC file <filename> <histoname> <obs2 label> <histolabel2> <min2> <max2>`

7. __A correlation between a model parameter and an observable__ can
   be obtained with:

  `%ModelParaVsObs <name> <par name> <par histolabel> <par min> <par max> <obs label> <obs histolabel> <obs min> <obs max>`


### %Monte Carlo configuration file:

The parameters and options of the Monte Carlo run are specified in
a configuration file, separated from the one for model parameters,
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

where you can place '#' at the beginnig of each line to comment it
out.


### Run:

After making the configuration files, run the command: 

~~~~~~~~~~~~~~~
  $ analysis <model conf> <Monte Carlo conf> <options>
~~~~~~~~~~~~~~~

where the available options are:

* <b>-\-rootfile=\<filename\></b> output root filename without extension (default: MCout)

* <b>-\-job_tag=\<arg\></b> job tag, appended to output files (default: none)

* <b>-\-thRange</b> output the minimun and maximum of theory values of observables to the file `Observables/HistoLog.txt`


MPI support.

Explain output files.


Single event mode
------------------

Using the model configuration file used in the Monte Carlo mode, you
can obtain predictions of observabes for the central values of the
model parameters with the command: 

~~~~~~~~~~~~~~~
  $ analysis <model conf> --noMC
~~~~~~~~~~~~~~~

where the results are printed on the standard output. 


Library mode 
------------- 

Write instructions for the library mode!


