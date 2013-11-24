/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MONTECARLO_H
#define	MONTECARLO_H

#include <InputParser.h>
#include "MonteCarloEngine.h"

/**
 * @addtogroup MonteCarlo
 * @brief A project for Markov Chain Monte Carlo based on
 * the <a href="https://www.mppmu.mpg.de/bat/" target=blank>Bayesian Analysis Tool (BAT)</a>
 * @{
 */
/**
 * @class MonteCarlo
 * @brief A class for Monte Carlo. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is responsible for the MCMC runs using the
 * MCMCEngine class together with input parameters from the InputParser
 * class. It allows or two modes of runs.
 * @li <b> Event Generator Mode </b>: If --noMC is specified as a run-time flag in the
 * argument, the run will be conducted as a randomly generated events as no Markov
 * Chain runs will be initiated. The observables that are specified in the
 * SomeModel.conf file will be claculated with the central values of the parameters
 * specified in the SomeModel.conf file as the first set of valuses. After which
 * observables will be generated with random sets of parameters. This can be run with
 * or without a seed. This mode is primarily to be used for debugging purposes
 * or for test runs.
 * @li <b> Monte Carlo Mode </b>: This is the default run mode if the --noMC flag is
 * not specified.
 */
class MonteCarlo {
public:
    /**
     * @brief The default constructor.
     * @details The default constructor sets the names of the configuration files
     * and the names of the output root file and specific output directory if a job ID is
     * specified as an argument to the executable. The constructor also sets the default 
     * values of some flags related to the Monte Carlo run. The boolean flags (FindModeWithMinuit,
     * PrintAllMarginalized, PrintCorrelationMatrix, PrintKnowledgeUpdatePlots, PrintParameterPlot)
     * attain their run time values after the MonteCarlo.conf file is parsed by 
     * MonteCarlo::Run().
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * the model parameters and observables to be calculated (Observables, Observables2D,
     * Model Parameters vs. Observables and Correlated Gaussian Observables)
     * @param[in] MonteCarloConf_i the name of the Monte Carlo configuration file that
     * specifies the parameters of the monte carlo run like no. of chains, no. of pre run 
     * iterations etc.
     * @param[in] OutFile_i the name of the root output file to be given without the .root
     * extention
     * @param[in] JobTag_i optional job tag that might be specified
     * @param[in] checkTheoryRange_i (default = false)
     */
    MonteCarlo(const std::string& ModelConf_i,
               const std::string& MonteCarloConf_i,
               const std::string& OutFile_i,
               const std::string& JobTag_i,
               const bool checkTheoryRange_i=false);
    /**
     * @brief The default destructor.
     */
    virtual ~MonteCarlo();
    
    /**
     * @brief This member reponsible for setting the Monte Carlo run parameters and conducting
     * the Monte Carlo run including initiating all output generation
     * @details The algorithm implemented by this member is as follows:
     *
     * \li Initiate InputParser.ReadParameters() which read the SomeModel.conf file for setting the
     * model parameters and the observables to be generated. This call also passes ont he name of the
     * model to the private member ModelName.
     *
     * \li Map the model parameter mean values with the map DP and calculate buffsize which can be used
     * for implementing MPI runs. The variable buffsize is incremented only for those model parameters
     * that are varied in the Monte Carlo run.
     *
     * @param[in] rank rank = MPI::COMM_WORLD.Get_rank(), specifies the rank of the process. This 
     * carries a non zero value only when the executable is compiled with the parallelalized version
     * of BAT and run as parallel processes with mpi.
     */
    void Run(const int rank);
    
    /**
     * @brief
     * @param[in]
     * @param[in]
     * @return
     */
    void generateEvent(const int rank, int unsigned nIteration, bool noMC, int i = 0);
private:
    InputParser myInputParser; /**< an oject of the InputParser() class */
    MonteCarloEngine MCEngine; /**< an object of the MonteCarloEngine() class*/
    std::vector<ModelParameter> ModPars; /**< vector for the model parameters defined in SomeModel.conf*/
    std::vector<Observable> Obs; /**< vector for the observables defined in SomeModel.conf*/
    std::vector<Observable2D> Obs2D; /**< vector for the Observables2D defined in SomeModel.conf*/
    std::vector<CorrelatedGaussianObservables> CGO; /**< vector for the Correlated Gaussian Observables defined in SomeModel.conf*/
    std::vector<ModelParaVsObs> ParaObs; /**< vector for the ModelParaVsObs defined in SomeModel.conf*/
    std::string ModelConf; /**< string for the name of the SomeModel.conf file*/
    std::string MCMCConf; /**< string for the name of the MonteCarlo.conf file*/
    std::string OutFile; /**< string for the name of the output root file without the .root extension*/
    std::string JobTag; /**< string for the optional JobTag argument to be passes to the executable*/
    std::string ObsDirName; /**< string for the output directory name*/
    bool noMC; /**< flag to specify the non Monte Carlo runs passed as an optional argument to the executable*/
    bool FindModeWithMinuit; /**< Flag for using Minuit libraries*/
    bool PrintAllMarginalized; /**< Flag for printing all Marginalized distributions to be passed on to the BAT routines*/
    bool PrintCorrelationMatrix; /**< Flag for printing the correlation matrix*/
    bool PrintKnowledgeUpdatePlots; /**< */
    bool PrintParameterPlot; /**< */
    int checkMode; /**< */
};

/** 
 * @}
 */

#endif	/* MONTECARLO_H */
