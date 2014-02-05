/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
 * @brief A module for Markov Chain Monte Carlo based on
 * the <a href="https://www.mppmu.mpg.de/bat/" target=blank>Bayesian Analysis Tool (BAT)</a>. 
 * @{
 */
/**
 * @class MonteCarlo
 * @brief A class for Monte Carlo. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is responsible for the MCMC runs using the
 * MCMCEngine class together with input parameters from the InputParser
 * class. The Monte Carlo mode is the default run mode if the --noMC flag is
 * not specified. The settings for the Monte Carlo run can be fixed in in the
 * format specified in the MonteCarlo.conf file, the name of this configuration
 * file should be specified as the second argument to the executable.
 */
class MonteCarlo {
public:
    /**
     * @brief Constructor.
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
     * \li Initiate InputParser::ReadParameters() which read the SomeModel.conf file for setting the
     * model parameters and the observables to be generated. This call also passes ont he name of the
     * model to the private member ModelName.
     *
     * \li Map the model parameter mean values with the map DP and calculate buffsize which can be used
     * for implementing MPI runs. The variable buffsize is incremented only for those model parameters
     * that are varied in the Monte Carlo run.
     *
     * \li Initiate the model using the map DP checking for consistency between required model parameters
     * and the ones supplied in the SomeModel.conf file
     *
     * \li Set the model name and initialize the MCMC engine with the model at hand.
     *
     * \li For MPI runs, send the evaluates of the loglikelihood for the prerun to the slave processes 
     * (rank != 0)
     *
     * \li The master in an MPI run or the process in a serial run then parses the MonterCarlo.conf file
     * to read the parameters for the Monte Carlo run.
     * 
     * \li The root ouput file is handed over to the object out of type BCModelOutput and some
     * <a href="https://www.mppmu.mpg.de/bat/" target=blank>BAT</a> options are set for the log and 
     * output histograms.
     *
     * \li The MCMC is run and all the out put is generated according to the options set in the 
     * MonteCarlo.conf file.
     *
     * \li For a MPI run, the final call to the MPI class prepares the processes for finalizing.
     *
     * The details for BCSummaryTool, BCLog, BCAux and BCModelOutput can be found in the
     * <a href="https://www.mppmu.mpg.de/bat/" target=blank>BAT website</a>. These are used mainly
     * to generate logs and output.
     *
     * The detials of the object MCEngine of type MonteCarloEngine which overloads the BCEngineMCMC
     * class can be found in our documentation of the former class.
     *
     * @param[in] rank = MPI::COMM_WORLD.Get_rank(), specifies the rank of the process. This
     * carries a non zero value only when the executable is compiled with the parallelalized version
     * of <a href="https://www.mppmu.mpg.de/bat/" target=blank>BAT</a> and run as parallel processes with mpi.
     */
    void Run(const int rank);
private:
    InputParser myInputParser; ///< An oject of the InputParser class.
    MonteCarloEngine MCEngine; ///< An object of the MonteCarloEngine class.
    std::vector<ModelParameter> ModPars; ///< Vector for the model parameters defined in SomeModel.conf.
    std::vector<Observable> Obs; ///< Vector for the observables defined in SomeModel.conf.
    std::vector<Observable2D> Obs2D; ///< Vector for the Observables2D defined in SomeModel.conf.
    std::vector<CorrelatedGaussianObservables> CGO; ///< Vector for the Correlated Gaussian Observables defined in SomeModel.conf.
    std::vector<ModelParaVsObs> ParaObs; ///< Vector for the ModelParaVsObs defined in SomeModel.conf.
    std::string ModelConf; ///< String for the name of the SomeModel.conf file.
    std::string MCMCConf; ///< String for the name of the MonteCarlo.conf file.
    std::string OutFile; ///< String for the name of the output root file without the .root extension.
    std::string JobTag; ///< String for the optional JobTag argument to be passes to the executable.
    std::string ObsDirName; ///< String for the output directory name.
    bool noMC; ///< flag to Specify the non Monte Carlo runs passed as an optional argument to the executable.
    bool FindModeWithMinuit; ///< Flag for using Minuit libraries.
    bool PrintAllMarginalized; ///< Flag for printing all Marginalized distributions to be passed on to the <a href="https://www.mppmu.mpg.de/bat/" target=blank>BAT</a> routines.
    bool PrintCorrelationMatrix; ///< Flag for printing the correlation matrix.
    bool PrintKnowledgeUpdatePlots; ///< Flag for printing plots to compare prior vs. posterior knowledge of parameters.
    bool PrintParameterPlot; ///< Flag for printing the overview parameter plots.
};

/** 
 * @}
 */

#endif	/* MONTECARLO_H */
