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
    MonteCarlo(const std::string& ModelConf_i,
               const std::string& MonteCarloConf_i,
               const std::string& OutFile_i,
               const std::string& JobTag_i,
               const bool checkTheoryRange_i=false);
    virtual ~MonteCarlo();
    void Run(const int rank);
    void generateEvent(const int rank, int unsigned nIteration, bool noMC, int i = 0);
private:
    InputParser myInputParser;
    MonteCarloEngine MCEngine;
    std::vector<ModelParameter> ModPars;
    std::vector<Observable> Obs;
    std::vector<Observable2D> Obs2D;   
    std::vector<CorrelatedGaussianObservables> CGO;
    std::vector<ModelParaVsObs> ParaObs;
    std::string ModelConf, MCMCConf, OutFile, JobTag, ObsDirName;
    bool noMC;
    bool FindModeWithMinuit, PrintAllMarginalized;
    bool PrintCorrelationMatrix, PrintKnowledgeUpdatePlots, PrintParameterPlot;
    int checkMode;
};

/** 
 * @}
 */

#endif	/* MONTECARLO_H */