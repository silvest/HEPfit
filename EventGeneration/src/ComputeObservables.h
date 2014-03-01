/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef COMPUTEOBSERVABLES_H
#define	COMPUTEOBSERVABLES_H

#include <InputParser.h>
#include <ThFactory.h>
#include <Observable.h>
#include <Observable2D.h>
#include <StandardModel.h>
#include <NPSTU.h>
#include <NPSTUVWXY.h>
#include <NPEpsilons.h>
#include <NPEpsilons_pureNP.h>
#include <NPHiggs.h>
#include <NPZbbbar.h>
#include <NPEffective1.h>
#include <NPEffective2.h>
#include <GeneralSUSY.h>
#include <pMSSM.h>
#include <SUSYMassInsertion.h>
#include <MFV.h>
#include <SUSY.h>
#include <THDM.h>
#include <CorrelatedGaussianObservables.h>
#include <ModelParaVsObs.h>
#include <ModelParameter.h>
#include <Model.h>

/**
 * @addtogroup EventGeneration
 * @brief A module for accessing the observables without a MCMC run.
 * @details This module is for using the implementations of the observables without
 * running a Markov Chain Monte Carlo. It contains code that allows for generations
 * of events using the random number generator from ROOT. It also allows for acessing
 * the observables in a library mode where the user can specify the parameters and 
 * compute the observables.
 * @{
 */

/**
 * @class ComputeObservables
 * @brief A class for providing access to the computation of observables without a Monte Carlo run.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class can be used to create an object that takes a map of model parameters that need to be varied 
 * in the analysis and passes out a map of the observables that it is told to compute.
 */
class ComputeObservables {
public:
    /**
     * @brief Constructor.
     * @details This constructor passes the name of the SomeModel.conf file.
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * the model parameters and observables to be calculated.
     * @param[in] DObs the map of observables to be computed
     */
    ComputeObservables(const std::string& ModelConf_i,
                       std::map<std::string, double> DObs_i);
    /**
     * @brief Constructor.
     * @details This constructor passes the model parameters, model name and model flags.
     * @param[in] ModelName_i the name of the model being used
     * @param[in] DParas_i the mandatory parameters of the model being used
     * @param[in] DObs_i the map of observables to be computed
     */

    ComputeObservables(const std::string& ModelName_i,
                       std::map<std::string, double> DPars_i,
                       std::map<std::string, double> DObs_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ComputeObservables();
    
    /**
     * @brief This method sets the necessary flag for the requested mode.
     * @param[in] DFlags_i the flags for the model being used
     */    
void setFlags(std::map<std::string, std::string> DFlags_i);
    
    /**
     * @brief The method used to compute observables
     * @param[in] DP the map of parameters being varied
     */
std::map<std::string, double> compute(std::map<std::string, double> DP);
    

    
private:
    
    StandardModel* ModelDictionary();
    
    ThFactory* thf;///< Pointer to an object of type ThFactory.
    std::string ModelName;///< Name of the Model to be used.
    InputParser myInputParser; ///< An oject of the InputParser class.
    std::map<std::string, double> DPars; ///< Map of parameters to be passed to Model.
    std::map<std::string, double> DObs; ///< Map of parameters to be passed to Model.
    std::map<std::string, std::string> DFlags; ///< Map of model flags to be passed to Model.
    StandardModel* Mod; ///< Name of the model as defined in SomeModel.conf
    std::vector<ModelParameter> ModPars; ///< Vector for the model parameters defined in SomeModel.conf.
    std::vector<ModelParameter> ModParsVar; ///< Vector for the model parameters varied in SomeModel.conf.
    std::vector<Observable> Obs; ///< Vector for the observables defined in SomeModel.conf. */
    std::vector<Observable2D> Obs2D; ///< Vector for the Observables2D defined in SomeModel.conf.
    std::vector<CorrelatedGaussianObservables> CGO; ///< vector for the Correlated Gaussian Observables defined in SomeModel.conf.
    std::vector<ModelParaVsObs> ParaObs; ///< Vector for the ModelParaVsObs defined in SomeModel.conf.
    std::string ModelConf; ///< String for the name of the SomeModel.conf file.
    std::vector<std::string> paraNames;///< The vector of allowed parameter names.
    bool checkPara; ///< The boolean cheack for consistency in parameter names.
};

/**
 * @}
 */

#endif	/* COMPUTEOBSERVABLES_H */



