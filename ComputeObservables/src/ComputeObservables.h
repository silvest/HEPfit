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
 * @addtogroup ComputeObservables
 * @brief A module for accessing the observables without a MCMC run.
 * @details This module is for using the implementations of the observables without
 * running a Markov Chain Monte Carlo. It allows for accessing
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
     * It is to be used to compute observables useing of a SomeModel.conf file to initialize
     * the mandatory parameters.
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * the model parameters and observables to be calculated
     */
    ComputeObservables(const std::string& ModelConf_i);
    
    /**
     * @brief Constructor.
     * @details This constructor passes the  model name, model parameters and model flags.
     * It is to be used to compute observables without the use of a SomeModel.conf file.
     * @param[in] ModelName_i the name of the model being used
     * @param[in] DPars_i the mandatory parameters of the model being used
     * @param[in] DObs_i the map of observables to be computed
     */
    ComputeObservables(const std::string& ModelName_i, std::map<std::string, double> DPars_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ComputeObservables();
    
    /**
     * @brief This method sets the necessary flag for the requested model
     * when SomeModel.conf is not used to pass the input values (c.f.
     * ComputeObservables(ModelName_i, DPars_i, DObs_i)).
     * @param[in] DFlags_i the flags for the model being used
     */
    void setFlags(std::map<std::string, std::string> DFlags_i);
    
    /**
     * @brief The method used to compute observables using an object built with either
     * ComputeObservables(ModelName_i, DPars_i) or ComputeObservables(ModelName_i, DPars_i, DObs_i).
     * @param[in] DP the map of parameters being varied
     */
    std::map<std::string, double> compute(std::map<std::string, double> DP);
    
    /**
     * @brief A method to add observables to the list of observables.
     * @param[in] ObsName the name of the observable to be added
     */
    void RemoveObservable(std::string ObsName);
    
    /**
     * @brief A method to remove observables to the list of observables.
     * @param[in] ObsName the name of the observable to be removed
     */
    void AddObservable(std::string ObsName);
    
    /**
     * @brief A method to get the map of observables.
     */
    std::map<std::string, double> getObservables()
    {
        return(DObs);
    };
    
    
private:
    
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
};

/**
 * @}
 */

#endif	/* COMPUTEOBSERVABLES_H */



