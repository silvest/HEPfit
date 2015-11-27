/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef COMPUTEOBSERVABLES_H
#define	COMPUTEOBSERVABLES_H

#include "InputParser.h"
#include "StandardModel.h"
#include "ModelParameter.h"
#include "ModelFactory.h"
#include "ThObsFactory.h"

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
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class can be used to create an object that takes a map of model parameters that need to be varied
 * in the analysis and passes out a map of the observables that it is told to compute.
 */
class ComputeObservables {
public:

    /**
     * @brief Constructor.
     * @details This constructor passes the name of the SomeModel.conf file.
     * It is to be used to compute observables using of a SomeModel.conf file to initialize
     * the mandatory parameters.
     * @param[in] ModelF
     * @param[in] ThObsF
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * @param[in] rank_i the rank of the process in a MPI run (set to 0 for serial run)
     * the model parameters and observables to be calculated
     */
    ComputeObservables(ModelFactory& ModelF, ThObsFactory& ThObsF,
            const std::string& ModelConf_i, const int rank_i = 0);

    /**
     * @brief Constructor.
     * @details This constructor passes the model name and model parameters.
     * It is to be used to compute observables without the use of a SomeModel.conf file.
     * @param[in] ModelF
     * @param[in] ThObsF
     * @param[in] ModelName_i the name of the model being used
     * @param[in] DPars_i the mandatory parameters of the model being used
     * @param[in] rank_i the rank of the process in a MPI run (set to 0 for serial run)
     */
    ComputeObservables(ModelFactory& ModelF, ThObsFactory& ThObsF,
            const std::string& ModelName_i, std::map<std::string, double> DPars_i,
            const int rank_i = 0);

    /**
     * @brief The default destructor.
     */
    virtual ~ComputeObservables();

    /**
     * @brief This method sets the necessary flag for the requested model. 
     * @param[in] DFlags_i the flags for the model being used
     */
    void setFlags(std::map<std::string, std::string> DFlags_i);

    /**
     * @brief The method used to compute observables.
     * @param[in] DP the map of parameters being varied
     */
    std::map<std::string, double> compute(std::map<std::string, double> DP);

    /**
     * @brief A method to add an observable to the list of observables.
     * @param[in] ObsName the name of the observable to be added
     */
    void RemoveObservable(std::string ObsName);

    /**
     * @brief A method to remove an observable from the list of observables.
     * @param[in] ObsName the name of the observable to be removed
     */
    void AddObservable(std::string ObsName);

    /**
     * @brief A method to get the map of observables.
     * @return the map of observables
     */
    std::map<std::string, double> getObservables()
    {
        return (DObs);
    };

    /**
     * @brief A method to get the map of parameters.
     * @return the map of parameters
     */
    std::map<std::string, double> getParameters()
    {
        return (DPars);
    };
    
   void addCustomObservableType(const std::string name, boost::function<Observable*() > funct);

private:

    std::string ModelName; ///< Name of the Model to be used.
    StandardModel* Mod; ///< Pointer to an object of the class StandardModel.
    InputParser myInputParser; ///< An object of the InputParser class.
    std::map<std::string, double> DPars; ///< Map of the parameters to be passed to Model.
    std::map<std::string, double> DObs; ///< Map of the observables to be computed. 
    std::map<std::string, std::string> DFlags; ///< Map of the model flags to be passed to Model.
    std::vector<std::string> paraNames; ///< The vector of allowed parameter names.
    std::map<std::string, ThObservable*> DThObs;
    const int rank; ///<< Rank of the MPI process. Set to 0 for serial run. 

};

/**
 * @}
 */

#endif	/* COMPUTEOBSERVABLES_H */



