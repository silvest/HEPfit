/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef COMPUTEOBSERVABLES_H
#define	COMPUTEOBSERVABLES_H

#include <InputParser.h>
#include <TF1.h>
#include <Observable.h>
#include <Observable2D.h>
#include <CorrelatedGaussianObservables.h>
#include <ModelParaVsObs.h>
#include <ModelParameter.h>
#include <Model.h>


/**
 * @class ComputeObservables
 * @ingroup EventGeneration
 * @brief A class for generating events.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class ComputeObservables {
public:
    /**
     * @brief The default constructor.
     * @details The default constructor passes the name of the SomeModel.conf file
     * @param[in] ModelConf_i the name of the input configuration file for the model name,
     * the model parameters and observables to be calculated.
     * @param[in] DObs the map of observables to be computed
     */
    ComputeObservables(const std::string& ModelConf_i,
                       std::map<std::string, double> DObs);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ComputeObservables();
    
    /**
     * @brief The method used to compute observables
     * @details 
     * @param[in] DPars the map of parameters being varied
     */
std::map<std::string, double> compute(std::map<std::string, double> DPars);
    

    
private:
    
    InputParser myInputParser; /**< An oject of the InputParser() class. */
    std::map<std::string, double> DPars; /**< Map of parameters to be passed to Model(). */
    std::map<std::string, double> DObs; /**< Map of parameters to be passed to Model(). */
    Model* Mod; /**< Name of the model as defined in SomeModel.conf*/
    std::vector<ModelParameter> ModPars; /**< Vector for the model parameters defined in SomeModel.conf. */
    std::vector<ModelParameter> ModParsVar; /**< Vector for the model parameters varied in SomeModel.conf. */
    std::vector<Observable> Obs; /**< Vector for the observables defined in SomeModel.conf. */
    std::vector<Observable2D> Obs2D; /**< Vector for the Observables2D defined in SomeModel.conf. */
    std::vector<CorrelatedGaussianObservables> CGO; /**< vector for the Correlated Gaussian Observables defined in SomeModel.conf. */
    std::vector<ModelParaVsObs> ParaObs; /**< Vector for the ModelParaVsObs defined in SomeModel.conf. */
    std::string ModelConf; /**< String for the name of the SomeModel.conf file. */
    std::vector<std::string> paraNames;/**< The vector of allowed parameter names. */
    bool checkPara; /**< The boolean cheack for consistency in parameter names.*/
};


#endif	/* COMPUTEOBSERVABLES_H */



