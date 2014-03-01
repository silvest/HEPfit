/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef INPUTPARSER_H
#define	INPUTPARSER_H

#include "ThFactory.h"
#include <Observable.h>
#include <Observable2D.h>
#include <CorrelatedGaussianObservables.h>
#include <ModelParaVsObs.h>
#include <ThObservable.h>
#include <ModelParameter.h>
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
#include <iostream>
#include <fstream>
#include <istream>
#include <boost/tokenizer.hpp>
#include <string>

/**
 * @addtogroup InputParser 
 * @brief A module for interfaces to model parameters and observables.
 * @details The InputParser class parses all the input given to the executable
 * through the SomeModel.conf file and the output directives in the same. \n
 * The ThFactory defines and initiates all available observables that can
 * be computed.
 * @{
 */

/**
 * @class InputParser
 * @brief A class for reading input parameters and output directives
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The InputParser parses the SomeModel.conf file for the name of the model, ModelParameters, 
 * %Observables, * %Observables2D, %ModelParaVsObs and %CorrelatedGaussianObservables. The observables 
 * (or set of observables) are parsed for their inclusion in the Monte Carlo run.
 */
class InputParser {
public:
    
    /**
     * @brief Constructor.
     */
    InputParser();
    
    /**
     * @brief The copy constructor.
     * @param[in] orig a reference to an object of type InputParser()
     */
    InputParser(const InputParser& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~InputParser();

    /**
     * @brief The member that parses the Observable directives from SomeModel.conf file
     * @param[in] beg an iterator over words in a line separated by a specific separator character
     */
    Observable ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg);
    
     /**
     * @brief The method used to create an object of the model being used.
     */
    StandardModel* ModelFactory(std::string& ModelName);
    
    /**
     * @brief Responsible for parsing the SomeModel.conf file.
     * @details This method parses the SomeModel.conf file for all input instructions. The algorithm
     * is as folows:
     * \li Search for the predefined names of the available models, initialize the necessary model
     * class and build a ThFactory for that model to define the observables.
     * \li Seach and read the ModelParameter list.
     * \li Search and read the Observable list using the ParseObservable() method.
     * \li Search and read the Observable2D list using the ParseObservable() method.
     * \li Search and read the CorrelatedGaussianObservables.
     * \li Search and read the ModelParaVsObs.
     * \li Return the name of the model initialized.
     * @param[in] filename name of the SomeModel.conf file passed as an argument to the executable
     * @param[out] ModelPars the vector of model parameters
     * @param[out] Observables the vector of observables
     * @param[out] Observables2D the vector of observable pairs
     * @param[out] CGO the vector of correlated Gaussian observables
     * @param[out] ParaObs The vector of model parameter - observable pairs
     * @return modname the name of the model initialized
     */
    std::string ReadParameters(const std::string filename,
            std::vector<ModelParameter>& ModelPars,
            std::vector<Observable>& Observables,
            std::vector<Observable2D>& Observables2D,
            std::vector<CorrelatedGaussianObservables>& CGO,
            std::vector<ModelParaVsObs>& ParaObs);

    /**
     * @brief A get method to access the pointer to the object of the StandardModel class.
     * @return myModel apointer to the object of the StandardModel class
     */
    StandardModel* getMyModel() const
    {
        return myModel;
    }
    
    ThFactory* getMyThFactory() const
    {
        return thf;
    }

private:
    StandardModel* myModel;///< Pointer to an object of type StandardModel.
    ThFactory* thf;///< Pointer to an object of type ThFactory.
    std::string modname;///< A string to store the model name in.
    bool modelset;
};

/** 
 * @}
 */

#endif	/* INPUTPARSER_H */
