/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
#include <StandardModelMatching.h>
#include <NPSTU.h>
#include <NPSTUVWXY.h>
#include <NPEpsilons.h>
#include <NPEpsilons_pureNP.h>
#include <NPHiggsST.h>
#include <NPZbbbar.h>
#include <NPEffective1.h>
#include <NPEffective2.h>
#include <GeneralSUSY.h>
#include <pMSSM.h>
#include <SUSYMassInsertion.h>
#include <SUSYMassInsertionMatching.h>
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
 * @brief A project for interfaces to model parameters and observables. 
 * @{
 */

/**
 * @class InputParser
 * @brief A class for reading input parameters. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class InputParser {
public:
    
    /**
     * @brief The default constructor.
     */
    InputParser();
    
    /**
     * @brief The copy constructor.
     * @param[in] orig
     */
    InputParser(const InputParser& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~InputParser();

    Observable ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg);
    
    /**
     * @brief responsible for parsing the SomeModel.conf file
     * @details 
     * @param[in] filename name of the SomeModel.conf file passed as an argument to the executable
     * @param[out] ModelPars the vector of model parameters
     * @param[out] Observables the vector of observables
     * @param[out] Observables2D the vector of observable pairs
     * @param[out] CGO the vector of correlated gaussian observables
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
     * @brief a get method to access the pointer to the object of the StandardModel() class
     * @return myModel apointer to the object of the StandardModel() class
     */
    StandardModel* getMyModel() const
    {
        return myModel;
    }

    /**
     * @brief a get method to access the pointer to the object of the StandardModelMatching() class
     * @return myModelMatching apointer to the object of the StandardModelMatching() class
     */
    StandardModelMatching* getMyModelMatching() const
    {
        return myModelMatching;
    }

private:
    StandardModel* myModel;/**< an object of the StandardModel() class*/
    StandardModelMatching* myModelMatching;/**< an object of the StandardModelMatching() class*/
    ThFactory* thf;/**< an object of the ThFactory() class*/
    std::string modname;/**< a string to store the model name in*/
};

/** 
 * @}
 */

#endif	/* INPUTPARSER_H */
