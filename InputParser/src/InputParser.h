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
#include <ThObservable.h>
#include <ModelParameter.h>
#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <NPSTU.h>
#include <NPSTUVWXY.h>
#include <NPEpsilons.h>
#include <NPHiggsST.h>
#include <NPZbbbar.h>
#include <NPEffective.h>
#include <GeneralSUSY.h>
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
    InputParser();
    InputParser(const InputParser& orig);
    virtual ~InputParser();

    Observable ParseObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg);
    
    std::string ReadParameters(const std::string filename,
            std::vector<ModelParameter>& ModelPars,
            std::vector<Observable>& Observables,
            std::vector<Observable2D>& Observables2D,
            std::vector<CorrelatedGaussianObservables>& CGO);

    StandardModel* getMyModel() const 
    {
        return myModel;
    }

    StandardModelMatching* getMyModelMatching() const 
    {
        return myModelMatching;
    }

private:
    StandardModel* myModel;
    StandardModelMatching* myModelMatching;

    ThFactory* thf;
    std::string modname;
};

/** 
 * @}
 */

#endif	/* INPUTPARSER_H */
