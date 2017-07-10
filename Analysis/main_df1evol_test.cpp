/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example libmode_config.cpp
 * This is an example of how to compute observables from the input parameters
 * defined in a model configuration file.
 *
 */

#include <iostream>
#include <ComputeObservables.h>
#include "HeffDF1.h"

int main(void) 
{
        
        /* Define the model configuration file.                        */
        /* Here it is passed as the first argument to the executable.  */
        /* The model configuration file provides the default values of */
        /* the mandatory model parameters.                             */
        std::string ModelConf = "StandardModel.conf";
        
        /* Define a map for the parameters to be varied. */
        std::map<std::string, double> DPars;
        
        /* Create objects of the classes ModelFactory and ThObsFactory */
        ModelFactory ModelF;
        ThObsFactory ThObsF;

        /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
        /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */
        
        /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
        /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/
        
        /* Create an object of the class ComputeObservables. */
        ComputeObservables CO(ModelF, ThObsF, ModelConf);
        StandardModel& mySM = *CO.getModel();
    
        std::cout << "%SUITE_STARTING% Evolutor" << std::endl;
        std::cout << "%SUITE_STARTED%" << std::endl;

        HeffDF1 Heff(6, "CP", mySM);

        std::cout << Heff.getEvol().DF1Evol(5.,90.,NLO) << std::endl;

}