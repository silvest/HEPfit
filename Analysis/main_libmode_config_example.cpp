/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
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

int main(int argc, char** argv) 
{
    try {
        
        if(argc != 2){
            /* Print usage and exit. */
            std::cout << "\nusage: " << argv[0] << " ModelConf.conf\n" << std::endl;
            return EXIT_SUCCESS;
        }
        
        /* Define the model configuration file.                        */
        /* Here it is passed as the first argument to the executable.  */
        /* The model configuration file provides the default values of */
        /* the mandatory model parameters.                             */
        std::string ModelConf = argv[1];
        
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
        
        /* Add the observables to be returned. */
        CO.AddObservable("Mw");
        CO.AddObservable("GammaZ");
        CO.AddObservable("AFBbottom");
        
        /* Remove a previously added observable if necessary. */
        //CO.RemoveObservable("AFBbottom");
        
        /* Set the flags for the model being used, if necessary.                         */
        /* The flags have to correspond to the model specified in the model config file. */
        std::map<std::string, std::string> DFlags;
        // DFlags["FLAG"] = "TRUE";
        CO.setFlags(DFlags);
        
        /* Get the map of observables if necessary. */
        std::map<std::string, double> DObs = CO.getObservables();
        
        for (int i = 0; i < 2; i++) {
            
            /* Vary the parameters that need to be varied in the analysis. */
            DPars["Mz"] = 91.1875 + 0.0001 * i;
            DPars["AlsMz"] = 0.1184 + 0.000001 * i;
            
            /* Get the map of observables with the parameter values defined above. */
            DObs = CO.compute(DPars);
        
            std::cout << "\nParameters[" << i + 1 << "]:"<< std::endl;
            for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) {
                std::cout << it->first << " = " << it->second << std::endl;
            }
            std::cout << "\nObservables[" << i + 1 << "]:" << std::endl;
            for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
                std::cout << it->first << " = " << it->second << std::endl;
            }
        }
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}