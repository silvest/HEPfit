/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <ComputeObservables.h>
#include <InputParameters.h>

int main(int argc, char** argv) 
{
    try {
        std::map<std::string, double> DObs; /* Map of observables */
        std::map<std::string, double> DPars_IN; /* Map of initialization parameters */
        std::map<std::string, double> DPars; /* Map of parameters being varied. */
        std::map<std::string, std::string> DFlags; /* Map of flags */
        
        /* Initialize  the parameters of the model */
        std::string ModelName = "StandardModel";
        InputParameters IP;
        DPars_IN = IP.getInputParameters(ModelName);
        
        /* Initialize the Observables to be returned */
        DObs["Mw"] = 0;
        //DObs["GammaW"] = 0.;
        //DObs["GammaZ"] = 0.;
        DObs["Mz"] = 0.;
        ComputeObservables CO(ModelName, DPars_IN, DObs);
        CO.setFlags(DFlags);
        
        for (int i = 0; i < 2; i++) {
        /* Pass values to the parameters being varied */
            DPars["Mz"] = 91.1875 + 0.0001 * i;
            DPars["AlsMz"] = 0.1184 + 0.000001 * i;
            
            DObs = CO.compute(DPars);
        
            std::cout << "\nParameters[" << i + 1 << "]:"<< std::endl;
            for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++){
                std::cout << it->first << " = " << it->second << std::endl;
            }
            std::cout << "\nObservables:[" << i + 1 << "]:" << std::endl;
            for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++){
                std::cout << it->first << " = " << it->second << std::endl;
            }
        }
        
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}
