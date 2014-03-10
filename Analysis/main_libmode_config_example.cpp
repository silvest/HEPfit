#include <iostream>
#include <ComputeObservables.h>

int main(int argc, char** argv) 
{
    try {
        
        /* Define the model configuration file.                        */
        /* Here it is passed as the first argument to the executable.  */
        /* The model configuration file provides the default values of */
        /* the mandatory model parameters.                             */
        std::string ModelConf = argv[1];
        
        /* Define a map for the parameters to be varied. */
        std::map<std::string, double> DPars;
        
        /* Create and object of the class ComputeObservables. */
        ComputeObservables CO(ModelConf);
        
        /* Add the Observables to be returned */
        CO.AddObservables("Mw");
        CO.AddObservables("GammaZ");
        CO.AddObservables("AFBbottom");
        
        /* Get the map of observables if necessary. */
        std::map<std::string, double> DObs = CO.getObservables();
        
        for (int i = 0; i < 2; i++) {
            
            /* Vary the parameters that need to be varied in the analysis. */
            DPars["Mz"] = 91.1875 + 0.0001 * i;
            DPars["AlsMz"] = 0.1184 + 0.000001 * i;
            
            /* Get the map of observables with the parameter values defined above. */
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