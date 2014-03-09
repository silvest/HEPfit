#include <iostream>
#include <ComputeObservables.h>

int main(int argc, char** argv) 
{
    try {
        
        /* Define the model configuration file.                        */
        /* Here it is passed as the first argument to the executable.  */
        /* The model configuration file provides the default values of */
        /* the mandatory model parameters, the list of observables and */
        /* corrlated Gaussian observables to be computed.              */
        std::string ModelConf = argv[1];
        
        /* Define a map for the observables and corellated Gaussian observables. */
        std::pair<std::map<std::string, double>, std::map<std::string, std::map<std::string, double> > > DObs_CGO;
        
        /* Define a map for the parameters to be varied. */
        std::map<std::string, double> DPars;
        
        /* Create and object of the class ComputeObservables. */
        ComputeObservables CO(ModelConf);
        
        for (int i = 0; i < 10; i++) {
            
        /* Vary the parameters that need to be varied in the analysis. */
            DPars["Mz"] = 91.1875 + 0.0001 * i;
            DPars["AlsMz"] = 0.1184 + 0.000001 * i;
            
            /* Get the map of observables and corellated Gaussian observables */
            /* with the parameter values defined above.                       */
            DObs_CGO = CO.compute_Obs_CGO(DPars);
        
            std::cout << "\nParameters[" << i + 1 << "]:"<< std::endl;
            for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++){
                std::cout << it->first << " = " << it->second << std::endl;
            }
            std::cout << "\nObservables:[" << i + 1 << "]:" << std::endl;
            for (std::map<std::string, double>::iterator it = DObs_CGO.first.begin(); it != DObs_CGO.first.end(); it++) {
                std::cout << it->first << " = " << it->second << std::endl;
            }
            std::cout << "\nCorellated Gaussian Observables:[" << i + 1 << "]:" << std::endl;
            for (std::map<std::string, std::map<std::string, double> >::iterator it1 = DObs_CGO.second.begin(); it1 != DObs_CGO.second.end(); it1++){
            std::cout << it1->first << " : " << std::endl;
            for (std::map<std::string, double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                std::cout << it2->first << " = " << it2->second << std::endl;
            }
            std::cout << std::endl;
            }
        }
        return EXIT_SUCCESS;
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}