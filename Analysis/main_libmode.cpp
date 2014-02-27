/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <ComputeObservables.h>
#include <boost/program_options.hpp>
#ifdef _MPI
#include <mpi.h>
#endif

using namespace boost::program_options;
using namespace std;

int main(int argc, char** argv) 
{

#ifdef _MPI
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    MPI::Status status;
#else
    int rank = 0;
#endif

    string ModelConf;
    if (rank == 0)
        cout << "\n *** SusyFit Markov Chain Montecarlo ***\n" << endl;

    try {
        options_description desc("Allowed (positional) options");
        desc.add_options()
                ("modconf", value<string > (), "model config filename (1st)")
                ("help", "help message")
                ;
        positional_options_description pd;
        variables_map vm;
            pd.add("modconf", 1);
        
        try {
            store(command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
            notify(vm);

            if (vm.count("help")) {
                cout << desc << endl;
                return EXIT_SUCCESS;
            }

            if (vm.count("modconf"))
                ModelConf = vm["modconf"].as<string > ();
            else
                throw runtime_error("Please specify mandatory model config filename");
            } catch(error& e) {
                cerr << "ERROR: " << e.what() << std::endl << std::endl
                     << desc << std::endl;
                return EXIT_FAILURE;
        }
        std::map<std::string, double> DObs;
        std::map<std::string, double> DPars;
        
        // Initialize the Observables to be returned
        DObs["Mw"] = 0;
        DObs["GammaW"] = 0.;
        DObs["GammaZ"] = 0.;
        DObs["Mz"] = 0.;
        ComputeObservables CO(ModelConf, DObs);
        
        for (int i = 0; i < 10000; i++) {
        /* Initilize/pass values to the parameters being varied */
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

#ifdef _MPI
        MPI::Finalize();
#endif
        
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
