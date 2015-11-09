/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <MonteCarlo.h>
#include <GenerateEvent.h>
#include <InputParser.h>
#include <Observable.h>
#include <boost/program_options.hpp>
//#include <mtbar.h>
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

    string ModelConf, JobTag;
    string FolderOut;
    int nIterations = 0;
    bool noMC;

    try {
        options_description desc("\nAllowed (positional) options");
        desc.add_options()
                ("modconf", value<string > (), "MANDATORY model config filename (1st)")
                ("mcconf", value<string > (), "MANDATORY Monte Carlo config filename (2nd) for MCMC run")
                ("test", "run in test mode to generate a single event with --test (no MCMC run but requires the Monte Carlo config file)")
                ("noMC", "run in generate event mode with --noMC (no MCMC run)")
                ("it", value<int > ()->default_value(0),
                 "no. of iterations in generate event mode, specify with --it=# or --it #")
                ("rootfile", value<string > ()->default_value("MCout"),
                "output root filename (without extension)")
                ("output_folder", value<string > ()->default_value(""),
                "output folder for Generate Event mode to be specified for printing to file, please specify with --output_folder <name>")
                ("job_tag", value<string > ()->default_value(""),
                "job tag, please specify with --job_tag <tag>")
                ("weight", "run in generate event mode with --weight for weights")
                ("help", "help message")
                ;
        string coderun = "\n *** HEPfit Routines ***\n"
                         "\nMonte Carlo mode: analysis Model.conf MonteCarlo.conf [--rootfile <name>] [--job_tag <tag>] [--thRange]"
                         "\nSingle Event mode: analysis Model.conf MonteCarlo.conf --test"
                         "\nGenerate Event mode: analysis Model.conf --noMC [--it #] [--weight] [--output_folder <name> [--job_tag <tag>]]\n";
        positional_options_description pd;
        variables_map vm;
        pd.add("modconf", 1);
        pd.add("mcconf", 1);

        
        try {
            store(command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
            notify(vm);

            if (vm.count("help")) {
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                }
                return EXIT_SUCCESS;
            }       
            
            if (vm.count("modconf"))
                ModelConf = vm["modconf"].as<string > ();
            else
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                    throw runtime_error("\nERROR: Please specify mandatory model config filename.\n");
                }

            if (vm.count("noMC") && vm["output_folder"].as<string > () == "" && !(vm["job_tag"].as<string > () == "")){
                if (rank == 0) cout << "\nWARNING: --job_tag not being used as --output_folder is not specified.\nOutput being sent to terminal only.\n";
            } else {
                JobTag = vm["job_tag"].as<string > ();
            }
            
            if (vm.count("noMC") && vm.count("mcconf"))
            {
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                    throw std::runtime_error("\nERROR: Please specify EITHER --noMC OR Monte Carlo configuration file as the second free argument.\nAll other arguments must be specified with option specifiers (--nI=#, --job_tag <tag>, etc.)\n");
                }
            }
            else if (vm.count("noMC") && !vm.count("mcconf"))
            {
                if (rank == 0)
                    cout << "\n *** HEPfit Event Generation ***\n" << endl;
                noMC = true;
                FolderOut = vm["output_folder"].as<string > ();
                bool weight = false;
                nIterations = vm["it"].as<int > ();
                if (vm.count("weight")) weight = true;
                ThObsFactory ThObsF;
                ModelFactory ModelF;
                GenerateEvent GE(ModelF, ThObsF, ModelConf, FolderOut, JobTag, noMC);
//                GE.addCustomParser("PS", boost::factory<InputParser*>());
//                GE.addCustomObservableType("Poisson", boost::factory<Observable*>());
//                GE.linkParserToObservable("Poisson", "PS");
                GE.generate(nIterations, 1, weight);
            }
            else if (!vm.count("noMC") && vm.count("mcconf"))
            {
                if (rank == 0)
                    cout << "\n *** HEPfit Markov Chain Monte Carlo ***\n" << endl;
                noMC = false;
                string MCMCConf = vm["mcconf"].as<string > ();
                string FileOut = vm["rootfile"].as<string > ();
                
                ThObsFactory ThObsF;
                ModelFactory ModelF;
//                ThObsF.addObsToFactory("mtbar", boost::factory<mtbar*>());
                
                MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag);
                if (vm.count("test")){
//                    MC.addCustomParser("PS", boost::factory<InputParser*>());
//                    MC.addCustomObservableType("Poisson", boost::factory<Observable*>());
//                    MC.linkParserToObservable("Poisson", "PS");
                    MC.TestRun(rank);
                } else {
                    MC.Run(rank);
                }
                
            } else {
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                throw runtime_error("\nERROR:  Please specify mandatory Monte Carlo config filename (required even with --test)\n\tOR specify --noMC for Generate Event Mode (no Monte Carlo run)\n");
                }
            }
            

        } catch(error& e) {
            if (rank == 0)
                cout << coderun << endl;
            cerr << "\nERROR: " << e.what() << std::endl << std::endl
                 << desc << std::endl;
            return EXIT_FAILURE;
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
