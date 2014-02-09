/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <MonteCarlo.h>
#include <GenerateEvent.h>
#include <boost/program_options.hpp>
#ifdef _MPI
#include <mpi.h>
#endif

using namespace boost::program_options;
using namespace std;

/**
 * @addtogroup Analysis
 * @brief A module for the main function of the SusyFit project.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @{
 */

int main(int argc, char** argv) 
{

#ifdef _MPI
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    MPI::Status status;
#else
    int rank = 0;
#endif

    string ModelConf, MCMCConf, FileOut, JobTag;
    string FolderOut;
    int nIterations = 0;
    bool noMC;
    bool checkTheoryRange = false;

    if (rank == 0)
        cout << "\n *** SusyFit Markov Chain Montecarlo ***\n" << endl;

    try {
        options_description desc("Allowed (positional) options");
        desc.add_options()
                ("modconf", value<string > (), "model config filename (1st)")
                ("mcconf", value<string > (), "montecarlo config filename (2nd)")
                ("rootfile", value<string > ()->default_value("MCout"),
                "output root filename (without extension)")
                ("output_folder", value<string > ()->default_value(""),
                "output folder for Generate Event mode to be specified for printing to file")
                ("job_tag", value<string > ()->default_value(""),
                "job tag")
                ("noMC", "run in generate event mode with --noMC")
                ("nI", value<int > ()->default_value(0),
                "no. of iterations in generate event mode, specify with --nI=#")
                ("thRange", "output the min and max of theory values to HistoLog.txt")
                ("help", "help message")
                ;
        positional_options_description pd;
        variables_map vm;
        if (!vm.count("noMC")) {
            pd.add("modconf", 1);
            pd.add("mcconf", 1);
            //pd.add("rootfile", 1);
            //pd.add("job_tag", 1);
        } else {
            pd.add("modconf", 1);
            //pd.add("output_folder", 1);
            //pd.add("job_tag", 1);
        }
        
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
            
            if (vm.count("noMC") && vm.count("mcconf"))
            {
                throw std::runtime_error("\nArgument Error:\nPlease specify EITHER --noMC OR Monte Carlo configuration file as the second argument.\n");
            }
            else if (vm.count("noMC") && !vm.count("mcconf"))
            {
                noMC = true;
                MCMCConf = "";
                FolderOut = vm["output_folder"].as<string > ();
                nIterations = vm["nI"].as<int > ();
            }
            else if (!vm.count("noMC") && vm.count("mcconf"))
            {
                noMC = false;
                MCMCConf = vm["mcconf"].as<string > ();
                FileOut = vm["rootfile"].as<string > ();
                    
            } else {
                throw runtime_error("Please specify mandatory Monte Carlo config filename\nor specify --noMC for Generate Event Mode (no Monte Carlo run)");
            }
            
            JobTag = vm["job_tag"].as<string > ();

            if (vm.count("thRange"))
                checkTheoryRange = true;

        } catch(error& e) {
            cerr << "ERROR: " << e.what() << std::endl << std::endl
                 << desc << std::endl;
            return EXIT_FAILURE;
        }
        
        if (!noMC) {
            MonteCarlo MC(ModelConf, MCMCConf, FileOut, JobTag, checkTheoryRange);
            MC.Run(rank);
        } else {
            GenerateEvent GE(ModelConf, FolderOut, JobTag, noMC);
            GE.generate(rank, nIterations);
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

/** 
 * @}
 */
