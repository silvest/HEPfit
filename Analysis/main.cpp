/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <MonteCarlo.h>
#include <boost/program_options.hpp>
#include <mpi.h>

using namespace boost::program_options;
using namespace std;

/**
 * @addtogroup Analysis
 * @brief A project for the main function of the SusyFit project. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 *  @{
 */

int main(int argc, char** argv) 
{

    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();
    MPI::Status status;

    string ModelConf, MCMCConf, FileOut, JobTag;
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
                ("job_tag", value<string > ()->default_value(""),
                "job tag")
                ("noMC", "run a single event")
                ("thRange", "output the min and max of theory values to HistoLog.txt")
                ("help", "help message")
                ;
        positional_options_description pd;
        pd.add("modconf", 1);
        pd.add("mcconf", 1);
        //pd.add("rootfile", 1);
        //pd.add("job_tag", 1);

        variables_map vm;
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
                throw runtime_error("missing mandatory model config filename");

            if (vm.count("mcconf"))
                MCMCConf = vm["mcconf"].as<string > ();
            else
                throw runtime_error("missing mandatory montecarlo config filename");

            FileOut = vm["rootfile"].as<string > ();

            JobTag = vm["job_tag"].as<string > ();

            if (vm.count("noMC"))
                FileOut = "";

            if (vm.count("thRange"))
                checkTheoryRange = true;

        } catch(error& e) {
            cerr << "ERROR: " << e.what() << std::endl << std::endl
                 << desc << std::endl;
            return EXIT_FAILURE;
        }

        MonteCarlo MC(ModelConf, MCMCConf, FileOut, JobTag, checkTheoryRange);
        
        MC.Run(rank);

        MPI::Finalize();

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

/** 
 * @}
 */
