/* 
 * File:   main.cpp
 * Author: silvest
 *
 * Created on March 30, 2011, 12:54 PM
 */

#include <iostream>
#include <MonteCarlo.h>
#include <boost/program_options.hpp>

using namespace boost::program_options;
using namespace std;

int main(int argc, char** argv) {

    string ModelConf, MCMCConf, FileOut, JobTag;

    cout << "\n *** SusyFit Markov Chain Montecarlo ***\n" << endl;

    try {
        options_description desc("Allowed (positional) options");
        desc.add_options()
                ("modconf", value<string > (), "model config filename (1st)")
                ("mcconf", value<string > (), "montecarlo config filename (2nd)")
                ("rootfile", value<string > ()->default_value("MCout"),
                "output root filename (without extension) (3rd)")
                ("job_tag", value<string > ()->default_value(""),
                "job tag (4th)")
                ("help", "help message")
                ;
        positional_options_description pd;
        pd.add("modconf", 1);
        pd.add("mcconf", 1);
        pd.add("rootfile", 1);
        pd.add("job_tag", 1);

        variables_map vm;
        store(command_line_parser(argc,
                argv).options(desc).positional(pd).run(), vm);
        notify(vm);

        if (vm.count("help")) {
            cout << desc << endl;
            return EXIT_SUCCESS;
        }

        if (vm.count("modconf"))
            ModelConf = vm["modconf"].as<string > ();
        else throw
            runtime_error("missing mandatory model config filename");

        if (vm.count("mcconf"))
            MCMCConf = vm["mcconf"].as<string > ();
        else throw
            runtime_error("missing mandatory montecarlo config filename");

        FileOut = vm["rootfile"].as<string > ();
        
        JobTag = vm["job_tag"].as<string> ();

        MonteCarlo MC(ModelConf, MCMCConf, FileOut, JobTag);
        MC.Run();

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
