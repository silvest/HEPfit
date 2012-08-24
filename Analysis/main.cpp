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

    string ModelConf, MCMCConf, FileOut;

    cout << "\n *** SusyFit Markov Chain Montecarlo ***\n" << endl;

    try {
        options_description desc("Allowed (positional) options");
        desc.add_options()
                ("modconf", value<string > (), "model config filename (1st)")
                ("mcconf", value<string > (), "montecarlo config filename (2nd)")
                ("rootfile", value<string > ()->default_value("MCout.root"),
                "output root filename (3rd)")
                ("help", "help message")
                ;
        positional_options_description pd;
        pd.add("modconf", 1);
        pd.add("mcconf", 1);
        pd.add("rootfile", 1);

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
            error("missing mandatory model config filename");

        if (vm.count("mcconf"))
            MCMCConf = vm["mcconf"].as<string > ();
        else throw
            error("missing mandatory montecarlo config filename");

        FileOut = vm["rootfile"].as<string > ();

        MonteCarlo MC(ModelConf, MCMCConf, FileOut);
        MC.Run();

        return EXIT_SUCCESS;
    } catch (const error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
