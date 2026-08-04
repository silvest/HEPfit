/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <MonteCarlo.h>
#include <GenerateEvent.h>
#include <InputParser.h>
#include <Observable.h>
//#include <mtbar.h>
#ifdef _MPI
#include <mpi.h>
#endif

using namespace std;

namespace {

struct ParsedOptions {
    string modconf;
    string mcconf;
    bool hasModconf = false;
    bool hasMcconf = false;
    bool test = false;
    bool noMC = false;
    int it = 0;
    string rootfile = "MCout";
    string output_folder;
    string job_tag;
    bool weight = false;
    bool loadRun = false;
    string inFile = "MCin";
    bool help = false;
};

static string requireValue(const char* optionName, int& index, int argc, char** argv)
{
    if (index + 1 >= argc) {
        throw runtime_error(string("Missing value for option ") + optionName + ".");
    }
    ++index;
    return argv[index];
}

static ParsedOptions parseOptions(int argc, char** argv)
{
    ParsedOptions parsed;
    std::vector<string> positional;

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];

        if (arg == "--help") {
            parsed.help = true;
        } else if (arg == "--test") {
            parsed.test = true;
        } else if (arg == "--noMC") {
            parsed.noMC = true;
        } else if (arg == "--weight") {
            parsed.weight = true;
        } else if (arg == "--loadRun") {
            parsed.loadRun = true;
        } else if (arg.rfind("--it=", 0) == 0) {
            parsed.it = stoi(arg.substr(5));
        } else if (arg == "--it") {
            parsed.it = stoi(requireValue("--it", i, argc, argv));
        } else if (arg.rfind("--rootfile=", 0) == 0) {
            parsed.rootfile = arg.substr(11);
        } else if (arg == "--rootfile") {
            parsed.rootfile = requireValue("--rootfile", i, argc, argv);
        } else if (arg.rfind("--output_folder=", 0) == 0) {
            parsed.output_folder = arg.substr(16);
        } else if (arg == "--output_folder") {
            parsed.output_folder = requireValue("--output_folder", i, argc, argv);
        } else if (arg.rfind("--job_tag=", 0) == 0) {
            parsed.job_tag = arg.substr(10);
        } else if (arg == "--job_tag") {
            parsed.job_tag = requireValue("--job_tag", i, argc, argv);
        } else if (arg.rfind("--inFile=", 0) == 0) {
            parsed.inFile = arg.substr(9);
        } else if (arg == "--inFile") {
            parsed.inFile = requireValue("--inFile", i, argc, argv);
        } else if (arg.rfind("--", 0) == 0) {
            throw runtime_error(string("Unknown option: ") + arg);
        } else {
            positional.push_back(arg);
        }
    }

    if (positional.size() > 2U) {
        throw runtime_error("Too many positional arguments. Expected at most: modconf mcconf.");
    }
    if (!positional.empty()) {
        parsed.modconf = positional[0];
        parsed.hasModconf = true;
    }
    if (positional.size() > 1U) {
        parsed.mcconf = positional[1];
        parsed.hasMcconf = true;
    }

    return parsed;
}

}

int main(int argc, char** argv) 
{

#ifdef _MPI
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //MPI_Status status;
#else
    int rank = 0;
#endif

    string ModelConf, JobTag;
    string FolderOut;
    int nIterations = 0;
    //bool noMC;

    try {
        const string desc =
                "\nAllowed (positional) options\n"
                "  modconf            MANDATORY model config filename (1st)\n"
                "  mcconf             MANDATORY Monte Carlo config filename (2nd) for MCMC run\n"
                "  --test             run in test mode to generate a single event\n"
                "  --noMC             run in generate event mode\n"
                "  --it <n>           no. of iterations in generate event mode\n"
                "  --rootfile <name>  output root filename (without extension)\n"
                "  --output_folder <name>\n"
                "                     output folder for Generate Event mode\n"
                "  --job_tag <tag>    job tag\n"
                "  --weight           run in generate event mode with weights\n"
                "  --loadRun          load previous run data from file\n"
                "  --inFile <name>    input root filename (without extension)\n"
                "  --help             help message\n";
        string coderun = "\n *** HEPfit Routines ***\n"
                         "\nMonte Carlo mode: analysis Model.conf MonteCarlo.conf [--rootfile <name>] [--job_tag <tag>] [--thRange] [--loadRun --inFile <name>]"
                         "\nSingle Event mode: analysis Model.conf MonteCarlo.conf --test"
                         "\nGenerate Event mode: analysis Model.conf --noMC [--it #] [--weight] [--output_folder <name> [--job_tag <tag>]]\n";
        ParsedOptions opts;

        
        try {
            opts = parseOptions(argc, argv);

            if (opts.help) {
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                }
                return EXIT_SUCCESS;
            }       
            
            if (opts.hasModconf)
                ModelConf = opts.modconf;
            else
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                    throw runtime_error("\nERROR: Please specify mandatory model config filename.\n");
                }

            if (opts.noMC && opts.output_folder == "" && !(opts.job_tag == "")){
                if (rank == 0) cout << "\nWARNING: --job_tag not being used as --output_folder is not specified.\nOutput being sent to terminal only.\n";
            } else {
                JobTag = opts.job_tag;
            }
            
            if (opts.noMC && opts.hasMcconf)
            {
                if (rank == 0){
                    cout << coderun << endl;
                    cout << desc << endl;
                    throw std::runtime_error("\nERROR: Please specify EITHER --noMC OR Monte Carlo configuration file as the second free argument.\nAll other arguments must be specified with option specifiers (--nI=#, --job_tag <tag>, etc.)\n");
                }
            }
            else if (opts.noMC && !opts.hasMcconf)
            {
                if (rank == 0)
                    cout << "\n *** HEPfit Event Generation ***\n" << endl;
//                noMC = true;
                FolderOut = opts.output_folder;
                bool weight = false;
                nIterations = opts.it;
                if (opts.weight) weight = true;
                ThObsFactory ThObsF;
                ModelFactory ModelF;
                GenerateEvent GE(ModelF, ThObsF, ModelConf, FolderOut, JobTag);
//                GE.addCustomParser("PS", boost::factory<InputParser*>());
//                GE.addCustomObservableType("Poisson", boost::factory<Observable*>());
//                GE.linkParserToObservable("Poisson", "PS");
                GE.generate(nIterations, 1, weight);
            }
            else if (!opts.noMC && opts.hasMcconf)
            {
                if (rank == 0)
                    cout << "\n *** HEPfit Markov Chain Monte Carlo ***\n" << endl;
//                noMC = false;
                string MCMCConf = opts.mcconf;
                string FileOut = opts.rootfile;
                if (opts.loadRun) cout << "\nLoading previous run data from file " << opts.inFile << ".root\n" << endl;
                if (opts.loadRun && opts.inFile == FileOut)
                    if (rank == 0) throw runtime_error("\nERROR:  Input and output root filenames cannot be the same when loading previous run data from file.\n");
                
                
                ThObsFactory ThObsF;
                ModelFactory ModelF;
//                ThObsF.addObsToFactory("mtbar", boost::factory<mtbar*>());
                std::string FileIn = opts.loadRun ? opts.inFile + ".root" : "";
                MonteCarlo MC(ModelF, ThObsF, ModelConf, MCMCConf, FileOut, JobTag, FileIn);
                if (opts.test){
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
            

        } catch(const exception& e) {
            if (rank == 0)
                cout << coderun << endl;
            cerr << "\nERROR: " << e.what() << std::endl << std::endl
                 << desc << std::endl;
            return EXIT_FAILURE;
        }
        
#ifdef _MPI
        MPI_Finalize();
#endif
        
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
