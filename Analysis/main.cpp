/* 
 * File:   main.cpp
 * Author: silvest
 *
 * Created on March 30, 2011, 12:54 PM
 */

#include <cstdlib>
#include <iostream>
#include <Observable.h>
#include <ModelParameter.h>
#include <MonteCarlo.h>
#include <InputParser.h>
#include <BCSummaryTool.h>
#include <BCAux.h>
#include <BCLog.h>
#include <BCH1D.h>

/*
 * 
 */
int main(int argc, char** argv) {
    char * ModelConf, * MCMCConf;
    switch (argc) {
        case 3:
            ModelConf = argv[1];
            MCMCConf = argv[2];
            break;
        default:
            std::cout << "Usage: " << argv[0] << " ModelConf MCMCConf\n";
            exit(EXIT_FAILURE);
    }
    InputParser myInputParser;
    std::vector<ModelParameter> ModPars;
    std::vector<Observable> Obs;
    myInputParser.ReadParameters(ModelConf,ModPars,Obs);
    std::map<std::string,double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) 
        DP[it->name] = it->ave;
    
    if(!myInputParser.GetMyModel()->init(DP)){
        std::cout << "parameter(s) missing in model initialization" <<std::endl;
        exit(EXIT_FAILURE);
    }
        
    std::cout << ModPars.size() << " parameters defined." << std::endl;
    std::cout << Obs.size() << " observables defined." << std::endl;
    MonteCarlo myMonteCarlo("MonteCarlo", myInputParser.GetMyModel(),ModPars,Obs);
    //ToDo Implement MonteCarlo configuration parser

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);

    // set precision
    myMonteCarlo.MCMCSetPrecision(BCEngineMCMC::kHigh);

    // run the MCMC and marginalize w.r.t. to all parameters
    myMonteCarlo.MarginalizeAll();

    // find mode using the best fit parameters as start values
    myMonteCarlo.FindMode(myMonteCarlo.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    myMonteCarlo.PrintAllMarginalized("MonteCarlo_plots.ps");

    // print results of the analysis into a text file
    myMonteCarlo.PrintResults("MonteCarlo_results.txt");

    // print ratio
    myMonteCarlo.PrintHistogram();

    // close log file
    BCLog::CloseLog();

    return 0;
}

