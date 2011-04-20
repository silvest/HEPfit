/* 
 * File:   main.cpp
 * Author: silvest
 *
 * Created on March 30, 2011, 12:54 PM
 */

#include <cstdlib>
#include <iostream>
#include <istream>
#include <Observable.h>
#include <ModelParameter.h>
#include <MonteCarlo.h>
#include <InputParser.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCModelOutput.h>
#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCH1D.h>
#include <root/TRandom3.h>

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
    std::vector<Observable2D> Obs2D;   
    myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D);
    std::map<std::string,double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) 
        DP[it->name] = it->ave;
    
    if(!myInputParser.GetMyModel()->init(DP)){
        std::cout << "parameter(s) missing in model initialization" <<std::endl;
        exit(EXIT_FAILURE);
    }
        
    std::cout << ModPars.size() << " parameters defined." << std::endl;
    std::cout << Obs.size() << " observables defined." << std::endl;
    MonteCarlo myMonteCarlo("MonteCarlo", myInputParser.GetMyModel(),ModPars,Obs,Obs2D);
    BCModelOutput out(&myMonteCarlo,"MCout.root");
    //MonteCarlo configuration parser
    std::ifstream ifile(MCMCConf);
    std::string line;
    while(!getline(ifile,line).eof()){
        if(line.at(0)=='#') continue;
        boost::char_separator<char> sep(" ");
        boost::tokenizer<boost::char_separator<char> > tok(line,sep);
        boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin();
        if(beg->compare("NChains")==0) {
            ++beg;
            myMonteCarlo.MCMCSetNChains(atoi((*beg).c_str()));
        }
        else if(beg->compare("PrerunMaxIter")==0){
            ++beg;
            myMonteCarlo.MCMCSetNIterationsMax(atoi((*beg).c_str()));
        }
        else if(beg->compare("Seed")==0){
            ++beg;
            int seed = atoi((*beg).c_str());
            if(seed != 0)
                myMonteCarlo.MCMCGetTRandom3()->SetSeed(seed);
        }
        else if(beg->compare("Iterations")==0){
            ++beg;
            myMonteCarlo.MCMCSetNIterationsRun(atoi((*beg).c_str()));
        }
        else if(beg->compare("WriteChain")==0){
            ++beg;
            if(beg->compare("true")==0)
            {
                out.WriteMarkovChain(true);
                myMonteCarlo.AddChains();
            }
        } else {
            std::cout << "wrong keyword in MonteCarlo config file: " << *beg << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog("log.txt");
    BCLog::SetLogLevel(BCLog::detail);
    
    // set precision
    //myMonteCarlo.MCMCSetPrecision(BCEngineMCMC::kLow);
       
   // run the MCMC and marginalize w.r.t. to all parameters
    myMonteCarlo.MarginalizeAll();

    // find mode using the best fit parameters as start values
//    myMonteCarlo.FindMode(myMonteCarlo.GetBestFitParameters());

    // draw all marginalized distributions into a PostScript file
    myMonteCarlo.PrintAllMarginalized("MonteCarlo_plots.ps");

    // print results of the analysis into a text file
    myMonteCarlo.PrintResults("MonteCarlo_results.txt");

    // print ratio
    myMonteCarlo.PrintHistogram(out);
    
    out.WriteMarginalizedDistributions();
    out.FillAnalysisTree();
    out.Close();

    // close log file
    BCLog::CloseLog();

    return 0;
}

