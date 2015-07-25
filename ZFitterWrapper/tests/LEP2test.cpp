/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <InputParser.h>
#include <ModelParameter.h>
#include <Observable.h>
#include <Observable2D.h>
#include "ZFitterWrapper.h"
#include "ZFEWObservables.h"

using namespace boost::program_options;
using namespace std;

int main(int argc, char** argv) {

    string ModelConf, FileOut;

    try {

        // get the name of the current directory
        char dir[255];
        getcwd(dir,255);
        
        // set filenames
        options_description desc("Allowed (positional) options");
        string inputFile(dir), outputFile(dir);
        inputFile += "/tests/StandardModel.conf";
        outputFile += "/tests/LEP2test.root";
        desc.add_options()
            ("modconf", value<string>()->default_value(inputFile), 
                "model config filename (1st)")
            ("rootfile", value<string>()->default_value(outputFile), 
                "output root filename (2nd)")
            ("help", "help message")
            ;
        positional_options_description pd;
        pd.add("modconf", 1);
        pd.add("rootfile", 1);
        variables_map vm;
        store(command_line_parser(argc, argv).options(desc).positional(pd).run(), vm);
        notify(vm);
        if (vm.count("help")) {
            cout << desc << endl;
            return EXIT_SUCCESS;
        }
        ModelConf = vm["modconf"].as<string > ();
        FileOut = vm["rootfile"].as<string>();
        cout << "# set " << ModelConf << " for the model config file" << endl;
        cout << "# set " << FileOut << " for the output file" << endl;
        
        // read parameters
        InputParser myInputParser;
        std::vector<ModelParameter> ModPars;
        std::vector<Observable> Obs;
        std::vector<Observable2D> Obs2D;
        std::vector<CorrelatedGaussianObservables> CGO;
        std::vector<ModelParaVsObs> ParaObs;
        std::map<string, double> DP;
        cout << "# Model: " 
             << myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs)
             << endl;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++)
            DP[it->name] = it->ave;
        if (!myInputParser.getMyModel()->Init(DP)) {
            cout << "parameter(s) missing in model initialization" << endl;
            exit(EXIT_FAILURE);
        }
        
        //////////////////////////////////////////////////////////////////////

        // TEST for inputs
        //cout << "Mz = " << myInputParser.getMyModel()->getMz() << endl;
        
        ZFitterWrapper* ZF = new ZFitterWrapper(*myInputParser.getMyModel());
        ZF->FlagInfo();
        
        // TESTs
        //ZF->printIntermediateResults();
        //ZF->printInputs();
        //ZF->printConstants();
        
        ZFMw* myZFMw = new ZFMw(*ZF);
        cout << "Mw = " << myZFMw->computeThValue() << " " << endl << endl;
        
        //////////////////////////////////////////////////////////////////////        
        // Differential cross sections

        const double cos_theta[11] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0,  
                                      0.2, 0.4, 0.6, 0.8, 1.0};
        ZFDsigmaQuarksLEP2* myZFDsigmaQuarks[11];
        ZFDsigmaMuLEP2* myZFDsigmaMu[11];
        ZFDsigmaTauLEP2* myZFDsigmaTau[11];
        cout << "sqrt{s} cos(theta) dsigma(q)/dcos_theta dsigma(mu)/dcos_theta dsigma(tau)/dcos_theta"
             << endl;
        for (int i=0; i<11; i++) {
            myZFDsigmaQuarks[i] = new ZFDsigmaQuarksLEP2(*ZF,200.0,cos_theta[i]);
            myZFDsigmaMu[i] = new ZFDsigmaMuLEP2(*ZF,200.0,cos_theta[i]);
            myZFDsigmaTau[i] = new ZFDsigmaTauLEP2(*ZF,200.0,cos_theta[i]);        
            cout << " 200.0  " << setw(6) <<  cos_theta[i]
                 << setw(20) << myZFDsigmaQuarks[i]->computeThValue() 
                 << setw(20) << myZFDsigmaMu[i]->computeThValue() 
                 << setw(20) << myZFDsigmaTau[i]->computeThValue()
                 << endl;
        }
        cout << endl;
        
        //////////////////////////////////////////////////////////////////////

        // LEP2 CM energies
        const double sqrt_s[12] = {130.0, 136.0, 161.0, 172.0, 183.0, 189.0, 
                                   192.0, 196.0, 200.0, 202.0, 205.0, 207.0};
        // LEP2 CM energies for heavy flavors
        const double sqrt_s_HF[10] = {133.0, 167.0, 183.0, 189.0, 192.0, 
                                      196.0, 200.0, 202.0, 205.0, 207.0};
        
        ZFsigmaQuarksLEP2* myZFsigmaQuarks[12];
        ZFsigmaMuLEP2* myZFsigmaMu[12];
        ZFsigmaTauLEP2* myZFsigmaTau[12];
        ZFAFBmuLEP2* myZFAFBmu[12];
        ZFAFBtauLEP2* myZFAFBtau[12];
        cout << "sqrt{s} sigma(q) sigma(mu) sigma(tau) A_FB(mu) A_FB(tau)"
             << endl;
        for (int i=0; i<12; i++) {
            myZFsigmaQuarks[i] = new ZFsigmaQuarksLEP2(*ZF,sqrt_s[i]);
            myZFsigmaMu[i] = new ZFsigmaMuLEP2(*ZF,sqrt_s[i]);
            myZFsigmaTau[i] = new ZFsigmaTauLEP2(*ZF,sqrt_s[i]);
            myZFAFBmu[i] = new ZFAFBmuLEP2(*ZF,sqrt_s[i]);
            myZFAFBtau[i] = new ZFAFBtauLEP2(*ZF,sqrt_s[i]);
            cout << setw(6) << sqrt_s[i]
                 << setw(10) << myZFsigmaQuarks[i]->computeThValue() 
                 << setw(10) << myZFsigmaMu[i]->computeThValue() 
                 << setw(10) << myZFsigmaTau[i]->computeThValue()
                 << setw(10) << myZFAFBmu[i]->computeThValue() 
                 << setw(10) << myZFAFBtau[i]->computeThValue()
                 << endl;
            //ZF->CutInfo();// TEST
        }
        cout << endl;

        //////////////////////////////////////////////////////////////////////
       
        const int INTF_NEW = 2;
        ZF->setFlag("INTF",INTF_NEW); // with ISR/FSR interference contribution
        cout << "Flag update: INTF=2" << endl << endl;
        
        ZFRbottomLEP2* myZFRbottom[10];
        ZFAFBbottomLEP2* myZFAFBbottom[10];
        ZFRcharmLEP2* myZFRcharm[10];
        ZFAFBcharmLEP2* myZFAFBcharm[10];        
        cout << "sqrt{s}    R_b     A_FB(b)     R_c     A_FB(c)" << endl;
        for (int i=0; i<10; i++) {
           myZFRbottom[i] = new ZFRbottomLEP2(*ZF,sqrt_s_HF[i]);
           myZFAFBbottom[i] = new ZFAFBbottomLEP2(*ZF,sqrt_s_HF[i]);
           myZFRcharm[i] = new ZFRcharmLEP2(*ZF,sqrt_s_HF[i]);
           myZFAFBcharm[i] = new ZFAFBcharmLEP2(*ZF,sqrt_s_HF[i]);
           cout << setw(6) << sqrt_s_HF[i]
                << setw(10) << myZFRbottom[i]->computeThValue() 
                << setw(10) << myZFAFBbottom[i]->computeThValue() 
                << setw(10) << myZFRcharm[i]->computeThValue() 
                << setw(10) << myZFAFBcharm[i]->computeThValue()
                << endl;
        }
        cout << endl;
        
        //////////////////////////////////////////////////////////////////////        
        cout << "Test finished" << endl;
        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}
