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
#include <string>
#include <vector>
#include <unistd.h>
#include <InputParser.h>
#include <ModelParameter.h>
#include <Observable.h>
#include <Observable2D.h>
#include "ZFitterWrapper.h"
#include "ZFEWObservables.h"

using namespace std;

namespace {

struct ParsedOptions {
    string modconf;
    string rootfile;
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

static ParsedOptions parseOptions(int argc, char** argv, const string& defaultModconf, const string& defaultRootfile)
{
    ParsedOptions parsed;
    parsed.modconf = defaultModconf;
    parsed.rootfile = defaultRootfile;
    std::vector<string> positional;

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        if (arg == "--help") {
            parsed.help = true;
        } else if (arg.rfind("--modconf=", 0) == 0) {
            parsed.modconf = arg.substr(10);
        } else if (arg == "--modconf") {
            parsed.modconf = requireValue("--modconf", i, argc, argv);
        } else if (arg.rfind("--rootfile=", 0) == 0) {
            parsed.rootfile = arg.substr(11);
        } else if (arg == "--rootfile") {
            parsed.rootfile = requireValue("--rootfile", i, argc, argv);
        } else if (arg.rfind("--", 0) == 0) {
            throw runtime_error(string("Unknown option: ") + arg);
        } else {
            positional.push_back(arg);
        }
    }

    if (positional.size() > 2U) {
        throw runtime_error("Too many positional arguments. Expected at most: modconf rootfile.");
    }
    if (!positional.empty()) {
        parsed.modconf = positional[0];
    }
    if (positional.size() > 1U) {
        parsed.rootfile = positional[1];
    }

    return parsed;
}

}

int main(int argc, char** argv) {

    string ModelConf, FileOut;

    try {

        // get the name of the current directory
        char dir[255];
        getcwd(dir,255);
        
        // set filenames
        const string desc =
            "Allowed (positional) options\n"
            "  modconf              model config filename (1st)\n"
            "  rootfile             output root filename (2nd)\n"
            "  --modconf <path>     model config filename\n"
            "  --rootfile <path>    output root filename\n"
            "  --help               help message\n";
        string inputFile(dir), outputFile(dir);
        inputFile += "/tests/StandardModel.conf";
        outputFile += "/tests/LEP2test.root";
        ParsedOptions opts = parseOptions(argc, argv, inputFile, outputFile);
        if (opts.help) {
            cout << desc << endl;
            return EXIT_SUCCESS;
        }
        ModelConf = opts.modconf;
        FileOut = opts.rootfile;
        cout << "# set " << ModelConf << " for the model config file" << endl;
        cout << "# set " << FileOut << " for the output file" << endl;
        
        ThObsFactory ThObsF;
        ModelFactory ModelF;
        
        // read parameters
        InputParser myInputParser(ModelF, ThObsF);
        std::vector<ModelParameter> ModPars;
        boost::ptr_vector<Observable> Obs;
        std::vector<Observable2D> Obs2D;
        std::vector<CorrelatedGaussianObservables> CGO;
        std::vector<CorrelatedGaussianParameters> CGP;
        //std::vector<ModelParaVsObs> ParaObs;
        int rank = 0;
        std::map<string, double> DP;
        cout << "# Model: " 
             << myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, CGP)
             << endl;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++)
            DP[it->getname()] = it->getave();
        if (!myInputParser.getModel()->Init(DP)) {
            cout << "parameter(s) missing in model initialization" << endl;
            exit(EXIT_FAILURE);
        }
        
        //////////////////////////////////////////////////////////////////////

        // TEST for inputs
        //cout << "Mz = " << myInputParser.getMyModel()->getMz() << endl;
        
        ZFitterWrapper* ZF = new ZFitterWrapper(*myInputParser.getModel());
        ZF->FlagInfo();
        
        // TESTs
        //ZF->printIntermediateResults();
        //ZF->printInputs();
        //ZF->printConstants();
        
        ZFMw* myZFMw = new ZFMw(*myInputParser.getModel());
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
            myZFDsigmaQuarks[i] = new ZFDsigmaQuarksLEP2(*myInputParser.getModel(),200.0,cos_theta[i]);
            myZFDsigmaMu[i] = new ZFDsigmaMuLEP2(*myInputParser.getModel(),200.0,cos_theta[i]);
            myZFDsigmaTau[i] = new ZFDsigmaTauLEP2(*myInputParser.getModel(),200.0,cos_theta[i]);        
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
            myZFsigmaQuarks[i] = new ZFsigmaQuarksLEP2(*myInputParser.getModel(),sqrt_s[i]);
            myZFsigmaMu[i] = new ZFsigmaMuLEP2(*myInputParser.getModel(),sqrt_s[i]);
            myZFsigmaTau[i] = new ZFsigmaTauLEP2(*myInputParser.getModel(),sqrt_s[i]);
            myZFAFBmu[i] = new ZFAFBmuLEP2(*myInputParser.getModel(),sqrt_s[i]);
            myZFAFBtau[i] = new ZFAFBtauLEP2(*myInputParser.getModel(),sqrt_s[i]);
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
           myZFRbottom[i] = new ZFRbottomLEP2(*myInputParser.getModel(),sqrt_s_HF[i]);
           myZFAFBbottom[i] = new ZFAFBbottomLEP2(*myInputParser.getModel(),sqrt_s_HF[i]);
           myZFRcharm[i] = new ZFRcharmLEP2(*myInputParser.getModel(),sqrt_s_HF[i]);
           myZFAFBcharm[i] = new ZFAFBcharmLEP2(*myInputParser.getModel(),sqrt_s_HF[i]);
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
