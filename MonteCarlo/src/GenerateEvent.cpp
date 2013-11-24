/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GenerateEvent.h"
#include <TSystem.h>
#ifdef _MPI
#include <mpi.h>
#endif

GenerateEvent::GenerateEvent(const std::string& ModelConf_i,
                             const std::string& OutDirName_i,
                             const std::string& JobTag_i,
                             const bool noMC_i,
                             const bool checkTheoryRange_i)
: myInputParser(), MCEngine(ModPars, Obs, Obs2D, CGO, ParaObs, checkTheoryRange_i)
{
    ModelConf = ModelConf_i;
    JobTag = JobTag_i;
    OutDirName = OutDirName_i + JobTag;
    ObsDirName = OutDirName + JobTag + "/Observables";
    ParsDirName = OutDirName + JobTag + "/Parameters";
    noMC = noMC_i;
}

GenerateEvent::~GenerateEvent() 
{
}

void GenerateEvent::generate(const int rank, int unsigned nIteration, int seed)
{
        if(!noMC)
        throw std::runtime_error("\nMonteCarlo::generateEvent():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    try {
        /* set model parameters */
        std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
        MCEngine.SetName(ModelName.c_str());
        MCEngine.MCMCSetRandomSeed(seed);
        int npar = 0;
        std::map<std::string, double> DP;
        std::vector<double> INI;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
            if (it->errg > 0. || it->errf > 0.){
                INI.push_back(it->ave);
                npar++;}
        }
        int buffsize = npar++;
        if (!myInputParser.getMyModel()->Init(DP)) {
            throw std::runtime_error("parameter(s) missing in model initialization");
        }
        MCEngine.DefineParameters();
        if (nIteration == 0) {
            std::cout.precision(10);
            std::cout << std::endl << "Running in Single Event mode...\n" << std::endl;
        }
        for (int unsigned i = 0; i < nIteration + 1; i++) {
            MCEngine.generateRandomEvent(myInputParser.getMyModel(),i);
            for (std::vector<Observable>::iterator it = Obs.begin();
                 it < Obs.end(); it++) {
                double th = it->computeTheoryValue();
                std::cout << it->getName() << " = " << th << std::endl;
            }
            std::cout << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin();
                 it < CGO.end(); it++) {
                std::vector<Observable> ObsInCGO = it->getObs();
                for (std::vector<Observable>::iterator it2 = ObsInCGO.begin();
                     it2 < ObsInCGO.end(); it2++) {
                    double th = it2->computeTheoryValue();
                    std::cout << it2->getName() << " = " << th << std::endl;
                }
                std::cout << std::endl;
            }
        }
        return;
        } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
        }
}

