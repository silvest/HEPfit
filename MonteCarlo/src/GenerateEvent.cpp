/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GenerateEvent.h"
#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#ifdef _MPI
#include <mpi.h>
#endif

GenerateEvent::GenerateEvent(const std::string& ModelConf_i,
                             const std::string& OutDirName_i,
                             const std::string& JobTag_i,
                             const bool noMC_i,
                             const bool checkTheoryRange_i)
: myInputParser()
{
    outputTerm = 0;
    ModelConf = ModelConf_i;
    JobTag = JobTag_i;
    if (OutDirName_i != ""){
        OutDirName = "GeneratedEvents/" + OutDirName_i + JobTag;
        OldOutDirName = ("GeneratedEvents/OLD/" + OutDirName_i + JobTag);
        ObsDirName = OutDirName + "/Observables";
        ParsDirName = OutDirName + "/Parameters";
        outputTerm = 1;
    }
    noMC = noMC_i;
    nPar = 0;
    nObs = 0;
}

GenerateEvent::~GenerateEvent() 
{
}

void GenerateEvent::generate(const int rank, int unsigned nIteration, int seed)
{
        if (!noMC)
        throw std::runtime_error("\nGenerateEvent::generate():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    if (rank == 0){
        FileStat_t info;
        if (gSystem->GetPathInfo("GeneratedEvents", info) != 0) {
            gSystem->MakeDirectory("GeneratedEvents");
        }
        if (gSystem->GetPathInfo("GeneratedEvents/OLD", info) != 0) {
            gSystem->MakeDirectory("GeneratedEvents/OLD");
        }
        if (gSystem->GetPathInfo(OutDirName.c_str(), info) == 0){
            if (gSystem->GetPathInfo(OldOutDirName.c_str(), info) == 0){
                gSystem->Exec(("rm -rf " + OldOutDirName).c_str());
                std::cout << "WARNING: Removed " << OldOutDirName << std::endl;
            }
            gSystem->Exec(("mv -f " + OutDirName + " " + "GeneratedEvents/OLD/").c_str());
            std::cout << "WARNING: " << OutDirName << " exists!\nWARNING: " << OutDirName <<" moved to " << OldOutDirName << "\n"<< std::endl;
        }
        gSystem->MakeDirectory(OutDirName.c_str());
        gSystem->MakeDirectory(ObsDirName.c_str());
        gSystem->MakeDirectory(ParsDirName.c_str());
    }
    try {
        /* set model parameters */
        std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
        //MCEngine.SetName(ModelName.c_str());
        //MCEngine.MCMCSetRandomSeed(seed);
        std::map<std::string, double> DP;
        std::vector<boost::shared_ptr<ofstream> > ObsOut;
        std::vector<boost::shared_ptr<ofstream> > ParsOut;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
            if (it->errg > 0. || it->errf > 0.){
                ModParsVar.push_back(*it);
                DPars[it->name] = 0.;
                ParsOut.push_back( boost::make_shared<ofstream>((ParsDirName + "/" + it->name + ".txt").c_str(), std::ios::out));
                ParsOut[nPar]->close();
                nPar++;
            }
        }
        int buffsize = ModParsVar.size() + 1;
        if (!myInputParser.getMyModel()->Init(DP)) {
            throw std::runtime_error("parameter(s) missing in model initialization");
        }
        defineParameters();
        if (nIteration == 0) {
            outputTerm = 0;
            std::cout.precision(10);
            std::cout << std::endl << "Running in Single Event mode...\nNo data will be written to disk!\n" << std::endl;
        } else {
        std::cout << std::endl << "Generating " << nIteration << " random events...\n" << std::endl;
        }
        for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++){
            ObsOut.push_back( boost::make_shared<ofstream>((ObsDirName + "/" + it->getName() + ".txt").c_str(), std::ios::out));
            
            nObs++;
        }
        for (int unsigned i = 0; i < nIteration + 1; i++) {
            int j = 0;
            generateRandomEvent(myInputParser.getMyModel(),i);
            for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                double th = it->computeTheoryValue();
                if (outputTerm == 0) std::cout << it->getName() << " = " << th << std::endl;
            }
            if (outputTerm == 0) std::cout << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++) {
                std::vector<Observable> ObsInCGO = it->getObs();
                for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                    double th = it2->computeTheoryValue();
                    if (outputTerm == 0) std::cout << it2->getName() << " = " << th << std::endl;
                }
                if (outputTerm == 0) std::cout << std::endl;
            }
        }
        for (int unsigned i = 0; i < nObs; i++) {
             ObsOut[i]->close();
        }
        for (int unsigned i = 0; i < nPar; i++) {
             ParsOut[i]->close();
        }
        return;
        } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
        }
        if (outputTerm == 1) std::cout << "Output written to: " << OutDirName << std::endl;
}

void GenerateEvent::defineParameters() 
{
   /* // Add parameters to your model here.
    // You can then use them in the methods below by calling the
    // parameters.at(i) or parameters[i], where i is the index
    // of the parameter. The indices increase from 0 according to the
    // order of adding the parameters.
    std::cout << "Parameters varied in Monte Carlo:" << std::endl;
    int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModPars.begin();
            it < ModPars.end(); it++) {
        if (it->errf == 0. && it->errg == 0.)
            continue;
        //AddParameter(it->name.c_str(), it->min, it->max);
        std::cout << "  " << it->name << " " << k << std::endl;
        DPars[it->name] = 0.;
        if (it->errf == 0.) SetPriorGauss(k, it->ave, it->errg);
        else if (it->errg == 0.) SetPriorConstant(k);
        else {
            TF1 * combined = new TF1(it->name.c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->min, it->max);
            combined->SetParameter(0, it->ave);
            combined->SetParameter(1, it->errg);
            combined->SetParameter(2, it->errf);
            SetPrior(k, combined);
            delete combined;
        }
        k++;
    }*/
}

void GenerateEvent::generateRandomEvent(Model* Mod_i, int iterationNo)
{
    Mod = Mod_i;
    //if (iterationNo == 0) { MCMCSetFlagInitialPosition(0); }
    //else { MCMCSetFlagInitialPosition(1); }
    
    //MCMCInitialize();
    if (iterationNo == 0){
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
        DPars[it->name] = it->ave;
        }
    } else {
        throw std::runtime_error("\nGenerateEvent::generateRandomEvent():\nRandom number generator not implemented.\n");
    }
    Mod->Update(DPars);
}
