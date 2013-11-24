/* 
 * File:   GenerateEvent.cpp
 * Author: ayanpaul
 * 
 * Created on November 24, 2013, 2:44 PM
 */

#include "GenerateEvent.h"
#include <TSystem.h>
#include <TF1.h>
#include <TMath.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <fstream>
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
    ModelConf = ModelConf_i;
    JobTag = JobTag_i;
    OutDirName = OutDirName_i + JobTag;
    ObsDirName = OutDirName + JobTag + "/Observables";
    ParsDirName = OutDirName + JobTag + "/Parameters";
    noMC = noMC_i;
    npar = 0;
}

GenerateEvent::~GenerateEvent() 
{
}

void GenerateEvent::generate(const int rank, int unsigned nIteration, int seed)
{
        if(!noMC)
        throw std::runtime_error("\nGenerateEvent::generate():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    try {
        /* set model parameters */
        std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
        //MCEngine.SetName(ModelName.c_str());
        //MCEngine.MCMCSetRandomSeed(seed);
        std::map<std::string, double> DP;
        //std::vector<double> INI;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
            if (it->errg > 0. || it->errf > 0.){
                ModParsVar.push_back(*it);
                DPars[it->name] = 0.;
                //INI.push_back(it->ave);
            }
            npar++;
        }
        int buffsize = ModParsVar.size() + 1;
        if (!myInputParser.getMyModel()->Init(DP)) {
            throw std::runtime_error("parameter(s) missing in model initialization");
        }
        defineParameters();
        if (nIteration == 0) {
            std::cout.precision(10);
            std::cout << std::endl << "Running in Single Event mode...\n" << std::endl;
        } else {
        std::cout << std::endl << "Generating " << nIteration << " random events...\n" << std::endl;
        }
        for (int unsigned i = 0; i < nIteration + 1; i++) {
            generateRandomEvent(myInputParser.getMyModel(),i);
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