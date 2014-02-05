/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GenerateEvent.h"
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
#ifdef _MPI
#include <mpi.h>
#endif

GenerateEvent::GenerateEvent(const std::string& ModelConf_i,
                             const std::string& OutDirName_i,
                             const std::string& JobTag_i,
                             const bool noMC_i)
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
        CGODirName = OutDirName + "/CGO";
        outputTerm = 1;
    }
    noMC = noMC_i;
}

GenerateEvent::~GenerateEvent() 
{
}

void GenerateEvent::generate(const int rank, int unsigned nIteration, int seed)
{
        if (!noMC)
        throw std::runtime_error("\nGenerateEvent::generate():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    gRandom->SetSeed(seed);
    if (nIteration == 0) outputTerm = 0;
    if (rank == 0 && nIteration > 0){
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
        gSystem->MakeDirectory(CGODirName.c_str());
        gSystem->MakeDirectory(ParsDirName.c_str());
    }
    try {
        /* set model parameters */
        std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
        if (Obs.size() == 0 && CGO.size() == 0) throw std::runtime_error("\nGenerateEvent::generate(): No observables or correlated Gaussian observables defined in SomeModel.conf file\n");
        std::ofstream summary;
        if (outputTerm == 1){
            summary.open((OutDirName + "/Summary.txt").c_str(), std::ios::out);
            summary << "Model\t" << ModelName << "\n";
        }
        std::map<std::string, double> DP;
        for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
            DP[it->name] = it->ave;
            if (it->errg > 0. || it->errf > 0.){
                ModParsVar.push_back(*it);
                DPars[it->name] = 0.;
                if (outputTerm == 1){
                    DDist[it->name] = NULL;
                    ParsOut[it->name] = boost::make_shared<std::ofstream>((ParsDirName + "/" + it->name + ".txt").c_str(), std::ios::out);
                    summary << "Parameter\t" << it->name << "\n";
                }
            }
        }
        //int buffsize = ModParsVar.size() + 1;
        if (!myInputParser.getMyModel()->Init(DP)) {
            throw std::runtime_error("parameter(s) missing in model initialization");
        }
        defineParameterDistributions();
        if (nIteration == 0) {
            outputTerm = 0;
            std::cout.precision(10);
            std::cout << std::endl << "Running in Single Event mode...\nNo data will be written to disk!\n" << std::endl;
        } else {
            std::cout << "\nRunning in Event Generation mode... \nWARNING: The output should not be used for any statistical analysis. \n         Neither randomness or completness of the sample is gauranteed!!\n" << std::endl;
        }
        for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++){
            if (outputTerm == 1) {
                ObsOut[it->getName()] = boost::make_shared<std::ofstream>((ObsDirName + "/" + it->getName() + ".txt").c_str(), std::ios::out);
                summary << "Observable\t" << it->getName() << "\n";
            }
        }
        for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++){
            if (outputTerm == 1) {
                CGOOut[it->getName()] = boost::make_shared<std::ofstream>((CGODirName + "/" + it->getName() + ".txt").c_str(), std::ios::out);
                summary << "CorrelatedGaussianObservables\t" << it->getName() << "\n";
            }
        }
        Mod = myInputParser.getMyModel();
        for (int unsigned i = 0; i < nIteration + 1; i++) {
            if (i == 1) std::cout << std::endl << "Generating " << nIteration << " random events...\n" << std::endl;
            generateRandomEvent(i);
            if (outputTerm == 0) std::cout << "\nOservables: \n" << std::endl;
            for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                double th = it->computeTheoryValue();
                if (outputTerm == 0) {
                    std::cout << it->getName() << " = " << th << std::endl;
                } else {
                    (*ObsOut[it->getName()]) << th << std::endl;
                    if (i == nIteration) ObsOut[it->getName()]->close();
                }
            }
            if (outputTerm == 0) std::cout << "\nCorrelated Gaussian Oservables: \n" << std::endl;
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                if (outputTerm == 0) std::cout << it1->getName() << ":" << std::endl;
                std::vector<Observable> ObsInCGO = it1->getObs();
                for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                    double th = it2->computeTheoryValue();
                    if (outputTerm == 0) {
                        std::cout << it2->getName() << " = " << th << std::endl;
                    } else {
                        (*CGOOut[it1->getName()]) << th << "\t";
                    }
                }
                if (outputTerm == 0) {
                    std::cout << std::endl;
                } else {
                    (*CGOOut[it1->getName()]) << std::endl;
                    if (i == nIteration) CGOOut[it1->getName()]->close();
                }
            }
            if (outputTerm == 0) std::cout << std::endl;
        }
        if (outputTerm == 1){
            for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
                ParsOut[it->name]->close();
            }
            summary << "Events\t" << nIteration << "\n";
            summary.close();
            std::cout << "Output written to: " << OutDirName << "\n" << std::endl;
        }
        return;
        } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
        }
}

void GenerateEvent::defineParameterDistributions()
{
    if (outputTerm == 0) std::cout << "\nParameters varied in Event Generation:" << std::endl;
    for (std::vector<ModelParameter>::const_iterator it = ModParsVar.begin();
            it < ModParsVar.end(); it++) {
        //if (it->errf == 0. && it->errg == 0.)
            //continue;
        if (outputTerm == 0) std::cout << it->name << ", "; //<< k << std::endl;
        if (it->errf == 0. && it->errg != 0.){
            DDist[it->name] = new TF1(it->name.c_str(),
            "1./sqrt(2.*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2./[1]/[1])",
            it->min,it->max);
            DDist[it->name]->SetParameter(0, it->ave);
            DDist[it->name]->SetParameter(1, it->errg);
        }
        else if (it->errg == 0. && it->errf != 0.){
            DDist[it->name] = new TF1(it->name.c_str(),
                    "1/([1]-[0])",
                    it->min, it->max);
            DDist[it->name]->SetParameter(0,it->min);
            DDist[it->name]->SetParameter(1,it->max);
        }
        else {
            DDist[it->name] = new TF1(it->name.c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->min, it->max);
            DDist[it->name]->SetParameter(0, it->ave);
            DDist[it->name]->SetParameter(1, it->errg);
            DDist[it->name]->SetParameter(2, it->errf);
        }
    }
    if (outputTerm == 0) std::cout << std::endl;
}

void GenerateEvent::generateRandomEvent(int iterationNo)
{
    if (iterationNo == 0){
        if (outputTerm == 0) std::cout << "\nThe central event is: " <<std::endl;
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
            DPars[it->name] = it->ave;
            if (outputTerm == 1) (*ParsOut[it->name]) << DPars[it->name] << std::endl;
        }
    } else {
        if (outputTerm == 0) std::cout << "\nEvent no.: " << iterationNo <<std::endl;
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
            DPars[it->name] = DDist[it->name]->GetRandom();
            if (outputTerm == 1) (*ParsOut[it->name]) << DPars[it->name] << std::endl;
        }
    }
    Mod->Update(DPars);
}
