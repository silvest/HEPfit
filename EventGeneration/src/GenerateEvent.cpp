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

GenerateEvent::GenerateEvent(ModelFactory& ModelF, ThObsFactory& ThObsF,
        const std::string& ModelConf_i,
                             const std::string& OutDirName_i,
                             const std::string& JobTag_i,
                             const bool noMC_i)
: myInputParser(ModelF, ThObsF)
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
        SMDebugDirName = OutDirName + "/SMDebug";
        SUSYDebugDirName = OutDirName + "/SUSYDebug";
        outputTerm = 1;
    }
    noMC = noMC_i;
    
#ifdef _MPI
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&procnum);
    MPI_Get_processor_name(processorName,&nameLen);
#else
    rank = 0;
    procnum = 1;
#endif
}

GenerateEvent::~GenerateEvent() 
{
}

void GenerateEvent::generate(int unsigned nIteration_i, int seed)
{
    if (!noMC)
        throw std::runtime_error("\nGenerateEvent::generate():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    nIteration = nIteration_i;
    if (nIteration == 0) outputTerm = 0;
    int it_comp = 0;
    int rem_iteration = nIteration + 1;
    int itno = 0;
#ifdef _MPI
    gRandom->SetSeed(seed * ((unsigned)time(NULL) + rank * procnum + nameLen));
#else
    gRandom->SetSeed(seed);
#endif

    try {
        /* initialize model and set model parameters */
        initModel();
        
        /* define parameter distributions */
        defineParameterDistributions();
        
        /* create directory structure and necessary files */
        if (rank == 0) createDirectories();
        
        /* memory allocation for MPI */
        sendbuff = new double[buffersize];
        sendbuff_int = new int[1];
        if( rank == 0 ){
            buff = new double*[procnum];
            buff[0]=new double[procnum * buffersize];
            for(int i = 1; i < procnum; i++) buff[i] = buff[i - 1] + buffersize;
        } else {
            buff = new double*[1];
            buff[0] = new double[1];
        }
        buff_int = new int*[procnum];
        buff_int[0]=new int[procnum];
        for(int i = 1; i < procnum; i++) buff_int[i] = buff_int[i - 1] + 1;
        
        /* computation of observables */
        while (rem_iteration > 0) {
            itno = rank + nIteration - rem_iteration + 1;
            if (rank < rem_iteration){
                generateRandomEvent(itno);
                for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                    sendbuff[positionID++] = it->computeTheoryValue();
                }
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                    std::vector<Observable> ObsInCGO = it1->getObs();
                    for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                        sendbuff[positionID++] = it2->computeTheoryValue();
                    }
                }
                sendbuff_int[0] = rank + 1;
            } else {
                for (int i = 0; i < buffersize; i++){
                    sendbuff[i] = 0.;
                }
                sendbuff_int[0] = 0;
            }
            
#ifdef _MPI
                MPI::COMM_WORLD.Gather(sendbuff, buffersize, MPI::DOUBLE, buff[0], buffersize, MPI::DOUBLE, 0);
                MPI::COMM_WORLD.Allgather(sendbuff_int, 1, MPI::INT, buff_int[0], 1, MPI::INT);
#else
                buff[0] = sendbuff;
                buff_int[0] = sendbuff_int;
#endif
                
            /* output generation */
            if ( rank == 0 ){
                itno = nIteration - rem_iteration + 1;
                for(int iproc = 0; iproc < procnum; iproc++){
                    positionID = 0;
                    if (buff_int[iproc][0] != 0) {
                        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
                            if (outputTerm == 1) (*ParsOut[it->name]) << buff[iproc][positionID] << std::endl;
                            positionID++;
                        }
                        if (outputTerm == 0 && itno != 0) std::cout << "\nEvent No.: " << itno++ << std::endl;
                        else if (outputTerm == 0 && itno++ == 0) std::cout << "\nThe central event is: " <<std::endl;
                        if (outputTerm == 0 && Obs.size() > 0) std::cout << "\nObservables: \n" << std::endl;
                        for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                            if (outputTerm == 0) {
                                std::cout << it->getName() << " = " << buff[iproc][positionID] << std::endl;
                                positionID++;
                            } else {
                                (*ObsOut[it->getName()]) << buff[iproc][positionID] << std::endl;
                                positionID++;
                            }
                        }
                        if (outputTerm == 0 && CGO.size() > 0) std::cout << "\nCorrelated Gaussian Observables: \n" << std::endl;
                        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                            if (outputTerm == 0) std::cout << it1->getName() << ":" << std::endl;
                            std::vector<Observable> ObsInCGO = it1->getObs();
                            for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                                if (outputTerm == 0) {
                                    std::cout << it2->getName() << " = " << buff[iproc][positionID] << std::endl;
                                    positionID++;
                                } else {
                                    (*CGOOut[it1->getName()]) << buff[iproc][positionID] << "\t";
                                    positionID++;
                                }
                            }
                            if (outputTerm == 0) {
                                std::cout << std::endl;
                            } else {
                                (*CGOOut[it1->getName()]) << std::endl;
                            }
                        }
                    }
                }
            }
                
            /* computation of remaining iterations */
            it_comp = 0;
            for(int iproc = 0; iproc < procnum; iproc++){
                it_comp = std::max(it_comp, buff_int[iproc][0]);
            }
            rem_iteration -= it_comp;
        }
        
        /* closing of output streams */
        if (outputTerm == 1 && rank == 0){
            for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
                ParsOut[it->name]->close();
            }
            for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++){
                ObsOut[it->getName()]->close();
            }
            for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++){
                CGOOut[it1->getName()]->close();
            }
            summary << "Events\t" << nIteration << "\n";
            summary.close();
            std::cout << "Output written to: " << OutDirName << "\n" << std::endl;
        } else {
            if (rank == 0) std::cout << std::endl;
        }
        return;
    } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
    }
}

void GenerateEvent::createDirectories()
{
    if (nIteration > 0 && outputTerm == 1){
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
                std::cout << "\nWARNING: Removed " << OldOutDirName << std::endl;
            }
            gSystem->Exec(("mv -f " + OutDirName + " " + "GeneratedEvents/OLD/").c_str());
            std::cout << "WARNING: " << OutDirName << " exists!\nWARNING: " << OutDirName <<" moved to " << OldOutDirName << "\n"<< std::endl;
        }
        gSystem->MakeDirectory(OutDirName.c_str());
        gSystem->MakeDirectory(ObsDirName.c_str());
        gSystem->MakeDirectory(CGODirName.c_str());
        gSystem->MakeDirectory(ParsDirName.c_str());
        
        if (outputTerm == 1){
            summary.open((OutDirName + "/Summary.txt").c_str(), std::ios::out);
            summary << "Model\t" << ModelName << "\n";
        }
        
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
            ParsOut[it->name] = boost::make_shared<std::ofstream>((ParsDirName + "/" + it->name + ".txt").c_str(), std::ios::out);
            summary << "Parameter\t" << it->name << "\n";
        }
        
        for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++){
            ObsOut[it->getName()] = boost::make_shared<std::ofstream>((ObsDirName + "/" + it->getName() + ".txt").c_str(), std::ios::out);
            summary << "Observable\t" << it->getName() << "\n";
        }
        for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++){
            CGOOut[it->getName()] = boost::make_shared<std::ofstream>((CGODirName + "/" + it->getName() + ".txt").c_str(), std::ios::out);
            summary << "CGO\t" << it->getName() << "\n";
        }
        
        std::cout << "\nRunning in Event Generation mode... \nWARNING: The output should not be used for any statistical analysis. \n         Neither randomness or completness of the sample is gauranteed!!\n" << std::endl;
    } else if(nIteration > 0 && outputTerm == 0) {
        std::cout << "\nRunning in Event Generation mode... \nWARNING: Output being sent to terminal, no data written to disk!! \nWARNING: The output should not be used for any statistical analysis. \n         Neither randomness or completness of the sample is gauranteed!!\n" << std::endl;
    } else {
        std::cout.precision(10);
        std::cout << std::endl << "Running in Single Event mode...\nNo data will be written to disk!\n" << std::endl;
    }
}

void GenerateEvent::initModel(){
    ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO);
    Mod = myInputParser.getModel();
    if (Obs.size() == 0 && CGO.size() == 0 && rank == 0) throw std::runtime_error("\nGenerateEvent::generate(): No observables or correlated Gaussian observables defined in " + ModelConf +" file\n");
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        DP[it->name] = it->ave;
        if (it->errg > 0. || it->errf > 0.){
            ModParsVar.push_back(*it);
        }
    }
    if (!myInputParser.getModel()->Init(DP)) {
        throw std::runtime_error("parameter(s) missing in model initialization");
    }
    
    buffersize = Obs.size() + ModParsVar.size();
    for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++){
        buffersize += it->getObs().size();
    }
}

void GenerateEvent::defineParameterDistributions()
{
    if (outputTerm == 0 && rank == 0) std::cout << "\nParameters varied in Event Generation:" << std::endl;
    for (std::vector<ModelParameter>::const_iterator it = ModParsVar.begin();
            it < ModParsVar.end(); it++) {
        //if (it->errf == 0. && it->errg == 0.)
            //continue;
        if (outputTerm == 0 && rank == 0) std::cout << it->name << ", "; //<< k << std::endl;
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
    positionID = 0;
    if (iterationNo == 0){
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
            DPars[it->name] = it->ave;
            sendbuff[positionID++] = DPars[it->name];
        }
    } else {
        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
            DPars[it->name] = DDist[it->name]->GetRandom();
            sendbuff[positionID++] = DPars[it->name];
        }
    }
    Mod->Update(DPars);
}
