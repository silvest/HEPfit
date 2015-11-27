/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GenerateEvent.h"
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
#include <iomanip>  

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
        OldOutDirName = ("GeneratedEvents/DELETED/" + OutDirName_i + JobTag);
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

void GenerateEvent::generate(int unsigned nIteration_i, int seed, bool weight)
{
    if (!noMC)
        throw std::runtime_error("\nGenerateEvent::generate():\n--noMC was not specified as an argument.\nPlease do so for running in Generate Event mode.\n");
    nIteration = nIteration_i;
    if (nIteration == 0) outputTerm = 0;
    int it_comp = 0;
    int rem_iteration = nIteration + 1;
    int itno = 0;
#ifdef _MPI
    gRandom->SetSeed(seed * (/*(unsigned)time(NULL)*/ seed + rank * procnum + nameLen));
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
        sendbuff_w = new double[buffersize];
        sendbuff_int = new int[1];
        if( rank == 0 ){
            buff = new double*[procnum];
            buff[0]=new double[procnum * buffersize];
            for(int i = 1; i < procnum; i++) buff[i] = buff[i - 1] + buffersize;
            
            buff_w = new double*[procnum];
            buff_w[0]=new double[procnum * buffersize];
            for(int i = 1; i < procnum; i++) buff_w[i] = buff_w[i - 1] + buffersize;
        } else {
            buff = new double*[1];
            buff[0] = new double[1];
            buff_w = new double*[1];
            buff_w[0] = new double[1];
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
                    sendbuff[positionID] = it->computeTheoryValue();
                    if (weight && it->getDistr().compare("noweight") != 0) sendbuff_w[positionID] = it->computeWeight();
                    else sendbuff_w[positionID] = 0.;
                    positionID++;
                }
                for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                    std::vector<Observable> ObsInCGO = it1->getObs();
                    for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                        sendbuff[positionID] = it2->computeTheoryValue();
                        if (weight && it2->getDistr().compare("noweight") != 0) sendbuff_w[positionID] = it1->computeWeight();
                        else sendbuff_w[positionID] = 0.;
                        positionID++;
                    }
                }
                sendbuff_int[0] = rank + 1;
            } else {
                for (int i = 0; i < buffersize; i++){
                    sendbuff[i] = 0.;
                    sendbuff_w[i] = 0.;
                }
                sendbuff_int[0] = 0;
            }
            
#ifdef _MPI
                MPI::COMM_WORLD.Gather(sendbuff, buffersize, MPI::DOUBLE, buff[0], buffersize, MPI::DOUBLE, 0);
                if (weight) MPI::COMM_WORLD.Gather(sendbuff_w, buffersize, MPI::DOUBLE, buff_w[0], buffersize, MPI::DOUBLE, 0);
                MPI::COMM_WORLD.Allgather(sendbuff_int, 1, MPI::INT, buff_int[0], 1, MPI::INT);
#else
                buff[0] = sendbuff;
                buff_w[0] = sendbuff_w;
                buff_int[0] = sendbuff_int;
#endif
                
            /* output generation */
            if ( rank == 0 ){
                itno = nIteration - rem_iteration + 1;
                for(int iproc = 0; iproc < procnum; iproc++){
                    positionID = 0;
                    if (buff_int[iproc][0] != 0) {
                        for (std::vector<ModelParameter>::iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
                            if (outputTerm == 1) (*ParsOut[it->getname()]) << buff[iproc][positionID] << std::endl;
                            positionID++;
                        }
                        if (outputTerm == 0 && itno != 0) std::cout << "\nEvent No.: " << itno++ << std::endl;
                        else if (outputTerm == 0 && itno++ == 0) std::cout << "\nThe central event is: " <<std::endl;
                        if (outputTerm == 0 && Obs.size() > 0) std::cout << "\nObservables: \n" << std::endl;
                        for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
                            if (outputTerm == 0) {
                                if (weight && it->getDistr().compare("noweight") != 0) std::cout << it->getName() << " = " << buff[iproc][positionID] << " (weight: " << buff_w[iproc][positionID] << ")"<< std::endl;
                                else std::cout << it->getName() << " = " << buff[iproc][positionID] << std::endl;
                                positionID++;
                            } else {
                                if (weight && it->getDistr().compare("noweight") != 0) (*ObsOut[it->getName()]) << buff[iproc][positionID] << "\t" << buff_w[iproc][positionID] <<std::endl;
                                else (*ObsOut[it->getName()]) << buff[iproc][positionID] << std::endl;
                                positionID++;
                            }
                        }
                        if (outputTerm == 0 && CGO.size() > 0) std::cout << "\nCorrelated Gaussian Observables: \n" << std::endl;
                        for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
                            if (outputTerm == 0) {
                                if (weight) std::cout << it1->getName() << ": (weight: " << buff_w[iproc][positionID] << ")" <<std::endl;
                                else std::cout << it1->getName() <<std::endl;
                            }
                            std::vector<Observable> ObsInCGO = it1->getObs();
                            for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
                                if (outputTerm == 0) {
                                    if (weight && it2->getDistr().compare("noweight") != 0) std::cout << it2->getName() << " = " << buff[iproc][positionID] << std::endl;
                                    else std::cout << it2->getName() << " = " << buff[iproc][positionID] << std::endl;
                                    positionID++;
                                } else {
                                    if (weight && it2->getDistr().compare("noweight") != 0) (*CGOOut[it1->getName()]) << buff[iproc][positionID] << "\t";
                                    else (*CGOOut[it1->getName()]) << buff[iproc][positionID] << "\t";
                                    positionID++;
                                }
                            }
                            if (outputTerm == 0) std::cout << std::endl;
                            else (*CGOOut[it1->getName()]) << buff_w[iproc][positionID - 1] << std::endl;
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
                ParsOut[it->getname()]->close();
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
    delete buff[0];
    delete [] buff;
    delete [] sendbuff;
    
    delete buff_w[0];
    delete [] buff_w;
    delete [] sendbuff_w;
    
    delete buff_int[0];
    delete [] buff_int;
    delete [] sendbuff_int;
}

void GenerateEvent::createDirectories()
{
    if (nIteration > 0 && outputTerm == 1){
        FileStat_t info;
        if (gSystem->GetPathInfo("GeneratedEvents", info) != 0) {
            gSystem->MakeDirectory("GeneratedEvents");
        }
        if (gSystem->GetPathInfo("GeneratedEvents/DELETED", info) != 0) {
            gSystem->MakeDirectory("GeneratedEvents/DELETED");
        }
        if (gSystem->GetPathInfo(OutDirName.c_str(), info) == 0){
            if (gSystem->GetPathInfo(OldOutDirName.c_str(), info) == 0){
                gSystem->Exec(("rm -rf " + OldOutDirName).c_str());
                std::cout << "\nWARNING: Removed " << OldOutDirName << std::endl;
            }
            gSystem->Exec(("mv -f " + OutDirName + " " + "GeneratedEvents/DELETED/").c_str());
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
            ParsOut[it->getname()] = boost::make_shared<std::ofstream>((ParsDirName + "/" + it->getname() + ".txt").c_str(), std::ios::out);
            summary << "Parameter\t" << it->getname() << "\n";
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
    ModelName = myInputParser.ReadParameters(ModelConf, rank, ModPars, Obs, Obs2D, CGO, CGP);
    Mod = myInputParser.getModel();
    if (Obs.size() == 0 && CGO.size() == 0) {
        if (rank == 0) throw std::runtime_error("\nGenerateEvent::generate(): No observables or correlated Gaussian observables defined in " + ModelConf +" file\n");
        else sleep (2);
    }
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        if (it->IsCorrelated()) {
                for (int i = 0; i < CGP.size(); i++) {
                    if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                        std::string index = it->getname().substr(CGP[i].getName().size());
                        long int lindex = strtol(index.c_str(), NULL, 10);
                        if (lindex > 0) DP[CGP[i].getPar(lindex - 1).getname()] = CGP[i].getPar(lindex - 1).getave();
                        else {
                            std::stringstream out;
                            out << it->getname();
                            throw std::runtime_error("MonteCarlo::Run(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                        }
                    }
                }
            } else DP[it->getname()] = it->getave();
        
        if (it->geterrg() > 0. || it->geterrf() > 0.) ModParsVar.push_back(*it);
    }
    if (!myInputParser.getModel()->Init(DP)) {
        if (rank == 0) throw std::runtime_error("\nERROR: Model cannot be initialization");
        else sleep(2);
    }
    
    buffersize = Obs.size() + ModParsVar.size();
    for (std::vector<CorrelatedGaussianObservables>::iterator it = CGO.begin(); it < CGO.end(); it++){
        buffersize += it->getObs().size();
    }
}

void GenerateEvent::defineParameterDistributions()
{
    if (outputTerm == 0 && rank == 0) std::cout << "\nParameters varied in Event Generation:" << std::endl;
    for (std::vector<ModelParameter>::const_iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++) {
        if (it->IsCorrelated()) {
                for (int i = 0; i < CGP.size(); i++) {
                    if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                        std::string index = it->getname().substr(CGP[i].getName().size());
                        long int lindex = strtol(index.c_str(), NULL, 10);
                        if (lindex > 0) {
                            if (outputTerm == 0 && rank == 0) std::cout << CGP[i].getPar(lindex - 1).getname() << ", ";
                        } else {
                            std::stringstream out;
                            out << it->getname();
                            throw std::runtime_error("MonteCarlo::Run(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                        }
                    }
                }
            } else if (outputTerm == 0 && rank == 0) std::cout << it->getname() << ", "; //<< k << std::endl;
        if (it->geterrf() == 0. && it->geterrg() != 0.){
            DDist[it->getname()] = new TF1(it->getname().c_str(),
            "1./sqrt(2.*TMath::Pi())/[1] * exp(-(x-[0])*(x-[0])/2./[1]/[1])",
            it->getmin(),it->getmax());
            DDist[it->getname()]->SetParameter(0, it->getave());
            DDist[it->getname()]->SetParameter(1, it->geterrg());
        }
        else if (it->geterrg() == 0. && it->geterrf() != 0.){
            DDist[it->getname()] = new TF1(it->getname().c_str(),
                    "1/([1]-[0])",
                    it->getmin(), it->getmax());
            DDist[it->getname()]->SetParameter(0,it->getmin());
            DDist[it->getname()]->SetParameter(1,it->getmax());
        }
        else {
            DDist[it->getname()] = new TF1(it->getname().c_str(),
                    "(TMath::Erf((x-[0]+[2])/sqrt(2.)/[1])-TMath::Erf((x-[0]-[2])/sqrt(2.)/[1]))/4./[2]",
                    it->getmin(), it->getmax());
            DDist[it->getname()]->SetParameter(0, it->getave());
            DDist[it->getname()]->SetParameter(1, it->geterrg());
            DDist[it->getname()]->SetParameter(2, it->geterrf());
        }
    }
    if (outputTerm == 0) std::cout << std::endl;
}

void GenerateEvent::generateRandomEvent(int iterationNo)
{
    positionID = 0;
    
    std::vector<double> vec(ModParsVar.size(),0.);
    for (int i=0; i< vec.size();i++){
        if (iterationNo == 0) vec[i] = ModParsVar[i].getave();
        else vec[i] = DDist[ModParsVar[i].getname()]->GetRandom();
    }
    
    setDParsFromParameters(vec,DPars);
    for (std::vector<ModelParameter>::const_iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++){
        if (it->IsCorrelated()) {
                for (int i = 0; i < CGP.size(); i++) {
                    if (CGP[i].getName().compare(it->getCgp_name()) == 0) {
                        std::string index = it->getname().substr(CGP[i].getName().size());
                        long int lindex = strtol(index.c_str(), NULL, 10);
                        if (lindex > 0) {
                            sendbuff[positionID] = DPars[CGP[i].getPar(lindex - 1).getname()];
                            sendbuff_w[positionID] = DPars[CGP[i].getPar(lindex - 1).getname()];
                            positionID++;
                        }
                        else {
                            std::stringstream out;
                            out << it->getname();
                            throw std::runtime_error("MonteCarlo::Run(): " + out.str() + "seems to be part of a CorrelatedGaussianParameters object, but I couldn't find the corresponding object");
                        }
                    }
                }
            } else {
                sendbuff[positionID] = DPars[it->getname()];  
                sendbuff_w[positionID] = DPars[it->getname()];
                positionID++;
            }
    }
    Mod->Update(DPars);
}

void GenerateEvent::setDParsFromParameters(const std::vector<double>& parameters, std::map<std::string,double>& DPars_i) 
{
    std::map<std::string, std::vector<double> > cgpmap;

   unsigned int k = 0;
    for (std::vector<ModelParameter>::const_iterator it = ModParsVar.begin(); it < ModParsVar.end(); it++){
        
        if (it->IsCorrelated()) {
            std::string index = it->getname().substr(it->getCgp_name().size());
            long int lindex = strtol(index.c_str(),NULL,10);
            if (lindex - 1 == cgpmap[it->getCgp_name()].size()) cgpmap[it->getCgp_name()].push_back(parameters[k]);
            else {
                std::stringstream out;
                out << it->getname() << " " << lindex;
                throw std::runtime_error("GenerateEvent::setDParsFromParameters(): " + out.str() + "seems to be a CorrelatedGaussianParameters object but the corresponding parameters are missing or not in the right order");
            }

        } else DPars_i[it->getname()] = parameters[k];
        k++;
    }
    
    for (int j = 0; j < CGP.size(); j++) {
        std::vector<double> current = cgpmap.at(CGP[j].getName());
        if (current.size() != CGP[j].getPars().size()) {
            std::stringstream out;
            out << CGP[j].getName();
            throw std::runtime_error("GenerateEvent::setDParsFromParameters(): " + out.str() + " appears to be represented in cgpmap with a wrong size");
        }
        
        std::vector<double> porig = CGP[j].getOrigParsValue(current);

        for(int l = 0; l < porig.size(); l++) {
            DPars_i[CGP[j].getPar(l).getname()] = porig[l];
        }
    }
}

void GenerateEvent::addCustomObservableType(const std::string name, boost::function<Observable*() > funct){
        myInputParser.addCustomObservableType(name, funct);
}
