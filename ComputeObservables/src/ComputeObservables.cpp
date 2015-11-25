/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Observable.h"
#include "Observable2D.h"
#include "CorrelatedGaussianObservables.h"
#include "ComputeObservables.h"
#include <boost/algorithm/string/predicate.hpp>

ComputeObservables::ComputeObservables(ModelFactory& ModelF, ThObsFactory& ThObsF, 
        const std::string& ModelConf_i, const int rank_i)
: myInputParser(ModelF, ThObsF),
  rank(rank_i)
{
    std::vector<ModelParameter> ModPars;
    boost::ptr_vector<Observable> Obs;
    std::vector<Observable2D> Obs2D;
    std::vector<CorrelatedGaussianObservables> CGO;
    std::vector<CorrelatedGaussianParameters> CGP;

    std::string ModelName = myInputParser.ReadParameters(ModelConf_i, rank, ModPars, Obs, Obs2D, CGO, CGP);
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        DP[it->getname()] = it->getave();
        paraNames.push_back(it->getname());
    }
    DPars = DP;
    
    for (boost::ptr_vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
        DObs.insert(std::pair<std::string, double> (it->getThname() , 0.));
        DThObs.insert(std::pair<std::string, ThObservable*> (it->getThname() , it->getTho()));
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
        std::vector<Observable> ObsInCGO = it1->getObs();
        for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
            DObs.insert(std::pair<std::string, double> (it2->getThname() , 0.));
            DThObs.insert(std::pair<std::string, ThObservable*> (it2->getThname() , it2->getTho()));
        }
    }
    
    Mod = myInputParser.getModel();
    if (!Mod->Init(DP)) 
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
}

ComputeObservables::ComputeObservables(ModelFactory& ModelF, ThObsFactory& ThObsF,
        const std::string& ModelName_i, std::map<std::string, double> DPars_i, const int rank_i)
: ModelName(ModelName_i),
  myInputParser(ModelF, ThObsF),
  rank(rank_i)
{
    for (std::map<std::string, double>::iterator it = DPars_i.begin(); it != DPars_i.end(); it++) {
        paraNames.push_back(it->first);
    }
    Mod = myInputParser.getModelFactory().CreateModel(ModelName);
    Mod->setModelName(ModelName);
    Mod->InitializeModel();
    if (Mod->IsModelInitialized()) {
        if (rank == 0) std::cout << "\nModel Initialized: " << ModelName << std::endl;
    } else {
        throw std::runtime_error("\nERROR: " + ModelName + " not initialized successfully.\n");
    }
    if (!Mod->Init(DPars_i))
        throw std::runtime_error("\nERROR: Model cannot be initialized initialization.\n");
}

ComputeObservables::~ComputeObservables()
{
}

void ComputeObservables::setFlags(std::map<std::string, std::string> DFlags_i)
{
    for (std::map<std::string, std::string>::iterator it = DFlags_i.begin(); it != DFlags_i.end(); it++) {
        if (!(boost::iequals(it->second, "true")) && !(boost::iequals(it->second, "false")) && !Mod->setFlagStr(it->first, it->second)) {
            throw std::runtime_error("ERROR: setFlagStr error for " + it->first);
        } else if (boost::iequals(it->second, "true") && !Mod->setFlag(it->first, 1)) {
            throw std::runtime_error("ERROR: setFlag error for " + it->first);
        } else if (boost::iequals(it->second, "false") && !Mod->setFlag(it->first, 0)) {
            throw std::runtime_error("ERROR: setFlag error for " + it->first);
        } else {
            if (rank == 0) std::cout << "set flag " << it->first << " = " << it->second << std::endl;
        }
    }
}

std::map<std::string, double> ComputeObservables::compute(std::map<std::string, double> DP)
{
    if (DP != DPars)
        for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++) {
            if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end()))
                throw std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
        }
    DPars = DP;
    Mod->Update(DPars);
    for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
        DObs[it->first] = DThObs[it->first]->computeThValue();
    }
    return(DObs);
}

void ComputeObservables::RemoveObservable(std::string ObsName)
{
    if(DObs.find(ObsName) == DObs.end()) 
        throw std::runtime_error("\nERROR: Observable cannot be removed since it has not been added.\n");
    DObs.erase(ObsName);
    if(DThObs.find(ObsName) == DThObs.end()) 
        throw std::runtime_error("\nERROR: ThObservable cannot be removed since it has not been added.\n");
    DThObs.erase(ObsName);
}

void ComputeObservables::AddObservable(std::string ObsName)
{
    DThObs[ObsName] = myInputParser.getObsFactory().CreateThMethod(ObsName, *Mod);
    DObs.insert(std::pair<std::string, double> (ObsName , 0.));
}
 
void ComputeObservables::addCustomObservableType(const std::string name, boost::function<Observable*() > funct){
        myInputParser.addCustomObservableType(name, funct);
}
