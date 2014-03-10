/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ComputeObservables.h"
#include <boost/algorithm/string/predicate.hpp>

ComputeObservables::ComputeObservables(const std::string& ModelConf_i)
: myInputParser(),
ModelConf(ModelConf_i)
{
    std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        DP[it->name] = it->ave;
        paraNames.push_back(it->name);
    }
    
    for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
        DObs.insert(std::pair<std::string, double> (it->getThname() , 0.));
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
        std::vector<Observable> ObsInCGO = it1->getObs();
        for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
            DObs.insert(std::pair<std::string, double> (it2->getThname() , 0.));
        }
    }
    
    thf = myInputParser.getMyThFactory();
    Mod = myInputParser.getMyModel();
    if (!Mod->Init(DP)) {
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
    }
}

ComputeObservables::ComputeObservables(const std::string& ModelName_i, std::map<std::string, double> DPars_i)
: ModelName(ModelName_i),
myInputParser()
{
    for (std::map<std::string, double>::iterator it = DPars_i.begin(); it != DPars_i.end(); it++) {
        paraNames.push_back(it->first);
    }
    Mod = myInputParser.ModelFactory(ModelName);
    Mod->setModelName(ModelName);
    Mod->InitializeModel();
    thf = new ThFactory(*Mod);
    if (!Mod->Init(DPars_i)) {
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
    }
}

ComputeObservables::~ComputeObservables()
{
}

void ComputeObservables::setFlags(std::map<std::string, std::string> DFlags_i)
{
    for (std::map<std::string, std::string>::iterator it = DFlags_i.begin(); it != DFlags_i.end(); it++) {
        if (!(boost::iequals(it->second, "true")) && !(boost::iequals(it->second, "false")) && !Mod->setFlagStr(it->first, it->second))
            throw std::runtime_error("ERROR: setFlagStr error for " + it->first);
        else if (boost::iequals(it->second, "true") && !Mod->setFlag(it->first, 1))
            throw std::runtime_error("ERROR: setFlag error for " + it->first);
        else if (boost::iequals(it->second, "false") && !Mod->setFlag(it->first, 0))
            throw std::runtime_error("ERROR: setFlag error for " + it->first);
        else
            std::cout << "set flag " << it->first << "=" << it->second << std::endl;
    }
}

std::map<std::string, double> ComputeObservables::compute(std::map<std::string, double> DP)
{
    if (DP != DPars)
        for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++) {
            if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end()))
                throw  std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
        }
    DPars = DP;
    Mod->Update(DPars);
    for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
        DObs[it->first] = thf->getThMethod(it->first)->computeThValue();
    }
    return(DObs);
}

void ComputeObservables::RemoveObservables(std::string ObsName)
{
    DObs.erase(ObsName);
}

void ComputeObservables::AddObservables(std::string ObsName)
{
    DObs.insert(std::pair<std::string, double> (ObsName , 0.));
}