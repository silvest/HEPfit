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
    if (Obs.size() == 0 && CGO.size() == 0) throw std::runtime_error("\nMonteCarlo::TestRun(): No observables or correlated Gaussian observables defined in SomeModel.conf file\n");
    
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
        DP[it->name] = it->ave;
        paraNames.push_back(it->name);
    }
    
    Mod = myInputParser.getMyModel();
    if (!Mod->Init(DP)) {
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
    }
    checkPara = 0;
    constructorIdx = 1;
}

ComputeObservables::ComputeObservables(const std::string& ModelConf_i,
                                       std::map<std::string, double> DObs_i)
: myInputParser(),
  DObs(DObs_i),
  ModelConf(ModelConf_i)
{
    std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
         DP[it->name] = it->ave;
         paraNames.push_back(it->name);
    }
    thf = myInputParser.getMyThFactory();
    Mod = myInputParser.getMyModel();
    if (!Mod->Init(DP)) {
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
    }
    checkPara = 0;
    constructorIdx = 2;
}

ComputeObservables::ComputeObservables(const std::string& ModelName_i,
                                       std::map<std::string, double> DPars_i,
                                       std::map<std::string, double> DObs_i)
: ModelName(ModelName_i),
  myInputParser(),
  DPars(DPars_i),
  DObs(DObs_i)
{
    for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) {
        paraNames.push_back(it->first);
    }
    Mod = myInputParser.ModelFactory(ModelName);
    Mod->setModelName(ModelName);
    Mod->InitializeModel();
    thf = new ThFactory(*Mod);
    if (!Mod->Init(DPars)) {
        throw std::runtime_error("ERROR: Parameter(s) missing in model initialization.\n");
    }
    checkPara = 0;
    constructorIdx = 3;
}

ComputeObservables::~ComputeObservables()
{
}

void ComputeObservables::setFlags(std::map<std::string, std::string> DFlags_i){
    
    if (constructorIdx != 3) {
        throw std::runtime_error("ERROR: ComputeObservables::setFlags() can only be used with ComputeObservables(ModelName, DPars, DObs) construction.\n");
    }
    
    DFlags = DFlags_i;
    for (std::map<std::string, std::string>::iterator it = DFlags.begin(); it != DFlags.end(); it++) {
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
    if (constructorIdx == 1) {
        throw std::runtime_error("ERROR: ComputeObservables::compute() can only be used with ComputeObservables(ModelName, DPars, DObs) or ComputeObservables(ModelConf, DObs) construction.\n");
    }
    if (checkPara == 0)
        for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++) {
            if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end()))
                throw  std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
            checkPara = 1;
        }
    
    Mod->Update(DP);
    for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
        DObs[it->first] = thf->getThMethod(it->first)->computeThValue();
    }
    return(DObs);
}

std::pair<std::map<std::string, double>, std::map<std::string, std::map<std::string, double> > > ComputeObservables::compute_Obs_CGO(std::map<std::string, double> DP){
    if (constructorIdx != 1) {
        throw std::runtime_error("ERROR: ComputeObservables::compute_Obs_CGO() can only be used with ComputeObservables(ModelConf) construction.\n");
    }
    if (checkPara == 0)
        for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++) {
            if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end()))
                throw  std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
            checkPara = 1;
        }
    
    Mod->Update(DP);
    std::map<std::string, double> DObs_temp;
    std::map<std::string, std::map<std::string, double> > DCGO_temp;
    
    for (std::vector<Observable>::iterator it = Obs.begin(); it < Obs.end(); it++) {
        DObs_temp[it->getName()] = it->computeTheoryValue();
    }
    for (std::vector<CorrelatedGaussianObservables>::iterator it1 = CGO.begin(); it1 < CGO.end(); it1++) {
        std::vector<Observable> ObsInCGO = it1->getObs();
        std::map<std::string, double> DObs_CGO_temp;
        for (std::vector<Observable>::iterator it2 = ObsInCGO.begin(); it2 < ObsInCGO.end(); it2++) {
            DObs_CGO_temp[it2->getName()] = it2->computeTheoryValue();
        }
        DCGO_temp[it1->getName()] = DObs_CGO_temp;
    }
    std::cout << std::endl;
    
    return (std::make_pair(DObs_temp, DCGO_temp));
}
