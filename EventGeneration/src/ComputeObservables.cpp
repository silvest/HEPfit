/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ComputeObservables.h"
#include <boost/algorithm/string/predicate.hpp>

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
        throw std::runtime_error("parameter(s) missing in model initialization");
    }
    checkPara = 0;
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
        throw std::runtime_error("parameter(s) missing in model initialization");
    }
    checkPara = 0;
}

ComputeObservables::~ComputeObservables()
{
}

void ComputeObservables::setFlags(std::map<std::string, std::string> DFlags_i){
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
    if (checkPara == 0)
        for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++) {
            if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end()))
                throw  std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
            checkPara = 1;
        }
    
    try {
        Mod->Update(DP);
        for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
            DObs[it->first] = thf->getThMethod(it->first)->computeThValue();
        }
        return(DObs);
    } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
    }
}
