/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ComputeObservables.h"

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
                                       std::map<std::string, double> DObs_i,
                                       std::map<std::string, std::string> DFlags_i)
:ModelName(ModelName_i),
DPars(DPars_i),
DObs(DObs_i),
DFlags(DFlags_i)
{
    for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++) {
        paraNames.push_back(it->first);
    }
    Mod = ModelDictionary();
    Mod->InitializeModel();
    thf = new ThFactory(*Mod);
    if (!Mod->Init(DPars)) {
        throw std::runtime_error("parameter(s) missing in model initialization");
    }
    checkPara = 0;
    
    for (std::map<std::string, std::string>::iterator it = DFlags.begin(); it != DFlags.end(); it++) {
            if (it->second.compare("true") == 0 && !Mod->setFlag(it->first, 1))
                    throw std::runtime_error("ERROR: setFlag error for " + it->first);
            if (it->second.compare("false") == 0 && !Mod->setFlag(it->first, 0))
                    throw std::runtime_error("ERROR: setFlag error for " + it->first);
            if (it->second.compare("true") != 0 && it->second.compare("false") != 0 && !Mod->setFlagStr(it->first, it->second))
                throw std::runtime_error("ERROR: setFlag error for " + it->first);
            else
                std::cout << "set flag " << it->first << "=" << it->second << std::endl;
        }
}

ComputeObservables::~ComputeObservables()
{
}

std::map<std::string, double> ComputeObservables::compute(std::map<std::string, double> DP)
{
    if (checkPara == 0)
    for (std::map<std::string, double>::iterator it = DP.begin(); it != DP.end(); it++){
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

StandardModel* ComputeObservables::ModelDictionary(){
    std::map<std::string, StandardModel* > DModel;
    
    DModel["StandardModel"] = new StandardModel();
    DModel["NPSTU"] = new NPSTU();
    DModel["NPSTUVWXY"] = new NPSTUVWXY();
    DModel["NPEpsilons"] =  new NPEpsilons();
    DModel["NPEpsilons_pureNP"] = new NPEpsilons_pureNP();
    DModel["NPHiggsST"] = new NPHiggsST();
    DModel["NPZbbbar"] = new NPZbbbar();
    DModel["NPEffective1"] = new NPEffective1();
    DModel["NPEffective2"] = new NPEffective2();
    DModel["MFV"] = new MFV();
    DModel["GeneralSUSY"] =  new GeneralSUSY();
    DModel["pMSSM"] = new pMSSM();
    DModel["SusyMI"] = new SUSYMassInsertion();
    DModel["THDM"] = new THDM();
    
    return DModel[ModelName];
}
