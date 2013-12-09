/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ComputeObservables.h"
#include <TSystem.h>
#include <TMath.h>
#include <TRandom3.h>
#ifdef _MPI
#include <mpi.h>
#endif

ComputeObservables::ComputeObservables(const std::string& ModelConf_i,
                                       std::map<std::string, double> DObs_i)
: DObs(DObs_i), 
  myInputParser(), 
  ModelConf(ModelConf_i)
{
    std::string ModelName = myInputParser.ReadParameters(ModelConf, ModPars, Obs, Obs2D, CGO, ParaObs);
    std::map<std::string, double> DP;
    for (std::vector<ModelParameter>::iterator it = ModPars.begin(); it < ModPars.end(); it++) {
         DP[it->name] = it->ave;
         paraNames.push_back(it->name);
    }
    Mod = myInputParser.getMyModel();
    if (!Mod->Init(DP)) {
        throw std::runtime_error("parameter(s) missing in model initialization");
    }
    checkPara = 0;
}

ComputeObservables::~ComputeObservables() 
{
}

std::map<std::string, double> ComputeObservables::compute(std::map<std::string, double> DPars)
{
    if (checkPara == 0)
    for (std::map<std::string, double>::iterator it = DPars.begin(); it != DPars.end(); it++){
        if(!(std::find(paraNames.begin(), paraNames.end(), it->first) != paraNames.end())) 
              throw  std::runtime_error("\nERROR: Incorrect parameter name passed to ComputeObservable");
    checkPara = 1;
    }
    
    try {
        Mod->Update(DPars);
        for (std::map<std::string, double>::iterator it = DObs.begin(); it != DObs.end(); it++) {
             DObs[it->first] = myInputParser.getMyThFactory()->getThMethod(it->first)->computeThValue();
        }
        return(DObs);
        } catch (std::string message) {
        std::cerr << message << std::endl;
        exit(EXIT_FAILURE);
        }
}
