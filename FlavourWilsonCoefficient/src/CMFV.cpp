/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CMFV.h"

const std::string CMFV::CMFVvars[NCMFVvars] = {"Ftt", "CMFV_dX", "CMFV_dY", "CMFV_dZ", "CMFV_dC7", "CMFV_muNP"};

CMFV::CMFV() : StandardModel(), CMFVM(*this) {
    setModelName("CMFV");
    SMM.setObj((StandardModelMatching&) CMFVM.getObj());
    ModelParamMap.insert(std::make_pair("Ftt", std::cref(Ftt)));
    // Delta F = 1 master functions
    ModelParamMap.insert(std::make_pair("CMFV_dX", std::cref(dX)));
    ModelParamMap.insert(std::make_pair("CMFV_dY", std::cref(dY)));
    ModelParamMap.insert(std::make_pair("CMFV_dZ", std::cref(dZ)));
    ModelParamMap.insert(std::make_pair("CMFV_dC7", std::cref(dC7)));
    ModelParamMap.insert(std::make_pair("CMFV_muNP", std::cref(muNP)));
}

///////////////////////////////////////////////////////////////////////////
// Initialization

void CMFV::setParameter(const std::string name, const double& value){
    if(name.compare("Ftt") == 0)
        Ftt = value;
    else if(name.compare("CMFV_dX") == 0)
        dX = value;
    else if(name.compare("CMFV_dY") == 0)
        dY = value;
    else if(name.compare("CMFV_dZ") == 0)
        dZ = value;
    else if(name.compare("CMFV_dC7") == 0)
        dC7 = value;
    else if(name.compare("CMFV_muNP") == 0)
        muNP = value;
    else
        StandardModel::setParameter(name,value);
}

bool CMFV::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NCMFVvars; i++) {
        if (DPars.find(CMFVvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory CMFV parameter " << CMFVvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(CMFVvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}
