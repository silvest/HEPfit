/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CMFV.h"

const std::string CMFV::CMFVvars[NCMFVvars] = {"Ftt"};

CMFV::CMFV() : StandardModel(), CMFVM(*this) {   

    SMM.setObj((StandardModelMatching&) CMFVM.getObj());
    ModelParamMap.insert(std::make_pair("Ftt", std::cref(Ftt)));
}

///////////////////////////////////////////////////////////////////////////
// Initialization

void CMFV::setParameter(const std::string name, const double& value){    
    if(name.compare("Ftt") == 0) 
        Ftt = value;
    else
        StandardModel::setParameter(name,value);
}

bool CMFV::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NCMFVvars; i++) {
        if (DPars.find(CMFVvars[i]) == DPars.end()) {
            std::cout << "missing mandatory CMFV parameter " << CMFVvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}
