/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsKvKf.h"

const std::string HiggsKvKf::HKvKfvars[NHKvKfvars] = {
    "Kv", "Kf"
};

// For the moment we just call the SM initializer, but this must be changed once we implement the EW class with Kv and Kf
bool HiggsKvKf::InitializeModel()
{
    return(StandardModel::InitializeModel());
}

void HiggsKvKf::setParameter(const std::string name, const double& value)
{
    if (name.compare("Kv") == 0)
        Kv = value;
    else if (name.compare("Kf") == 0)
        Kf = value;    
    else
        StandardModel::setParameter(name, value);
}

bool HiggsKvKf::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NHKvKfvars; i++) {
        if (DPars.find(HKvKfvars[i]) == DPars.end()) {
            std::cout << "missing mandatory HiggsKvKf parameter " << HKvKfvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}


