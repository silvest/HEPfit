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

HiggsKvKf::HiggsKvKf() : HiggsExtensionModel()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kv", boost::cref(Kv)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Kf", boost::cref(Kf)));
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


