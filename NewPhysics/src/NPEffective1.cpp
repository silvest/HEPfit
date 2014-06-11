/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective1.h"


const std::string NPEffective1::NPEffectiveVars[NNPEffectiveVars]
        = {"cWB_NP", "cH_NP", "cLL_NP", "cHLp_NP", "cHQp_NP",
    "cHL_NP", "cHQ_NP", "cHE_NP", "cHU_NP", "cHD_NP", "Lambda_NP"};

NPEffective1::NPEffective1()
: NPEffective()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cWB_NP", boost::cref(cWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cH_NP", boost::cref(cH)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cLL_NP", boost::cref(cL1L1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHLp_NP", boost::cref(cHL1p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQp_NP", boost::cref(cHQ1p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL_NP", boost::cref(cHL1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ_NP", boost::cref(cHQ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE_NP", boost::cref(cHE1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU_NP", boost::cref(cHU1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD_NP", boost::cref(cHD1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(LambdaNP)));
}

void NPEffective1::setParameter(const std::string name, const double& value)
{
    if (name.compare("cWB_NP") == 0)
        cWB = value;
    else if (name.compare("cH_NP") == 0)
        cH = value;
    else if (name.compare("cLL_NP") == 0) {
        cL1L1 = value;
        cL1L2 = value;
        cL1L3 = value;
        cL2L2 = value;
        cL2L3 = value;
        cL3L3 = value;
    } else if (name.compare("cHLp_NP") == 0) {
        cHL1p = value;
        cHL2p = value;
        cHL3p = value;
    } else if (name.compare("cHQp_NP") == 0) {
        cHQ1p = value;
        cHQ2p = value;
        cHQ3p = value;
    } else if (name.compare("cHL_NP") == 0) {
        cHL1 = value;
        cHL2 = value;
        cHL3 = value;
    } else if (name.compare("cHQ_NP") == 0) {
        cHQ1 = value;
        cHQ2 = value;
        cHQ3 = value;
    } else if (name.compare("cHE_NP") == 0) {
        cHE1 = value;
        cHE2 = value;
        cHE3 = value;
    } else if (name.compare("cHU_NP") == 0) {
        cHU1 = value;
        cHU2 = value;
        cHU3 = value;
    } else if (name.compare("cHD_NP") == 0) {
        cHD1 = value;
        cHD2 = value;
        cHD3 = value;
    } else if (name.compare("Lambda_NP") == 0)
        LambdaNP = value;
    else
        NPEffective::setParameter(name, value);
}

bool NPEffective1::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective1 parameter "
                    << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return (NPEffective::CheckParameters(DPars));
}


