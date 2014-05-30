/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPEffective2.h"


const std::string NPEffective2::NPEffectiveVars[NNPEffectiveVars]
        = {"cWB_NP", "cH_NP", "cLL_NP",
    "cHLp_NP", "cHL_NP",
    "cHQ1p_NP", "cHQ2p_NP", "cHQ3p_NP",
    "cHQ1_NP", "cHQ2_NP", "cHQ3_NP",
    "cHU1_NP", "cHU2_NP", "cHU3_NP",
    "cHD1_NP", "cHD2_NP", "cHD3_NP",
    "cHE_NP",
    "Lambda_NP"};

NPEffective2::NPEffective2()
: NPEffective()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cWB_NP", boost::cref(cWB)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cH_NP", boost::cref(cH)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cLL_NP", boost::cref(cL1L1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHLp_NP", boost::cref(cHL1p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHL_NP", boost::cref(cHL1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ1p_NP", boost::cref(cHQ1p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ2p_NP", boost::cref(cHQ2p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ3p_NP", boost::cref(cHQ3p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ1_NP", boost::cref(cHQ1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ2_NP", boost::cref(cHQ2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHQ3_NP", boost::cref(cHQ3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU1_NP", boost::cref(cHU1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU2_NP", boost::cref(cHU2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHU3_NP", boost::cref(cHU3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD1_NP", boost::cref(cHD1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD2_NP", boost::cref(cHD2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHD3_NP", boost::cref(cHD3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("cHE_NP", boost::cref(cHE1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Lambda_NP", boost::cref(LambdaNP)));
}

void NPEffective2::setParameter(const std::string name, const double& value)
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
    } else if (name.compare("cHL_NP") == 0) {
        cHL1 = value;
        cHL2 = value;
        cHL3 = value;
    } else if (name.compare("cHQ1p_NP") == 0) {
        cHQ1p = value;
    } else if (name.compare("cHQ2p_NP") == 0) {
        cHQ2p = value;
    } else if (name.compare("cHQ3p_NP") == 0) {
        cHQ3p = value;
    } else if (name.compare("cHQ1_NP") == 0) {
        cHQ1 = value;
    } else if (name.compare("cHQ2_NP") == 0) {
        cHQ2 = value;
    } else if (name.compare("cHQ3_NP") == 0) {
        cHQ3 = value;
    } else if (name.compare("cHU1_NP") == 0) {
        cHU1 = value;
    } else if (name.compare("cHU2_NP") == 0) {
        cHU2 = value;
    } else if (name.compare("cHU3_NP") == 0) {
        cHU3 = value;
    } else if (name.compare("cHD1_NP") == 0) {
        cHD1 = value;
    } else if (name.compare("cHD2_NP") == 0) {
        cHD2 = value;
    } else if (name.compare("cHD3_NP") == 0) {
        cHD3 = value;
    } else if (name.compare("cHE_NP") == 0) {
        cHE1 = value;
        cHE2 = value;
        cHE3 = value;
    } else if (name.compare("Lambda_NP") == 0)
        LambdaNP = value;
    else
        NPEffective::setParameter(name, value);
}

bool NPEffective2::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective2 parameter "
                    << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return (NPEffective::CheckParameters(DPars));
}


