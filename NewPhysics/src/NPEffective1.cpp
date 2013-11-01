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


bool NPEffective1::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective1 parameter"
                      << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}


void NPEffective1::setParameters(const std::string name, const double& value)
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
        StandardModel::setParameters(name, value);
}





