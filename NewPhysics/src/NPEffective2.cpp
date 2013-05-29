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
   "cHLp_NP", "cHL_NP", "cHE_NP",
   "cHQ2pMINUScHQ2_NP",
   "cHQ3pPLUScHQ3_NP",
   "cHU2_NP", "cHD3_NP",
   "l_NP", 
   "Lambda_NP"};


bool NPEffective2::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NNPEffectiveVars; i++) {
        if (DPars.find(NPEffectiveVars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEffective2 parameter"
                      << NPEffectiveVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}


void NPEffective2::SetParameter(const std::string name, const double& value)
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
    } else if (name.compare("cHE_NP") == 0) {
        cHE1 = value;
        cHE2 = value;
        cHE3 = value;
    } else if (name.compare("cHQ2pMINUScHQ2_NP") == 0) {
        /* cHQ2p + cHQ2p = 0 */
        cHQ2p = value/2.0;
        cHQ2 = - cHQ2p;
    } else if (name.compare("cHQ3pPLUScHQ3_NP") == 0) {
        /* cHQ3p - cHQ3p = 0 */
        cHQ3p = value/2.0;
        cHQ3 = cHQ3p;
    } else if (name.compare("cHU2_NP") == 0) {
        cHU2 = value;
    } else if (name.compare("cHD3_NP") == 0) {
        cHD3 = value;
    } else if (name.compare("l_NP") == 0) {

        throw std::runtime_error("NPEffective2::SetParameter(): Not implemented yet!");


        l_NP = value;
        cHQ1p = 0.0;
        cHQ1 = 0.0;
        cHU1 = 0.0;
        cHU3 = 0.0;
        cHD1 = 0.0;
        cHD2 = 0.0;
    } else if (name.compare("Lambda_NP") == 0)
        LambdaNP = value;
    else
        StandardModel::SetParameter(name, value);
}





