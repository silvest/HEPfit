/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <math.h>
#include "pMSSM.h"

const std::string pMSSM::pMSSMvars[NpMSSMvars] = {
    "msQ12", "msQ3", "msU12", "msU3", "msD12", "msD3",
    "msL12", "msL3", "msE12", "msE3",
    "AU", "AD", "AE"
};

pMSSM::pMSSM()
: SUSY()
{
}

bool pMSSM::InitializeModel()
{
    SetModelInitialized(SUSY::InitializeModel());
    return (IsModelInitialized());
}

bool pMSSM::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return (CheckParameters(DPars));
}

bool pMSSM::PreUpdate()
{    
    if(!SUSY::PreUpdate()) return (false);
    return (true);
}

bool pMSSM::Update(const std::map<std::string, double>& DPars)
{    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate()) return (false);
    
    return (true);
}

bool pMSSM::PostUpdate()
{
    if (!SUSY::PostUpdate()) return (false);
    return (true);
}

void pMSSM::SetParameter(const std::string name, const double& value)
{
    if(name.compare("msQ12") == 0)
        msQ12 = value;
    else if(name.compare("msQ3") == 0)
        msQ3 = value;
    else if(name.compare("msU12") == 0)
        msU12 = value;
    else if(name.compare("msU3") == 0)
        msU3 = value;
    else if(name.compare("msD12") == 0)
        msD12 = value;
    else if(name.compare("msD3") == 0)
        msD3 = value;
    else if(name.compare("msL12") == 0)
        msL12 = value;
    else if(name.compare("msL3") == 0)
        msL3 = value;
    else if(name.compare("msE12") == 0)
        msE12 = value;
    else if(name.compare("msE3") == 0)
        msE3 = value;
    else if(name.compare("AU") == 0)
        AU = value;
    else if(name.compare("AD") == 0)
        AD = value;
    else if(name.compare("AE") == 0)
        AE = value;
    else
        SUSY::SetParameter(name, value);
}

bool pMSSM::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NpMSSMvars; i++) {
        if (DPars.find(pMSSMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory pMSSM parameter " << pMSSMvars[i] << std::endl;
            return false;
        }
    }
    return(SUSY::CheckParameters(DPars));
}

void pMSSM::SetSoftTerms()
{
    MsQ2.assign(0,0, msQ12*msQ12);
    MsQ2.assign(1,1, msQ12*msQ12);
    MsQ2.assign(2,2, msQ3*msQ3);
    
    MsU2.assign(0,0, msU12*msU12);
    MsU2.assign(1,1, msU12*msU12);
    MsU2.assign(2,2, msU3*msU3);
    
    MsD2.assign(0,0, msD12*msD12);
    MsD2.assign(1,1, msD12*msD12);
    MsD2.assign(2,2, msD3*msD3);
    
    MsL2.assign(0,0, msL12*msL12);
    MsL2.assign(1,1, msL12*msL12);
    MsL2.assign(2,2, msL3*msL3);
    
    MsE2.assign(0,0, msE12*msE12);
    MsE2.assign(1,1, msE12*msE12);
    MsE2.assign(2,2, msE3*msE3);
    
    MsN2.assign(0,0, 1.e20);
    MsN2.assign(1,1, 1.e20);
    MsN2.assign(2,2, 1.e20);
    
    TU.assign(2,2, AU*Yu(2,2));
    
    TD.assign(2,2, AD*Yd(2,2));
    
    TE.assign(2,2, AE*Ye(2,2));
}

