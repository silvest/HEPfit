/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <math.h>
#include "GeneralSUSY.h"

const std::string GeneralSUSY::GeneralSUSYvars[NGeneralSUSYvars] = {
    "msQ2_11r","msQ2_12r","msQ2_12i","msQ2_13r","msQ2_13i","msQ2_22r","msQ2_23r","msQ2_23i","msQ2_33r",
    "msU2_11r","msU2_12r","msU2_12i","msU2_13r","msU2_13i","msU2_22r","msU2_23r","msU2_23i","msU2_33r",
    "msD2_11r","msD2_12r","msD2_12i","msD2_13r","msD2_13i","msD2_22r","msD2_23r","msD2_23i","msD2_33r",
    "msL2_11r","msL2_12r","msL2_12i","msL2_13r","msL2_13i","msL2_22r","msL2_23r","msL2_23i","msL2_33r",
    "msE2_11r","msE2_12r","msE2_12i","msE2_13r","msE2_13i","msE2_22r","msE2_23r","msE2_23i","msE2_33r",
    "msN2_11r","msN2_12r","msN2_12i","msN2_13r","msN2_13i","msN2_22r","msN2_23r","msN2_23i","msN2_33r",
    "TU_11r","TU_12r","TU_13r","TU_21r","TU_22r","TU_23r","TU_31r","TU_32r","TU_33r",
    "TU_11i","TU_12i","TU_13i","TU_21i","TU_22i","TU_23i","TU_31i","TU_32i","TU_33i",
    "TD_11r","TD_12r","TD_13r","TD_21r","TD_22r","TD_23r","TD_31r","TD_32r","TD_33r",
    "TD_11i","TD_12i","TD_13i","TD_21i","TD_22i","TD_23i","TD_31i","TD_32i","TD_33i",
    "TE_11r","TE_12r","TE_13r","TE_21r","TE_22r","TE_23r","TE_31r","TE_32r","TE_33r",
    "TE_11i","TE_12i","TE_13i","TE_21i","TE_22i","TE_23i","TE_31i","TE_32i","TE_33i",
    "TN_11r","TN_12r","TN_13r","TN_21r","TN_22r","TN_23r","TN_31r","TN_32r","TN_33r",
    "TN_11i","TN_12i","TN_13i","TN_21i","TN_22i","TN_23i","TN_31i","TN_32i","TN_33i"
};

GeneralSUSY::GeneralSUSY()
: SUSY()
{
}

bool GeneralSUSY::InitializeModel()
{
    SetModelInitialized(SUSY::InitializeModel());
    return (IsModelInitialized());
}

bool GeneralSUSY::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return (CheckParameters(DPars));
}

bool GeneralSUSY::PreUpdate()
{    
    if(!SUSY::PreUpdate()) return (false);
    return (true);
}

bool GeneralSUSY::Update(const std::map<std::string, double>& DPars)
{    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate()) return (false);
    
    return (true);
}

bool GeneralSUSY::PostUpdate()
{
    if (!SUSY::PostUpdate()) return (false);
    return (true);
}

void GeneralSUSY::SetParameter(const std::string name, const double& value)
{
    if(name.compare("msQ2_11r") == 0)
        msQ2_11r = value;
    else if(name.compare("msQ2_12r") == 0)
        msQ2_12r = value;
    else if(name.compare("msQ2_12i") == 0)
        msQ2_12i = value;
    else if(name.compare("msQ2_13r") == 0)
        msQ2_13r = value;
    else if(name.compare("msQ2_13i") == 0)
        msQ2_13i = value;
    else if(name.compare("msQ2_22r") == 0)
        msQ2_22r = value;
    else if(name.compare("msQ2_23r") == 0)
        msQ2_23r = value;
    else if(name.compare("msQ2_23i") == 0)
        msQ2_23i = value;
    else if(name.compare("msQ2_33r") == 0)
        msQ2_33r = value;
    else if(name.compare("msU2_11r") == 0)
        msU2_11r = value;
    else if(name.compare("msU2_12r") == 0)
        msU2_12r = value;
    else if(name.compare("msU2_12i") == 0)
        msU2_12i = value;
    else if(name.compare("msU2_13r") == 0)
        msU2_13r = value;
    else if(name.compare("msU2_13i") == 0)
        msU2_13i = value;
    else if(name.compare("msU2_22r") == 0)
        msU2_22r = value;
    else if(name.compare("msU2_23r") == 0)
        msU2_23r = value;
    else if(name.compare("msU2_23i") == 0)
        msU2_23i = value;
    else if(name.compare("msU2_33r") == 0)
        msU2_33r = value;
    else if(name.compare("msD2_11r") == 0)
        msD2_11r = value;
    else if(name.compare("msD2_12r") == 0)
        msD2_12r = value;
    else if(name.compare("msD2_12i") == 0)
        msD2_12i = value;
    else if(name.compare("msD2_13r") == 0)
        msD2_13r = value;
    else if(name.compare("msD2_13i") == 0)
        msD2_13i = value;
    else if(name.compare("msD2_22r") == 0)
        msD2_22r = value;
    else if(name.compare("msD2_23r") == 0)
        msD2_23r = value;
    else if(name.compare("msD2_23i") == 0)
        msD2_23i = value;
    else if(name.compare("msD2_33r") == 0)
        msD2_33r = value;
    else if(name.compare("msL2_11r") == 0)
        msL2_11r = value;
    else if(name.compare("msL2_12r") == 0)
        msL2_12r = value;
    else if(name.compare("msL2_12i") == 0)
        msL2_12i = value;
    else if(name.compare("msL2_13r") == 0)
        msL2_13r = value;
    else if(name.compare("msL2_13i") == 0)
        msL2_13i = value;
    else if(name.compare("msL2_22r") == 0)
        msL2_22r = value;
    else if(name.compare("msL2_23r") == 0)
        msL2_23r = value;
    else if(name.compare("msL2_23i") == 0)
        msL2_23i = value;
    else if(name.compare("msL2_33r") == 0)
        msL2_33r = value;
    else if(name.compare("msE2_11r") == 0)
        msE2_11r = value;
    else if(name.compare("msE2_12r") == 0)
        msE2_12r = value;
    else if(name.compare("msE2_12i") == 0)
        msE2_12i = value;
    else if(name.compare("msE2_13r") == 0)
        msE2_13r = value;
    else if(name.compare("msE2_13i") == 0)
        msE2_13i = value;
    else if(name.compare("msE2_22r") == 0)
        msE2_22r = value;
    else if(name.compare("msE2_23r") == 0)
        msE2_23r = value;
    else if(name.compare("msE2_23i") == 0)
        msE2_23i = value;
    else if(name.compare("msE2_33r") == 0)
        msE2_33r = value;
    else if(name.compare("msN2_11r") == 0)
        msN2_11r = value;
    else if(name.compare("msN2_12r") == 0)
        msN2_12r = value;
    else if(name.compare("msN2_12i") == 0)
        msN2_12i = value;
    else if(name.compare("msN2_13r") == 0)
        msN2_13r = value;
    else if(name.compare("msN2_13i") == 0)
        msN2_13i = value;
    else if(name.compare("msN2_22r") == 0)
        msN2_22r = value;
    else if(name.compare("msN2_23r") == 0)
        msN2_23r = value;
    else if(name.compare("msN2_23i") == 0)
        msN2_23i = value;
    else if(name.compare("msN2_33r") == 0)
        msN2_33r = value;
    else if(name.compare("TU_11r") == 0)
        TU_11r = value;
    else if(name.compare("TU_11i") == 0)
        TU_11i = value;
    else if(name.compare("TU_12r") == 0)
        TU_12r = value;
    else if(name.compare("TU_12i") == 0)
        TU_12i = value;
    else if(name.compare("TU_13r") == 0)
        TU_13r = value;
    else if(name.compare("TU_13i") == 0)
        TU_13i = value;
    else if(name.compare("TU_21r") == 0)
        TU_21r = value;
    else if(name.compare("TU_21i") == 0)
        TU_21i = value;
    else if(name.compare("TU_22r") == 0)
        TU_22r = value;
    else if(name.compare("TU_22i") == 0)
        TU_22i = value;
    else if(name.compare("TU_23r") == 0)
        TU_23r = value;
    else if(name.compare("TU_23i") == 0)
        TU_23i = value;
    else if(name.compare("TU_31r") == 0)
        TU_31r = value;
    else if(name.compare("TU_31i") == 0)
        TU_31i = value;
    else if(name.compare("TU_32r") == 0)
        TU_32r = value;
    else if(name.compare("TU_32i") == 0)
        TU_32i = value;
    else if(name.compare("TU_33r") == 0)
        TU_33r = value;
    else if(name.compare("TU_33i") == 0)
        TU_33i = value;
    else if(name.compare("TD_11r") == 0)
        TD_11r = value;
    else if(name.compare("TD_11i") == 0)
        TD_11i = value;
    else if(name.compare("TD_12r") == 0)
        TD_12r = value;
    else if(name.compare("TD_12i") == 0)
        TD_12i = value;
    else if(name.compare("TD_13r") == 0)
        TD_13r = value;
    else if(name.compare("TD_13i") == 0)
        TD_13i = value;
    else if(name.compare("TD_21r") == 0)
        TD_21r = value;
    else if(name.compare("TD_21i") == 0)
        TD_21i = value;
    else if(name.compare("TD_22r") == 0)
        TD_22r = value;
    else if(name.compare("TD_22i") == 0)
        TD_22i = value;
    else if(name.compare("TD_23r") == 0)
        TD_23r = value;
    else if(name.compare("TD_23i") == 0)
        TD_23i = value;
    else if(name.compare("TD_31r") == 0)
        TD_31r = value;
    else if(name.compare("TD_31i") == 0)
        TD_31i = value;
    else if(name.compare("TD_32r") == 0)
        TD_32r = value;
    else if(name.compare("TD_32i") == 0)
        TD_32i = value;
    else if(name.compare("TD_33r") == 0)
        TD_33r = value;
    else if(name.compare("TD_33i") == 0)
        TD_33i = value;
    else if(name.compare("TE_11r") == 0)
        TE_11r = value;
    else if(name.compare("TE_11i") == 0)
        TE_11i = value;
    else if(name.compare("TE_12r") == 0)
        TE_12r = value;
    else if(name.compare("TE_12i") == 0)
        TE_12i = value;
    else if(name.compare("TE_13r") == 0)
        TE_13r = value;
    else if(name.compare("TE_13i") == 0)
        TE_13i = value;
    else if(name.compare("TE_21r") == 0)
        TE_21r = value;
    else if(name.compare("TE_21i") == 0)
        TE_21i = value;
    else if(name.compare("TE_22r") == 0)
        TE_22r = value;
    else if(name.compare("TE_22i") == 0)
        TE_22i = value;
    else if(name.compare("TE_23r") == 0)
        TE_23r = value;
    else if(name.compare("TE_23i") == 0)
        TE_23i = value;
    else if(name.compare("TE_31r") == 0)
        TE_31r = value;
    else if(name.compare("TE_31i") == 0)
        TE_31i = value;
    else if(name.compare("TE_32r") == 0)
        TE_32r = value;
    else if(name.compare("TE_32i") == 0)
        TE_32i = value;
    else if(name.compare("TE_33r") == 0)
        TE_33r = value;
    else if(name.compare("TE_33i") == 0)
        TE_33i = value;
    else if(name.compare("TN_11r") == 0)
        TN_11r = value;
    else if(name.compare("TN_11i") == 0)
        TN_11i = value;
    else if(name.compare("TN_12r") == 0)
        TN_12r = value;
    else if(name.compare("TN_12i") == 0)
        TN_12i = value;
    else if(name.compare("TN_13r") == 0)
        TN_13r = value;
    else if(name.compare("TN_13i") == 0)
        TN_13i = value;
    else if(name.compare("TN_21r") == 0)
        TN_21r = value;
    else if(name.compare("TN_21i") == 0)
        TN_21i = value;
    else if(name.compare("TN_22r") == 0)
        TN_22r = value;
    else if(name.compare("TN_22i") == 0)
        TN_22i = value;
    else if(name.compare("TN_23r") == 0)
        TN_23r = value;
    else if(name.compare("TN_23i") == 0)
        TN_23i = value;
    else if(name.compare("TN_31r") == 0)
        TN_31r = value;
    else if(name.compare("TN_31i") == 0)
        TN_31i = value;
    else if(name.compare("TN_32r") == 0)
        TN_32r = value;
    else if(name.compare("TN_32i") == 0)
        TN_32i = value;
    else if(name.compare("TN_33r") == 0)
        TN_33r = value;
    else if(name.compare("TN_33i") == 0)
        TN_33i = value;
    else
        SUSY::SetParameter(name, value);
}

bool GeneralSUSY::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NGeneralSUSYvars; i++) {
        if (DPars.find(GeneralSUSYvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeneralSUSY parameter " << GeneralSUSYvars[i] << std::endl;
            return false;
        }
    }
    return(SUSY::CheckParameters(DPars));
}

void GeneralSUSY::SetSoftTerms()
{
    MsQ2.assign(0,0, msQ2_11r);
    MsQ2.assign(0,1, gslpp::complex(msQ2_12r, msQ2_12i));
    MsQ2.assign(0,2, gslpp::complex(msQ2_13r, msQ2_13i));
    MsQ2.assign(1,1, msQ2_22r);
    MsQ2.assign(1,2, gslpp::complex(msQ2_23r, msQ2_23i));
    MsQ2.assign(1,0, MsQ2(0,1).conjugate());
    MsQ2.assign(2,0, MsQ2(0,2).conjugate());
    MsQ2.assign(2,1, MsQ2(1,2).conjugate());
    MsQ2.assign(2,2, msQ2_33r);
    
    MsU2.assign(0,0, msU2_11r);
    MsU2.assign(0,1, gslpp::complex(msU2_12r, msU2_12i));
    MsU2.assign(0,2, gslpp::complex(msU2_13r, msU2_13i));
    MsU2.assign(1,1, msU2_22r);
    MsU2.assign(1,2, gslpp::complex(msU2_23r, msU2_23i));
    MsU2.assign(2,2, msU2_33r);
    MsU2.assign(1,0, MsU2(0,1).conjugate());
    MsU2.assign(2,0, MsU2(0,2).conjugate());
    MsU2.assign(2,1, MsU2(1,2).conjugate());
    
    MsD2.assign(0,0, msD2_11r);
    MsD2.assign(0,1, gslpp::complex(msD2_12r, msD2_12i));
    MsD2.assign(0,2, gslpp::complex(msD2_13r, msD2_13i));
    MsD2.assign(1,1, msD2_22r);
    MsD2.assign(1,2, gslpp::complex(msD2_23r, msD2_23i));
    MsD2.assign(2,2, msD2_33r);
    MsD2.assign(1,0, MsD2(0,1).conjugate());
    MsD2.assign(2,0, MsD2(0,2).conjugate());
    MsD2.assign(2,1, MsD2(1,2).conjugate());
    
    MsL2.assign(0,0, msL2_11r);
    MsL2.assign(0,1, gslpp::complex(msL2_12r, msL2_12i));
    MsL2.assign(0,2, gslpp::complex(msL2_13r, msL2_13i));
    MsL2.assign(1,1, msL2_22r);
    MsL2.assign(1,2, gslpp::complex(msL2_23r, msL2_23i));
    MsL2.assign(2,2, msL2_33r);
    MsL2.assign(1,0, MsL2(0,1).conjugate());
    MsL2.assign(2,0, MsL2(0,2).conjugate());
    MsL2.assign(2,1, MsL2(1,2).conjugate());
    
    MsE2.assign(0,0, msE2_11r);
    MsE2.assign(0,1, gslpp::complex(msE2_12r, msE2_12i));
    MsE2.assign(0,2, gslpp::complex(msE2_13r, msE2_13i));
    MsE2.assign(1,1, msE2_22r);
    MsE2.assign(1,2, gslpp::complex(msE2_23r, msE2_23i));
    MsE2.assign(2,2, msE2_33r);
    MsE2.assign(1,0, MsE2(0,1).conjugate());
    MsE2.assign(2,0, MsE2(0,2).conjugate());
    MsE2.assign(2,1, MsE2(1,2).conjugate());
    
    MsN2.assign(0,0, msN2_11r);
    MsN2.assign(0,1, gslpp::complex(msN2_12r, msN2_12i));
    MsN2.assign(0,2, gslpp::complex(msN2_13r, msN2_13i));
    MsN2.assign(1,1, msN2_22r);
    MsN2.assign(1,2, gslpp::complex(msN2_23r, msN2_23i));
    MsN2.assign(2,2, msN2_33r);
    MsN2.assign(1,0, MsN2(0,1).conjugate());
    MsN2.assign(2,0, MsN2(0,2).conjugate());
    MsN2.assign(2,1, MsN2(1,2).conjugate());
    
    TU.assign(0,0, gslpp::complex(TU_11r, TU_11i));
    TU.assign(0,1, gslpp::complex(TU_12r, TU_12i));
    TU.assign(0,2, gslpp::complex(TU_13r, TU_13i));
    TU.assign(1,0, gslpp::complex(TU_21r, TU_21i));
    TU.assign(1,1, gslpp::complex(TU_22r, TU_22i));
    TU.assign(1,2, gslpp::complex(TU_23r, TU_23i));
    TU.assign(2,0, gslpp::complex(TU_31r, TU_31i));
    TU.assign(2,1, gslpp::complex(TU_32r, TU_32i));
    TU.assign(2,2, gslpp::complex(TU_33r, TU_33i));
    
    TD.assign(0,0, gslpp::complex(TD_11r, TD_11i));
    TD.assign(0,1, gslpp::complex(TD_12r, TD_12i));
    TD.assign(0,2, gslpp::complex(TD_13r, TD_13i));
    TD.assign(1,0, gslpp::complex(TD_21r, TD_21i));
    TD.assign(1,1, gslpp::complex(TD_22r, TD_22i));
    TD.assign(1,2, gslpp::complex(TD_23r, TD_23i));
    TD.assign(2,0, gslpp::complex(TD_31r, TD_31i));
    TD.assign(2,1, gslpp::complex(TD_32r, TD_32i));
    TD.assign(2,2, gslpp::complex(TD_33r, TD_33i));

    TE.assign(0,0, gslpp::complex(TE_11r, TE_11i));
    TE.assign(0,1, gslpp::complex(TE_12r, TE_12i));
    TE.assign(0,2, gslpp::complex(TE_13r, TE_13i));
    TE.assign(1,0, gslpp::complex(TE_21r, TE_21i));
    TE.assign(1,1, gslpp::complex(TE_22r, TE_22i));
    TE.assign(1,2, gslpp::complex(TE_23r, TE_23i));
    TE.assign(2,0, gslpp::complex(TE_31r, TE_31i));
    TE.assign(2,1, gslpp::complex(TE_32r, TE_32i));
    TE.assign(2,2, gslpp::complex(TE_33r, TE_33i));
    
    TN.assign(0,0, gslpp::complex(TN_11r, TN_11i));
    TN.assign(0,1, gslpp::complex(TN_12r, TN_12i));
    TN.assign(0,2, gslpp::complex(TN_13r, TN_13i));
    TN.assign(1,0, gslpp::complex(TN_21r, TN_21i));
    TN.assign(1,1, gslpp::complex(TN_22r, TN_22i));
    TN.assign(1,2, gslpp::complex(TN_23r, TN_23i));
    TN.assign(2,0, gslpp::complex(TN_31r, TN_31i));
    TN.assign(2,1, gslpp::complex(TN_32r, TN_32i));
    TN.assign(2,2, gslpp::complex(TN_33r, TN_33i));
}

