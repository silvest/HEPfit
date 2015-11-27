/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
// */

#include "GeneralSUSY.h"
#include <math.h>

const std::string GeneralSUSY::GeneralSUSYvars[NGeneralSUSYvars] = {
    "msQhat2_11r","msQhat2_12r","msQhat2_12i","msQhat2_13r","msQhat2_13i","msQhat2_22r","msQhat2_23r","msQhat2_23i","msQhat2_33r",
    "msUhat2_11r","msUhat2_12r","msUhat2_12i","msUhat2_13r","msUhat2_13i","msUhat2_22r","msUhat2_23r","msUhat2_23i","msUhat2_33r",
    "msDhat2_11r","msDhat2_12r","msDhat2_12i","msDhat2_13r","msDhat2_13i","msDhat2_22r","msDhat2_23r","msDhat2_23i","msDhat2_33r",
    "msLhat2_11r","msLhat2_12r","msLhat2_12i","msLhat2_13r","msLhat2_13i","msLhat2_22r","msLhat2_23r","msLhat2_23i","msLhat2_33r",
    "msEhat2_11r","msEhat2_12r","msEhat2_12i","msEhat2_13r","msEhat2_13i","msEhat2_22r","msEhat2_23r","msEhat2_23i","msEhat2_33r",
    "msNhat2_11r","msNhat2_12r","msNhat2_12i","msNhat2_13r","msNhat2_13i","msNhat2_22r","msNhat2_23r","msNhat2_23i","msNhat2_33r",
    "TUhat_11r","TUhat_12r","TUhat_13r","TUhat_21r","TUhat_22r","TUhat_23r","TUhat_31r","TUhat_32r","TUhat_33r",
    "TUhat_11i","TUhat_12i","TUhat_13i","TUhat_21i","TUhat_22i","TUhat_23i","TUhat_31i","TUhat_32i","TUhat_33i",
    "TDhat_11r","TDhat_12r","TDhat_13r","TDhat_21r","TDhat_22r","TDhat_23r","TDhat_31r","TDhat_32r","TDhat_33r",
    "TDhat_11i","TDhat_12i","TDhat_13i","TDhat_21i","TDhat_22i","TDhat_23i","TDhat_31i","TDhat_32i","TDhat_33i",
    "TEhat_11r","TEhat_12r","TEhat_13r","TEhat_21r","TEhat_22r","TEhat_23r","TEhat_31r","TEhat_32r","TEhat_33r",
    "TEhat_11i","TEhat_12i","TEhat_13i","TEhat_21i","TEhat_22i","TEhat_23i","TEhat_31i","TEhat_32i","TEhat_33i",
    "TNhat_11r","TNhat_12r","TNhat_13r","TNhat_21r","TNhat_22r","TNhat_23r","TNhat_31r","TNhat_32r","TNhat_33r",
    "TNhat_11i","TNhat_12i","TNhat_13i","TNhat_21i","TNhat_22i","TNhat_23i","TNhat_31i","TNhat_32i","TNhat_33i"
};

GeneralSUSY::GeneralSUSY()
: SUSY()
{
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_11r", boost::cref(msQhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_12r", boost::cref(msQhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_12i", boost::cref(msQhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_13r", boost::cref(msQhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_13i", boost::cref(msQhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_22r", boost::cref(msQhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_23r", boost::cref(msQhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_23i", boost::cref(msQhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msQhat2_33r", boost::cref(msQhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_11r", boost::cref(msUhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_12r", boost::cref(msUhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_12i", boost::cref(msUhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_13r", boost::cref(msUhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_13i", boost::cref(msUhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_22r", boost::cref(msUhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_23r", boost::cref(msUhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_23i", boost::cref(msUhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msUhat2_33r", boost::cref(msUhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_11r", boost::cref(msDhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_12r", boost::cref(msDhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_12i", boost::cref(msDhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_13r", boost::cref(msDhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_13i", boost::cref(msDhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_22r", boost::cref(msDhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_23r", boost::cref(msDhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_23i", boost::cref(msDhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msDhat2_33r", boost::cref(msDhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_11r", boost::cref(msLhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_12r", boost::cref(msLhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_12i", boost::cref(msLhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_13r", boost::cref(msLhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_13i", boost::cref(msLhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_22r", boost::cref(msLhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_23r", boost::cref(msLhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_23i", boost::cref(msLhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msLhat2_33r", boost::cref(msLhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_11r", boost::cref(msEhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_12r", boost::cref(msEhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_12i", boost::cref(msEhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_13r", boost::cref(msEhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_13i", boost::cref(msEhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_22r", boost::cref(msEhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_23r", boost::cref(msEhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_23i", boost::cref(msEhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msEhat2_33r", boost::cref(msEhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_11r", boost::cref(msNhat2_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_12r", boost::cref(msNhat2_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_12i", boost::cref(msNhat2_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_13r", boost::cref(msNhat2_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_13i", boost::cref(msNhat2_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_22r", boost::cref(msNhat2_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_23r", boost::cref(msNhat2_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_23i", boost::cref(msNhat2_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("msNhat2_33r", boost::cref(msNhat2_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_11r", boost::cref(TUhat_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_12r", boost::cref(TUhat_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_13r", boost::cref(TUhat_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_21r", boost::cref(TUhat_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_22r", boost::cref(TUhat_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_23r", boost::cref(TUhat_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_31r", boost::cref(TUhat_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_32r", boost::cref(TUhat_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_33r", boost::cref(TUhat_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_11i", boost::cref(TUhat_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_12i", boost::cref(TUhat_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_13i", boost::cref(TUhat_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_21i", boost::cref(TUhat_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_22i", boost::cref(TUhat_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_23i", boost::cref(TUhat_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_31i", boost::cref(TUhat_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_32i", boost::cref(TUhat_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TUhat_33i", boost::cref(TUhat_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_11r", boost::cref(TDhat_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_12r", boost::cref(TDhat_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_13r", boost::cref(TDhat_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_21r", boost::cref(TDhat_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_22r", boost::cref(TDhat_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_23r", boost::cref(TDhat_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_31r", boost::cref(TDhat_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_32r", boost::cref(TDhat_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_33r", boost::cref(TDhat_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_11i", boost::cref(TDhat_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_12i", boost::cref(TDhat_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_13i", boost::cref(TDhat_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_21i", boost::cref(TDhat_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_22i", boost::cref(TDhat_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_23i", boost::cref(TDhat_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_31i", boost::cref(TDhat_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_32i", boost::cref(TDhat_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TDhat_33i", boost::cref(TDhat_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_11r", boost::cref(TEhat_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_12r", boost::cref(TEhat_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_13r", boost::cref(TEhat_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_21r", boost::cref(TEhat_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_22r", boost::cref(TEhat_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_23r", boost::cref(TEhat_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_31r", boost::cref(TEhat_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_32r", boost::cref(TEhat_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_33r", boost::cref(TEhat_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_11i", boost::cref(TEhat_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_12i", boost::cref(TEhat_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_13i", boost::cref(TEhat_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_21i", boost::cref(TEhat_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_22i", boost::cref(TEhat_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_23i", boost::cref(TEhat_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_31i", boost::cref(TEhat_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_32i", boost::cref(TEhat_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TEhat_33i", boost::cref(TEhat_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_11r", boost::cref(TNhat_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_12r", boost::cref(TNhat_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_13r", boost::cref(TNhat_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_21r", boost::cref(TNhat_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_22r", boost::cref(TNhat_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_23r", boost::cref(TNhat_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_31r", boost::cref(TNhat_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_32r", boost::cref(TNhat_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_33r", boost::cref(TNhat_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_11i", boost::cref(TNhat_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_12i", boost::cref(TNhat_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_13i", boost::cref(TNhat_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_21i", boost::cref(TNhat_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_22i", boost::cref(TNhat_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_23i", boost::cref(TNhat_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_31i", boost::cref(TNhat_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_32i", boost::cref(TNhat_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("TNhat_33i", boost::cref(TNhat_33i)));
}

bool GeneralSUSY::InitializeModel()
{
    setModelInitialized(SUSY::InitializeModel());
    return (IsModelInitialized());
}

bool GeneralSUSY::Init(const std::map<std::string, double>& DPars)
{
    return(SUSY::Init(DPars));
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
        setParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate()) return (false);
    
    return (true);
}

bool GeneralSUSY::PostUpdate()
{
    if (!SUSY::PostUpdate()) return (false);
    return (true);
}

void GeneralSUSY::setParameter(const std::string name, const double& value)
{
    if(name.compare("msQhat2_11r") == 0)
        msQhat2_11r = value;
    else if(name.compare("msQhat2_12r") == 0)
        msQhat2_12r = value;
    else if(name.compare("msQhat2_12i") == 0)
        msQhat2_12i = value;
    else if(name.compare("msQhat2_13r") == 0)
        msQhat2_13r = value;
    else if(name.compare("msQhat2_13i") == 0)
        msQhat2_13i = value;
    else if(name.compare("msQhat2_22r") == 0)
        msQhat2_22r = value;
    else if(name.compare("msQhat2_23r") == 0)
        msQhat2_23r = value;
    else if(name.compare("msQhat2_23i") == 0)
        msQhat2_23i = value;
    else if(name.compare("msQhat2_33r") == 0)
        msQhat2_33r = value;
    else if(name.compare("msUhat2_11r") == 0)
        msUhat2_11r = value;
    else if(name.compare("msUhat2_12r") == 0)
        msUhat2_12r = value;
    else if(name.compare("msUhat2_12i") == 0)
        msUhat2_12i = value;
    else if(name.compare("msUhat2_13r") == 0)
        msUhat2_13r = value;
    else if(name.compare("msUhat2_13i") == 0)
        msUhat2_13i = value;
    else if(name.compare("msUhat2_22r") == 0)
        msUhat2_22r = value;
    else if(name.compare("msUhat2_23r") == 0)
        msUhat2_23r = value;
    else if(name.compare("msUhat2_23i") == 0)
        msUhat2_23i = value;
    else if(name.compare("msUhat2_33r") == 0)
        msUhat2_33r = value;
    else if(name.compare("msDhat2_11r") == 0)
        msDhat2_11r = value;
    else if(name.compare("msDhat2_12r") == 0)
        msDhat2_12r = value;
    else if(name.compare("msDhat2_12i") == 0)
        msDhat2_12i = value;
    else if(name.compare("msDhat2_13r") == 0)
        msDhat2_13r = value;
    else if(name.compare("msDhat2_13i") == 0)
        msDhat2_13i = value;
    else if(name.compare("msDhat2_22r") == 0)
        msDhat2_22r = value;
    else if(name.compare("msDhat2_23r") == 0)
        msDhat2_23r = value;
    else if(name.compare("msDhat2_23i") == 0)
        msDhat2_23i = value;
    else if(name.compare("msDhat2_33r") == 0)
        msDhat2_33r = value;
    else if(name.compare("msLhat2_11r") == 0)
        msLhat2_11r = value;
    else if(name.compare("msLhat2_12r") == 0)
        msLhat2_12r = value;
    else if(name.compare("msLhat2_12i") == 0)
        msLhat2_12i = value;
    else if(name.compare("msLhat2_13r") == 0)
        msLhat2_13r = value;
    else if(name.compare("msLhat2_13i") == 0)
        msLhat2_13i = value;
    else if(name.compare("msLhat2_22r") == 0)
        msLhat2_22r = value;
    else if(name.compare("msLhat2_23r") == 0)
        msLhat2_23r = value;
    else if(name.compare("msLhat2_23i") == 0)
        msLhat2_23i = value;
    else if(name.compare("msLhat2_33r") == 0)
        msLhat2_33r = value;
    else if(name.compare("msEhat2_11r") == 0)
        msEhat2_11r = value;
    else if(name.compare("msEhat2_12r") == 0)
        msEhat2_12r = value;
    else if(name.compare("msEhat2_12i") == 0)
        msEhat2_12i = value;
    else if(name.compare("msEhat2_13r") == 0)
        msEhat2_13r = value;
    else if(name.compare("msEhat2_13i") == 0)
        msEhat2_13i = value;
    else if(name.compare("msEhat2_22r") == 0)
        msEhat2_22r = value;
    else if(name.compare("msEhat2_23r") == 0)
        msEhat2_23r = value;
    else if(name.compare("msEhat2_23i") == 0)
        msEhat2_23i = value;
    else if(name.compare("msEhat2_33r") == 0)
        msEhat2_33r = value;
    else if(name.compare("msNhat2_11r") == 0)
        msNhat2_11r = value;
    else if(name.compare("msNhat2_12r") == 0)
        msNhat2_12r = value;
    else if(name.compare("msNhat2_12i") == 0)
        msNhat2_12i = value;
    else if(name.compare("msNhat2_13r") == 0)
        msNhat2_13r = value;
    else if(name.compare("msNhat2_13i") == 0)
        msNhat2_13i = value;
    else if(name.compare("msNhat2_22r") == 0)
        msNhat2_22r = value;
    else if(name.compare("msNhat2_23r") == 0)
        msNhat2_23r = value;
    else if(name.compare("msNhat2_23i") == 0)
        msNhat2_23i = value;
    else if(name.compare("msNhat2_33r") == 0)
        msNhat2_33r = value;
    else if(name.compare("TUhat_11r") == 0)
        TUhat_11r = value;
    else if(name.compare("TUhat_11i") == 0)
        TUhat_11i = value;
    else if(name.compare("TUhat_12r") == 0)
        TUhat_12r = value;
    else if(name.compare("TUhat_12i") == 0)
        TUhat_12i = value;
    else if(name.compare("TUhat_13r") == 0)
        TUhat_13r = value;
    else if(name.compare("TUhat_13i") == 0)
        TUhat_13i = value;
    else if(name.compare("TUhat_21r") == 0)
        TUhat_21r = value;
    else if(name.compare("TUhat_21i") == 0)
        TUhat_21i = value;
    else if(name.compare("TUhat_22r") == 0)
        TUhat_22r = value;
    else if(name.compare("TUhat_22i") == 0)
        TUhat_22i = value;
    else if(name.compare("TUhat_23r") == 0)
        TUhat_23r = value;
    else if(name.compare("TUhat_23i") == 0)
        TUhat_23i = value;
    else if(name.compare("TUhat_31r") == 0)
        TUhat_31r = value;
    else if(name.compare("TUhat_31i") == 0)
        TUhat_31i = value;
    else if(name.compare("TUhat_32r") == 0)
        TUhat_32r = value;
    else if(name.compare("TUhat_32i") == 0)
        TUhat_32i = value;
    else if(name.compare("TUhat_33r") == 0)
        TUhat_33r = value;
    else if(name.compare("TUhat_33i") == 0)
        TUhat_33i = value;
    else if(name.compare("TDhat_11r") == 0)
        TDhat_11r = value;
    else if(name.compare("TDhat_11i") == 0)
        TDhat_11i = value;
    else if(name.compare("TDhat_12r") == 0)
        TDhat_12r = value;
    else if(name.compare("TDhat_12i") == 0)
        TDhat_12i = value;
    else if(name.compare("TDhat_13r") == 0)
        TDhat_13r = value;
    else if(name.compare("TDhat_13i") == 0)
        TDhat_13i = value;
    else if(name.compare("TDhat_21r") == 0)
        TDhat_21r = value;
    else if(name.compare("TDhat_21i") == 0)
        TDhat_21i = value;
    else if(name.compare("TDhat_22r") == 0)
        TDhat_22r = value;
    else if(name.compare("TDhat_22i") == 0)
        TDhat_22i = value;
    else if(name.compare("TDhat_23r") == 0)
        TDhat_23r = value;
    else if(name.compare("TDhat_23i") == 0)
        TDhat_23i = value;
    else if(name.compare("TDhat_31r") == 0)
        TDhat_31r = value;
    else if(name.compare("TDhat_31i") == 0)
        TDhat_31i = value;
    else if(name.compare("TDhat_32r") == 0)
        TDhat_32r = value;
    else if(name.compare("TDhat_32i") == 0)
        TDhat_32i = value;
    else if(name.compare("TDhat_33r") == 0)
        TDhat_33r = value;
    else if(name.compare("TDhat_33i") == 0)
        TDhat_33i = value;
    else if(name.compare("TEhat_11r") == 0)
        TEhat_11r = value;
    else if(name.compare("TEhat_11i") == 0)
        TEhat_11i = value;
    else if(name.compare("TEhat_12r") == 0)
        TEhat_12r = value;
    else if(name.compare("TEhat_12i") == 0)
        TEhat_12i = value;
    else if(name.compare("TEhat_13r") == 0)
        TEhat_13r = value;
    else if(name.compare("TEhat_13i") == 0)
        TEhat_13i = value;
    else if(name.compare("TEhat_21r") == 0)
        TEhat_21r = value;
    else if(name.compare("TEhat_21i") == 0)
        TEhat_21i = value;
    else if(name.compare("TEhat_22r") == 0)
        TEhat_22r = value;
    else if(name.compare("TEhat_22i") == 0)
        TEhat_22i = value;
    else if(name.compare("TEhat_23r") == 0)
        TEhat_23r = value;
    else if(name.compare("TEhat_23i") == 0)
        TEhat_23i = value;
    else if(name.compare("TEhat_31r") == 0)
        TEhat_31r = value;
    else if(name.compare("TEhat_31i") == 0)
        TEhat_31i = value;
    else if(name.compare("TEhat_32r") == 0)
        TEhat_32r = value;
    else if(name.compare("TEhat_32i") == 0)
        TEhat_32i = value;
    else if(name.compare("TEhat_33r") == 0)
        TEhat_33r = value;
    else if(name.compare("TEhat_33i") == 0)
        TEhat_33i = value;
    else if(name.compare("TNhat_11r") == 0)
        TNhat_11r = value;
    else if(name.compare("TNhat_11i") == 0)
        TNhat_11i = value;
    else if(name.compare("TNhat_12r") == 0)
        TNhat_12r = value;
    else if(name.compare("TNhat_12i") == 0)
        TNhat_12i = value;
    else if(name.compare("TNhat_13r") == 0)
        TNhat_13r = value;
    else if(name.compare("TNhat_13i") == 0)
        TNhat_13i = value;
    else if(name.compare("TNhat_21r") == 0)
        TNhat_21r = value;
    else if(name.compare("TNhat_21i") == 0)
        TNhat_21i = value;
    else if(name.compare("TNhat_22r") == 0)
        TNhat_22r = value;
    else if(name.compare("TNhat_22i") == 0)
        TNhat_22i = value;
    else if(name.compare("TNhat_23r") == 0)
        TNhat_23r = value;
    else if(name.compare("TNhat_23i") == 0)
        TNhat_23i = value;
    else if(name.compare("TNhat_31r") == 0)
        TNhat_31r = value;
    else if(name.compare("TNhat_31i") == 0)
        TNhat_31i = value;
    else if(name.compare("TNhat_32r") == 0)
        TNhat_32r = value;
    else if(name.compare("TNhat_32i") == 0)
        TNhat_32i = value;
    else if(name.compare("TNhat_33r") == 0)
        TNhat_33r = value;
    else if(name.compare("TNhat_33i") == 0)
        TNhat_33i = value;
    else
        SUSY::setParameter(name, value);
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
    msQhat2.assign(0,0, msQhat2_11r);
    msQhat2.assign(0,1, gslpp::complex(msQhat2_12r, msQhat2_12i));
    msQhat2.assign(0,2, gslpp::complex(msQhat2_13r, msQhat2_13i));
    msQhat2.assign(1,1, msQhat2_22r);
    msQhat2.assign(1,2, gslpp::complex(msQhat2_23r, msQhat2_23i));
    msQhat2.assign(1,0, msQhat2(0,1).conjugate());
    msQhat2.assign(2,0, msQhat2(0,2).conjugate());
    msQhat2.assign(2,1, msQhat2(1,2).conjugate());
    msQhat2.assign(2,2, msQhat2_33r);
    
    msUhat2.assign(0,0, msUhat2_11r);
    msUhat2.assign(0,1, gslpp::complex(msUhat2_12r, msUhat2_12i));
    msUhat2.assign(0,2, gslpp::complex(msUhat2_13r, msUhat2_13i));
    msUhat2.assign(1,1, msUhat2_22r);
    msUhat2.assign(1,2, gslpp::complex(msUhat2_23r, msUhat2_23i));
    msUhat2.assign(2,2, msUhat2_33r);
    msUhat2.assign(1,0, msUhat2(0,1).conjugate());
    msUhat2.assign(2,0, msUhat2(0,2).conjugate());
    msUhat2.assign(2,1, msUhat2(1,2).conjugate());
    
    msDhat2.assign(0,0, msDhat2_11r);
    msDhat2.assign(0,1, gslpp::complex(msDhat2_12r, msDhat2_12i));
    msDhat2.assign(0,2, gslpp::complex(msDhat2_13r, msDhat2_13i));
    msDhat2.assign(1,1, msDhat2_22r);
    msDhat2.assign(1,2, gslpp::complex(msDhat2_23r, msDhat2_23i));
    msDhat2.assign(2,2, msDhat2_33r);
    msDhat2.assign(1,0, msDhat2(0,1).conjugate());
    msDhat2.assign(2,0, msDhat2(0,2).conjugate());
    msDhat2.assign(2,1, msDhat2(1,2).conjugate());
    
    msLhat2.assign(0,0, msLhat2_11r);
    msLhat2.assign(0,1, gslpp::complex(msLhat2_12r, msLhat2_12i));
    msLhat2.assign(0,2, gslpp::complex(msLhat2_13r, msLhat2_13i));
    msLhat2.assign(1,1, msLhat2_22r);
    msLhat2.assign(1,2, gslpp::complex(msLhat2_23r, msLhat2_23i));
    msLhat2.assign(2,2, msLhat2_33r);
    msLhat2.assign(1,0, msLhat2(0,1).conjugate());
    msLhat2.assign(2,0, msLhat2(0,2).conjugate());
    msLhat2.assign(2,1, msLhat2(1,2).conjugate());
    
    msEhat2.assign(0,0, msEhat2_11r);
    msEhat2.assign(0,1, gslpp::complex(msEhat2_12r, msEhat2_12i));
    msEhat2.assign(0,2, gslpp::complex(msEhat2_13r, msEhat2_13i));
    msEhat2.assign(1,1, msEhat2_22r);
    msEhat2.assign(1,2, gslpp::complex(msEhat2_23r, msEhat2_23i));
    msEhat2.assign(2,2, msEhat2_33r);
    msEhat2.assign(1,0, msEhat2(0,1).conjugate());
    msEhat2.assign(2,0, msEhat2(0,2).conjugate());
    msEhat2.assign(2,1, msEhat2(1,2).conjugate());
    
    msNhat2.assign(0,0, msNhat2_11r);
    msNhat2.assign(0,1, gslpp::complex(msNhat2_12r, msNhat2_12i));
    msNhat2.assign(0,2, gslpp::complex(msNhat2_13r, msNhat2_13i));
    msNhat2.assign(1,1, msNhat2_22r);
    msNhat2.assign(1,2, gslpp::complex(msNhat2_23r, msNhat2_23i));
    msNhat2.assign(2,2, msNhat2_33r);
    msNhat2.assign(1,0, msNhat2(0,1).conjugate());
    msNhat2.assign(2,0, msNhat2(0,2).conjugate());
    msNhat2.assign(2,1, msNhat2(1,2).conjugate());
    
    TUhat.assign(0,0, gslpp::complex(TUhat_11r, TUhat_11i));
    TUhat.assign(0,1, gslpp::complex(TUhat_12r, TUhat_12i));
    TUhat.assign(0,2, gslpp::complex(TUhat_13r, TUhat_13i));
    TUhat.assign(1,0, gslpp::complex(TUhat_21r, TUhat_21i));
    TUhat.assign(1,1, gslpp::complex(TUhat_22r, TUhat_22i));
    TUhat.assign(1,2, gslpp::complex(TUhat_23r, TUhat_23i));
    TUhat.assign(2,0, gslpp::complex(TUhat_31r, TUhat_31i));
    TUhat.assign(2,1, gslpp::complex(TUhat_32r, TUhat_32i));
    TUhat.assign(2,2, gslpp::complex(TUhat_33r, TUhat_33i));
    
    TDhat.assign(0,0, gslpp::complex(TDhat_11r, TDhat_11i));
    TDhat.assign(0,1, gslpp::complex(TDhat_12r, TDhat_12i));
    TDhat.assign(0,2, gslpp::complex(TDhat_13r, TDhat_13i));
    TDhat.assign(1,0, gslpp::complex(TDhat_21r, TDhat_21i));
    TDhat.assign(1,1, gslpp::complex(TDhat_22r, TDhat_22i));
    TDhat.assign(1,2, gslpp::complex(TDhat_23r, TDhat_23i));
    TDhat.assign(2,0, gslpp::complex(TDhat_31r, TDhat_31i));
    TDhat.assign(2,1, gslpp::complex(TDhat_32r, TDhat_32i));
    TDhat.assign(2,2, gslpp::complex(TDhat_33r, TDhat_33i));

    TEhat.assign(0,0, gslpp::complex(TEhat_11r, TEhat_11i));
    TEhat.assign(0,1, gslpp::complex(TEhat_12r, TEhat_12i));
    TEhat.assign(0,2, gslpp::complex(TEhat_13r, TEhat_13i));
    TEhat.assign(1,0, gslpp::complex(TEhat_21r, TEhat_21i));
    TEhat.assign(1,1, gslpp::complex(TEhat_22r, TEhat_22i));
    TEhat.assign(1,2, gslpp::complex(TEhat_23r, TEhat_23i));
    TEhat.assign(2,0, gslpp::complex(TEhat_31r, TEhat_31i));
    TEhat.assign(2,1, gslpp::complex(TEhat_32r, TEhat_32i));
    TEhat.assign(2,2, gslpp::complex(TEhat_33r, TEhat_33i));
    
    TNhat.assign(0,0, gslpp::complex(TNhat_11r, TNhat_11i));
    TNhat.assign(0,1, gslpp::complex(TNhat_12r, TNhat_12i));
    TNhat.assign(0,2, gslpp::complex(TNhat_13r, TNhat_13i));
    TNhat.assign(1,0, gslpp::complex(TNhat_21r, TNhat_21i));
    TNhat.assign(1,1, gslpp::complex(TNhat_22r, TNhat_22i));
    TNhat.assign(1,2, gslpp::complex(TNhat_23r, TNhat_23i));
    TNhat.assign(2,0, gslpp::complex(TNhat_31r, TNhat_31i));
    TNhat.assign(2,1, gslpp::complex(TNhat_32r, TNhat_32i));
    TNhat.assign(2,2, gslpp::complex(TNhat_33r, TNhat_33i));
}

