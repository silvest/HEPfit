/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RealWeakEFTLFV.h"

const std::string RealWeakEFTLFV::RealWeakEFTLFVvars[NRealWeakEFTLFVvars] = {"C7", "C7p", "C8", "C8p", 
                                                                             "C9_11", "C9p_11", "C10_11", "C10p_11", "CS_11", "CSp_11", "CP_11", "CPp_11", 
                                                                             "C9_22", "C9p_22", "C10_22", "C10p_22", "CS_22", "CSp_22", "CP_22", "CPp_22",
                                                                             "WCscale"};

RealWeakEFTLFV::RealWeakEFTLFV() : StandardModel(), ReWEFTM(*this) {   

    SMM.setObj((StandardModelMatching&) ReWEFTM.getObj());
    ModelParamMap.insert(std::make_pair("C7", std::cref(C7)));
    ModelParamMap.insert(std::make_pair("C7p", std::cref(C7p)));
    ModelParamMap.insert(std::make_pair("C8", std::cref(C8)));
    ModelParamMap.insert(std::make_pair("C8p", std::cref(C8p)));
    
    ModelParamMap.insert(std::make_pair("C9_11", std::cref(C9_11)));
    ModelParamMap.insert(std::make_pair("C9p_11", std::cref(C9p_11)));
    ModelParamMap.insert(std::make_pair("C10_11", std::cref(C10_11)));
    ModelParamMap.insert(std::make_pair("C10p_11", std::cref(C10p_11)));
    ModelParamMap.insert(std::make_pair("CS_11", std::cref(CS_11)));
    ModelParamMap.insert(std::make_pair("CSp_11", std::cref(CSp_11)));
    ModelParamMap.insert(std::make_pair("CP_11", std::cref(CP_11)));
    ModelParamMap.insert(std::make_pair("CPp_11", std::cref(CPp_11)));
    
    ModelParamMap.insert(std::make_pair("C9_22", std::cref(C9_22)));
    ModelParamMap.insert(std::make_pair("C9p_22", std::cref(C9p_22)));
    ModelParamMap.insert(std::make_pair("C10_22", std::cref(C10_22)));
    ModelParamMap.insert(std::make_pair("C10p_22", std::cref(C10p_22)));
    ModelParamMap.insert(std::make_pair("CS_22", std::cref(CS_22)));
    ModelParamMap.insert(std::make_pair("CSp_22", std::cref(CSp_22)));
    ModelParamMap.insert(std::make_pair("CP_22", std::cref(CP_22)));
    ModelParamMap.insert(std::make_pair("CPp_22", std::cref(CPp_22)));
    
    ModelParamMap.insert(std::make_pair("WCscale", std::cref(WCscale)));
}

RealWeakEFTLFV::~RealWeakEFTLFV(){
    if (IsModelInitialized()) {
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool RealWeakEFTLFV::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return(true);
}
    
bool RealWeakEFTLFV::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool RealWeakEFTLFV::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool RealWeakEFTLFV::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool RealWeakEFTLFV::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and FlavourWC and FlavourWC-derived parameters in FlavourWCMatching */
    ReWEFTM.getObj().updateRealWeakEFTLFVParameters();

    return (true);
}

void RealWeakEFTLFV::setParameter(const std::string name, const double& value){    
    if(name.compare("C7") == 0) 
        C7 = value;
    else if(name.compare("C7p") == 0) 
        C7p = value;
    else if(name.compare("C8") == 0)
        C8 = value;
    else if(name.compare("C8p") == 0)
        C8p = value;
    else if(name.compare("C9_11") == 0)
        C9_11 = value;
    else if(name.compare("C9p_11") == 0)
        C9p_11 = value;
    else if(name.compare("C10_11") == 0)
        C10_11 = value;
    else if(name.compare("C10p_11") == 0)
        C10p_11 = value;
    else if(name.compare("CS_11") == 0)
        CS_11 = value;
    else if(name.compare("CSp_11") == 0)
        CSp_11 = value;
    else if(name.compare("CP_11") == 0)
        CP_11 = value;
    else if(name.compare("CPp_11") == 0)
        CPp_11 = value;
    else if(name.compare("C9_22") == 0)
        C9_22 = value;
    else if(name.compare("C9p_22") == 0)
        C9p_22 = value;
    else if(name.compare("C10_22") == 0)
        C10_22 = value;
    else if(name.compare("C10p_22") == 0)
        C10p_22 = value;
    else if(name.compare("CS_22") == 0)
        CS_22 = value;
    else if(name.compare("CSp_22") == 0)
        CSp_22 = value;
    else if(name.compare("CP_22") == 0)
        CP_22 = value;
    else if(name.compare("CPp_22") == 0)
        CPp_22 = value;
    else if(name.compare("WCscale") == 0)
        WCscale = value;
    else
        StandardModel::setParameter(name,value);
}

bool RealWeakEFTLFV::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NRealWeakEFTLFVvars; i++) {
        if (DPars.find(RealWeakEFTLFVvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory RealWeakEFTLFV parameter " << RealWeakEFTLFVvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(RealWeakEFTLFVvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool RealWeakEFTLFV::setFlag(const std::string name, const bool value)
{
    return StandardModel::setFlag(name,value);
}

