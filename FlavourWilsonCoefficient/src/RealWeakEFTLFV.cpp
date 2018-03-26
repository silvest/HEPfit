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
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C7", boost::cref(C7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C7p", boost::cref(C7p)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C8", boost::cref(C8)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C8p", boost::cref(C8p)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C9_11", boost::cref(C9_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C9p_11", boost::cref(C9p_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C10_11", boost::cref(C10_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C10p_11", boost::cref(C10p_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CS_11", boost::cref(CS_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CSp_11", boost::cref(CSp_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CP_11", boost::cref(CP_11)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CPp_11", boost::cref(CPp_11)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C9_22", boost::cref(C9_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C9p_22", boost::cref(C9p_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C10_22", boost::cref(C10_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("C10p_22", boost::cref(C10p_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CS_22", boost::cref(CS_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CSp_22", boost::cref(CSp_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CP_22", boost::cref(CP_22)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("CPp_22", boost::cref(CPp_22)));
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("WCscale", boost::cref(WCscale)));
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
            std::cout << "missing mandatory FlavourWilsonCoefficient parameter " << RealWeakEFTLFVvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool RealWeakEFTLFV::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}

