/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RealWeakEFTCC.h"

const std::string RealWeakEFTCC::RealWeakEFTCCvars[NRealWeakEFTCCvars] = {"CS", "CP", "CV", "CA", "CT"};
const std::string RealWeakEFTCC::RealWeakEFTCCPMvars[NRealWeakEFTCCvars] = {"CSL", "CSR", "CVL", "CVR", "CT"};

RealWeakEFTCC::RealWeakEFTCC(const unsigned int basis) : StandardModel(), basis(basis)
{
    if (basis == 1) {
        ModelParamMap.insert(std::make_pair("CSL", std::cref(C1)));
        ModelParamMap.insert(std::make_pair("CSR", std::cref(C2)));
        ModelParamMap.insert(std::make_pair("CVL", std::cref(C3)));
        ModelParamMap.insert(std::make_pair("CVR", std::cref(C4)));
        ModelParamMap.insert(std::make_pair("CT", std::cref(C5)));
    } else if (basis == 0) {
        ModelParamMap.insert(std::make_pair("CS", std::cref(C1)));
        ModelParamMap.insert(std::make_pair("CP", std::cref(C2)));
        ModelParamMap.insert(std::make_pair("CV", std::cref(C3)));
        ModelParamMap.insert(std::make_pair("CA", std::cref(C4)));
        ModelParamMap.insert(std::make_pair("CT", std::cref(C5)));
    } else
        throw std::runtime_error("Basis not defined in RealWeakEFTCC::RealWeakEFTCC()");
}

RealWeakEFTCC::~RealWeakEFTCC()
{
    if (IsModelInitialized()) {
    }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool RealWeakEFTCC::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return (true);
}

bool RealWeakEFTCC::Init(const std::map<std::string, double>& DPars)
{
    return (StandardModel::Init(DPars));
}

bool RealWeakEFTCC::PreUpdate()
{
    if (!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool RealWeakEFTCC::Update(const std::map<std::string, double>& DPars)
{

    if (!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if (!PostUpdate()) return (false);

    return (true);
}

bool RealWeakEFTCC::PostUpdate()
{
    if (!StandardModel::PostUpdate()) return (false);

    return (true);
}

void RealWeakEFTCC::setParameter(const std::string name, const double& value)
{
    if (basis == 1) {
        if (name.compare("CSL") == 0)
            C1 = value;
        else if (name.compare("CSR") == 0)
            C2 = value;
        else if (name.compare("CVL") == 0)
            C3 = value;
        else if (name.compare("CVR") == 0)
            C4 = value;
        else if (name.compare("CT") == 0)
            C5 = value;
        else
            StandardModel::setParameter(name, value);
    } else if (basis == 0) {
        if (name.compare("CS") == 0)
            C1 = value;
        else if (name.compare("CP") == 0)
            C2 = value;
        else if (name.compare("CV") == 0)
            C3 = value;
        else if (name.compare("CA") == 0)
            C4 = value;
        else if (name.compare("CT") == 0)
            C5 = value;
        else
            StandardModel::setParameter(name, value);
    } else 
        throw std::runtime_error("Basis not defined in RealWeakEFTCC::RealWeakEFTCC()");
}

bool RealWeakEFTCC::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NRealWeakEFTCCvars; i++) {
        if (basis == 1) {
            if (DPars.find(RealWeakEFTCCPMvars[i]) == DPars.end()) {
                std::cout << "ERROR: missing mandatory RealWeakEFTCCPM parameter " << RealWeakEFTCCPMvars[i] << std::endl;
                raiseMissingModelParameterCount();
                addMissingModelParameter(RealWeakEFTCCPMvars[i]);
            }
        } else if (basis == 0) {
            if (DPars.find(RealWeakEFTCCvars[i]) == DPars.end()) {
                std::cout << "ERROR: missing mandatory RealWeakEFTCC parameter " << RealWeakEFTCCvars[i] << std::endl;
                raiseMissingModelParameterCount();
                addMissingModelParameter(RealWeakEFTCCvars[i]);
            }
        }
    }
    return (StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool RealWeakEFTCC::setFlag(const std::string name, const bool value)
{
    return StandardModel::setFlag(name, value);
}

