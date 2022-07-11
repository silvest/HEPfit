/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RealWETHeffDF1.h"

const std::string RealWETHeffDF1::RealWETHeffDF1vars[NRealWETHeffDF1vars] = {"C9_11",  "C9_22",  "C10_11",  "C10_22",
                                                                             "C9p_11", "C9p_22", "C10p_11", "C10p_22"};

RealWETHeffDF1::RealWETHeffDF1() : StandardModel()
{
    ModelParamMap.insert(std::make_pair("C9_11",   std::cref(C9_11)));
    ModelParamMap.insert(std::make_pair("C9_22",   std::cref(C9_22)));
    ModelParamMap.insert(std::make_pair("C10_11",  std::cref(C10_11)));
    ModelParamMap.insert(std::make_pair("C10_22",  std::cref(C10_22)));
    ModelParamMap.insert(std::make_pair("C9p_11",  std::cref(C9p_11)));
    ModelParamMap.insert(std::make_pair("C9p_22",  std::cref(C9p_22)));
    ModelParamMap.insert(std::make_pair("C10p_11", std::cref(C10p_11)));
    ModelParamMap.insert(std::make_pair("C10p_22", std::cref(C10p_22)));
}

RealWETHeffDF1::~RealWETHeffDF1()
{
    if (IsModelInitialized()) {
    }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool RealWETHeffDF1::InitializeModel()
{
    setModelInitialized(StandardModel::InitializeModel());
    return (true);
}

bool RealWETHeffDF1::Init(const std::map<std::string, double>& DPars)
{
    return (StandardModel::Init(DPars));
}

bool RealWETHeffDF1::PreUpdate()
{
    if (!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool RealWETHeffDF1::Update(const std::map<std::string, double>& DPars)
{

    if (!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if (!PostUpdate()) return (false);

    return (true);
}

bool RealWETHeffDF1::PostUpdate()
{
    if (!StandardModel::PostUpdate()) return (false);

    return (true);
}

void RealWETHeffDF1::setParameter(const std::string name, const double& value)
{
    if (name.compare("C9_11") == 0)
        C9_11 = value;
    else if (name.compare("C9_22") == 0)
        C9_22 = value;
    else if (name.compare("C10_11") == 0)
        C10_11 = value;
    else if (name.compare("C10_22") == 0)
        C10_22 = value;
    else if (name.compare("C9p_11") == 0)
        C9p_11 = value;
    else if (name.compare("C9p_22") == 0)
        C9p_22 = value;
    else if (name.compare("C10p_11") == 0)
        C10p_11 = value;
    else if (name.compare("C10p_22") == 0)
        C10p_22 = value;
    else
        StandardModel::setParameter(name, value);
}

bool RealWETHeffDF1::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NRealWETHeffDF1vars; i++) {
        if (DPars.find(RealWETHeffDF1vars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory RealWETHeffDF1 parameter " << RealWETHeffDF1vars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(RealWETHeffDF1vars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool RealWETHeffDF1::setFlag(const std::string name, const bool value)
{
    return StandardModel::setFlag(name, value);
}
