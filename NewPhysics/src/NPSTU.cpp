/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPSTU.h"


const std::string NPSTU::STUvars[NSTUvars] 
= {"obliqueS", "obliqueT", "obliqueU"};


NPSTU::NPSTU() 
: NPZbbbar() 
{
}


bool NPSTU::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!NPZbbbar::Update(DPars)) return (false);

    return (true);
}


bool NPSTU::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPSTU::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NSTUvars; i++) {
        if (DPars.find(STUvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTU parameter" 
                      << STUvars[i] << std::endl;
            return false;
        }
    }
    return(NPZbbbar::CheckParameters(DPars));
}

    
void NPSTU::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("obliqueS") == 0)
        myObliqueS = value;
    else if (name.compare("obliqueT") == 0)
        myObliqueT = value;
    else if (name.compare("obliqueU") == 0)
        myObliqueU = value;    
    else
        NPZbbbar::SetParameter(name, value);       
}


bool NPSTU::InitializeModel() 
{
    SetModelInitialized(NPZbbbar::InitializeModel());
    return (IsModelInitialized());
}


void NPSTU::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPSTU::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("EWABC") == 0)
        throw std::runtime_error("ERROR: Flag EWABC is not applicable to NPSTU"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("ERROR: Flag EWABC2 is not applicable to NPSTU"); 
    else
        res = NPZbbbar::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////

double NPSTU::epsilon1() const
{
    throw std::runtime_error("ERROR: NPSTU::epsilon1() is not implemented");
}


double NPSTU::epsilon2() const
{
    throw std::runtime_error("ERROR: NPSTU::epsilon2() is not implemented");
}


double NPSTU::epsilon3() const
{
    throw std::runtime_error("ERROR: NPSTU::epsilon3() is not implemented");
}


double NPSTU::epsilonb() const
{
    throw std::runtime_error("ERROR: NPSTU::epsilonb() is not implemented");
}


