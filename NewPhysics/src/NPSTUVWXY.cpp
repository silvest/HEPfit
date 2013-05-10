/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPSTUVWXY.h"


const std::string NPSTUVWXY::STUVWXYvars[NSTUVWXYvars] 
= {"obliqueShat", "obliqueThat", "obliqueUhat", 
   "obliqueV", "obliqueW", "obliqueX", "obliqueY"};


NPSTUVWXY::NPSTUVWXY()
: NPZbbbar() 
{
}


bool NPSTUVWXY::Update(const std::map<std::string,double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!NPZbbbar::Update(DPars)) return (false);

    return (true);
}


bool NPSTUVWXY::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPSTUVWXY::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NSTUVWXYvars; i++) {
        if (DPars.find(STUVWXYvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPSTUVWXY parameter" 
                      << STUVWXYvars[i] << std::endl;
            return false;
        }
    }
    return(NPZbbbar::CheckParameters(DPars));
}

    
void NPSTUVWXY::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("obliqueShat") == 0)
        myObliqueShat = value;
    else if (name.compare("obliqueThat") == 0)
        myObliqueThat = value;
    else if (name.compare("obliqueUhat") == 0)
        myObliqueUhat = value;    
    else if (name.compare("obliqueV") == 0)
        myObliqueV = value;    
    else if (name.compare("obliqueW") == 0)
        myObliqueW = value;    
    else if (name.compare("obliqueX") == 0)
        myObliqueX = value;    
    else if (name.compare("obliqueY") == 0)
        myObliqueY = value;    
    else
        NPZbbbar::SetParameter(name, value);       
}


bool NPSTUVWXY::InitializeModel() 
{
    SetModelInitialized(NPZbbbar::InitializeModel());
    return (IsModelInitialized());
}


void NPSTUVWXY::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPSTUVWXY::SetFlag(const std::string name, const bool& value)
{
    bool res = false;
    if (name.compare("EWABC") == 0)
        throw std::runtime_error("ERROR: Flag EWABC is not applicable to NPSTUVWXY"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("ERROR: Flag EWABC2 is not applicable to NPSTUVWXY"); 
    else
        res = NPZbbbar::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////

double NPSTUVWXY::epsilon1() const
{
    throw std::runtime_error("ERROR: NPSTUVWXY::epsilon1() is not implemented");
}


double NPSTUVWXY::epsilon2() const
{
    throw std::runtime_error("ERROR: NPSTUVWXY::epsilon2() is not implemented");
}


double NPSTUVWXY::epsilon3() const
{
    throw std::runtime_error("ERROR: NPSTUVWXY::epsilon3() is not implemented");
}


double NPSTUVWXY::epsilonb() const
{
    throw std::runtime_error("ERROR: NPSTUVWXY::epsilonb() is not implemented");
}



