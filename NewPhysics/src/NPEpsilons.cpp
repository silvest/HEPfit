/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPEpsilons.h"
#include "EWNPEpsilons.h"

const std::string NPEpsilons::EPSILONvars[NEPSILONvars] 
= {"epsilon_1", "epsilon_2", "epsilon_3", "epsilon_b"};

const std::string NPEpsilons::EPSILONflags[NEPSILONflags] 
= {"epsilon1SM", "epsilon2SM", "epsilon3SM", "epsilonbSM"};


NPEpsilons::NPEpsilons() 
: NPbase()
{
    FlagEpsilon1SM = false;
    FlagEpsilon2SM = false;
    FlagEpsilon3SM = false;
    FlagEpsilonbSM = false;
}


bool NPEpsilons::InitializeModel()
{
    /* do not use setModelInitialized(NPbase::InitializeModel()); */
    myEWSM = new EWNPEpsilons(*this);
    this->setEWSMflags(*myEWSM);
    setModelInitialized(true);
    return(IsModelInitialized());
}


void NPEpsilons::setEWSMflags(EWSM& myEWSM)
{
    NPbase::setEWSMflags(myEWSM);
}


bool NPEpsilons::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}


bool NPEpsilons::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}

    
void NPEpsilons::setParameter(const std::string name, const double& value) 
{
    if (name.compare("epsilon_1") == 0)
        myEpsilon_1 = value;
    else if (name.compare("epsilon_2") == 0)
        myEpsilon_2 = value;
    else if (name.compare("epsilon_3") == 0)
        myEpsilon_3 = value;    
    else if (name.compare("epsilon_b") == 0)
        myEpsilon_b = value;    
    else
        NPbase::setParameter(name, value);
}


bool NPEpsilons::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NEPSILONvars; i++) {
        if (DPars.find(EPSILONvars[i]) == DPars.end()) {
            std::cout << "ERROR: Missing mandatory NPEpsilons parameter "
                      << EPSILONvars[i] << std::endl;
            return false;
        }
    }
    return(NPbase::CheckParameters(DPars));
}


bool NPEpsilons::setFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("epsilon1SM") == 0) {
        FlagEpsilon1SM = value;
        res = true;
    } else if (name.compare("epsilon2SM") == 0) {
        FlagEpsilon2SM = value;
        res = true;
    } else if (name.compare("epsilon3SM") == 0) {
        FlagEpsilon3SM = value;
        res = true;
    } else if (name.compare("epsilonbSM") == 0) {
        FlagEpsilonbSM = value;
        res = true;
    } else
        res = NPbase::setFlag(name,value);

    return(res);
}


bool NPEpsilons::CheckFlags() const
{
    if ( !IsFlagNoApproximateGqOverGb()
            && !IsFlagRhoZbFromGuOverGb() && !IsFlagRhoZbFromGdOverGb()
            && !IsFlagTestSubleadingTwoLoopEW())
        throw std::runtime_error("ERROR: The current flags are not compatible with NPEpsilons model.");
    
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////

double NPEpsilons::Mw() const 
{
    return (static_cast<const EWNPEpsilons*> (myEWSM))->Mw_NPEpsilons();
}


double NPEpsilons::cW2() const 
{
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NPEpsilons::sW2() const 
{
    return ( 1.0 - cW2() );
}


double NPEpsilons::GammaW() const
{
    throw std::runtime_error("NPEpsilons::GammaW() is not implemented.");
}


