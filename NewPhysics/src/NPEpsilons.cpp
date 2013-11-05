/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPEpsilons.h"

const std::string NPEpsilons::EPSILONvars[NEPSILONvars] 
= {"epsilon_1", "epsilon_2", "epsilon_3", "epsilon_b"};

const std::string NPEpsilons::EPSILONflags[NEPSILONflags] 
= {"epsilon1SM", "epsilon2SM", "epsilon3SM", "epsilonbSM", "EWABC", "EWABC2"};


NPEpsilons::NPEpsilons() 
: NPbase()
{
    FlagEpsilon1SM = false;
    FlagEpsilon2SM = false;
    FlagEpsilon3SM = false;
    FlagEpsilonbSM = false;
    FlagEWABC = false;
    FlagEWABC2 = false;
}


bool NPEpsilons::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if(!NPbase::Update(DPars)) return (false);
    return (true);
}


bool NPEpsilons::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
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


bool NPEpsilons::InitializeModel() 
{
    setModelInitialized(NPbase::InitializeModel());
    myEWepsilons = new EWepsilons(*this);
    return (IsModelInitialized());
}


void NPEpsilons::setEWSMflags(EWSM& myEWSM) 
{
    NPbase::setEWSMflags(myEWSM);
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
    } else if (name.compare("EWABC") == 0) {
        FlagEWABC = value;
        res = true;
    } else if (name.compare("EWABC2") == 0) {
        FlagEWABC2 = value;
        res = true;
    } else
        res = NPbase::setFlag(name,value);

    return(res);
}

bool NPEpsilons::CheckFlags() const
{
    if ( IsFlagApproximateGqOverGb()
            && !IsFlagRhoZbFromGuOverGb() && !IsFlagRhoZbFromGdOverGb()
            && !IsFlagTestSubleadingTwoLoopEW())
        throw std::runtime_error("ERROR: The current flags are not compatible with NPEpsilons model.");

    if ( FlagEWABC && FlagEWABC2 )
        throw std::runtime_error("ERROR: Flags EWABC and EWABC2 are incompatible with each other.");
    
    return(NPbase::CheckFlags());
}


////////////////////////////////////////////////////////////////////////     

double NPEpsilons::Mw() const 
{
    return myEWepsilons->Mw(epsilon1(), epsilon2(), epsilon3());
}


double NPEpsilons::cW2() const 
{
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NPEpsilons::sW2() const 
{
    return ( 1.0 - cW2() );
}

    
complex NPEpsilons::rhoZ_l(const StandardModel::lepton l) const 
{
    return myEWepsilons->rhoZ_l(l, epsilon1());
}

    
complex NPEpsilons::rhoZ_q(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->rhoZ_q(q, epsilon1());
        case StandardModel::BOTTOM:
            return myEWepsilons->rhoZ_b(epsilon1(), epsilonb());
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPEpsilons::rhoZ_q()");        
    }
}


complex NPEpsilons::kappaZ_l(const StandardModel::lepton l) const 
{
    return myEWepsilons->kappaZ_l(l, epsilon1(), epsilon3());
}


complex NPEpsilons::kappaZ_q(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->kappaZ_q(q, epsilon1(), epsilon3());
        case StandardModel::BOTTOM:
            return myEWepsilons->kappaZ_b(epsilon1(), epsilon3(), epsilonb());
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPEpsilons::kappaZ_q()");        
    }
}
      
    
complex NPEpsilons::gVl(const StandardModel::lepton l) const 
{
    return myEWepsilons->gVl(l, epsilon1(), epsilon3());
}


complex NPEpsilons::gVq(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->gVq(q, epsilon1(), epsilon3());
        case StandardModel::BOTTOM:
            return myEWepsilons->gVb(epsilon1(), epsilon3(), epsilonb());
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPEpsilons::gVq()");        
    }
}


complex NPEpsilons::gAl(const StandardModel::lepton l) const 
{
    return myEWepsilons->gAl(l, epsilon1());
}


complex NPEpsilons::gAq(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->gAq(q, epsilon1());
        case StandardModel::BOTTOM:
            return myEWepsilons->gAb(epsilon1(), epsilonb());
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPEpsilons::gAq()");        
    }
}

    
double NPEpsilons::GammaW() const 
{
    throw std::runtime_error("NPEpsilons::GammaW() is not implemented.");         
}
 
