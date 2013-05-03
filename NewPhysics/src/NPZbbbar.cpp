/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "NPZbbbar.h"


const std::string NPZbbbar::ZbbbarVars[NZbbbarVars] 
                  = {"deltaGVb", "deltaGAb"};


NPZbbbar::NPZbbbar() 
: StandardModel() 
{
}


bool NPZbbbar::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NPZbbbar::Init(const std::map<std::string, double>& DPars) 
{
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPZbbbar::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NZbbbarVars; i++) {
        if (DPars.find(ZbbbarVars[i]) == DPars.end()) {
            std::cout << "missing mandatory NPZbbbar parameter " 
                      << ZbbbarVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NPZbbbar::SetParameter(const std::string name, const double& value) 
{
    if (name.compare("deltaGVb") == 0)
        myDeltaGVb = value;
    else if (name.compare("deltaGAb") == 0)
        myDeltaGAb = value;
    else
        StandardModel::SetParameter(name, value);       
}


bool NPZbbbar::InitializeModel() 
{
    SetModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


void NPZbbbar::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


////////////////////////////////////////////////////////////////////////     

bool NPZbbbar::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("EWABC") == 0) 
        throw std::runtime_error("Flag EWABC is not applicable to NPSTU"); 
    else if (name.compare("EWABC2") == 0)
        throw std::runtime_error("Flag EWABC2 is not applicable to NPSTU"); 
    else if (name.compare("EWBURGESS") == 0) 
        throw std::runtime_error("Flag EWBURGESS is not applicable to NPZbbbar"); 
    else if (name.compare("EWCHMN") == 0) 
        throw std::runtime_error("Flag EWCHMN is not applicable to NPZbbbar"); 
    else if (name.compare("epsilon1SM") == 0) 
        throw std::runtime_error("Flag epsilon1SM is not applicable to NPZbbbar"); 
    else if (name.compare("epsilon2SM") == 0) 
        throw std::runtime_error("Flag epsilon2SM is not applicable to NPZbbbar"); 
    else if (name.compare("epsilon3SM") == 0) 
        throw std::runtime_error("Flag epsilon3SM is not applicable to NPZbbbar"); 
    else if (name.compare("epsilonbSM") == 0) 
        throw std::runtime_error("Flag epsilonbSM is not applicable to NPZbbbar"); 
    else {
        res = StandardModel::SetFlag(name,value);
    }
    return(res);
}


////////////////////////////////////////////////////////////////////////     

    
complex NPZbbbar::rhoZ_l(const StandardModel::lepton l) const 
{
    return StandardModel::rhoZ_l(l);
}

    
complex NPZbbbar::rhoZ_q(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return StandardModel::rhoZ_q(q);
        case StandardModel::BOTTOM:
            return ( gAq(q)*gAq(q)/getQuarks(q).getIsospin()/getQuarks(q).getIsospin() );
        default:
            throw std::runtime_error("Error in NPZbbbar::rhoZ_q()");        
    }
}


complex NPZbbbar::kappaZ_l(const StandardModel::lepton l) const 
{
    return StandardModel::kappaZ_l(l);
}


complex NPZbbbar::kappaZ_q(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return StandardModel::kappaZ_q(q);
        case StandardModel::BOTTOM:
            return ( (1.0 - gVq(q)/gAq(q))/(4.0*fabs(getQuarks(q).getCharge())*sW2()) );
        default:
            throw std::runtime_error("Error in NPZbbbar::kappaZ_q()");        
    }
}
      
    
complex NPZbbbar::gVl(const StandardModel::lepton l) const 
{
    return StandardModel::gVl(l);
}


complex NPZbbbar::gVq(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return StandardModel::gVq(q);
        case StandardModel::BOTTOM:
            if (IsFlagNPZbbbarLinearize())
                return StandardModel::gVq(q);            
            else
                return ( StandardModel::gVq(q) + myDeltaGVb );
        default:
            throw std::runtime_error("Error in NPZbbbar::gVq()");        
    }
}


complex NPZbbbar::gAl(const StandardModel::lepton l) const 
{
    return StandardModel::gAl(l);
}


complex NPZbbbar::gAq(const StandardModel::quark q) const 
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return StandardModel::gAq(q);
        case StandardModel::BOTTOM:
            if (IsFlagNPZbbbarLinearize())
                return StandardModel::gAq(q);
            else
                return ( StandardModel::gAq(q) + myDeltaGAb );
        default:
            throw std::runtime_error("Error in NPZbbbar::gAq()");        
    }
}


double NPZbbbar::epsilon1() const
{ 
    return epsilon1_SM();
}


double NPZbbbar::epsilon2() const 
{
    return epsilon2_SM();    
}


double NPZbbbar::epsilon3() const 
{
    return epsilon3_SM();
}


double NPZbbbar::epsilonb() const 
{
    complex kappaZe = kappaZ_l(ELECTRON);
    complex kappaZb = kappaZ_q(BOTTOM);
    if (IsFlagWithoutNonUniversalVC()) 
        return ( kappaZe.real()/kappaZb.real() - 1.0 ); 
    else 
        return ( (kappaZe.real() + myEWSM->kappaZ_q_SM_FlavorDep(BOTTOM).real())
                 /kappaZb.real() - 1.0 );   
}



