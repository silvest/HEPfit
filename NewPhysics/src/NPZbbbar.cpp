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


const std::string NPZbbbar::Zbbbarflags[NZbbbarflags]
= {"NPZbbbarLR",  "NotLinearizedNP"};


NPZbbbar::NPZbbbar() 
: StandardModel() 
{
    FlagNPZbbbarLR = false;
    FlagNotLinearizedNP = false;
}


bool NPZbbbar::Update(const std::map<std::string,double>& DPars) 
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameters(it->first, it->second);
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
            std::cout << "ERROR: Missing mandatory NPZbbbar parameter" 
                      << ZbbbarVars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NPZbbbar::setParameters(const std::string name, const double& value) 
{
    if (name.compare("deltaGVb") == 0)
        myDeltaGVb = value;
    else if (name.compare("deltaGAb") == 0)
        myDeltaGAb = value;
    else
        StandardModel::setParameters(name, value);       
}


bool NPZbbbar::InitializeModel() 
{
    setModelInitialized(StandardModel::InitializeModel());
    return (IsModelInitialized());
}


void NPZbbbar::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPZbbbar::SetFlag(const std::string name, const bool& value) 
{
    bool res = false;
    if (name.compare("NPZbbbarLR") == 0) {
        FlagNPZbbbarLR = value;
        res = true;
    } else if (name.compare("NotLinearizedNP") == 0) {
        FlagNotLinearizedNP = value;
        res = true;
    } else
        res = StandardModel::SetFlag(name,value);

    return(res);
}


////////////////////////////////////////////////////////////////////////     


double NPZbbbar::deltaGVl(StandardModel::lepton l) const
{
    return StandardModel::deltaGVl(l);
}


double NPZbbbar::deltaGVq(StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
        case StandardModel::TOP:
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            return StandardModel::deltaGVq(q);
        case StandardModel::BOTTOM:
            if (FlagNPZbbbarLR)
                // delta g_L^b + delta g_R^b
                return ( myDeltaGVb + StandardModel::deltaGVq(q)
                         + myDeltaGAb + StandardModel::deltaGAq(q)); 
            else
                return ( myDeltaGVb + StandardModel::deltaGVq(q) );
        default:
            throw std::runtime_error("Error in NPZbbbar::deltaGVq()");
    }
}


double NPZbbbar::deltaGAl(StandardModel::lepton l) const
{
    return StandardModel::deltaGAl(l);
}


 double NPZbbbar::deltaGAq(StandardModel::quark q) const
 {
     switch (q) {
         case StandardModel::UP:
         case StandardModel::CHARM:
         case StandardModel::TOP:
         case StandardModel::DOWN:
         case StandardModel::STRANGE:
             return StandardModel::deltaGAq(q);
         case StandardModel::BOTTOM:
             if (FlagNPZbbbarLR)
                // delta g_L^b - delta g_R^b
                return ( myDeltaGVb + StandardModel::deltaGVq(q)
                         - myDeltaGAb - StandardModel::deltaGAq(q));
             else
                 return ( myDeltaGAb + StandardModel::deltaGAq(q) );
         default:
             throw std::runtime_error("Error in NPZbbbar::deltaGAq()");
     }
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
            if (IsFlagNotLinearizedNP())
                return ( StandardModel::gVq(q) + deltaGVq(q) );
            else
                return StandardModel::gVq(q);
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
            if (IsFlagNotLinearizedNP())
                return ( StandardModel::gAq(q) + deltaGAq(q) );
            else
                return StandardModel::gAq(q);
        default:
            throw std::runtime_error("Error in NPZbbbar::gAq()");        
    }
}


double NPZbbbar::epsilonb() const 
{
    complex kappaZe = kappaZ_l(ELECTRON);
    complex kappaZb;
    if (IsFlagNotLinearizedNP())
        kappaZb = kappaZ_q(BOTTOM); /* In this case, kappaZ_q(BOTTOM) includes deltaGVb and deltaGAb. */
    else {
        complex gVb = StandardModel::gVq(BOTTOM) + deltaGVq(BOTTOM);
        complex gAb = StandardModel::gAq(BOTTOM) + deltaGAq(BOTTOM);
        kappaZb = (1.0 - gVb/gAb)/(4.0*fabs(getQuarks(BOTTOM).getCharge())*sW2());
    }
    
    if (IsFlagWithoutNonUniversalVC()) 
        return ( kappaZe.real()/kappaZb.real() - 1.0 ); 
    else 
        return ( (kappaZe.real() + myEWSM->kappaZ_q_SM_FlavorDep(BOTTOM).real())
                /kappaZb.real() - 1.0 );   
}



