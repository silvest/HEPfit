/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NewPhysicsEpsilons.h"
#include "EWSM.h"

const std::string NewPhysicsEpsilons::EPSILONvars[NEPSILONvars] 
                  = {"epsilon_1", "epsilon_2", "epsilon_3", "epsilon_b"};

const std::string NewPhysicsEpsilons::EPSILONflags[NEPSILONflags] 
    = {"epsilon1SM", "epsilon2SM", "epsilon3SM", "epsilonbSM", 
       "withoutNonUniversalVCinEpsilons"};


NewPhysicsEpsilons::NewPhysicsEpsilons() : StandardModel() {
   FlagEpsilon1SM = false;
   FlagEpsilon2SM = false;
   FlagEpsilon3SM = false;
   FlagEpsilonbSM = false;
}


bool NewPhysicsEpsilons::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NewPhysicsEpsilons::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NewPhysicsEpsilons::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NEPSILONvars; i++) {
        if (DPars.find(EPSILONvars[i]) == DPars.end()) {
            std::cout << "missing mandatory NewPhysicsEpsilons parameter " 
                      << EPSILONvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NewPhysicsEpsilons::SetParameter(const std::string name, const double& value) {
    if (name.compare("epsilon_1") == 0)
        myEpsilon_1 = value;
    else if (name.compare("epsilon_2") == 0)
        myEpsilon_2 = value;
    else if (name.compare("epsilon_3") == 0)
        myEpsilon_3 = value;    
    else if (name.compare("epsilon_b") == 0)
        myEpsilon_b = value;    
    else
        StandardModel::SetParameter(name, value);       
}


bool NewPhysicsEpsilons::InitializeModel() {
    SetModelInitialized(StandardModel::InitializeModel());
    myEWepsilons = new EWepsilons(*this);
    return (IsModelInitialized());
}

void NewPhysicsEpsilons::SetEWSMflags(EWSM& myEWSM) {
    /* The flags below are used to compute the SM values of the epsilons. */
    myEWSM.setSchemeMw(EWSM::APPROXIMATEFORMULA);
    myEWSM.setSchemeRhoZ(EWSM::OMSI);
    myEWSM.setSchemeKappaZ(EWSM::APPROXIMATEFORMULA);
}


////////////////////////////////////////////////////////////////////////     

bool NewPhysicsEpsilons::SetFlag(const std::string name, const bool& value) {
    bool res = false;
    if (name.compare("EWBURGESS") == 0){
        throw std::runtime_error("Flag EWBURGESS is not applicable to NewPhysicsEpsilons"); 
    } else if (name.compare("EWCHMN") == 0) {
        throw std::runtime_error("Flag EWCHMN is not applicable to NewPhysicsEpsilons"); 
    } else if (name.compare("epsilon1SM") == 0) {
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
    } else if (name.compare("withoutNonUniversalVCinEpsilons") == 0) {
        myEWepsilons->setFlagWithoutNonUniversalVC(value);
        res = true;
    } else if (name.compare("R0bApproximate") == 0) {
        if (value) 
            throw std::runtime_error("R0bApproximate=true is not applicable to NewPhysicsEpsilons"); 
    } else {
        res = StandardModel::SetFlag(name,value);
    }
    return(res);
}


////////////////////////////////////////////////////////////////////////     

double NewPhysicsEpsilons::Mw() const {
    return myEWepsilons->Mw(epsilon1(), epsilon2(), epsilon3());
}


double NewPhysicsEpsilons::cW2() const {
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NewPhysicsEpsilons::sW2() const {
    return ( 1.0 - cW2() );
}

    
complex NewPhysicsEpsilons::rhoZ_l(const StandardModel::lepton l) const {
    return myEWepsilons->rhoZ_l(l, epsilon1());
}

    
complex NewPhysicsEpsilons::rhoZ_q(const StandardModel::quark q) const {
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
            throw std::runtime_error("Error in NewPhysicsEpsilons::rhoZ_q()");        
    }
}


complex NewPhysicsEpsilons::kappaZ_l(const StandardModel::lepton l) const {
    return myEWepsilons->kappaZ_l(l, epsilon1(), epsilon3());
}


complex NewPhysicsEpsilons::kappaZ_q(const StandardModel::quark q) const {
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
            throw std::runtime_error("Error in NewPhysicsEpsilons::kappaZ_q()");        
    }
}
      
    
complex NewPhysicsEpsilons::gVl(const StandardModel::lepton l) const {
    return myEWepsilons->gVl(l, epsilon1(), epsilon3());
}


complex NewPhysicsEpsilons::gVq(const StandardModel::quark q) const {
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
            throw std::runtime_error("Error in NewPhysicsEpsilons::gVq()");        
    }
}


complex NewPhysicsEpsilons::gAl(const StandardModel::lepton l) const {
    return myEWepsilons->gAl(l, epsilon1());
}


complex NewPhysicsEpsilons::gAq(const StandardModel::quark q) const {
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
            throw std::runtime_error("Error in NewPhysicsEpsilons::gAq()");        
    }
}

    
double NewPhysicsEpsilons::GammaW() const {
    throw std::runtime_error("NewPhysicsEpsilons::GammaW() is not implemented.");         
}
 
