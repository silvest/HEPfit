/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "NPHiggs.h"
#include <EWSM.h>


const std::string NPHiggs::NPHIGGSvars[NNPHIGGSvars] 
                  = {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4"};


NPHiggs::NPHiggs() 
: NPZbbbar() 
{
}


bool NPHiggs::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!NPZbbbar::Update(DPars)) return (false);

    return (true);
}


bool NPHiggs::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NPHiggs::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NNPHIGGSvars; i++) {
        if (DPars.find(NPHIGGSvars[i]) == DPars.end()) {
            std::cout << "missing mandatory NPHiggs parameter " 
                      << NPHIGGSvars[i] << std::endl;
            return false;
        }
    }
    return(NPZbbbar::CheckParameters(DPars));
}

    
void NPHiggs::SetParameter(const std::string name, const double& value) {
    if (name.compare("a") == 0)
        a = value;
    else if (name.compare("b") == 0)
        b = value;
    else if (name.compare("c_u") == 0)
        c_u = value;
    else if (name.compare("c_d") == 0)
        c_d = value;
    else if (name.compare("c_e") == 0)
        c_e = value;
    else if (name.compare("d_3") == 0)
        d_3 = value;
    else if (name.compare("d_4") == 0)
        d_4 = value;
    else
        NPZbbbar::SetParameter(name, value);       
}


bool NPHiggs::InitializeModel() {
    SetModelInitialized(NPZbbbar::InitializeModel());
    myEWepsilons = new EWepsilons(*this);
    return (IsModelInitialized());
}


void NPHiggs::SetEWSMflags(EWSM& myEWSM) 
{
    StandardModel::SetEWSMflags(myEWSM);
}


bool NPHiggs::SetFlag(const std::string name, const bool& value) {
    bool res = false;
    if (name.compare("EWBURGESS") == 0)
        throw std::runtime_error("Flag EWBURGESS is not applicable to NPHiggs"); 
    else if (name.compare("EWCHMN") == 0) 
        throw std::runtime_error("Flag EWCHMN is not applicable to NPHiggs"); 
    else if (name.compare("epsilon1SM") == 0) 
        throw std::runtime_error("Flag epsilon1SM is not applicable to NPHiggs"); 
    else if (name.compare("epsilon2SM") == 0) 
        throw std::runtime_error("Flag epsilon2SM is not applicable to NPHiggs"); 
    else if (name.compare("epsilon3SM") == 0) 
        throw std::runtime_error("Flag epsilon3SM is not applicable to NPHiggs"); 
    else if (name.compare("epsilonbSM") == 0) 
        throw std::runtime_error("Flag epsilonbSM is not applicable to NPHiggs"); 
    else {
        res = NPZbbbar::SetFlag(name,value);
    }
    return(res);
}


double NPHiggs::epsilon1() const{ 
    double Lambda;
    if (fabs(1.0-a*a) < pow(10.0, -32.0) ) 
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    double DeltaEps1 = - 3.0/16.0/M_PI*alphaMz()/c02()*(1.0 - a*a)
                       *log(Lambda*Lambda/mHl/mHl);
    return ( epsilon1_SM() + DeltaEps1 );
}


double NPHiggs::epsilon2() const {
    return epsilon2_SM();    
}
    

double NPHiggs::epsilon3() const {
    double Lambda;
    if (fabs(1.0-a*a) < pow(10.0, -32.0) ) 
        Lambda = pow(10.0, 19.0);
    else
        Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));

    double DeltaEps3 = 1.0/48.0/M_PI*alphaMz()/s02()*(1.0 - a*a)
                       *log(Lambda*Lambda/mHl/mHl);
    return ( epsilon3_SM() + DeltaEps3 );
}


double NPHiggs::epsilonb() const {
    //return epsilonb_SM();    
    return NPZbbbar::epsilonb();
}


////////////////////////////////////////////////////////////////////////     

double NPHiggs::Mw() const {
    return myEWepsilons->Mw(epsilon1(), epsilon2(), epsilon3());
}


double NPHiggs::cW2() const {
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NPHiggs::sW2() const {
    return ( 1.0 - cW2() );
}

    
complex NPHiggs::rhoZ_l(const StandardModel::lepton l) const {
    return myEWepsilons->rhoZ_l(l, epsilon1());
}

    
complex NPHiggs::rhoZ_q(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->rhoZ_q(q, epsilon1());
        case StandardModel::BOTTOM:
            return ( gAq(q)*gAq(q)/getQuarks(q).getIsospin()/getQuarks(q).getIsospin() );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPHiggs::rhoZ_q()");        
    }
}


complex NPHiggs::kappaZ_l(const StandardModel::lepton l) const {
    return myEWepsilons->kappaZ_l(l, epsilon1(), epsilon3());
}


complex NPHiggs::kappaZ_q(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->kappaZ_q(q, epsilon1(), epsilon3());
        case StandardModel::BOTTOM:
            return ( (1.0 - gVq(q)/gAq(q))/(4.0*fabs(getQuarks(q).getCharge())*sW2()) );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPHiggs::kappaZ_q()");        
    }
}
      
    
complex NPHiggs::gVl(const StandardModel::lepton l) const {
    return myEWepsilons->gVl(l, epsilon1(), epsilon3());
}


complex NPHiggs::gVq(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->gVq(q, epsilon1(), epsilon3());
        case StandardModel::BOTTOM:
            return ( myEWepsilons->gVb(epsilon1(), epsilon3(), epsilonb_SM()) + myDeltaGVb );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPHiggs::gVq()");        
    }
}


complex NPHiggs::gAl(const StandardModel::lepton l) const {
    return myEWepsilons->gAl(l, epsilon1());
}


complex NPHiggs::gAq(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return myEWepsilons->gAq(q, epsilon1());
        case StandardModel::BOTTOM:
            return ( myEWepsilons->gAb(epsilon1(), epsilonb_SM()) + myDeltaGAb );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NPHiggs::gAq()");        
    }
}

    
double NPHiggs::GammaW() const {
    throw std::runtime_error("NPHiggs::GammaW() is not implemented.");         
}
 


