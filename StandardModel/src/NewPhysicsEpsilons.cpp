/* 
 * File:   NewPhysicsEpsilons.cpp
 * Author: mishima
 */

#include <stdexcept>
#include "NewPhysicsEpsilons.h"
#include "EWSM.h"


const std::string NewPhysicsEpsilons::EPSILONvars[NEPSILONvars] 
                  = {"epsilon_1", "epsilon_2", "epsilon_3", "epsilon_b"};


NewPhysicsEpsilons::NewPhysicsEpsilons() : StandardModel() {
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


////////////////////////////////////////////////////////////////////////     

double NewPhysicsEpsilons::Mw() const {
    double Delta_r = 1.0 - (1.0 - DeltaAlpha())*(1.0 - Delta_rW());

    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    if (tmp/(1.0 - Delta_r) > 1.0) 
        throw std::runtime_error("Error in NewPhysicsEpsilons::Mw()"); 
    
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp/(1.0 - Delta_r))) );
}


double NewPhysicsEpsilons::cW2() const {
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NewPhysicsEpsilons::sW2() const {
    return ( 1.0 - cW2() );
}

    
complex NewPhysicsEpsilons::rhoZ_l(const StandardModel::lepton l) const {
//    return ( rhoZ_e() );
    return ( rhoZ_e() + myEWSM->rhoZ_l_SM_FlavorDep(l) );
}

    
complex NewPhysicsEpsilons::rhoZ_q(const StandardModel::quark q) const {
//    return ( rhoZ_e() );
    return ( rhoZ_e() + myEWSM->rhoZ_q_SM_FlavorDep(q) );
}


complex NewPhysicsEpsilons::kappaZ_l(const StandardModel::lepton l) const {
//    return ( kappaZ_e() );
    return ( kappaZ_e() + myEWSM->kappaZ_l_SM_FlavorDep(l) );
}


complex NewPhysicsEpsilons::kappaZ_q(const StandardModel::quark q) const {
//    return ( kappaZ_e() );
    return ( kappaZ_e() + myEWSM->kappaZ_q_SM_FlavorDep(q) );
}
      
    
complex NewPhysicsEpsilons::gVl(const StandardModel::lepton l) const {
    double I3f = getLeptons(l).getIsospin();
    double Qf = getLeptons(l).getCharge();
    return ( sqrt(rhoZ_l(l).abs())*I3f*(1.0 - 4.0*fabs(Qf)*kappaZ_l(l)*s02()) );
}


complex NewPhysicsEpsilons::gVq(const StandardModel::quark q) const {
    double I3f = getQuarks(q).getIsospin();
    double Qf = getQuarks(q).getCharge();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( sqrt(rhoZ_q(q).abs())*I3f*(1.0 - 4.0*fabs(Qf)*kappaZ_q(q)*s02()) );
        case StandardModel::BOTTOM:
            return ( (1.0 - 4.0/3.0*(1.0 + Delta_kappaPrime())*s02() + myEpsilon_b)
                     /(1.0 + myEpsilon_b)*gAq(BOTTOM) );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsEpsilons::gVq()");        
    }
}


complex NewPhysicsEpsilons::gAl(const StandardModel::lepton l) const {
    double I3f = getLeptons(l).getIsospin();
    return ( complex(sqrt(rhoZ_l(l).abs())*I3f, 0.0, false) );
}


complex NewPhysicsEpsilons::gAq(const StandardModel::quark q) const {
    double I3f = getQuarks(q).getIsospin();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( sqrt(rhoZ_q(q).abs())*I3f );
        case StandardModel::BOTTOM:
            return ( gAe()*(1.0 + myEpsilon_b) );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsEpsilons::gAq()");        
    }
}

    
double NewPhysicsEpsilons::GammaW() const {
    throw std::runtime_error("NewPhysicsEpsilons::GammaW() is not implemented.");         
}
 

////////////////////////////////////////////////////////////////////////     

double NewPhysicsEpsilons::Delta_rW() const {
    return (  (c02() - s02())/s02()
              *(myEpsilon_2 - c02()*myEpsilon_1 + 2.0*s02()*Delta_kappaPrime()) );
}


double NewPhysicsEpsilons::Delta_kappaPrime() const {
    return ( (myEpsilon_3 - c02()*myEpsilon_1)/(c02() - s02()) );
}


complex NewPhysicsEpsilons::rhoZ_e() const {
    return ( 4.0*gAe()*gAe() );
}


complex NewPhysicsEpsilons::kappaZ_e() const {
    //return ( (1.0 - gVe()/gAe())/(4.0*sW2()) ); // wrong!  
    return ( (1.0 - gVe()/gAe())/(4.0*s02()) ); // corrected on Nov.2, 2012
}


complex NewPhysicsEpsilons::gVe() const {
    return ( (1.0 - 4.0*(1.0 + Delta_kappaPrime())*s02())*gAe() );
}


complex NewPhysicsEpsilons::gAe() const {
    return complex( -(1.0 + myEpsilon_1/2.0)/2.0, 0.0, false);    
}


