/* 
 * File:   NewPhysicsEpsilons.cpp
 * Author: mishima
 */

#include <stdexcept>
#include "NewPhysicsEpsilons.h"


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

double NewPhysicsEpsilons::Delta_rW() const {
    return (  (c02() - s02())/s02()
              *(myEpsilon_2 - c02()*myEpsilon_1 + 2.0*s02()*Delta_kappaPrime()) );
}


double NewPhysicsEpsilons::Delta_kappaPrime() const {
    return ( (myEpsilon_3 - c02()*myEpsilon_1)/(c02() - s02()) );
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
    double I3f = getLeptons(l).getIsospin();
    return ( gAl(l)*gAl(l)/I3f/I3f );
}

    
complex NewPhysicsEpsilons::rhoZ_q(const StandardModel::quark q) const {
    double I3f = getQuarks(q).getIsospin();
    return ( gAq(q)*gAq(q)/I3f/I3f );    
}


complex NewPhysicsEpsilons::kappaZ_l(const StandardModel::lepton l) const {
    double Qf = getLeptons(l).getCharge();
    return ( (1.0 - gVl(l)/gAl(l))/(4.0*fabs(Qf)*sW2()) );
}


complex NewPhysicsEpsilons::kappaZ_q(const StandardModel::quark q) const {
    double Qf = getQuarks(q).getCharge();
    return ( (1.0 - gVq(q)/gAq(q))/(4.0*fabs(Qf)*sW2()) );    
}
      
    
complex NewPhysicsEpsilons::gVl(const StandardModel::lepton l) const {
    double Qf = getLeptons(l).getCharge();
    return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02())*gAl(l) );    
}


complex NewPhysicsEpsilons::gVq(const StandardModel::quark q) const {
    double Qf = getQuarks(q).getCharge();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02())*gAq(q) );
        case StandardModel::BOTTOM:
            return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02() + myEpsilon_b)
                     /(1.0 + myEpsilon_b)*gAq(q) );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsEpsilons::gVq()");        
    }
}


complex NewPhysicsEpsilons::gAl(const StandardModel::lepton l) const {
    return complex( -(1.0 + myEpsilon_1/2.0)/2.0, 0.0, false);
}


complex NewPhysicsEpsilons::gAq(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return complex( -(1.0 + myEpsilon_1/2.0)/2.0, 0.0, false);
        case StandardModel::BOTTOM:
            return complex( -(1.0 + myEpsilon_1/2.0)/2.0*(1.0 + myEpsilon_b), 0.0, false);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsEpsilons::gAq()");        
    }
}

    
//double NewPhysicsEpsilons::GammaW() const {
//    
//}
 



