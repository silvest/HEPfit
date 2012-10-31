/* 
 * File:   NewPhysicsHiggs.cpp
 * Author: mishima
 */

#include <stdexcept>
#include "NewPhysicsHiggs.h"


const std::string NewPhysicsHiggs::NPHIGGSvars[NNPHIGGSvars] 
                  = {"a", "b", "c_u", "c_d", "c_e", "d_3", "d_4"};


NewPhysicsHiggs::NewPhysicsHiggs() : StandardModel() {
}


bool NewPhysicsHiggs::Update(const std::map<std::string,double>& DPars) {
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    if(!StandardModel::Update(DPars)) return (false);

    return (true);
}


bool NewPhysicsHiggs::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars)); 
}


bool NewPhysicsHiggs::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NNPHIGGSvars; i++) {
        if (DPars.find(NPHIGGSvars[i]) == DPars.end()) {
            std::cout << "missing mandatory NewPhysicsHiggs parameter " 
                      << NPHIGGSvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

    
void NewPhysicsHiggs::SetParameter(const std::string name, const double& value) {
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
        StandardModel::SetParameter(name, value);       
}


double NewPhysicsHiggs::epsilon1() const{ 
    double Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));
    double DeltaEps1 = - 3.0/16.0/M_PI*alphaMz()/c02()*(1.0 - a*a)
                       *log(Lambda*Lambda/mHl/mHl);

    return ( epsilon1_SM() + DeltaEps1 );
}


double NewPhysicsHiggs::epsilon2() const {
    return epsilon2_SM();    
}
    

double NewPhysicsHiggs::epsilon3() const {
    double Lambda = 4.0*M_PI*v()/sqrt(fabs(1.0 - a*a));
    double DeltaEps3 = 1.0/48.0/M_PI*alphaMz()/s02()*(1.0 - a*a)
                       *log(Lambda*Lambda/mHl/mHl);

    return ( epsilon3_SM() + DeltaEps3 );
}


double NewPhysicsHiggs::epsilonb() const {
    return epsilonb_SM();    
}


////////////////////////////////////////////////////////////////////////     

double NewPhysicsHiggs::Delta_rW() const {
    return (  (c02() - s02())/s02()
              *(epsilon2() - c02()*epsilon1() + 2.0*s02()*Delta_kappaPrime()) );
}


double NewPhysicsHiggs::Delta_kappaPrime() const {
    return ( (epsilon3() - c02()*epsilon1())/(c02() - s02()) );
}


////////////////////////////////////////////////////////////////////////     

double NewPhysicsHiggs::Mw() const {
    double Delta_r = 1.0 - (1.0 - DeltaAlpha())*(1.0 - Delta_rW());

    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    if (tmp/(1.0 - Delta_r) > 1.0) 
        throw std::runtime_error("Error in NewPhysicsHiggs::Mw()"); 
    
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp/(1.0 - Delta_r))) );
}


double NewPhysicsHiggs::cW2() const {
    return ( Mw()*Mw()/Mz/Mz );
}

    
double NewPhysicsHiggs::sW2() const {
    return ( 1.0 - cW2() );
}

    
complex NewPhysicsHiggs::rhoZ_l(const StandardModel::lepton l) const {
    double I3f = getLeptons(l).getIsospin();
    return ( gAl(l)*gAl(l)/I3f/I3f );
}

    
complex NewPhysicsHiggs::rhoZ_q(const StandardModel::quark q) const {
    double I3f = getQuarks(q).getIsospin();
    return ( gAq(q)*gAq(q)/I3f/I3f );    
}


complex NewPhysicsHiggs::kappaZ_l(const StandardModel::lepton l) const {
    double Qf = getLeptons(l).getCharge();
    return ( (1.0 - gVl(l)/gAl(l))/(4.0*fabs(Qf)*sW2()) );
}


complex NewPhysicsHiggs::kappaZ_q(const StandardModel::quark q) const {
    double Qf = getQuarks(q).getCharge();
    return ( (1.0 - gVq(q)/gAq(q))/(4.0*fabs(Qf)*sW2()) );    
}
      
    
complex NewPhysicsHiggs::gVl(const StandardModel::lepton l) const {
    double Qf = getLeptons(l).getCharge();
    return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02())*gAl(l) );    
}


complex NewPhysicsHiggs::gVq(const StandardModel::quark q) const {
    double Qf = getQuarks(q).getCharge();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02())*gAq(q) );
        case StandardModel::BOTTOM:
            return ( (1.0 - 4.0*fabs(Qf)*(1.0 + Delta_kappaPrime())*s02() + epsilonb())
                     /(1.0 + epsilonb())*gAq(q) );
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsHiggs::gVq()");        
    }
}


complex NewPhysicsHiggs::gAl(const StandardModel::lepton l) const {
    return complex( -(1.0 + epsilon1()/2.0)/2.0, 0.0, false);
}


complex NewPhysicsHiggs::gAq(const StandardModel::quark q) const {
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return complex( -(1.0 + epsilon1()/2.0)/2.0, 0.0, false);
        case StandardModel::BOTTOM:
            return complex( -(1.0 + epsilon1()/2.0)/2.0*(1.0 + epsilonb()), 0.0, false);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in NewPhysicsHiggs::gAq()");        
    }
}

    
//double NewPhysicsHiggs::GammaW() const {
//    
//}
 





