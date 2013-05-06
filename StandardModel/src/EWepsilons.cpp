/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "EWSM.h"
#include "EWepsilons.h"


////////////////////////////////////////////////////////////////////////     

double EWepsilons::Mw(const double eps1, const double eps2, const double eps3) const 
{
    double Delta_r = 1.0 - (1.0 - SM.DeltaAlpha())*(1.0 - Delta_rW(eps1,eps2,eps3));
    
    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    if (tmp/(1.0 - Delta_r) > 1.0) 
        throw std::runtime_error("Error in EWepsilons::Mw()"); 
    
    return ( SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp/(1.0 - Delta_r))) );
}


complex EWepsilons::rhoZ_l(const StandardModel::lepton l, const double eps1) const 
{
    if (SM.IsFlagWithoutNonUniversalVC())    
        return ( rhoZ_e(eps1) );
    else
        return ( rhoZ_e(eps1) + SM.getEWSM()->rhoZ_l_SM_FlavorDep(l).real() 
                 + SM.getEWSM()->delRhoZ_l(l) - SM.getEWSM()->delRhoZ_l(SM.ELECTRON) );
}

    
complex EWepsilons::rhoZ_q(const StandardModel::quark q, const double eps1) const 
{
    if(q==SM.BOTTOM || q==SM.TOP)
        throw std::runtime_error("Error in EWepsilons::rhoZ_q()"); 
    if (SM.IsFlagWithoutNonUniversalVC()) 
        return ( rhoZ_e(eps1) );
    else
        return ( rhoZ_e(eps1) + SM.getEWSM()->rhoZ_q_SM_FlavorDep(q).real() 
                 + SM.getEWSM()->delRhoZ_q(q) - SM.getEWSM()->delRhoZ_l(SM.ELECTRON) );
}


complex EWepsilons::kappaZ_l(const StandardModel::lepton l, 
                             const double eps1, const double eps3) const 
{
    if (SM.IsFlagWithoutNonUniversalVC()) 
        return ( kappaZ_e(eps1,eps3) );
    else
        return ( kappaZ_e(eps1,eps3) + SM.getEWSM()->kappaZ_l_SM_FlavorDep(l).real() );
}


complex EWepsilons::kappaZ_q(const StandardModel::quark q, 
                             const double eps1, const double eps3) const 
{
    if(q==SM.BOTTOM || q==SM.TOP)
        throw std::runtime_error("Error in EWepsilons::kappaZ_q()"); 
    if (SM.IsFlagWithoutNonUniversalVC()) 
        return ( kappaZ_e(eps1,eps3) );
    else
        return ( kappaZ_e(eps1,eps3) + SM.getEWSM()->kappaZ_q_SM_FlavorDep(q).real() );
}
      
    
complex EWepsilons::gVl(const StandardModel::lepton l, const double eps1, 
                        const double eps3) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();
    return ( sqrt(rhoZ_l(l,eps1).real())*I3f
             *(1.0 - 4.0*fabs(Qf)*kappaZ_l(l,eps1,eps3)*SM.sW2()) );
}


complex EWepsilons::gVq(const StandardModel::quark q, const double eps1, 
                        const double eps3) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    double Qf = SM.getQuarks(q).getCharge();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( sqrt(rhoZ_q(q,eps1).real())*I3f
                     *(1.0 - 4.0*fabs(Qf)*kappaZ_q(q,eps1,eps3)*SM.sW2()) );
        case StandardModel::BOTTOM:
        case StandardModel::TOP:
        default:
            throw std::runtime_error("Error in EWepsilons::gVq()");        
    }
}


complex EWepsilons::gAl(const StandardModel::lepton l, const double eps1) const 
{
    double I3f = SM.getLeptons(l).getIsospin();
    return ( complex(sqrt(rhoZ_l(l,eps1).real())*I3f, 0.0, false) );
}


complex EWepsilons::gAq(const StandardModel::quark q, const double eps1) const 
{
    double I3f = SM.getQuarks(q).getIsospin();
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            return ( sqrt(rhoZ_q(q,eps1).real())*I3f );
        case StandardModel::BOTTOM:
        case StandardModel::TOP:
        default:
            throw std::runtime_error("Error in EWepsilons::gAq()");        
    }
}


complex EWepsilons::rhoZ_b(const double eps1, const double epsb) const 
{
    complex rhoZe = rhoZ_l(SM.ELECTRON, eps1);
    if (SM.IsFlagWithoutNonUniversalVC()) 
        return ( rhoZe*(1.0 + epsb)*(1.0 + epsb) );
    else {            
        double DeltaRhoZb = SM.getEWSM()->rhoZ_q_SM_FlavorDep(SM.BOTTOM).real();
        return ( (rhoZe + DeltaRhoZb)*(1.0 + epsb)*(1.0 + epsb) 
                 + SM.getEWSM()->delRhoZ_q(SM.BOTTOM) 
                 - SM.getEWSM()->delRhoZ_l(SM.ELECTRON) );
    }
}


complex EWepsilons::kappaZ_b(const double eps1, const double eps3, 
                             const double epsb) const 
{
    complex kappaZe = kappaZ_l(SM.ELECTRON, eps1, eps3);
    if (SM.IsFlagWithoutNonUniversalVC()) 
        return ( kappaZe/(1.0 + epsb) );
    else {
        double DeltaKappaZb = SM.getEWSM()->kappaZ_q_SM_FlavorDep(SM.BOTTOM).real();
        return ( (kappaZe + DeltaKappaZb)/(1.0 + epsb) );
    }
}


complex EWepsilons::gVb(const double eps1, const double eps3, 
                        const double epsb) const 
{
    double Qb = SM.getQuarks(SM.BOTTOM).getCharge();
    double I3b = SM.getQuarks(SM.BOTTOM).getIsospin();
    return ( sqrt(rhoZ_b(eps1, epsb).real())*I3b
            *(1.0 - 4.0*fabs(Qb)*kappaZ_b(eps1,eps3,epsb)*SM.sW2()) );
}


complex EWepsilons::gAb(const double eps1, const double epsb) const 
{
    double I3b = SM.getQuarks(SM.BOTTOM).getIsospin();
    return ( sqrt(rhoZ_b(eps1,epsb).real())*I3b );
}


////////////////////////////////////////////////////////////////////////     

double EWepsilons::Delta_rW(const double eps1, const double eps2, 
                            const double eps3) const
{
    double s02 = SM.s02(), c02 = SM.c02();
    return (  (c02 - s02)/s02
              *(eps2 - c02*eps1 + 2.0*s02*Delta_kappaPrime(eps1,eps3)) );
}


double EWepsilons::Delta_kappaPrime(const double eps1, const double eps3) const 
{
    double s02 = SM.s02(), c02 = SM.c02();
    return ( (eps3 - c02*eps1)/(c02 - s02) );
}


complex EWepsilons::rhoZ_e(const double eps1) const
{
    return ( 4.0*gAe(eps1)*gAe(eps1) );
}


complex EWepsilons::kappaZ_e(const double eps1, const double eps3) const 
{
    return ( (1.0 - gVe(eps1,eps3)/gAe(eps1))/(4.0*SM.sW2()) );
}


complex EWepsilons::gVe(const double eps1, const double eps3) const 
{
    double s02 = SM.s02();
    return ( (1.0 - 4.0*(1.0 + Delta_kappaPrime(eps1,eps3))*s02)*gAe(eps1) );
}


complex EWepsilons::gAe(const double eps1) const 
{
    return complex( -(1.0 + eps1/2.0)/2.0, 0.0, false);    
}





