/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWNPEpsilons.h"
#include "NPEpsilons.h"

EWNPEpsilons::EWNPEpsilons(const StandardModel& SM_i)
: EWSM(SM_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWNPEpsilons::Mw_NPEpsilons() const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    double eps2 = (static_cast<const NPEpsilons*> (&SM))->epsilon2();
    double eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
    return Mw(eps1, eps2, eps3);
}


////////////////////////////////////////////////////////////////////////

complex EWNPEpsilons::rhoZ_l(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    return rhoZ_l(l, eps1);
}


complex EWNPEpsilons::rhoZ_q(const StandardModel::quark q) const
{
    double eps1, epsb;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            return rhoZ_q(q, eps1);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return rhoZ_b(eps1, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::rhoZ_q()");
    }
}


complex EWNPEpsilons::kappaZ_l(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    double eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
    return kappaZ_l(l, eps1, eps3);
}


complex EWNPEpsilons::kappaZ_q(const StandardModel::quark q) const
{
    double eps1, eps3, epsb;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            return kappaZ_q(q, eps1, eps3);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return kappaZ_b(eps1, eps3, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::kappaZ_q()");
    }
}


complex EWNPEpsilons::gVl(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    double eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
    return gVl(l, eps1, eps3);
}


complex EWNPEpsilons::gVq(const StandardModel::quark q) const
{
    double eps1, eps3, epsb;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            return gVq(q, eps1, eps3);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return gVb(eps1, eps3, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::gVq()");
    }
}


complex EWNPEpsilons::gAl(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    return gAl(l, eps1);
}


complex EWNPEpsilons::gAq(const StandardModel::quark q) const
{
    double eps1, epsb;
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            return gAq(q, eps1);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return gAb(eps1, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::gAq()");
    }
}


////////////////////////////////////////////////////////////////////////

double EWNPEpsilons::Mw(const double eps1, const double eps2, const double eps3) const
{
    double Delta_r = 1.0 - (1.0 - SM.DeltaAlpha())*(1.0 - Delta_rW(eps1,eps2,eps3));

    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    if (tmp/(1.0 - Delta_r) > 1.0)
        throw std::runtime_error("Error in EWNPEpsilons::Mw()");

    return ( SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp/(1.0 - Delta_r))) );
}


complex EWNPEpsilons::rhoZ_l(const StandardModel::lepton l, const double eps1) const
{
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( rhoZ_e(eps1) );
    else
        return ( rhoZ_e(eps1) + SM.getEWSM()->rhoZ_l_SM_FlavorDep(l).real() );
}


complex EWNPEpsilons::rhoZ_q(const StandardModel::quark q, const double eps1) const
{
    if(q==SM.BOTTOM || q==SM.TOP)
        throw std::runtime_error("Error in EWNPEpsilons::rhoZ_q()");
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( rhoZ_e(eps1) );
    else
        return ( rhoZ_e(eps1) + SM.getEWSM()->rhoZ_q_SM_FlavorDep(q).real() );
}


complex EWNPEpsilons::kappaZ_l(const StandardModel::lepton l,
                             const double eps1, const double eps3) const
{
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( kappaZ_e(eps1,eps3) );
    else
        return ( kappaZ_e(eps1,eps3) + SM.getEWSM()->kappaZ_l_SM_FlavorDep(l).real() );
}


complex EWNPEpsilons::kappaZ_q(const StandardModel::quark q,
                             const double eps1, const double eps3) const
{
    if(q==SM.BOTTOM || q==SM.TOP)
        throw std::runtime_error("Error in EWNPEpsilons::kappaZ_q()");
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( kappaZ_e(eps1,eps3) );
    else
        return ( kappaZ_e(eps1,eps3) + SM.getEWSM()->kappaZ_q_SM_FlavorDep(q).real() );
}


complex EWNPEpsilons::gVl(const StandardModel::lepton l, const double eps1,
                        const double eps3) const
{
    double I3f = SM.getLeptons(l).getIsospin();
    double Qf = SM.getLeptons(l).getCharge();
    return ( sqrt(rhoZ_l(l,eps1).real())*I3f
             *(1.0 - 4.0*fabs(Qf)*kappaZ_l(l,eps1,eps3)*SM.sW2()) );
}


complex EWNPEpsilons::gVq(const StandardModel::quark q, const double eps1,
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
            throw std::runtime_error("Error in EWNPEpsilons::gVq()");
    }
}


complex EWNPEpsilons::gAl(const StandardModel::lepton l, const double eps1) const
{
    double I3f = SM.getLeptons(l).getIsospin();
    return ( complex(sqrt(rhoZ_l(l,eps1).real())*I3f, 0.0, false) );
}


complex EWNPEpsilons::gAq(const StandardModel::quark q, const double eps1) const
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
            throw std::runtime_error("Error in EWNPEpsilons::gAq()");
    }
}


complex EWNPEpsilons::rhoZ_b(const double eps1, const double epsb) const
{
    complex rhoZe = rhoZ_l(SM.ELECTRON, eps1);
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( rhoZe*(1.0 + epsb)*(1.0 + epsb) );
    else {
        double DeltaRhoZb = SM.getEWSM()->rhoZ_q_SM_FlavorDep(SM.BOTTOM).real();
        return ( (rhoZe + DeltaRhoZb)*(1.0 + epsb)*(1.0 + epsb) );
    }
}


complex EWNPEpsilons::kappaZ_b(const double eps1, const double eps3,
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


complex EWNPEpsilons::gVb(const double eps1, const double eps3,
                        const double epsb) const
{
    double Qb = SM.getQuarks(SM.BOTTOM).getCharge();
    double I3b = SM.getQuarks(SM.BOTTOM).getIsospin();
    return ( sqrt(rhoZ_b(eps1, epsb).real())*I3b
            *(1.0 - 4.0*fabs(Qb)*kappaZ_b(eps1,eps3,epsb)*SM.sW2()) );
}


complex EWNPEpsilons::gAb(const double eps1, const double epsb) const
{
    double I3b = SM.getQuarks(SM.BOTTOM).getIsospin();
    return ( sqrt(rhoZ_b(eps1,epsb).real())*I3b );
}


////////////////////////////////////////////////////////////////////////

double EWNPEpsilons::Delta_rW(const double eps1, const double eps2,
                            const double eps3) const
{
    double s02 = SM.getEWSM()->s02(), c02 = SM.getEWSM()->c02();
    return (  (c02 - s02)/s02
              *(eps2 - c02*eps1 + 2.0*s02*Delta_kappaPrime(eps1,eps3)) );
}


double EWNPEpsilons::Delta_kappaPrime(const double eps1, const double eps3) const
{
    double s02 = SM.getEWSM()->s02(), c02 = SM.getEWSM()->c02();
    return ( (eps3 - c02*eps1)/(c02 - s02) );
}


complex EWNPEpsilons::rhoZ_e(const double eps1) const
{
    return ( 4.0*gAe(eps1)*gAe(eps1) );
}


complex EWNPEpsilons::kappaZ_e(const double eps1, const double eps3) const
{
    return ( (1.0 - gVe(eps1,eps3)/gAe(eps1))/(4.0*SM.sW2()) );
}


complex EWNPEpsilons::gVe(const double eps1, const double eps3) const
{
    double s02 = SM.getEWSM()->s02();
    return ( (1.0 - 4.0*(1.0 + Delta_kappaPrime(eps1,eps3))*s02)*gAe(eps1) );
}


complex EWNPEpsilons::gAe(const double eps1) const
{
    return complex( -(1.0 + eps1/2.0)/2.0, 0.0, false);
}









