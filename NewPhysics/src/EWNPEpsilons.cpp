/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWNPEpsilons.h"
#include "NPEpsilons.h"

EWNPEpsilons::EWNPEpsilons(const StandardModel& SM_i)
: SM(SM_i), EWSM(SM_i)
{
    myEWepsilons = new EWepsilons(SM_i);
}


////////////////////////////////////////////////////////////////////////

double EWNPEpsilons::Mw_NPEpsilons() const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    double eps2 = (static_cast<const NPEpsilons*> (&SM))->epsilon2();
    double eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
    return myEWepsilons->Mw(eps1, eps2, eps3);
}


////////////////////////////////////////////////////////////////////////

complex EWNPEpsilons::rhoZ_l(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    return myEWepsilons->rhoZ_l(l, eps1);
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
            return myEWepsilons->rhoZ_q(q, eps1);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return myEWepsilons->rhoZ_b(eps1, epsb);
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
    return myEWepsilons->kappaZ_l(l, eps1, eps3);
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
            return myEWepsilons->kappaZ_q(q, eps1, eps3);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return myEWepsilons->kappaZ_b(eps1, eps3, epsb);
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
    return myEWepsilons->gVl(l, eps1, eps3);
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
            return myEWepsilons->gVq(q, eps1, eps3);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            eps3 = (static_cast<const NPEpsilons*> (&SM))->epsilon3();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return myEWepsilons->gVb(eps1, eps3, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::gVq()");
    }
}


complex EWNPEpsilons::gAl(const StandardModel::lepton l) const
{
    double eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
    return myEWepsilons->gAl(l, eps1);
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
            return myEWepsilons->gAq(q, eps1);
        case StandardModel::BOTTOM:
            eps1 = (static_cast<const NPEpsilons*> (&SM))->epsilon1();
            epsb = (static_cast<const NPEpsilons*> (&SM))->epsilonb();
            return myEWepsilons->gAb(eps1, epsb);
        case StandardModel::TOP:
            return complex(0.0, 0.0, false);
        default:
            throw std::runtime_error("Error in EWNPEpsilons::gAq()");
    }
}





