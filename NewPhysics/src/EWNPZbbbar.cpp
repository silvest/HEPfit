/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPZbbbar.h"
#include "EWNPZbbbar.h"

EWNPZbbbar::EWNPZbbbar(const StandardModel& SM_i)
: SM(SM_i), EWSM(SM_i)
{
}


////////////////////////////////////////////////////////////////////////

complex EWNPZbbbar::rhoZ_l(const StandardModel::lepton l) const
{
    return rhoZ_l_SM(l);
}


complex EWNPZbbbar::rhoZ_q(const StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return rhoZ_q_SM(q);
        case StandardModel::BOTTOM:
            return ( gAq(q)*gAq(q)
                     /SM.getQuarks(q).getIsospin()/SM.getQuarks(q).getIsospin() );
        default:
            throw std::runtime_error("Error in EWNPZbbbar::rhoZ_q()");
    }
}


complex EWNPZbbbar::kappaZ_l(const StandardModel::lepton l) const
{
    return kappaZ_l_SM(l);
}


complex EWNPZbbbar::kappaZ_q(const StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return kappaZ_q_SM(q);
        case StandardModel::BOTTOM:
            return ( (1.0 - gVq(q)/gAq(q))
                     /(4.0*fabs(SM.getQuarks(q).getCharge())*SM.sW2()) );
        default:
            throw std::runtime_error("Error in EWNPZbbbar::kappaZ_q()");
    }
}


complex EWNPZbbbar::gVl(const StandardModel::lepton l) const
{
    return gVl_SM(l);
}


complex EWNPZbbbar::gVq(const StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return gVq_SM(q);
        case StandardModel::BOTTOM:
            if ((static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                return ( gVq_SM(q)
                         + (static_cast<const NPZbbbar*> (&SM))->deltaGVq(q) );
            else {
                double Qq = SM.getQuarks(q).getCharge();
                return ( gAq_SM(q)*(1.0 - 4.0*fabs(Qq)*(kappaZ_q_SM(q))*SM.sW2()) );
            }
        default:
            throw std::runtime_error("Error in EWNPZbbbar::gVq()");
    }
}


complex EWNPZbbbar::gAl(const StandardModel::lepton l) const
{
    return gAl_SM(l);
}


complex EWNPZbbbar::gAq(const StandardModel::quark q) const
{
    switch (q) {
        case StandardModel::UP:
        case StandardModel::DOWN:
        case StandardModel::CHARM:
        case StandardModel::STRANGE:
        case StandardModel::TOP:
            return gAq_SM(q);
        case StandardModel::BOTTOM:
            if ((static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
                return ( gAq_SM(q)
                         + (static_cast<const NPZbbbar*> (&SM))->deltaGAq(q) );
            else {
                double I3q = SM.getQuarks(q).getIsospin();
                return ( sqrt(rhoZ_q_SM(q))*I3q );
            }
        default:
            throw std::runtime_error("Error in EWNPZbbbar::gAq()");
    }
}


