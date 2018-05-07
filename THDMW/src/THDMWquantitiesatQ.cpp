/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWquantitiesatQ.h"

Q_stTHDMW::Q_stTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double Q_stTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->Q_cutoff;
}


DeltaQ_THDMW::DeltaQ_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double DeltaQ_THDMW::computeThValue()
{
    return myTHDMW.getQ_THDMW() - myTHDMW.getMyTHDMWCache()->Q_cutoff;
}


lambda1atQTHDMW::lambda1atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double lambda1atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->lambda1_at_Q;
}

lambda2atQTHDMW::lambda2atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double lambda2atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->lambda2_at_Q;
}

lambda3atQTHDMW::lambda3atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double lambda3atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->lambda3_at_Q;
}

lambda4atQTHDMW::lambda4atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double lambda4atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->lambda4_at_Q;
}

mu1atQTHDMW::mu1atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu1atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu1_at_Q;
}

mu2atQTHDMW::mu2atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu2atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu2_at_Q;
}

mu3atQTHDMW::mu3atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu3atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu3_at_Q;
}

mu4atQTHDMW::mu4atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu4atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu4_at_Q;
}

mu5atQTHDMW::mu5atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu5atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu5_at_Q;
}

mu6atQTHDMW::mu6atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mu6atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mu6_at_Q;
}

nu1atQTHDMW::nu1atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double nu1atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->nu1_at_Q;
}

omega1atQTHDMW::omega1atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double omega1atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->omega1_at_Q;
}

kappa1atQTHDMW::kappa1atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double kappa1atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->kappa1_at_Q;
}

nu2atQTHDMW::nu2atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double nu2atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->nu2_at_Q;
}

omega2atQTHDMW::omega2atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double omega2atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->omega2_at_Q;
}

kappa2atQTHDMW::kappa2atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double kappa2atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->kappa2_at_Q;
}

nu4atQTHDMW::nu4atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double nu4atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->nu4_at_Q;
}

omega4atQTHDMW::omega4atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double omega4atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->omega4_at_Q;
}

nu3atQTHDMW::nu3atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double nu3atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->nu3_at_Q;
}

nu5atQTHDMW::nu5atQTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double nu5atQTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->nu5_at_Q;
}
