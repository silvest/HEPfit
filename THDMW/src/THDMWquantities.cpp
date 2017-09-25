/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWquantities.h"

m12sqTHDMW::m12sqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double m12sqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->m12sq;
}

m11sqTHDMW::m11sqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double m11sqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->m11sq;
}

m22sqTHDMW::m22sqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double m22sqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->m22sq;
}

mhsqTHDMW::mhsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mhsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mhsq;
}

mHHsqTHDMW::mHHsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mHsq;
}

mAsqTHDMW::mAsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mAsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mAsq;
}

mSRsqTHDMW::mSRsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSRsq;
}

mSIsqTHDMW::mSIsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSIsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSIsq;
}

mHpsqTHDMW::mHpsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHpsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mHpsq;
}

mSpsqTHDMW::mSpsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSpsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSpsq;
}
