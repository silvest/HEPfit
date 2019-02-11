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

mhTHDMW::mhTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mhTHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mhsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mHHsqTHDMW::mHHsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mHsq;
}

mHHTHDMW::mHHTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHTHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mHsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mAsqTHDMW::mAsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mAsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mAsq;
}

mATHDMW::mATHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mATHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mAsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mSRsqTHDMW::mSRsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSRsq;
}

mSRTHDMW::mSRTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRTHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mSRsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mSIsqTHDMW::mSIsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSIsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSIsq;
}

mSITHDMW::mSITHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSITHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mSIsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mHpsqTHDMW::mHpsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHpsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mHpsq;
}

mHpTHDMW::mHpTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHpTHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mHpsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mSpsqTHDMW::mSpsqTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSpsqTHDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->mSpsq;
}

mSpTHDMW::mSpTHDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSpTHDMW::computeThValue()
{
    double msq=myTHDMW.getMyTHDMWCache()->mSpsq;
    if (msq<0.)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    else
    {
        return sqrt(msq);
    }
}

mAmmHH_THDMW::mAmmHH_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mAmmHH_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mAsq) - sqrt(myTHDMW.getMyTHDMWCache()->mHsq);
}

mHHmmA_THDMW::mHHmmA_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHmmA_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mHsq) - sqrt(myTHDMW.getMyTHDMWCache()->mAsq);
}

mAmmSR_THDMW::mAmmSR_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mAmmSR_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mAsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSRsq);
}

mSRmmA_THDMW::mSRmmA_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRmmA_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSRsq) - sqrt(myTHDMW.getMyTHDMWCache()->mAsq);
}

mAmmSI_THDMW::mAmmSI_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mAmmSI_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mAsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSIsq);
}

mSImmA_THDMW::mSImmA_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSImmA_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSIsq) - sqrt(myTHDMW.getMyTHDMWCache()->mAsq);
}

mHHmmSR_THDMW::mHHmmSR_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHmmSR_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mHsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSRsq);
}

mSRmmHH_THDMW::mSRmmHH_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRmmHH_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSRsq) - sqrt(myTHDMW.getMyTHDMWCache()->mHsq);
}

mHHmmSI_THDMW::mHHmmSI_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mHHmmSI_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mHsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSIsq);
}

mSImmHH_THDMW::mSImmHH_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSImmHH_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSIsq) - sqrt(myTHDMW.getMyTHDMWCache()->mHsq);
}

mSRmmSI_THDMW::mSRmmSI_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRmmSI_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSRsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSIsq);
}

mSImmSR_THDMW::mSImmSR_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSImmSR_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSIsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSRsq);
}


mSpmmSI_THDMW::mSpmmSI_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSpmmSI_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSpsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSIsq);
}

mSRmmSp_THDMW::mSRmmSp_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSRmmSp_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSRsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSpsq);
}



mSImmSp_THDMW::mSImmSp_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSImmSp_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSIsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSpsq);
}


mSpmmSR_THDMW::mSpmmSR_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double mSpmmSR_THDMW::computeThValue()
{
    return sqrt(myTHDMW.getMyTHDMWCache()->mSpsq) - sqrt(myTHDMW.getMyTHDMWCache()->mSRsq);
}




rh_gg_THDMW::rh_gg_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double rh_gg_THDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->rh_gg;
}

rh_gaga_THDMW::rh_gaga_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double rh_gaga_THDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->rh_gaga;
}

rh_Zga_THDMW::rh_Zga_THDMW(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double rh_Zga_THDMW::computeThValue()
{
    return myTHDMW.getMyTHDMWCache()->rh_Zga;
}
