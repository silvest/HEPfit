/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMquantitiesatQ.h"

Q_st::Q_st(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double Q_st::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Q_cutoff;
}


DeltaQ_THDM::DeltaQ_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double DeltaQ_THDM::computeThValue()
{
    return myTHDM.getQ_THDM() - myTHDM.getMyTHDMCache()->Q_cutoff;
}


g1atQ::g1atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double g1atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->g1_at_Q;
}



g2atQ::g2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double g2atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->g2_at_Q;
}



g3atQ::g3atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double g3atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->g3_at_Q;
}



YtopatQ::YtopatQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double YtopatQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Ytop_at_Q;
}



YbottomatQ::YbottomatQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double YbottomatQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Ybottom1_at_Q + myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
}



YtauatQ::YtauatQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double YtauatQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Ytau1_at_Q + myTHDM.getMyTHDMCache()->Ytau2_at_Q;
}



m11_2atQ::m11_2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double m11_2atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->m11_2_at_Q;
}



m22_2atQ::m22_2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double m22_2atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->m22_2_at_Q;
}



m12_2atQ::m12_2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double m12_2atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->m12_2_at_Q;
}


lambda1atQ::lambda1atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double lambda1atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->lambda1_at_Q;
}


lambda2atQ::lambda2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double lambda2atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->lambda2_at_Q;
}


lambda3atQ::lambda3atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double lambda3atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->lambda3_at_Q;
}


lambda4atQ::lambda4atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double lambda4atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->lambda4_at_Q;
}


lambda5atQ::lambda5atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM&> (SM_i))
{}

double lambda5atQ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->lambda5_at_Q;
}
