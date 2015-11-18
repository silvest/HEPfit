/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMquantities.h"

THDMquantities::THDMquantities(const StandardModel& SM_i):
        ThObservable(SM_i),
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
{}

THDMquantities::~THDMquantities()
{}

double THDMquantities::computeThValue()
{
    return 0.0;
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

mass_mHh::mass_mHh(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double mass_mHh::computeThValue()
{
    return myTHDM->getmHh();
}


mass_mA::mass_mA(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double mass_mA::computeThValue()
{
    return myTHDM->getmA();
}


mass_mHp::mass_mHp(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double mass_mHp::computeThValue()
{
    return myTHDM->getmHp();
}


masssquare_mA::masssquare_mA(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double masssquare_mA::computeThValue()
{
    return myTHDM->getmA2();
}


masssquare_mHp::masssquare_mHp(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double masssquare_mHp::computeThValue()
{
    return myTHDM->getmHp2();
}


massdifference_mHhmmA::massdifference_mHhmmA(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mHhmmA::computeThValue()
{
    return myTHDM->getmHh() - myTHDM->getmA();
}


massdifference_mAmmHh::massdifference_mAmmHh(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mAmmHh::computeThValue()
{
    return myTHDM->getmA() - myTHDM->getmHh();
}


massdifference_mHhmmHp::massdifference_mHhmmHp(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mHhmmHp::computeThValue()
{
    return myTHDM->getmHh() - myTHDM->getmHp();
}


massdifference_mHpmmHh::massdifference_mHpmmHh(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mHpmmHh::computeThValue()
{
    return myTHDM->getmHp() - myTHDM->getmHh();
}


massdifference_mAmmHp::massdifference_mAmmHp(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mAmmHp::computeThValue()
{
    return myTHDM->getmA() - myTHDM->getmHp();
}


massdifference_mHpmmA::massdifference_mHpmmA(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double massdifference_mHpmmA::computeThValue()
{
    return myTHDM->getmHp() - myTHDM->getmA();
}


lambda1::lambda1(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double lambda1::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double cosb=myTHDM->getcosb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();
    double lambda6=myTHDM->getlambda6();
    double lambda7=myTHDM->getlambda7();

    return (mHh2*cosa*cosa+mHl*mHl*sina*sina-m12_2*tanb)/(vev*vev*cosb*cosb)
           -3./2.*lambda6*tanb +lambda7*tanb*tanb*tanb/2.;
}


lambda2::lambda2(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double lambda2::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();
    double lambda6=myTHDM->getlambda6();
    double lambda7=myTHDM->getlambda7();

    return (mHh2*sina*sina+mHl*mHl*cosa*cosa-m12_2/tanb)/(vev*vev*sinb*sinb)
           +1./2.*lambda6/(tanb*tanb*tanb) -3./2.*lambda7/tanb;
}


lambda3::lambda3(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double lambda3::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();
    double lambda6=myTHDM->getlambda6();
    double lambda7=myTHDM->getlambda7();

    return ((mHh2-mHl*mHl)*cosa*sina+2*mHp2*sinb*cosb-m12_2)/(vev*vev*sinb*cosb)
           -1./2.*lambda6/tanb-1./2.*lambda7*tanb;
}


lambda4::lambda4(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double lambda4::computeThValue()
{
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();
    double lambda6=myTHDM->getlambda6();
    double lambda7=myTHDM->getlambda7();

    return ((mA2-2*mHp2)*sinb*cosb+m12_2)/(vev*vev*sinb*cosb)
           -1./2.*lambda6/tanb-1./2.*lambda7*tanb;
}


lambda5::lambda5(const StandardModel& SM_i)
: THDMquantities(SM_i)
{}

double lambda5::computeThValue()
{
    double mA2=myTHDM->getmA2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();
    double lambda6=myTHDM->getlambda6();
    double lambda7=myTHDM->getlambda7();

    return (m12_2-mA2*sinb*cosb)/(vev*vev*sinb*cosb)
           -1./2.*lambda6/tanb-1./2.*lambda7*tanb;
}
