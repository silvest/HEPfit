/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "positivity.h"
#include "StandardModel.h"

positivity::positivity(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    mylambda1=new lambda1(SM_i);
    mylambda2=new lambda2(SM_i);
};

positivity::~positivity()
{
  delete mylambda1;
  delete mylambda2;
};

double positivity::computeThValue()
{
    return 0.0;

}

positivity1::positivity1(const StandardModel& SM_i)
: positivity(SM_i)
{}

double positivity1::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();
    double lambda1=mylambda1->computeThValue();
    double lambda2=mylambda2->computeThValue();

    if (lambda1>0 && lambda2>0)
    {
        return 0.0001*(2.0*mHp2*cosb*sinb -m12_2 +(mHh2-mHl*mHl)*cosa*sina
                       +sqrt(std::abs((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)
                             *(mHh2*cosa*cosa +mHl*mHl*sina*sina -m12_2*tanb))));
    }
    else
    {
        return -1.;
    }
}

positivity2::positivity2(const StandardModel& SM_i)
: positivity(SM_i)
{}

double positivity2::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();
    double lambda1=mylambda1->computeThValue();
    double lambda2=mylambda2->computeThValue();

    if (lambda1>0 && lambda2>0)
    {
        return 0.0001*(mA2 -std::abs(mA2-m12_2/cosb/sinb)
                       +((mHh2-mHl*mHl)*cosa*sina
                         +sqrt(std::abs((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)
                               *(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb))))
                        /cosb/sinb);
    }
    else
    {
        return -1.;
    }
}
