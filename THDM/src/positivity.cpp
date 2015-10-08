/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "positivity.h"
#include "lambda1.h"
#include "lambda2.h"
#include "StandardModel.h"

positivity::positivity(const StandardModel& SM_i, int obsFlag) 
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    mylambda1=new lambda1(SM_i);
    mylambda2=new lambda2(SM_i);
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in positivity() called from ThFactory"
            "::ThFactory() can only be 1 (positivity1) or 2 (positivity2)");
};

positivity::~positivity()
{
  delete mylambda1;
  delete mylambda2;
};

double positivity::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA=myTHDM->getMA();
    double mHh=myTHDM->getMHh();
    double mHp=myTHDM->getMHp();
    double vev=myTHDM->v();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double tanb=myTHDM->getTanb();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getCosb();
    double m12_2=myTHDM->getM12_2();
    double lambda1=mylambda1->computeThValue();
    double lambda2=mylambda2->computeThValue();

    double positivity1=0.0001*(2*mHp*mHp*cosb*sinb -m12_2 +(mHh*mHh-mHl*mHl)*cosa*sina
                       +sqrt(std::abs((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)
                             *(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina -m12_2*tanb))));
    double positivity2=0.0001*(mA*mA -std::abs(mA*mA-m12_2/cosb/sinb)
                       +((mHh*mHh-mHl*mHl)*cosa*sina
                         +sqrt(std::abs((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)
                               *(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb))))
                        /cosb/sinb);
    if (lambda1 > 0 && lambda2>0)
    {
     if (obs == 1) return( positivity1);
     if (obs == 2) return( positivity2);
     throw std::runtime_error("positivity::computeThValue(): Observable type not "
            "defined. Can be only 1 or 2.");
    }
    else
    {
     return -1.;
    }
    return (EXIT_FAILURE);

}
