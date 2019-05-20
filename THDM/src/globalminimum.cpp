/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "globalminimum.h"
#include "THDM.h"

globalminimum::globalminimum(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
};

double globalminimum::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double tanb=myTHDM->gettanb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return 0.0000001*m12_2*(m12_2 +(mHl2-mHh2)*cosa*sina
                                   +sqrt((mHl2*cosa*cosa +mHh2*sina*sina -m12_2/tanb)
                                         *(mHh2*cosa*cosa +mHl2*sina*sina -m12_2*tanb)));
}
