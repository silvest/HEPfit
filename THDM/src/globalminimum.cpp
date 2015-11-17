/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "globalminimum.h"
#include "StandardModel.h"

double globalminimum::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double tanb=myTHDM->gettanb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return 0.0000001*m12_2*(m12_2 +(mHl*mHl-mHh2)*cosa*sina
                                   +sqrt((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)
                                         *(mHh2*cosa*cosa +mHl*mHl*sina*sina -m12_2*tanb)));
}
