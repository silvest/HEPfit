/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "globalminimum.h"
#include "StandardModel.h"

double globalminimum::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double tanb=myTHDM->getTanb();
    double mHh=myTHDM->getMHh();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double m12_2=myTHDM->getM12_2();

    double globalminimum = 0.0000001*m12_2*(m12_2 +(mHl*mHl-mHh*mHh)*cosa*sina
                                   +sqrt((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)
                                         *(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina -m12_2*tanb)));

    return globalminimum;
}
