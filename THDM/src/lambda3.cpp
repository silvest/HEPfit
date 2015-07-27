/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lambda3.h"
#include "StandardModel.h"

double lambda3::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->getTanb();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getSinb();
    double mHh=myTHDM->getMHh();
    double mHp=myTHDM->getMHp();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double m12_2=myTHDM->getM12_2();
    double lambda6=myTHDM->getLambda6();
    double lambda7=myTHDM->getLambda7();
    double lambda3 = ((mHh*mHh-mHl*mHl)*cosa*sina+2*mHp*mHp*sinb*cosb-m12_2)/(vev*vev*sinb*cosb)
                   -1./2.*lambda6/tanb-1./2.*lambda7*tanb;
    
    return lambda3;
}

