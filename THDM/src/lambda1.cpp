/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lambda1.h"
#include "StandardModel.h"

double lambda1::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->getTanb();
    double cosb=myTHDM->getCosb();
    double mHh=myTHDM->getMHh();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double m12_2=myTHDM->getM12_2();
    double lambda6=myTHDM->getLambda6();
    double lambda7=myTHDM->getLambda7();
    double lambda1 = (mHh*mHh*cosa*cosa+mHl*mHl*sina*sina-m12_2*tanb)/(vev*vev*cosb*cosb)
                   -3./2.*lambda6*tanb +lambda7*tanb*tanb*tanb/2.;

    return lambda1;
}

