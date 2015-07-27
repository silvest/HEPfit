/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lambda2.h"
#include "StandardModel.h"

double lambda2::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double vev=myTHDM->v();
    double tanb=myTHDM->getTanb();
    double sinb=myTHDM->getSinb();
    double mHh=myTHDM->getMHh();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double m12_2=myTHDM->getM12_2();
    double lambda6=myTHDM->getLambda6();
    double lambda7=myTHDM->getLambda7();
    double lambda2 = (mHh*mHh*sina*sina+mHl*mHl*cosa*cosa-m12_2/tanb)/(vev*vev*sinb*sinb)
                   +1./2.*lambda6*pow(tanb,-3.) -3./2.*lambda7/tanb;
    
    return lambda2;
}
