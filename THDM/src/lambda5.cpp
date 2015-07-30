/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lambda5.h"
#include "StandardModel.h"

double lambda5::computeThValue()
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
    double lambda6=myTHDM->getLambda6();
    double lambda7=myTHDM->getLambda7();
    double lambda5=(m12_2-mA*mA*sinb*cosb)/(vev*vev*sinb*cosb)
                   -1./2.*lambda6/tanb-1./2.*lambda7*tanb;
    
    return lambda5;
}

