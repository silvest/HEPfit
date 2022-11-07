/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FBdSqrtBBd.h"
#include "StandardModel.h"

FBdSqrtBBd::FBdSqrtBBd(const StandardModel& SM_i) : ThObservable(SM_i) 
{
    SM.initializeBParameter("BBd");
};

double FBdSqrtBBd::computeThValue() 
{
    return(SM.getMesons(QCD::B_D).getDecayconst()*sqrt(SM.getBBd().getBpars()(0)));
}
