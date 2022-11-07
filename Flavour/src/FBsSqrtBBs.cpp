/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FBsSqrtBBs.h"
#include "StandardModel.h"

FBsSqrtBBs::FBsSqrtBBs(const StandardModel& SM_i) : ThObservable(SM_i) 
{
    SM.initializeBParameter("BBs");
};

double FBsSqrtBBs::computeThValue() 
{
    return(SM.getMesons(QCD::B_S).getDecayconst()*sqrt(SM.getBBs().getBpars()(0)));
}
