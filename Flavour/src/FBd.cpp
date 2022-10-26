/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FBd.h"
#include "StandardModel.h"

FBd::FBd(const StandardModel& SM_i) : ThObservable(SM_i) 
{
};

double FBd::computeThValue() 
{
    return(SM.getMesons(QCD::B_D).getDecayconst());
}
