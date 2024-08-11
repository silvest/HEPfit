/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBd.h"
#include "StandardModel.h"
#include "AmpDB2.h"
 
DmBd::DmBd(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(0);
    };

double  DmBd::computeThValue() 
{
    return(2. * SM.getCBd() * SM.getFlavour().getDB2(0).getM21(FULLNLO).abs());
}