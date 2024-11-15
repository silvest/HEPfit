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

double  RmBd::computeThValue()
{
    return SM.getFlavour().getDB2(0).getRB(FULLNLO).abs()-1.;
}

CBd::CBd(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(0);
    };

double  CBd::computeThValue() 
{
    return SM.getCBd();
}

PhiBd::PhiBd(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(0);
    };

double  PhiBd::computeThValue() 
{
    return SM.getPhiBd();
}