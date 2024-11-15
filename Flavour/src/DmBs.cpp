/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBs.h"
#include "StandardModel.h"
#include "AmpDB2.h"

DmBs::DmBs(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(1);
    };

double  DmBs::computeThValue()
{
    return(2. * SM.getCBs() * SM.getFlavour().getDB2(1).getM21(FULLNLO).abs());
}

double  RmBs::computeThValue()
{
    return SM.getFlavour().getDB2(1).getRB(FULLNLO).abs()-1.;
}

CBs::CBs(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(1);
    };

double  CBs::computeThValue()
{
    return SM.getCBs();
}

PhiBs::PhiBs(const StandardModel& SM_i) : ThObservable(SM_i){
        SM.getFlavour().getDB2(1);
    };

double  PhiBs::computeThValue()
{
    return SM.getPhiBs();
}