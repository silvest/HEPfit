/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "xi.h"
#include "StandardModel.h"

xi::xi(const StandardModel& SM_i) : ThObservable(SM_i) 
{
    SM.initializeBParameter("BBd");
    SM.initializeBParameter("BBs");
};

double xi::computeThValue() 
{
    return(SM.getMesons(QCD::B_S).getDecayconst()*sqrt(SM.getBBs().getBpars()(0))/(SM.getMesons(QCD::B_D).getDecayconst()*sqrt(SM.getBBd().getBpars()(0))));
}
