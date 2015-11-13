/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "bsgamma.h"
#include "StandardModel.h"

bsgammaTHDM::bsgammaTHDM(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i)),
        mySM (SM_i)
{
    mycache = new THDMcache();
}

bsgammaTHDM::~bsgammaTHDM()
{}

double bsgammaTHDM::computeThValue()
{
    std::string modelflag=myTHDM->getModelTypeflag();
    if (modelflag == "type2") {
        double logtb=myTHDM->getLogtb();
        double logmHp=log10(myTHDM->getMHp());
        double bsgtable=mycache->ex_bsgamma(logtb,logmHp);
        double B_BXsgamma_ctrl=bsgtable;
        double B_BXsgamma_err=0.07;
        double bsgamma_theoryerror=myTHDM->getbsgamma_theoryerror();

        return B_BXsgamma_ctrl * (1.0 + B_BXsgamma_err * bsgamma_theoryerror);
    }
    else {
        throw std::runtime_error("flag_model in bsgammaTHDM::computeThValue() can only be \"type2\" at the moment");
        return 0.;
    }
}