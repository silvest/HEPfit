/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "bsgammaTHDM.h"
#include "StandardModel.h"

bsgammaTHDM::bsgammaTHDM(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i))
{
    mycache = new THDMcache(SM_i);
}

bsgammaTHDM::~bsgammaTHDM()
{
    if (mycache != NULL) delete(mycache);
}

double bsgammaTHDM::computeThValue()
{
    std::string modelflag=myTHDM->getModelTypeflag();
    if (modelflag == "type2" || modelflag == "typeY" ) {
        double logtb=myTHDM->getlogtb();
        double logmHp=log10(myTHDM->getmHp());
        double B_BXsgamma_ctrl=mycache->ip_ex_bsgamma(logtb,logmHp);
        double B_BXsgamma_err=0.07;
        double bsgamma_theoryerror=myTHDM->getbsgamma_theoryerror();

        return B_BXsgamma_ctrl * (1.0 + B_BXsgamma_err * bsgamma_theoryerror);
    }
    else {
        throw std::runtime_error("flag_model in bsgammaTHDM::computeThValue() can only be \"type2\" or \"typeY\" at the moment");
        return 0.;
    }
}
