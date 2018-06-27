/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPlnuObservables.h"
#include "MPlnu.h"
#include "StandardModel.h"

RD_MPlnu::RD_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).initializeMPlnuParameters());
}

double RD_MPlnu::computeThValue() 
{
    double MM = SM.getMesons(meson).getMass();
    double MP = SM.getMesons(pseudoscalarM).getMass();
    double q2_min_lep1 = SM.getLeptons(lep1).getMass()*SM.getLeptons(lep1).getMass();
    double q2_min_lep2 = SM.getLeptons(lep2).getMass()*SM.getLeptons(lep2).getMass();
    double q2_min_lep3 = SM.getLeptons(lep3).getMass()*SM.getLeptons(lep3).getMass();
    double w_max_lep1 = (MM*MM+MP*MP)/(2.*MM*MP)-q2_min_lep1/(2.*MM*MP);
    double w_max_lep2 = (MM*MM+MP*MP)/(2.*MM*MP)-q2_min_lep2/(2.*MM*MP);
    double w_max_lep3 = (MM*MM+MP*MP)/(2.*MM*MP)-q2_min_lep3/(2.*MM*MP);
    double q2_max = (MM-MP)*(MM-MP);
    double w_min = (MM*MM+MP*MP)/(2.*MM*MP)-q2_max/(2.*MM*MP);
    double deltaGammadeltaw_lep1 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).getDeltaGammaDeltaw(w_min,w_max_lep1);
    double deltaGammadeltaw_lep2 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep2).getDeltaGammaDeltaw(w_min,w_max_lep2);
    double deltaGammadeltaw_lep3 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep3).getDeltaGammaDeltaw(w_min,w_max_lep3);
    
    return deltaGammadeltaw_lep1/((deltaGammadeltaw_lep2+deltaGammadeltaw_lep3)/2.);
}