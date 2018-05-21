/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVlnuObservables.h"
#include "MVlnu.h"
#include "StandardModel.h"

BR_MVlnu::BR_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
}

double BR_MVlnu::computeThValue() 
{
    double w_min = getBinMin();
    double w_max = getBinMax();
    double deltaBRdeltaw_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaBRDeltaw(w_min,w_max);
    double deltaBRdeltaw_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaBRDeltaw(w_min,w_max);
    double tauBinv = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getMwidth();
    
    return (deltaBRdeltaw_lep1+deltaBRdeltaw_lep2)/2./tauBinv;
}
