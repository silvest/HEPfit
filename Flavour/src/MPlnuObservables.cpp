/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MPlnuObservables.h"
#include "MPlnu.h"
#include "StandardModel.h"

Gammaw_MPlnu::Gammaw_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).initializeMPlnuParameters());
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep2).initializeMPlnuParameters());
}

double Gammaw_MPlnu::computeThValue() 
{
    double w_min = getBinMin();
    double w_max = getBinMax();
    double deltaGammadeltaw_lep1 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).getDeltaGammaDeltaw(w_min,w_max);
    double deltaGammadeltaw_lep2 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep2).getDeltaGammaDeltaw(w_min,w_max);
    
    return (deltaGammadeltaw_lep1+deltaGammadeltaw_lep2)/2./(w_max-w_min);
}

RD_MPlnu::RD_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).initializeMPlnuParameters());
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep2).initializeMPlnuParameters());
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep3).initializeMPlnuParameters());
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

UnitarityV_MPlnu::UnitarityV_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).initializeMPlnuParameters());
}

double UnitarityV_MPlnu::computeThValue() 
{
    return SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).get_unitarity_0plus_BGL();
}

UnitarityA_MPlnu::UnitarityA_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).initializeMPlnuParameters());
}

double UnitarityA_MPlnu::computeThValue() 
{
    return SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).get_unitarity_1min_BGL();
}

FFplus_MPlnu::FFplus_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep1 = lep_i;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).initializeMPlnuParameters());
}

double FFplus_MPlnu::computeThValue() 
{
    double MM = SM.getMesons(meson).getMass();
    double MP = SM.getMesons(pseudoscalarM).getMass();
    double w0 = (MM*MM+MP*MP)/2./MM/MP;
    double w_val = getBinMax();
    double q2_val = 2.*MM*MP*(w0-w_val);
    double res_fplus = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).get_fplus(q2_val);
    
    return res_fplus;
}

FF0_MPlnu::FF0_MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep1 = lep_i;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).initializeMPlnuParameters());
}

double FF0_MPlnu::computeThValue() 
{
    double MM = SM.getMesons(meson).getMass();
    double MP = SM.getMesons(pseudoscalarM).getMass();
    double w0 = (MM*MM+MP*MP)/2./MM/MP; 
    double w_val = getBinMax();
    double q2_val = 2.*MM*MP*(w0-w_val);
    double res_f0 = SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep1).get_f0(q2_val);
    
    return res_f0;
}
