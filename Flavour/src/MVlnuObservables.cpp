/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVlnuObservables.h"
#include "MVlnu.h"
#include "StandardModel.h"

Gammaw_MVlnu::Gammaw_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
}

double Gammaw_MVlnu::computeThValue() 
{
    double w_min = getBinMin();
    double w_max = getBinMax();
    double deltaGammadeltaw_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltaw(w_min,w_max);
    double deltaGammadeltaw_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltaw(w_min,w_max);
    
    return (deltaGammadeltaw_lep1+deltaGammadeltaw_lep2)/2.;
}

Gammacl_MVlnu::Gammacl_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
}

double Gammacl_MVlnu::computeThValue() 
{
    /* sign change for theta_l to go from exp to theory conventions */
    double cl_min = -getBinMax();
    double cl_max = -getBinMin();
    double deltaGammadeltacl_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltacl(cl_min,cl_max);
    double deltaGammadeltacl_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltacl(cl_min,cl_max);
    
    return (deltaGammadeltacl_lep1+deltaGammadeltacl_lep2)/2.;
}

GammacV_MVlnu::GammacV_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
}

double GammacV_MVlnu::computeThValue() 
{
    double cV_min = getBinMin();
    double cV_max = getBinMax();
    double deltaGammadeltacV_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltacV(cV_min,cV_max);
    double deltaGammadeltacV_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltacV(cV_min,cV_max);
    
    return (deltaGammadeltacV_lep1+deltaGammadeltacV_lep2)/2.;
}

Gammachi_MVlnu::Gammachi_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
}

double Gammachi_MVlnu::computeThValue() 
{
    double chi_min = getBinMin();
    double chi_max = getBinMax();
    double deltaGammadeltachi_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltachi(chi_min,chi_max);
    double deltaGammadeltachi_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltachi(chi_min,chi_max);
    
    return (deltaGammadeltachi_lep1+deltaGammadeltachi_lep2)/2.;
}

UnitarityV_MVlnu::UnitarityV_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}

double UnitarityV_MVlnu::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_unitarity_V_BGL();
}

UnitarityA_MVlnu::UnitarityA_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}

double UnitarityA_MVlnu::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_unitarity_A_BGL();
}

FF_hA1atw1::FF_hA1atw1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_hA1atw1::computeThValue()  
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1w1();
}
