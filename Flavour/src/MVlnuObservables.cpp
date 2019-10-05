/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVlnuObservables.h"
#include "MVlnu.h"
#include "MPlnu.h"
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
    double MM = SM.getMesons(meson).getMass();
    double MV = SM.getMesons(vectorM).getMass();
    double q2_min_lep1 = SM.getLeptons(lep1).getMass()*SM.getLeptons(lep1).getMass();
    double w_max_lep1 = (MM*MM+MV*MV)/(2.*MM*MV)-q2_min_lep1/(2.*MM*MV);
    double deltaGammadeltaw_lep2_extra = 0.;
    
    w_max_lep1 = std::min(w_max_lep1, w_max);
    
    if (w_max >  w_max_lep1) deltaGammadeltaw_lep2_extra = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltaw(w_max_lep1, w_max);
    
    double deltaGammadeltaw_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltaw(w_min,w_max_lep1);
    double deltaGammadeltaw_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltaw(w_min,w_max_lep1);
    
    return (deltaGammadeltaw_lep2_extra+0.5*(deltaGammadeltaw_lep1+deltaGammadeltaw_lep2)); // NOTE: DO NOT divide by the bin size because experimentalists do not.
}

RDstar_MVlnu::RDstar_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2, QCD::lepton lep_3)
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    lep3 = lep_3;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep1).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep2).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep3).initializeMVlnuParameters());
}

double RDstar_MVlnu::computeThValue() 
{
    double MM = SM.getMesons(meson).getMass();
    double MV = SM.getMesons(vectorM).getMass();
    double q2_min_lep1 = SM.getLeptons(lep1).getMass()*SM.getLeptons(lep1).getMass();
    double q2_min_lep2 = SM.getLeptons(lep2).getMass()*SM.getLeptons(lep2).getMass();
    double q2_min_lep3 = SM.getLeptons(lep3).getMass()*SM.getLeptons(lep3).getMass();
    double w_max_lep1 = (MM*MM+MV*MV)/(2.*MM*MV)-q2_min_lep1/(2.*MM*MV);
    double w_max_lep2 = (MM*MM+MV*MV)/(2.*MM*MV)-q2_min_lep2/(2.*MM*MV);
    double w_max_lep3 = (MM*MM+MV*MV)/(2.*MM*MV)-q2_min_lep3/(2.*MM*MV);
    double q2_max = (MM-MV)*(MM-MV);
    double w_min = (MM*MM+MV*MV)/(2.*MM*MV)-q2_max/(2.*MM*MV);
    double deltaGammadeltaw_lep1 = SM.getFlavour().getMVlnu(meson, vectorM, lep1).getDeltaGammaDeltaw(w_min,w_max_lep1);
    double deltaGammadeltaw_lep2 = SM.getFlavour().getMVlnu(meson, vectorM, lep2).getDeltaGammaDeltaw(w_min,w_max_lep2);
    double deltaGammadeltaw_lep3 = SM.getFlavour().getMVlnu(meson, vectorM, lep3).getDeltaGammaDeltaw(w_min,w_max_lep3);
    
    return deltaGammadeltaw_lep1/((deltaGammadeltaw_lep2+deltaGammadeltaw_lep3)/2.);
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

FL_MVlnu::FL_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i): ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}

double FL_MVlnu::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).getFL();
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

UnitarityP_MVlnu::UnitarityP_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}

double UnitarityP_MVlnu::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_unitarity_P_BGL();
}

FF_hA1atw1::FF_hA1atw1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
FF_hV::FF_hV(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_hV::computeThValue()  
{
    double w = getBinMin();
    
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1(w) * SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R1(w);
}

double FF_hA1atw1::computeThValue()  
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1w1();
}

FF_hA1::FF_hA1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_hA1::computeThValue()  
{
    double w = getBinMin();
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1(w);
}

FF_hA2::FF_hA2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_hA2::computeThValue()  
{
    double w = getBinMin();
//    double rD = SM.getMesons(vectorM).getMass()/SM.getMesons(meson).getMass();
//    double hA1 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1(w);
//    double R2 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R2(w);
//    double R0 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R0(w);
//    double R3 = (R2 * (1. - rD) + rD * (R0 * (1. + rD) - 2.))/(1. - rD)/(1. - rD);
    
//    return (R2 - R3) * hA1/2./rD;
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA2(w);
}

FF_hA3::FF_hA3(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_hA3::computeThValue()  
{
    double w = getBinMin();
//    double rD = SM.getMesons(vectorM).getMass()/SM.getMesons(meson).getMass();
//    double hA1 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA1(w);
//    double R2 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R2(w);
//    double R0 = SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R0(w);
//    double R3 = (R2 * (1. - rD) + rD * (R0 * (1. + rD) - 2.))/(1. - rD)/(1. - rD);
//    
//    return (R2 + R3) * hA1/2.;
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_hA3(w);
}

FF_R1::FF_R1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_R1::computeThValue()  
{
    double w = getBinMin();
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R1(w);
}

FF_R2::FF_R2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_R2::computeThValue()  
{
    double w = getBinMin();
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R2(w);
}

FF_R0::FF_R0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}
 
double FF_R0::computeThValue()  
{
    double w = getBinMin();
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_R0(w);
}

UnitarityV_D_Dst::UnitarityV_D_Dst(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::meson pseudoscalar_i,  QCD::lepton lep_i)
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    pseudoscalarM = pseudoscalar_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
    setParametersForObservable(SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).initializeMPlnuParameters());
}

double UnitarityV_D_Dst::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).get_unitarity_V_BGL() + SM.getFlavour().getMPlnu(meson, pseudoscalarM, lep).get_unitarity_1min_BGL();
}

Plep_MVlnu::Plep_MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i): ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVlnu(meson, vectorM, lep).initializeMVlnuParameters());
}

double Plep_MVlnu::computeThValue() 
{
    return SM.getFlavour().getMVlnu(meson, vectorM, lep).getPlep();
}