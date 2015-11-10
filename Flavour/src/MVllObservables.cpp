/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVllObservables.h"
#include "Flavour.h"
#include <gslpp.h>

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



P_1::P_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_1::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(4,q_min,q_max)/(2.* SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max));
}


P_2::P_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_2::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(7,q_min,q_max)/(8.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max));
}


P_3::P_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_3::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(11,q_min,q_max)/(4.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max));
}


P_4Prime::P_4Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_4Prime::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(5,q_min,q_max)/sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(2,q_min,q_max)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max));
}


P_5Prime::P_5Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_5Prime::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(6,q_min,q_max)/(2.*sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(2,q_min,q_max)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max)));
}


P_6Prime::P_6Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_6Prime::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(9,q_min,q_max)/(2.*sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(2,q_min,q_max)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max)));
}


P_8Prime::P_8Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_8Prime::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(10,q_min,q_max)/(sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(2,q_min,q_max)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max)));
}


GammaPrime::GammaPrime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double GammaPrime::computeGammaPrime(double qmin, double qmax, StandardModel::lepton lep)
{
    double q_min = qmin;
    double q_max = qmax;
    StandardModel::lepton lep_i = lep;
    
    return ((3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(0,q_min,q_max) - SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(2,q_min,q_max)) + 2.*(3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(1,q_min,q_max) - SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(3,q_min,q_max)))/4.;
}

double GammaPrime::computeThValue()
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return computeGammaPrime(q_min, q_max, lep);
}


A_FB::A_FB(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double A_FB::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -3. * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(7,q_min,q_max) / 4. / computeGammaPrime(q_min, q_max, lep);
}


BR_MVll::BR_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double BR_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeGammaPrime(q_min, q_max, lep)/SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getwidth() / ( q_max - q_min );
}


R_MVll::R_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2) 
: GammaPrime(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
}

double R_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeGammaPrime(q_min, q_max, lep1) / computeGammaPrime(q_min, q_max, lep2);
}


RL_MVll::RL_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2) 
: F_L(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
}

double RL_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (computeGammaPrime(q_min, q_max, lep1) * computeFL(q_min, q_max, lep1)) / (computeGammaPrime(q_min, q_max, lep2) * computeFL(q_min, q_max, lep2));
}


RT_MVll::RT_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2) 
: F_L(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
}

double RT_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (computeGammaPrime(q_min, q_max, lep1) * (1 - computeFL(q_min, q_max, lep1))) / (computeGammaPrime(q_min, q_max, lep2) * ( 1 - computeFL(q_min, q_max, lep2)));
}


R_6::R_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2) 
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
}

double R_6::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep1)->integrateSigma(7,q_min,q_max) / SM.getMyFlavour()->getMVll(meson, vectorM, lep2)->integrateSigma(7,q_min,q_max);
}


ACP_MVll::ACP_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double ACP_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
          
    return (3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(0,q_min,q_max) - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(2,q_min,q_max) + 2. * ( 3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(1,q_min,q_max) - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(3,q_min,q_max) ) )/(4.*computeGammaPrime(q_min, q_max, lep));
}


P3CP::P3CP(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P3CP::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(11,q_min,q_max)/(4.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(3,q_min,q_max));
}


F_L::F_L(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double F_L::computeFL(double qmin, double qmax, StandardModel::lepton lep) 
{
    double q_min = qmin;
    double q_max = qmax;
    StandardModel::lepton lep_i = lep;
    
    double sigma0 = SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(0,q_min,q_max);
    double sigma2 = SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->integrateSigma(2,q_min,q_max);
    
    return (3.*sigma0 - sigma2) / (4. * computeGammaPrime(q_min, q_max, lep_i)) ;
}

double F_L::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeFL(q_min, q_max, lep);
}


M_1Prime::M_1Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double M_1Prime::computeThValue() 
{
    double q_min = getBinMin();
    
    return ( SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_p(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_m(q_min).abs2() - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_p(q_min).abs2() - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_m(q_min).abs2() )/
            ( 2.*( SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_p(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_m(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_p(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_m(q_min).abs2() ) );
}


M_2Prime::M_2Prime(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double M_2Prime::computeThValue() 
{
    double q_min = getBinMin();
    
    return ( q_min/(2.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getMlep()*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getMlep())*( SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_P(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->beta(q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->beta(q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_S(q_min).abs2() ) + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_0(q_min).abs2() - SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_0(q_min).abs2() )/
            ( SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_V_0(q_min).abs2() + SM.getMyFlavour()->getMVll(meson, vectorM, lep)->H_A_0(q_min).abs2() );  
}


S_3::S_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_3::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(4,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_4::S_4(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_4::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(5,q_min,q_max) / computeGammaPrime(q_min, q_max, lep) /2.;
}


S_5::S_5(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_5::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(6,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_7::S_7(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_7::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(9,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_8::S_8(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_8::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(10,q_min,q_max) / computeGammaPrime(q_min, q_max, lep) / 2.;
}


S_9::S_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_9::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateSigma(11,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_6::A_6(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double A_6::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(7,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_9::A_9(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double A_9::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->integrateDelta(11,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}

V0::V0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double V0::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getV0(getBinMin());
}

Vp::Vp(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double Vp::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getVp(getBinMin());
}

Vm::Vm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double Vm::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getVm(getBinMin());
}

T0::T0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double T0::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getT0(getBinMin());
}

Tp::Tp(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double Tp::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getTp(getBinMin());
}

Tm::Tm(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double Tm::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getTm(getBinMin());
}

S::S(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double S::computeThValue() 
{
   return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getS(getBinMin());
}


gtilde_1::gtilde_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double gtilde_1::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_1_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_1: incorrect type");
}


gtilde_2::gtilde_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double gtilde_2::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_2_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_2: incorrect type");
}


gtilde_3::gtilde_3(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double gtilde_3::computeThValue() 
{
    double q_min = getBinMin();
    
    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getgtilde_3_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_3: incorrect type");
}

h_0::h_0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double h_0::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_0_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_0: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}


h_p::h_p(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double h_p::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_p_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_p: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}


h_m::h_m(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
}

double h_m::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_re(q_min));
    else if (typ == 2) return (SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_im(q_min));
    else if (typ == 3) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_im(q_min)).abs());
    else if (typ == 4) return ((SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_re(q_min) + gslpp::complex::i() * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->geth_m_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_m: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}

hp0_hm0::hp0_hm0(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i)
: ThObservable(SM_i)
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double hp0_hm0::computeThValue()
{
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->gethp0_hm0_abs();
}

hm0_h00::hm0_h00(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i)
: ThObservable(SM_i)
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
}

double hm0_h00::computeThValue()
{
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->gethm0_h00_abs();
}

/***********************************************************************************************************************************
FUNCTIONAL
***********************************************************************************************************************************/

P_1f::P_1f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_1f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(4,q_min)/(2.* SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min));
}


P_2f::P_2f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_2f::computeThValue() 
{   
    double q_min = getBinMin();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(7,q_min)/(8.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min));
}


P_3f::P_3f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_3f::computeThValue() 
{   
    double q_min = getBinMin();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(11,q_min)/(4.*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min));
}


P_4Primef::P_4Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_4Primef::computeThValue() 
{   
    double q_min = getBinMin();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(5,q_min)/sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(2,q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min));
}


P_5Primef::P_5Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_5Primef::computeThValue() 
{   
    double q_min = getBinMin();
    
    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(6,q_min)/(2.*sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(2,q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min)));
}


P_6Primef::P_6Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_6Primef::computeThValue() 
{
    double q_min = getBinMin();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(9,q_min)/(2.*sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(2,q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min)));
}


P_8Primef::P_8Primef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double P_8Primef::computeThValue() 
{
    double q_min = getBinMin();
    
    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(10,q_min)/(sqrt(-SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(2,q_min)*SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(3,q_min)));
}


GammaPrimef::GammaPrimef(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double GammaPrimef::computeGammaPrimef(double qmin, StandardModel::lepton lep)
{
    double q_min = qmin;
    StandardModel::lepton lep_i = lep;
    
    return ((3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(0,q_min) - SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(2,q_min)) + 2.*(3.*SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(1,q_min) - SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(3,q_min)))/4.;
}

double GammaPrimef::computeThValue()
{
    double q_min = getBinMin();

    return computeGammaPrimef(q_min, lep);
}


A_FBf::A_FBf(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double A_FBf::computeThValue() 
{
    double q_min = getBinMin();

    return -3. * SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(7,q_min) / 4. / computeGammaPrimef(q_min, lep);
}

F_Lf::F_Lf(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double F_Lf::computeFLf(double qmin, StandardModel::lepton lep) 
{
    double q_min = qmin;
    StandardModel::lepton lep_i = lep;
    
    double sigma0 = SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(0,q_min);
    double sigma2 = SM.getMyFlavour()->getMVll(meson, vectorM, lep_i)->getSigma(2,q_min);
    
    return (3.*sigma0 - sigma2) / (4. * computeGammaPrimef(q_min, lep_i)) ;
}


double F_Lf::computeThValue() 
{
    double q_min = getBinMin();
    
    return computeFLf(q_min, lep);
}

S_3f::S_3f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_3f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(4,q_min) / computeGammaPrimef(q_min, lep);
}


S_4f::S_4f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_4f::computeThValue() 
{
    double q_min = getBinMin();

    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(5,q_min) / computeGammaPrimef(q_min, lep) /2.;
}


S_5f::S_5f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_5f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(6,q_min) / computeGammaPrimef(q_min, lep);
}


S_7f::S_7f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_7f::computeThValue() 
{
    double q_min = getBinMin();

    return -SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(9,q_min) / computeGammaPrimef(q_min, lep);
}


S_8f::S_8f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

/*Returns experimental value, defined according to 1308.1707*/
double S_8f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(10,q_min) / computeGammaPrimef(q_min, lep) / 2.;
}


S_9f::S_9f(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}


double S_9f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getSigma(11,q_min) / computeGammaPrimef(q_min, lep);
}

BRf_MVll::BRf_MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
}

double BRf_MVll::computeThValue() 
{
    double q_min = getBinMin();
    
    return computeGammaPrimef(q_min, lep)/SM.getMyFlavour()->getMVll(meson, vectorM, lep)->getwidth();
}