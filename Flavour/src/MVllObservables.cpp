/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "MVllObservables.h"
#include "MVll.h"
#include "DGamma.h"
#include "StandardModel.h"
//#include "gslpp.h"

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/



P_1::P_1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_1::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(4,q_min,q_max)/(2.* SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
}

/*Returns experimental value, defined according to 1510.04239.*/
P_2::P_2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_2::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double P_2_theory = SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(7,q_min,q_max)/(8.*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
    
    return -P_2_theory;
}


P_3::P_3(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_3::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double P_3_theory = -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(11,q_min,q_max)/(4.*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
    
    return -P_3_theory;
}


P_4Prime::P_4Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_4Prime::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double P_4p_theory = SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(5,q_min,q_max)/sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
    
    return -1./2.*P_4p_theory;
}


P_5Prime::P_5Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_5Prime::computeThValue() 
{   
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(6,q_min,q_max)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)));
}


P_6Prime::P_6Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_6Prime::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(9,q_min,q_max)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)));
}


P_8Prime::P_8Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_8Prime::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double P_8p_theory = -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(10,q_min,q_max)/(sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(2,q_min,q_max)*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max)));
    
    return -1./2.*P_8p_theory;
}


GammaPrime::GammaPrime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double GammaPrime::computeGammaPrime(double qmin, double qmax, QCD::lepton lep)
{
    double q_min = qmin;
    double q_max = qmax;
    QCD::lepton lep_i = lep;
    
    return ((3.*SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(0,q_min,q_max) - SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(2,q_min,q_max)) + 2.*(3.*SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(1,q_min,q_max) - SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(3,q_min,q_max)))/4.;
}

double GammaPrime::computeThValue()
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return computeGammaPrime(q_min, q_max, lep);
}


A_FB::A_FB(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double A_FB::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -3. * SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(7,q_min,q_max) / 4. / computeGammaPrime(q_min, q_max, lep);
}


BR_MVll::BR_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double BR_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    double ys = DGamma_s_MSbar(SM).computeThValue()*SM.getMesons(QCD::B_S).getLifetime()/2.;
    
    switch(vectorM){
            case StandardModel::K_star:
            case StandardModel::K_star_P:
                return computeGammaPrime(q_min, q_max, lep)/SM.getFlavour().getMVll(meson, vectorM, lep).getwidth() / ( q_max - q_min );
                break;
            case StandardModel::PHI:
                return computeGammaPrime(q_min, q_max, lep)/SM.getFlavour().getMVll(meson, vectorM, lep).getwidth()/(1. - ys*ys) / ( q_max - q_min );
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("BR_MVll: vector " + out.str() + " not implemented");
        }
}


R_MVll::R_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: GammaPrime(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep1).initializeMVllParameters());
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep2).initializeMVllParameters());
}

double R_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeGammaPrime(q_min, q_max, lep1) / computeGammaPrime(q_min, q_max, lep2);
}


RL_MVll::RL_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: F_L(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep1).initializeMVllParameters());
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep2).initializeMVllParameters());
}

double RL_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (computeGammaPrime(q_min, q_max, lep1) * computeFL(q_min, q_max, lep1)) / (computeGammaPrime(q_min, q_max, lep2) * computeFL(q_min, q_max, lep2));
}


RT_MVll::RT_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: F_L(SM_i, meson_i, vector_i, lep_1) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep1).initializeMVllParameters());
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep2).initializeMVllParameters());
}

double RT_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return (computeGammaPrime(q_min, q_max, lep1) * (1 - computeFL(q_min, q_max, lep1))) / (computeGammaPrime(q_min, q_max, lep2) * ( 1 - computeFL(q_min, q_max, lep2)));
}


R_6::R_6(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_1, QCD::lepton lep_2) 
: ThObservable(SM_i) 
{  
    lep1 = lep_1;
    lep2 = lep_2;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep1).initializeMVllParameters());
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep2).initializeMVllParameters());
}

double R_6::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return SM.getFlavour().getMVll(meson, vectorM, lep1).integrateSigma(7,q_min,q_max) / SM.getFlavour().getMVll(meson, vectorM, lep2).integrateSigma(7,q_min,q_max);
}


ACP_MVll::ACP_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double ACP_MVll::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
          
    return (3.*SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(0,q_min,q_max) - SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(2,q_min,q_max) + 2. * ( 3.*SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(1,q_min,q_max) - SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(3,q_min,q_max) ) )/(4.*computeGammaPrime(q_min, q_max, lep));
}


P3CP::P3CP(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P3CP::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return - SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(11,q_min,q_max)/(4.*SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(3,q_min,q_max));
}


F_L::F_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double F_L::computeFL(double qmin, double qmax, QCD::lepton lep) 
{
    double q_min = qmin;
    double q_max = qmax;
    QCD::lepton lep_i = lep;
    
    double sigma0 = SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(0,q_min,q_max);
    double sigma2 = SM.getFlavour().getMVll(meson, vectorM, lep_i).integrateSigma(2,q_min,q_max);
    
    return (3.*sigma0 - sigma2) / (4. * computeGammaPrime(q_min, q_max, lep_i)) ;
}

double F_L::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    
    return computeFL(q_min, q_max, lep);
}


M_1Prime::M_1Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double M_1Prime::computeThValue() 
{
    double q_min = getBinMin();
    
    return ( SM.getFlavour().getMVll(meson, vectorM, lep).H_V_p(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).H_V_m(q_min, 0).abs2() - SM.getFlavour().getMVll(meson, vectorM, lep).H_A_p(q_min, 0).abs2() - SM.getFlavour().getMVll(meson, vectorM, lep).H_A_m(q_min, 0).abs2() )/
            ( 2.*( SM.getFlavour().getMVll(meson, vectorM, lep).H_V_p(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).H_V_m(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).H_A_p(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).H_A_m(q_min, 0).abs2() ) );
}


M_2Prime::M_2Prime(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double M_2Prime::computeThValue() 
{
    double q_min = getBinMin();
    
    return ( q_min/(2.*SM.getFlavour().getMVll(meson, vectorM, lep).getMlep()*SM.getFlavour().getMVll(meson, vectorM, lep).getMlep())*( SM.getFlavour().getMVll(meson, vectorM, lep).H_P(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).beta(q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).beta(q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).H_S(q_min, 0).abs2() ) + SM.getFlavour().getMVll(meson, vectorM, lep).H_V_0(q_min, 0).abs2() - SM.getFlavour().getMVll(meson, vectorM, lep).H_A_0(q_min, 0).abs2() )/
            ( SM.getFlavour().getMVll(meson, vectorM, lep).H_V_0(q_min, 0).abs2() + SM.getFlavour().getMVll(meson, vectorM, lep).H_A_0(q_min, 0).abs2() );  
}


S_3::S_3(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double S_3::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(4,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_4::S_4(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_4::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(5,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_5::S_5(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double S_5::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(6,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_7::S_7(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_7::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(9,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_8::S_8(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double S_8::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(10,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


S_9::S_9(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_9::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).integrateSigma(11,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_5::A_5(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double A_5::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(6,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_6::A_6(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double A_6::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(7,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_6c::A_6c(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double A_6c::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(8,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_8::A_8(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double A_8::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(10,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}


A_9::A_9(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrime(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double A_9::computeThValue() 
{
    double q_min = getBinMin();
    double q_max = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).integrateDelta(11,q_min,q_max) / computeGammaPrime(q_min, q_max, lep);
}

V0::V0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double V0::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getV0(getBinMin());
}

Vp::Vp(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double Vp::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getVp(getBinMin());
}

Vm::Vm(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double Vm::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getVm(getBinMin());
}

T0::T0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double T0::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getT0(getBinMin());
}

Tp::Tp(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double Tp::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getTp(getBinMin());
}

Tm::Tm(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double Tm::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getTm(getBinMin());
}

S::S(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
   lep = lep_i;
   meson = meson_i;
   vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double S::computeThValue() 
{
   return SM.getFlavour().getMVll(meson, vectorM, lep).getS(getBinMin());
}


gtilde_1::gtilde_1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double gtilde_1::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_1_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_1: incorrect type");
}


gtilde_2::gtilde_2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double gtilde_2::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_2_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_2: incorrect type");
}


gtilde_3::gtilde_3(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double gtilde_3::computeThValue() 
{
    double q_min = getBinMin();
    
    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).getgtilde_3_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::gtilde_3: incorrect type");
}

h_0::h_0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double h_0::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_0_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_0: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}


h_p::h_p(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double h_p::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_p_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_p: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}


h_m::h_m(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i, unsigned int typ_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    typ = typ_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double h_m::computeThValue() 
{
    double q_min = getBinMin();

    if (typ == 1) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_re(q_min));
    else if (typ == 2) return (SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_im(q_min));
    else if (typ == 3) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_im(q_min)).abs());
    else if (typ == 4) return ((SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_re(q_min) + gslpp::complex::i() * SM.getFlavour().getMVll(meson, vectorM, lep).geth_m_im(q_min)).arg());
    else throw std::runtime_error("MVllObservables::h_m: type can only be 1:real, 2:imaginary, 3:absolute and 4:argument");
}

/***********************************************************************************************************************************
FUNCTIONAL
***********************************************************************************************************************************/

P_1f::P_1f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_1f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(4,q_min)/(2.* SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
}


P_2f::P_2f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_2f::computeThValue() 
{   
    double q_min = getBinMin();
    
    double P_2f_theory = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(7,q_min)/(8.*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    
    return -P_2f_theory;
}


P_3f::P_3f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_3f::computeThValue() 
{   
    double q_min = getBinMin();
    
    double P_3f_theory = -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(11,q_min)/(4.*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    
    return -P_3f_theory;
}


P_4Primef::P_4Primef(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, defined according to 1510.04239.*/
double P_4Primef::computeThValue() 
{   
    double q_min = getBinMin();
    
    double P_4pf_theory = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(5,q_min)/sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    
    return -1./2.*P_4pf_theory;
}


P_5Primef::P_5Primef(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_5Primef::computeThValue() 
{   
    double q_min = getBinMin();
    
    return SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(6,q_min)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
}


P_6Primef::P_6Primef(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_6Primef::computeThValue() 
{
    double q_min = getBinMin();
    
    return -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(9,q_min)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
}

/*Returns experimental value, defined according to 1510.04239.*/
P_8Primef::P_8Primef(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_8Primef::computeThValue() 
{
    double q_min = getBinMin();
    
    double P_8pf_theory = -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(10,q_min)/(sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
    
    return -1./2.*P_8pf_theory;
}


GammaPrimef::GammaPrimef(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double GammaPrimef::computeGammaPrimef(double qmin, QCD::lepton lep)
{
    double q_min = qmin;
    QCD::lepton lep_i = lep;
    
    return ((3.*SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(0,q_min) - SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(2,q_min)) + 2.*(3.*SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(1,q_min) - SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(3,q_min)))/4.;
}

double GammaPrimef::computeThValue()
{
    double q_min = getBinMin();

    return computeGammaPrimef(q_min, lep);
}


A_FBf::A_FBf(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double A_FBf::computeThValue() 
{
    double q_min = getBinMin();

    return -3. * SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(7,q_min) / 4. / computeGammaPrimef(q_min, lep);
}

F_Lf::F_Lf(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double F_Lf::computeFLf(double qmin, QCD::lepton lep) 
{
    double q_min = qmin;
    QCD::lepton lep_i = lep;
    
    double sigma0 = SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(0,q_min);
    double sigma2 = SM.getFlavour().getMVll(meson, vectorM, lep_i).getSigma(2,q_min);
    
    return (3.*sigma0 - sigma2) / (4. * computeGammaPrimef(q_min, lep_i)) ;
}


double F_Lf::computeThValue() 
{
    double q_min = getBinMin();
    
    return computeFLf(q_min, lep);
}

S_3f::S_3f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double S_3f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(4,q_min) / computeGammaPrimef(q_min, lep);
}


S_4f::S_4f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_4f::computeThValue() 
{
    double q_min = getBinMin();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(5,q_min) / computeGammaPrimef(q_min, lep);
}


S_5f::S_5f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}


double S_5f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(6,q_min) / computeGammaPrimef(q_min, lep);
}


S_7f::S_7f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_7f::computeThValue() 
{
    double q_min = getBinMin();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(9,q_min) / computeGammaPrimef(q_min, lep);
}


S_8f::S_8f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double S_8f::computeThValue() 
{
    double q_min = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(10,q_min) / computeGammaPrimef(q_min, lep);
}


S_9f::S_9f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

/*Returns experimental value, due to difference in \theta_l and \phi definitions*/
double S_9f::computeThValue() 
{
    double q_min = getBinMin();

    return -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(11,q_min) / computeGammaPrimef(q_min, lep);
}

BRf_MVll::BRf_MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: GammaPrimef(SM_i, meson_i, vector_i, lep_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double BRf_MVll::computeThValue() 
{
    double q_min = getBinMin();
    
    return computeGammaPrimef(q_min, lep)/SM.getFlavour().getMVll(meson, vectorM, lep).getwidth();
}

P_relationf::P_relationf(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_relationf::computeThValue() 
{   
    double q_min = getBinMin();
    double P1 = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(4,q_min)/(2.* SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P2 = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(7,q_min)/(8.*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P4p = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(5,q_min)/sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P5p = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(6,q_min)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
    double beta = sqrt(1. - 4.*SM.getLeptons(lep).getMass()*SM.getLeptons(lep).getMass()/q_min/q_min);
    
    return 1./2.*(P4p*P5p + 1./beta *sqrt(std::abs((-1. + P1 + P4p*P4p)*(-1. - P1 + beta*beta*P5p*P5p)))) - P2;
}

P_relation_exactf::P_relation_exactf(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{  
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double P_relation_exactf::computeThValue() 
{   
    double q_min = getBinMin();
    double P1 = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(4,q_min)/(2.* SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P2 = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(7,q_min)/(8.*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P3 = -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(11,q_min)/(4.*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P4p = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(5,q_min)/sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min));
    double P5p = SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(6,q_min)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
    double P6p = -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(9,q_min)/(2.*sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
    double P8p = -SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(10,q_min)/(sqrt(-SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(2,q_min)*SM.getFlavour().getMVll(meson, vectorM, lep).getSigma(3,q_min)));
    
    double beta = sqrt(1. - 4.*SM.getLeptons(lep).getMass()*SM.getLeptons(lep).getMass()/q_min/q_min);
    
    /* k1 and k2 not fully implemented. arXiv:1402.6855*/
    double k1 = 1.;
    double k2 = 1.;
    
    double delta_1 = P6p*P8p;
    double delta_2 = -1. + k1*k1*k2*k2 + (1. - k1*k2)*(P4p*P4p + beta*beta*P5p*P5p) - 4.*k1*k1*P3*P3 + beta*beta*P6p*P8p*(2.*P4p*P5p + P6p*P8p) + k1*(beta*beta*P6p*(4.*P3*P5p - k2*P6p) - P8p*(4.*P3*P4p + k2*P8p));
    double delta_3 = (1. - k1)*P4p*P4p + beta*beta*((-1. + k1)*P5p*P5p - k1*P6p*P6p) + k1*P8p*P8p;
    double delta_4 = 1. - k1*k1;

    return 1./2./k1*((P4p*P5p + delta_1) + 1./beta *sqrt(std::abs((-1. + P1 + P4p*P4p)*(-1. - P1 + beta*beta*P5p*P5p) + delta_2 + delta_3*P1 + delta_4*P1*P1))) - P2;
}

QCDfC9_1f::QCDfC9_1f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9_1f::computeThValue() 
{
    double q2 = getBinMin();
    double cutoff = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9_1(q2, cutoff);
}

QCDfC9_2f::QCDfC9_2f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9_2f::computeThValue() 
{
    double q2 = getBinMin();
    double cutoff = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9_2(q2, cutoff);
}

QCDfC9_3f::QCDfC9_3f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9_3f::computeThValue() 
{
    double q2 = getBinMin();
    double cutoff = getBinMax();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9_3(q2, cutoff);
}

QCDfC9p_1f::QCDfC9p_1f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9p_1f::computeThValue() 
{
    double cutoff = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9p_1(cutoff);
}

QCDfC9p_2f::QCDfC9p_2f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9p_2f::computeThValue() 
{
    double cutoff = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9p_2(cutoff);
}

QCDfC9p_3f::QCDfC9p_3f(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) 
: ThObservable(SM_i) 
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVll(meson, vectorM, lep).initializeMVllParameters());
}

double QCDfC9p_3f::computeThValue() 
{
    double cutoff = getBinMin();

    return SM.getFlavour().getMVll(meson, vectorM, lep).getQCDfC9p_3(cutoff);
}