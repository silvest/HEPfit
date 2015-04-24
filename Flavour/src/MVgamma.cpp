/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVgamma.h"


MVgamma::MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
}


MVgamma::~MVgamma() 
{}

void MVgamma::updateParameters()
{
    GF = SM.getGF();
    ale = SM.getAle();
    MM = SM.getMesons(meson).getMass();
    MM2 = MM * MM;
    MV = SM.getMesons(vectorM).getMass();
    Mb = SM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Ms = SM.getQuarks(QCD::STRANGE).getMass();
    MW = SM.Mw();
    lambda_t = SM.computelamt_s();
    mu_b = SM.getMub();
    width = SM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV,2.);
    
    switch(vectorM){
        case StandardModel::K_star :
            r_1T1=SM.getr_1T1();
            r_2T1=SM.getr_1T3t() + SM.getr_2T3t() - r_1T1;//SM.getr_2T1();
            m_RT1=SM.getm_RT1();
            m_fit2T1=SM.getm_fit2T1();
            
            break;
        case StandardModel::PHI :
            r_1T1=SM.getr_1T1phi();
            r_2T1=SM.getr_1T3tphi() + SM.getr_2T3tphi() - r_1T1;//SM.getr_2T1();
            m_RT1=SM.getm_RT1phi();
            m_fit2T1=SM.getm_fit2T1phi();
            
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
    
    
    h[0]=SM.geth_p();    //h_plus
    h[1]=SM.geth_m();    //h_minus
    
    allcoeff = SM.getMyFlavour()->ComputeCoeffBMll(mu_b);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = SM.getMyFlavour()->ComputeCoeffprimeBMll(mu_b);   //check the mass scale, scheme fixed to NDR
    
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    
}


/*******************************************************************************
 * Form Factor                                                     *
 * ****************************************************************************/
double MVgamma::T_1()
{
    return SM.getMyFlavour()->getMVll(meson, vectorM, StandardModel::MU)->LCSR_fit1(0., r_1T1, r_2T1, pow(m_RT1, 2.), m_fit2T1);
}



/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
complex MVgamma::H_V_m() 
{
    return lambda_t * (C_7 *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*h[1]);
}

complex MVgamma::H_V_p() {
    return lambda_t * (- C_7p *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*h[0]);
}

complex MVgamma::H_V_m_bar() {
    return lambda_t.conjugate() * (C_7 *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*h[1]);
}

complex MVgamma::H_V_p_bar() {
    return lambda_t.conjugate() * (- C_7p *T_1() * lambda / MM2 - MM/(2*Mb)*16*M_PI*M_PI*h[0]);
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MVgamma::BR_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) 
: MVgamma(SM_i, meson_i, vector_i) 
{
    meson = meson_i;
    vectorM = vector_i;
}

double BR_MVgamma::computeThValue()
{
    updateParameters();
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda /(4. * width) * (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

ACP_MVgamma::ACP_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : MVgamma(SM_i, meson_i, vector_i) 
{
    meson = meson_i;
    vectorM = vector_i;
}

double ACP_MVgamma::computeThValue()
{
    updateParameters();
    
    return ((H_V_p().abs2() + H_V_m().abs2() - H_V_p_bar().abs2() - H_V_m_bar().abs2())) / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}