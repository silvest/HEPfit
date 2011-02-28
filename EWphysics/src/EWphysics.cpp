/* 
 * File:   EWphysics.cpp
 * Author: mishima
 * 
 * Created on February 28, 2011, 2:43 PM
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "EWphysics.h"


EWphysics::EWphysics(gslpp::complex gZf_i[10], gslpp::complex rhoZf_i[10],
                     double Delta_r_i,
                     double mu_i, double mc_i, double mt_i,
                     double md_i, double ms_i, double mb_i,
                     double mnu1_i, double mnu2_i, double mnu3_i,
                     double me_i, double mmu_i, double mtau_i,
                     double mZ_i, double mHl_i, double alsMz_i, double GF_i,
                     double ale_i, double aleMz_i){
    int INDF;
    for(INDF=0; INDF<10; INDF++) {
        gZf[INDF] = gZf_i[INDF];
        rhoZf[INDF] = rhoZf_i[INDF];
    }
    Delta_r = Delta_r_i;

    mu = mu_i;
    mc = mc_i;
    mt = mt_i;
    md = md_i;
    ms = ms_i;
    mb = mb_i;
    mnu1 = mnu1_i;
    mnu2 = mnu2_i;
    mnu3 = mnu3_i;
    me = me_i;
    mmu = mmu_i;
    mtau = mtau_i;

    mZ = mZ_i;
    mHl = mHl_i;
    alsMz = alsMz_i;
    GF = GF_i;
    ale = ale_i;
    aleMz = aleMz_i;
}

EWphysics::EWphysics(const EWphysics& orig) {
    int INDF;
    for(INDF=0; INDF<10; INDF++) {
        gZf[INDF] = orig.getGZf(INDF);
        rhoZf[INDF] = orig.getRhoZf(INDF);
    }
    Delta_r = orig.getDelta_r();

    /*
    mu = 
    mc = 
    mt = 
    md = 
    ms = 
    mb = 
    mnu1 = 
    mnu2 = 
    mnu3 = 
    me = 
    mmu = 
    mtau = 

    mZ = 
    mHl = 
    alsMz = 
    GF = 
    ale = 
    aleMz = 
     */
}

EWphysics::~EWphysics() {
}


///////////////////////////////////////////////////////////////////////////

int EWphysics::flavour_st_to_int(const std::string flavour) {
    if (flavour == "nu") {
        return 0;
    } else if (flavour == "e") {
        return 1;
    } else if (flavour == "mu") {
        return 2;
    } else if (flavour == "tau") {
        return 3;
    } else if (flavour == "u") {
        return 4;
    } else if (flavour == "d") {
        return 5;
    } else if (flavour == "c") {
        return 6;
    } else if (flavour == "s") {
        return 7;
    } else if (flavour == "t") {
        return 8;
    } else if (flavour == "b") {
        return 9;
    } else {
        std::cout << "flavour = nu, e, mu, tau, u, d, c, s, t, b" << std::endl;
        exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////////

double EWphysics::mW() {
    return ( mZ/2.0*(1.0
             + sqrt(1.0 - 4.0*M_PI*ale/sqrt(2.0)/mZ/mZ/GF*(1.0+Delta_r))) );
}

double EWphysics::Gamma_W() {
    return ( 0.0 ); //!!!!! 
}

double EWphysics::sw2() {
    return ( 1.0 - mW()*mW()/mZ/mZ );
}

double EWphysics::s2teff_f(const std::string flavour) {
    int INDF = flavour_st_to_int(flavour);
    double Qf[10] = {0., -1., -1., -1., 2./3., -1./3., 2./3., -1./3., 2./3., -1./3.};
    double Re_gZf = gZf[INDF].real();

    return ( (1.0-Re_gZf)/4.0/abs(Qf[INDF]) );
}

double EWphysics::Gamma_f(const std::string flavour) {
    int INDF = flavour_st_to_int(flavour);
    double m[10] = {0., me, mmu, mtau, mu, md, mc, ms, mt, mb};
    double Nc_f[10] = {1., 1., 1., 1., 3., 3., 3., 3., 3., 3.};
    double Qf[10] = {0., -1., -1., -1., 2./3., -1./3., 2./3., -1./3., 2./3., -1./3.};

    /* Radiator functions from the final-state QED and QCD corrections
     * to the vector and axial-vector currents */
    double RVf, RAf;
    /* non-factorizable EW-QCD corrections */
    double Delta_EWQCD;

    if (INDF < 4) {                  // leptons
        RVf = sqrt(1.0-4.0*m[INDF]*m[INDF]/mZ/mZ)
                * (1.0+2.0*m[INDF]*m[INDF]/mZ/mZ)
                * (1.0+3.0/4.0*aleMz/M_PI*Qf[INDF]*Qf[INDF]);
        RAf = RVf - sqrt(1.0-4.0*m[INDF]*m[INDF]/mZ/mZ)
                    * 6.0*m[INDF]*m[INDF]/mZ/mZ
                    * (1.0+3.0/4.0*aleMz/M_PI*Qf[INDF]*Qf[INDF]);
        Delta_EWQCD = 0.0;
    } else if (INDF==4 || INDF==6) { // u and c quarks
        RVf = 0.0; //!!!!! 
        RAf = 0.0; //!!!!!
        Delta_EWQCD = -0.000113;     // in GeV
    } else if (INDF==5 || INDF==7) { // d and s quarks
        RVf = 0.0; //!!!!!
        RAf = 0.0; //!!!!!
        Delta_EWQCD = -0.000160;     // in GeV
    } else if (INDF==9) {            // b quark
        RVf = 0.0; //!!!!!
        RAf = 0.0; //!!!!!
        Delta_EWQCD = -0.000040;     // in GeV
    } else if (INDF==8) {            // t quark
        return (0.0);
    } else {
        std::cout << "flavour = nu, e, mu, tau, u, d, c, s, t, b" << std::endl;
        exit(EXIT_FAILURE);
    }

    // If2 is already incorporated into gZf in the EWPOSM class
    //double If2 = 35.0/18.0*aleMz*aleMz*( 1 - 8.0/3.0*s2teff_f(flavour) );
    //gZf -= 4.0 * abs(Qf) * If2;

    double G0 = GF*mZ*mZ*mZ/24.0/sqrt(2.0)/M_PI;
    double gamma = Nc_f[INDF] * G0 * rhoZf[INDF].abs()
                   * ( gZf[INDF].abs2()*RVf + RAf ) + Delta_EWQCD;

    return gamma;
}

double EWphysics::Gamma_inv() {
    return ( 3.0*Gamma_f("nu") );
}

double EWphysics::Gamma_had() {
    return ( Gamma_f("u") + Gamma_f("c")
             + Gamma_f("d") + Gamma_f("s") + Gamma_f("b") );
}

double EWphysics::Gamma_Z() {
    return ( Gamma_f("e") + Gamma_f("mu") + Gamma_f("tau")
             + Gamma_inv() + Gamma_had() );
}

double EWphysics::sigma0_l(const std::string flavour_l) {
    if (flavour_l!="e" && flavour_l!="mu" && flavour_l!="tau") {
        std::cout << "flavour_l = e, mu, tau" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( 12.0*M_PI*Gamma_f("e")*Gamma_f(flavour_l)
            /mZ/mZ/Gamma_Z()/Gamma_Z() );
}

double EWphysics::sigma0_had() {
    return ( 12.0*M_PI*Gamma_f("e")*Gamma_had()
             /mZ/mZ/Gamma_Z()/Gamma_Z() );
}

double EWphysics::R0_l(const std::string flavour_l) {
    if (flavour_l!="e" && flavour_l!="mu" && flavour_l!="tau") {
        std::cout << "flavour_l = e, mu, tau" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( Gamma_had()/Gamma_f(flavour_l) );
}

double EWphysics::R0_q(const std::string flavour_q) {
    if (flavour_q!="b" && flavour_q!="c" && flavour_q!="s") {
        std::cout << "flavour_q = b, c, s" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( Gamma_f(flavour_q)/Gamma_had() );
}

double EWphysics::A_f(const std::string flavour) {
    double Re_gZf = gZf[flavour_st_to_int(flavour)].real();

    //!! alternatively use the two-loop result for sin^2\theta_eff^f
    return ( 2.0*Re_gZf/(Re_gZf*Re_gZf + 1.0) ); //!!!!!
}

double EWphysics::AFB0_f(const std::string flavour) {
    return ( 3.0/4.0*A_f("e")*A_f(flavour) );
}

double EWphysics::obliqueEpsilon1() {
    double Delta_rho = 2.0*(sqrt(rhoZf[1].real()) - 1.0);

    return Delta_rho;
}

double EWphysics::obliqueEpsilon2() {
    double Delta_rW = 1.0 - M_PI*aleMz/sqrt(2.0)/GF
                            /( 1.0 - mW()*mW()/mZ/mZ )/mW()/mW();
    double Delta_rho = 2.0*(rhoZf[1].real() - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f("e")/s02 - 1.0;

    //std::cout << "  Delta_rW = " << Delta_rW << std::endl; // TEST
    //std::cout << "  Delta_rho = " << Delta_rho << std::endl; // TEST
    //std::cout << "  Delta_k = " << Delta_k << std::endl; // TEST
    //std::cout << "  s02 = " << s02 << std::endl << std::endl; // TEST

    return ( c02*Delta_rho + s02*Delta_rW/(c02-s02) - 2.0*s02*Delta_k );
}

double EWphysics::obliqueEpsilon3() {
    double Delta_rho = 2.0*(sqrt(rhoZf[1].real()) - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f("e")/s02 - 1.0;

    return ( c02*Delta_rho + (c02-s02)*Delta_k );
}

double EWphysics::obliqueS() {
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);

    return ( obliqueEpsilon3()/ale*4.0*s02 );
}

double EWphysics::obliqueT() {
    return ( obliqueEpsilon1()/ale );
}

double EWphysics::obliqueU() {
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);

    return ( - obliqueEpsilon2()/ale*4.0*s02 );
}