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
#include <gsl/gsl_sf_zeta.h>
#include "EWphysics.h"


EWphysics::EWphysics(gslpp::complex gZf_i[10], gslpp::complex rhoZf_i[10],
                     double Delta_r_i,
                     double mcMz_i, double mbMz_i, double mt_i, 
                     double me_i, double mmu_i, double mtau_i,
                     double mZ_i, double mHl_i, double alsMz_i, double GF_i,
                     double ale_i, double aleMz_i){
    int INDF;
    for(INDF=0; INDF<10; INDF++) {
        gZf[INDF] = gZf_i[INDF];
        rhoZf[INDF] = rhoZf_i[INDF];
    }
    Delta_r = Delta_r_i;

    mcMz = mcMz_i;
    mbMz = mbMz_i;

    mt = mt_i;
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

//EWphysics::EWphysics(StandardModel& StandardModel_i) {
//}

//EWphysics::EWphysics(SUSY& SUSY_i) {
//}

//EWphysics::EWphysics(MFV& MFV_i) {
//}

EWphysics::EWphysics(const EWphysics& orig) {
    int INDF;
    for(INDF=0; INDF<10; INDF++) {
        gZf[INDF] = orig.getGZf(INDF);
        rhoZf[INDF] = orig.getRhoZf(INDF);
    }
    Delta_r = orig.getDelta_r();

    mcMz = orig.getMcMz();
    mbMz = orig.getMbMz();
    mt = orig.getMt();
    me = orig.getMe();
    mmu = orig.getMmu();
    mtau = orig.getMtau();

    mZ = orig.getMZ();
    mHl = orig.getMHl();
    alsMz = orig.getAlsMz();
    GF = orig.getGF();
    ale = orig.getAle();
    aleMz = orig.getAleMz();
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

std::string EWphysics::flavour_int_to_st(const int INDF) {
    std::string channel;
    if (INDF==0) channel="nu,nubar";
    if (INDF==1) channel="e+,e-";
    if (INDF==2) channel="mu+,mu-";
    if (INDF==3) channel="tau+,tau-";
    if (INDF==4) channel="u,ubar";
    if (INDF==5) channel="d,dbar";
    if (INDF==6) channel="c,cbar";
    if (INDF==7) channel="s,sbar";
    if (INDF==8) channel="t,tbar";
    if (INDF==9) channel="b,bbar";
    if (INDF==10) channel="hadron";
    if (INDF==11) channel="total";
    return channel;
}

double EWphysics::Qf(const int INDF) {
    if (INDF == 0) {
        return (0.0);
    } else if (INDF == 1 || INDF == 2 || INDF == 3) {
        return (-1.0);
    } else if (INDF == 4 || INDF == 6 || INDF == 8) {
        return (2.0/3.0);
    } else if (INDF == 5 || INDF == 7 || INDF == 9) {
        return (-1.0/3.0);
    } else {
        std::cout << "Qf(INDF = 0-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
}


///////////////////////////////////////////////////////////////////////////

double EWphysics::mW() {
    return ( mZ/sqrt(2.0)
             * sqrt(1.0 + sqrt(1.0 - 4.0*M_PI*ale/sqrt(2.0)/mZ/mZ/GF*(1.0+Delta_r))) );
}

double EWphysics::Gamma_W() {
    std::cout << "No code for Gamma_W()" << std::endl;
    return ( 0.0 ); //!!!!! 
}

double EWphysics::sw2() {
    return ( 1.0 - mW()*mW()/mZ/mZ );
}

double EWphysics::s2teff_f(const int INDF) {
    if (INDF < 0 || INDF > 9) {
        std::cout << "s2teff_f(INDF = 0-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (INDF==0 || INDF==8) return 0.0; /* neutrinos, top quark */
    double Re_gZf = gZf[INDF].real();
    //std::cout << abs(-1.0/3.0) << std::endl;      // TEST --> 0
    //std::cout << std::abs(-1.0/3.0) << std::endl; // TEST --> 0.33333
    return ( (1.0-Re_gZf)/4.0/std::abs(Qf(INDF)) );
}

double EWphysics::Gamma_l(const int INDF_l) {
    if (INDF_l < 0 || INDF_l > 3) {
        std::cout << "Gamma_l(INDF_l = 0-4)" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    double m_l[4] = {0., me, mmu, mtau};

    double xl = m_l[INDF_l]*m_l[INDF_l]/mZ/mZ;
    double G0 = GF*mZ*mZ*mZ/24.0/sqrt(2.0)/M_PI;
    double Gamma = G0 * rhoZf[INDF_l].abs() * sqrt(1.0 - 4.0*xl)
                   * ( (1.0 + 2.0*xl)*(gZf[INDF_l].abs2() + 1.0) - 6.0*xl )
                   * ( 1.0 + 3.0/4.0*aleMz/M_PI*Qf(INDF_l)*Qf(INDF_l) );
    return Gamma;
}

double EWphysics::Gamma_q(const int INDF_q) {
     if (INDF_q < 4 || INDF_q > 9) {
        std::cout << "Gamma_q(INDF_q = 4-9)" << std::endl;
        exit(EXIT_FAILURE);
    }

    /* Radiator functions from the final-state QED and QCD corrections
     * to the vector and axial-vector currents */
    double RVf, RAf;
    /* non-factorizable EW-QCD corrections in GeV */
    double Delta_EWQCD;
    /* z-component of isospin */
    double I3q;
    /* electric charge squared */
    double Qf2 = Qf(INDF_q)*Qf(INDF_q);

    if (INDF_q==4 || INDF_q==6) { /* up and charm */
        Delta_EWQCD = -0.000113;
        I3q = 1.0/2.0;
    } else if (INDF_q==5 || INDF_q==7) { /* down and strange */
        Delta_EWQCD = -0.000160;
        I3q = -1.0/2.0;
    } else if (INDF_q==9) { /* bottom */
        Delta_EWQCD = -0.000040;
        I3q = -1.0/2.0;
    } else if (INDF_q==8) { /* top */
        return (0.0);
    }

    /* s = mZ^2 */
    double s = mZ*mZ;

    /* products of the charm and bottom masses at mZ */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    if (INDF_q==9) {
        mqMz2 = mbMz*mbMz;
        mqdash4 = mcMz2*mcMz2;
    } else if (INDF_q==6) {
        mqMz2 = mcMz*mcMz;
        mqdash4 = mbMz2*mbMz2;
    } else {
        mqMz2 = 0.0;
        mqdash4 = 0.0;
    }

    /* logs */
    double log_t = log(mt*mt/s);
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    if (INDF_q==9 || INDF_q==6) {
        log_q = log(mqMz2/s);
    } else {
        log_q = 0.0;
    }

    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = gsl_sf_zeta_int(2);
    double zeta3 = gsl_sf_zeta_int(3);
    double zeta4 = gsl_sf_zeta_int(4);
    double zeta5 = gsl_sf_zeta_int(5);

    /* massless non-singlet corrections */
    double C02 = 365.0/24.0 - 11.0*zeta3 + (-11.0/12.0 + 2.0/3.0*zeta3)*nf;
    double C03 = 87029.0/288.0 - 121.0/8.0*zeta2 - 1103.0/4.0*zeta3
                 + 275.0/6.0*zeta5 
                 + (-7847.0/216.0 + 11.0/6.0*zeta2 + 262.0/9.0*zeta3
                    - 25.0/9.0*zeta5)*nf
                 + (151.0/162.0 - zeta2/18.0 - 19.0/27.0*zeta3)*nf*nf;
    double C04 = -156.61 + 18.77*nf - 0.7974*nf*nf + 0.0215*nf*nf*nf;
    //std::cout << "TEST: C02 = " << C02 << std::endl;// TEST (should be 1.40923)
    //std::cout << "TEST: C03 = " << C03 << std::endl;// TEST (should be -12.7671)
    //std::cout << "TEST: C04 = " << C04 << std::endl;// TEST (should be -80.0075)

    /* quadratic massive corrections */
    double C23  = -80.0 + 60.0*zeta3 + (32.0/9.0 - 8.0/3.0*zeta3)*nf;
    double C21V = 12.0;
    double C22V = 253.0/2.0 - 13.0/3.0*nf;
    double C23V = 2522.0 - 855.0/2.0*zeta2 + 310.0/3.0*zeta3 - 5225.0/6.0*zeta5
                  + (-4942.0/27.0 + 34.0*zeta2 - 394.0/27.0*zeta3
                     + 1045.0/27.0*zeta5)*nf
                  + (125.0/54.0 - 2.0/3.0*zeta2)*nf*nf;
    double C20A = -6.0;
    double C21A = -22.0;
    double C22A = -8221.0/24.0 + 57.0*zeta2 + 117.0*zeta3
                  + (151.0/12.0 - 2.0*zeta2 - 4.0*zeta3)*nf;
    double C23A = -4544045.0/864.0 + 1340.0*zeta2 + 118915.0/36.0*zeta3
                  - 127.0*zeta5
                  + (71621.0/162.0 - 209.0/2.0*zeta2 - 216.0*zeta3
                     + 5.0*zeta4 + 55.0*zeta5)*nf
                  + (-13171.0/1944.0 + 16.0/9.0*zeta2 + 26.0/9.0*zeta3)*nf*nf;

    /* quartic massive corrections */
    double C42  = 13.0/3.0 - 4.0*zeta3;
    //double C40V = -6.0; /* not used */
    double C41V = -22.0;
    double C42V = -3029.0/12.0 + 162.0*zeta2 + 112.0*zeta3
                  + (143.0/18.0 - 4.0*zeta2 - 8.0/3.0*zeta3)*nf;
    double C42VL= -11.0/2.0 + nf/3.0;
    double C40A = 6.0;
    double C41A = 10.0;
    double C42A = 3389.0/12.0 - 162.0*zeta2 - 220.0*zeta3
                  + (-41.0/6.0 + 4.0*zeta2 + 16.0/3.0*zeta3)*nf;
    double C42AL= 77.0/2.0 - 7.0/3.0*nf;

    /* power suppressed top-mass correction */
    double xt = s/mt/mt;
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* singlet axial-vector corrections */
    double I2 = -37.0/12.0 + (-log_t) + 7.0/81.0*xt + 0.0132*xt*xt;
    double I3 = -5075.0/216.0 + 23.0/6.0*zeta2 + zeta3 + 67.0/18.0*(-log_t)
                + 23.0/12.0*log_t*log_t;

    /* singlet vector corrections */
    //double RVh; /* not used */

    /* rescaled strong coupling constant */
    double alsMzPi  = alsMz/M_PI;
    double alsMzPi2 = alsMzPi*alsMzPi;
    double alsMzPi3 = alsMzPi2*alsMzPi;
    double alsMzPi4 = alsMzPi3*alsMzPi;
    double alsMzPi5 = alsMzPi4*alsMzPi;
    double alsMzPi6 = alsMzPi5*alsMzPi;

    /* radiator function to the vector current */
    RVf = 1.0 + 3.0/4.0*Qf2*aleMz/M_PI + alsMzPi - Qf2/4.0*aleMz/M_PI*alsMzPi
            + (C02 + C2t)*alsMzPi2 + C03*alsMzPi3 + C04*alsMzPi4
            //+ deltaC05*alsMzPi5 /* theoretical uncertainty (absent in ZFITTER) */
            + (mcMz2 + mbMz2)/s*C23*alsMzPi3
            + mqMz2/s*(C21V*alsMzPi + C22V*alsMzPi2 + C23V*alsMzPi3)
            + mcMz2*mcMz2/s/s*(C42 - log_c)*alsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*alsMzPi2
            + mqMz2*mqMz2/s/s*(C41V*alsMzPi + (C42V + C42VL*log_q)*alsMzPi2)
            + 12.0*mqdash4/s/s*alsMzPi2
            - mqMz2*mqMz2*mqMz2/s/s/s
              *(8.0+16.0/27.0*(155.0 + 6.0*log_q)*alsMzPi);

    /* radiator function to the axial-vector current */
    RAf = 1.0 + 3.0/4.0*Qf2*aleMz/M_PI + alsMzPi - Qf2/4.0*aleMz/M_PI*alsMzPi
            + (C02 + C2t - 2.0*I3q*I2)*alsMzPi2
            + (C03 - 2.0*I3q*I3)*alsMzPi3
            + C04*alsMzPi4 /* (absent in ZFITTER) */
            //- 2.0*I3q*deltaI4*alsMzPi4 /* theoretical uncertainty */
            //+ deltaC05*alsMzPi5 /* theoretical uncertainty (absent in ZFITTER) */
            + (mcMz2 + mbMz2)/s*C23*alsMzPi3
            + mqMz2/s*(C20A + C21A*alsMzPi + C22A*alsMzPi2
                       + 6.0*(3.0 + log_t)*alsMzPi2 + C23A*alsMzPi3)
            - 10.0*mqMz2/mt/mt*(8.0/81.0 + log_t/54.0)*alsMzPi2
            + mcMz2*mcMz2/s/s*(C42 - log_c)*alsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*alsMzPi2
            + mqMz2*mqMz2/s/s*(C40A + C41A*alsMzPi
                               + (C42A + C42AL*log_q)*alsMzPi2)
            - 12.0*mqdash4/s/s*alsMzPi2 ;

    double G0 = GF*mZ*mZ*mZ/24.0/sqrt(2.0)/M_PI;
    double Gamma = 3.0 * G0 * rhoZf[INDF_q].abs()
                   * ( gZf[INDF_q].abs2()*RVf + RAf ) + Delta_EWQCD;
    return Gamma;
}

double EWphysics::Gamma_f(const int INDF) {
    if (INDF >= 0 && INDF <= 3) {
        return ( Gamma_l(INDF) );
    } else if (INDF >= 4 && INDF <= 9) {
        return ( Gamma_q(INDF) );
    } else {
        std::cout << "Gamma_f(INDF = 0-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double EWphysics::Gamma_inv() {
    return ( 3.0*Gamma_f(0) );
}

double EWphysics::Gamma_had() {
    return ( Gamma_f(4) + Gamma_f(5) + Gamma_f(6) + Gamma_f(7) + Gamma_f(9) );
}

double EWphysics::Gamma_Z() {
    return ( Gamma_f(1) + Gamma_f(2) + Gamma_f(3) + Gamma_inv() + Gamma_had() );
}

double EWphysics::sigma0_l(const int INDF_l) {
    if (INDF_l < 0 || INDF_l > 3) {
        std::cout << "sigma0_l(INDF_l = 0-3)" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( 12.0*M_PI*Gamma_f(1)*Gamma_f(INDF_l)
            /mZ/mZ/Gamma_Z()/Gamma_Z() );
}

double EWphysics::sigma0_had() {
    return ( 12.0*M_PI*Gamma_f(1)*Gamma_had()
             /mZ/mZ/Gamma_Z()/Gamma_Z() );
}

double EWphysics::R0_l(const int INDF_l) {
    if (INDF_l < 0 || INDF_l > 3) {
        std::cout << "R0_l(INDF_l = 0-3)" << std::endl;
        exit(EXIT_FAILURE);
    }  
    return ( Gamma_had()/Gamma_f(INDF_l) );
}

double EWphysics::R0_q(const int INDF_q) {
    if (INDF_q < 4 || INDF_q > 9) {
        std::cout << "R0_q(INDF_q = 4-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( Gamma_f(INDF_q)/Gamma_had() );
}

double EWphysics::A_f(const int INDF) {
    if (INDF < 0 || INDF > 9) {
        std::cout << "A_f(INDF = 0-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
    double Re_gZf = gZf[INDF].real();
    return ( 2.0*Re_gZf/(Re_gZf*Re_gZf + 1.0) ); //!!!!!
}

double EWphysics::AFB0_f(const int INDF) {
    if (INDF < 0 || INDF > 9) {
        std::cout << "AFB0_f(INDF = 0-9)" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( 3.0/4.0*A_f(1)*A_f(INDF) );
}

double EWphysics::obliqueEpsilon1() {
    double Delta_rho = 2.0*(sqrt(rhoZf[1].real()) - 1.0);
    return Delta_rho;
}

double EWphysics::obliqueEpsilon2() {
    double Delta_rW = 1.0 - M_PI*aleMz/sqrt(2.0)/GF
                            /( 1.0 - mW()*mW()/mZ/mZ )/mW()/mW();
    double Delta_rho = 2.0*(sqrt(rhoZf[1].real()) - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f(1)/s02 - 1.0;
    return ( c02*Delta_rho + s02*Delta_rW/(c02-s02) - 2.0*s02*Delta_k );
}

double EWphysics::obliqueEpsilon3() {
    double Delta_rho = 2.0*(sqrt(rhoZf[1].real()) - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*aleMz/sqrt(2.0)/GF/mZ/mZ);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f(1)/s02 - 1.0;
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


///////////////////////////////////////////////////////////////////////////

double EWphysics::s2teff_f(const std::string flavour) {
    return ( s2teff_f(flavour_st_to_int(flavour)) );
}

double EWphysics::Gamma_l(const std::string flavour_l) {
    return ( Gamma_l(flavour_st_to_int(flavour_l)) );
}

double EWphysics::Gamma_q(const std::string flavour_q) {
    return ( Gamma_q(flavour_st_to_int(flavour_q)) );
}

double EWphysics::Gamma_f(const std::string flavour) {
    return ( Gamma_f(flavour_st_to_int(flavour)) );
}

double EWphysics::sigma0_l(const std::string flavour_l) {
    return ( sigma0_l(flavour_st_to_int(flavour_l)) );
}

double EWphysics::R0_l(const std::string flavour_l) {
    return ( R0_l(flavour_st_to_int(flavour_l)) );
}

double EWphysics::R0_q(const std::string flavour_q) {
    return ( R0_q(flavour_st_to_int(flavour_q)) );
}

double EWphysics::A_f(const std::string flavour) {
    return ( A_f(flavour_st_to_int(flavour)) );
}

double EWphysics::AFB0_f(const std::string flavour) {
    return ( AFB0_f(flavour_st_to_int(flavour)) );
}

