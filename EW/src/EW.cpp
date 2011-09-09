/* 
 * File:   EW.cpp
 * Author: mishima
 */

#include "EW.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_sf_zeta.h>


EW::EW(const StandardModel& SM_i) : ThObsType(SM_i), myEWSM(SM_i), myZFitter(SM_i) {

    //mcMz = 0.580624; // TEST!!!
    //mbMz = 2.84386; // TEST!!!

    //mcMz = 0.55696435; 
    //mbMz = 2.8095955;


    mcMz = 0.55381685; 
    mbMz = 2.8194352;


}

//EW::EW(const EW& orig) : ThObsType(orig.SM), EWSM(orig.SM) {
//}

EW::~EW() {
}


////////////////////////////////////////////////////////////////////////

void EW::ComputeEWSM(const schemes_EW schemeMw, 
                     const schemes_EW schemeRhoZ,
                     const schemes_EW schemeKappaZ,
                     const bool flag_order[EWSM::orders_EW_size]) {

    myEWSM.ComputeDeltaAlpha(flag_order);
    
    DeltaAlpha_l5q = myEWSM.getDeltaAlpha_l(EWSM::orders_EW_size) + SM.getDAle5Mz();
    DeltaAlpha = DeltaAlpha_l5q + myEWSM.getDeltaAlpha_t(EWSM::orders_EW_size);
    
#define TEST_DEBUG
#ifdef TEST_DEBUG
    DeltaAlpha = 0.10;
#endif    
    
    alphaMz = SM.getAle()/(1.0 - DeltaAlpha);
    
    if (schemeMw==APPROXIMATEFORMULA) {
        myApproximateFormulae = new ApproximateFormulae(SM, DeltaAlpha);
        Mw = myApproximateFormulae->Mw();        
        delete myApproximateFormulae;
    } else {
        myEWSM.ComputeCC(SM.Mw_tree(), flag_order);
        Mw = resumMw(schemeMw);
        
        /* Mw from iterations */
        double Mw_org = 0.0;
        while (fabs(Mw - Mw_org) > 0.0000001) {
            Mw_org = Mw;
            myEWSM.ComputeCC(Mw, flag_order);
            Mw = resumMw(schemeMw);
            /* TEST */
            int prec_def = std::cout.precision();
            std::cout << std::setprecision(12) << "TEST: Mw_org = " << Mw_org 
                      << "  Mw_new = " << Mw << std::endl;
            std::cout.precision(prec_def);
        }
    }

    myEWSM.ComputeNC(Mw, flag_order);

    /* Resummations */    
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;   
        double deltaRho_rem_l_real[EWSM::orders_EW_size+1],
               deltaRho_rem_q_real[EWSM::orders_EW_size+1],
               deltaKappa_rem_l_real[EWSM::orders_EW_size+1],
               deltaKappa_rem_q_real[EWSM::orders_EW_size+1];
        for (int j=0; j<=EWSM::orders_EW_size; j++) {
            EWSM::orders_EW j_order = (EWSM::orders_EW) j;
            deltaRho_rem_l_real[j] = myEWSM.getDeltaRho_rem_l(i_l,j_order).real();
            deltaRho_rem_q_real[j] = myEWSM.getDeltaRho_rem_q(i_q,j_order).real(); 
            deltaKappa_rem_l_real[j] = myEWSM.getDeltaRho_rem_l(i_l,j_order).real();
            deltaKappa_rem_q_real[j] = myEWSM.getDeltaRho_rem_q(i_q,j_order).real();
        }
        rhoZ_l[i].real() = resumRhoZ(schemeRhoZ, deltaRho_rem_l_real);
        rhoZ_q[i].real() = resumRhoZ(schemeRhoZ, deltaRho_rem_q_real);
        kappaZ_l[i].real() = resumKappaZ(schemeKappaZ, deltaKappa_rem_l_real);
        kappaZ_q[i].real() = resumKappaZ(schemeKappaZ, deltaKappa_rem_q_real);
    }

    /* Imaginary parts */
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;   
        rhoZ_l[i].imag() = myEWSM.getDeltaRho_rem_l(i_l,EWSM::orders_EW_size).imag();    
        rhoZ_q[i].imag() = myEWSM.getDeltaRho_rem_q(i_q,EWSM::orders_EW_size).imag();    
        kappaZ_l[i].imag() = myEWSM.getDeltaKappa_rem_l(i_l,EWSM::orders_EW_size).imag();    
        kappaZ_q[i].imag() = myEWSM.getDeltaKappa_rem_q(i_q,EWSM::orders_EW_size).imag();            
    }
    
    /* Other contributions to Im[kappa_Z^f] */
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double cW2 = 1.0 - sW2;
    for (int i=0; i<6; i++) {
//        kappaZ_l[i].imag() -= SM.getAle()*SM.getAlsMz()/24.0/M_PI*(cW2-sW2)/sW2/sW2; //???
//        kappaZ_q[i].imag() -= SM.getAle()*SM.getAlsMz()/24.0/M_PI*(cW2-sW2)/sW2/sW2; //???
    
//        kappaZ_l[i].imag() -= 3.0*Xt/f_AlphaToGF*pow()*
//        kappaZ_q[i].imag() -=    
    }
    
    /* Using the approximate formula for the real parts of kappa_Z^f*/
    if (schemeKappaZ==EW::APPROXIMATEFORMULA) {
        myApproximateFormulae = new ApproximateFormulae(SM, DeltaAlpha);
        double sin2thetaEff_l[6], sin2thetaEff_q[6];
        double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
        for (int i=0; i<6; i++) {
            StandardModel::lepton i_l = (StandardModel::lepton) i;
            StandardModel::quark i_q = (StandardModel::quark) i;        
            sin2thetaEff_l[i] = myApproximateFormulae->sin2thetaEff(i_l);
            sin2thetaEff_q[i] = myApproximateFormulae->sin2thetaEff(i_q);
            kappaZ_l[i].real() = sin2thetaEff_l[i]/sW2; 
            kappaZ_q[i].real() = sin2thetaEff_q[i]/sW2;
        }
        delete myApproximateFormulae;
    }    
}

void EW::ComputeZFitter(const schemes_EW schemeMw, 
                        const schemes_EW schemeRhoZ,
                        const schemes_EW schemeKappaZ,
                        const bool flag_order[EWSM::orders_EW_size]) {

    SetZFitterFlags(schemeMw, schemeRhoZ, schemeKappaZ, flag_order);
    myZFitter.calcCommonBlocks();
    
    /* Outputs */
    alphaMz = myZFitter.getAlphaMZ();
    DeltaAlpha = 1.0 - 1.0/alphaMz*SM.getAle();
    DeltaAlpha_l5q = DeltaAlpha - SM.getDAle5Mz();
    Mw = myZFitter.getMw();
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        rhoZ_l[i] = myZFitter.getRhoZ_l(i_l);
        rhoZ_q[i] = myZFitter.getRhoZ_q(i_q);
        kappaZ_l[i] = myZFitter.getKappaZ_l(i_l);
        kappaZ_q[i] = myZFitter.getKappaZ_q(i_q);
    }
}


////////////////////////////////////////////////////////////////////////

double EW::Qf(const StandardModel::lepton l) const {
    switch(l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return (0.0);
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return (-1.0);
        default:
            throw "Error in EW::Qf()";  
    }
}

double EW::Qf(const StandardModel::quark q) const {
    switch(q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
        case StandardModel::TOP:
            return (2.0/3.0);
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
        case StandardModel::BOTTOM:
            return (-1.0/3.0);
        default:
            throw "Error in EW::Qf()";  
    }
}

double EW::sin2thetaEff(const StandardModel::lepton l) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    return (  kappaZ_l[l].real() * sW2 );
}

double EW::sin2thetaEff(const StandardModel::quark q) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    return (  kappaZ_q[q].real() * sW2 );
}

double EW::Gamma_l(const StandardModel::lepton l) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    complex gV_over_gA = 1.0 - 4.0*fabs(Qf(l))*kappaZ_l[l]*sW2;
    double xl = pow(SM.getLeptons(l).getMass()/SM.getMz(), 2.0);
    double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;
    double Gamma = G0*rhoZ_l[l].abs()*sqrt(1.0 - 4.0*xl)
                   * ( (1.0 + 2.0*xl)*(gV_over_gA.abs2() + 1.0) - 6.0*xl )
                   * ( 1.0 + 3.0/4.0*alphaMz/M_PI*Qf(l)*Qf(l) );
    return Gamma;
}

double EW::Gamma_q(const StandardModel::quark q) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    complex gV_over_gA = 1.0 - 4.0*fabs(Qf(q))*kappaZ_q[q]*sW2;

    /* Radiator functions from the final-state QED and QCD corrections
     * to the vector and axial-vector currents */
    double RVf, RAf;
    /* non-factorizable EW-QCD corrections in GeV*/
    double Delta_EWQCD;
    /* z-component of isospin */
    double I3q;
    /* electric charge squared */
    double Qf2 = Qf(q)*Qf(q);

    switch(q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
            Delta_EWQCD = -0.000113;
            I3q = 1.0/2.0;
            break;
        case StandardModel::TOP:
            return 0.0;
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            Delta_EWQCD = -0.000160;
            I3q = -1.0/2.0;
            break;
        case StandardModel::BOTTOM:
            Delta_EWQCD = -0.000040;
            I3q = -1.0/2.0;
            break;
        default:
            throw "Error in EW::Gamma_q()";  
    }

    /* s = Mz^2 */
    double s = SM.getMz()*SM.getMz();

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch(q) {
        case StandardModel::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case StandardModel::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* logs */
    double log_t = log(pow(SM.getQuarks(SM.TOP).getMass(),2.0)/s);
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    switch(q) {
        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            log_q = log(mqMz2/s);
            break;
        default:
            log_q = 0.0;
            break;
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
    double xt = s/pow(SM.getQuarks(SM.TOP).getMass(),2.0);
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* singlet axial-vector corrections */
    double I2 = -37.0/12.0 + (-log_t) + 7.0/81.0*xt + 0.0132*xt*xt;
    double I3 = -5075.0/216.0 + 23.0/6.0*zeta2 + zeta3 + 67.0/18.0*(-log_t)
                + 23.0/12.0*log_t*log_t;

    /* singlet vector corrections */
    //double RVh; /* not used */
    
    /* rescaled strong coupling constant */
    double AlsMzPi  = SM.getAlsMz()/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    /* radiator function to the vector current */
    RVf = 1.0 + 3.0/4.0*Qf2*alphaMz/M_PI + AlsMzPi - Qf2/4.0*alphaMz/M_PI*AlsMzPi
            + (C02 + C2t)*AlsMzPi2 + C03*AlsMzPi3 + C04*AlsMzPi4
            //+ deltaC05*AlsMzPi5 /* theoretical uncertainty (absent in ZFITTER) */
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C21V*AlsMzPi + C22V*AlsMzPi2 + C23V*AlsMzPi3)
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C41V*AlsMzPi + (C42V + C42VL*log_q)*AlsMzPi2)
            + 12.0*mqdash4/s/s*AlsMzPi2
            - mqMz2*mqMz2*mqMz2/s/s/s
              *(8.0+16.0/27.0*(155.0 + 6.0*log_q)*AlsMzPi);

    /* radiator function to the axial-vector current */
    RAf = 1.0 + 3.0/4.0*Qf2*alphaMz/M_PI + AlsMzPi - Qf2/4.0*alphaMz/M_PI*AlsMzPi
            + (C02 + C2t - 2.0*I3q*I2)*AlsMzPi2
            + (C03 - 2.0*I3q*I3)*AlsMzPi3
            + C04*AlsMzPi4 /* (absent in ZFITTER) */
            //- 2.0*I3q*deltaI4*AlsMzPi4 /* theoretical uncertainty */
            //+ deltaC05*AlsMzPi5 /* theoretical uncertainty (absent in ZFITTER) */
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C20A + C21A*AlsMzPi + C22A*AlsMzPi2
                       + 6.0*(3.0 + log_t)*AlsMzPi2 + C23A*AlsMzPi3)
            - 10.0*mqMz2/pow(SM.getQuarks(SM.TOP).getMass(),2.0)
              *(8.0/81.0 + log_t/54.0)*AlsMzPi2
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C40A + C41A*AlsMzPi
                               + (C42A + C42AL*log_q)*AlsMzPi2)
            - 12.0*mqdash4/s/s*AlsMzPi2 ;    
    
    double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;    
    double Gamma = 3.0* G0*rhoZ_q[q].abs()
                   * ( gV_over_gA.abs2()*RVf + RAf ) + Delta_EWQCD;
    return Gamma;
}

double EW::Gamma_inv() const {
    return ( Gamma_l(SM.NEUTRINO_1) + Gamma_l(SM.NEUTRINO_2) 
             + Gamma_l(SM.NEUTRINO_3) );
}

double EW::Gamma_had() const {
    return ( Gamma_q(SM.UP) + Gamma_q(SM.DOWN) + Gamma_q(SM.CHARM)
             + Gamma_q(SM.STRANGE) + Gamma_q(SM.BOTTOM) );
}

double EW::Gamma_Z() const {
    return ( Gamma_l(SM.ELECTRON) + Gamma_l(SM.MU) + Gamma_l(SM.TAU) 
             + Gamma_inv() + Gamma_had() );
}

double EW::sigma0_l(const StandardModel::lepton l) const {
    return ( 12.0*M_PI*Gamma_l(SM.ELECTRON)*Gamma_l(l)
             /SM.getMz()/SM.getMz()/Gamma_Z()/Gamma_Z() );
}

double EW::sigma0_had() const {
    return ( 12.0*M_PI*Gamma_l(SM.ELECTRON)*Gamma_had()
             /SM.getMz()/SM.getMz()/Gamma_Z()/Gamma_Z() );
}

double EW::A_l(const StandardModel::lepton l) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double Re_gV_over_gA = 1.0 - 4.0*fabs(Qf(l))*kappaZ_l[l].real()*sW2;
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}

double EW::A_q(const StandardModel::quark q) const {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double Re_gV_over_gA = 1.0 - 4.0*fabs(Qf(q))*kappaZ_q[q].real()*sW2;
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}


////////////////////////////////////////////////////////////////////////

double EW::resumMw(const schemes_EW schemeMw) {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double cW2 = 1.0 - sW2;
    double cW2_to_sW2 = cW2/sW2;
    double R;

    if (myEWSM.getDeltaR_rem(EWSM::EW1QCD2)!=0.0) throw "Error in EW::resumMw()";
    if (myEWSM.getDeltaR_rem(EWSM::EW2QCD1)!=0.0) throw "Error in EW::resumMw()";
    if (myEWSM.getDeltaR_rem(EWSM::EW3)!=0.0) throw "Error in EW::resumMw()";
    
    switch (schemeMw) {
        case NORESUM: 
            // R = 1 + DeltaRho
            R = 1.0 + DeltaAlpha_l5q 
                - cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size)
                + myEWSM.getDeltaR_rem(EWSM::orders_EW_size);
            break;
        case OMSI:
            // R = 1/(1 - DeltaRho)
            R = 1.0/(1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size))
                /(1.0 - DeltaAlpha_l5q 
                  - myEWSM.getDeltaR_rem(EWSM::EW1)
                  - myEWSM.getDeltaR_rem(EWSM::EW1QCD1) 
                  - myEWSM.getDeltaR_rem(EWSM::EW2));
            break;
        case INTERMEDIATE:
            // R = 1/(1 - DeltaRho)
            R = 1.0/( (1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size))
                      *(1.0 - DeltaAlpha_l5q - myEWSM.getDeltaR_rem(EWSM::EW1)) 
                      - myEWSM.getDeltaR_rem(EWSM::EW1QCD1) 
                      - myEWSM.getDeltaR_rem(EWSM::EW2) );
            break;        
        case OMSII:
            // R = 1/(1 - DeltaRho)
            R = 1.0/( (1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size))
                      *(1.0 - DeltaAlpha_l5q)
                      - (1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::EW1))
                        *myEWSM.getDeltaR_rem(EWSM::EW1)
                      - myEWSM.getDeltaR_rem(EWSM::EW1QCD1)
                      - myEWSM.getDeltaR_rem(EWSM::EW2) );
            break;
        default:
            throw "Error in EW::resumMw()";            
            break;
    }   

    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    if (tmp*R > 1.0) throw "Negative (1-tmp*R) in EW::resumMw()";
    
    return (SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R)));
}

double EW::resumRhoZ(const schemes_EW schemeRhoZ, 
                     const double deltaRho_rem[EWSM::orders_EW_size]) {
    double rhoZ;
    double deltaRho_rem_G = deltaRho_rem[EWSM::EW1]
                            + deltaRho_rem[EWSM::EW1QCD1];
    if (deltaRho_rem[EWSM::EW1QCD2]!=0.0) throw "Error in EW::resumRhoZ()";
    if (deltaRho_rem[EWSM::EW2QCD1]!=0.0) throw "Error in EW::resumRhoZ()";    
    if (deltaRho_rem[EWSM::EW3]!=0.0) throw "Error in EW::resumRhoZ()";    
    
    /* Real parts */
    switch (schemeRhoZ) {
        case NORESUM: 
            rhoZ = 1.0 + myEWSM.getDeltaRho(EWSM::orders_EW_size)
                   + deltaRho_rem[EWSM::orders_EW_size];
            break;
        case OMSI:
            rhoZ = (1.0 + deltaRho_rem_G + deltaRho_rem[EWSM::EW2])
                   /(1.0 - myEWSM.getDeltaRho(EWSM::orders_EW_size)
                           *(1.0 - myEWSM.getDeltaRbar_rem()));
            break;
        case INTERMEDIATE:
            rhoZ = (1.0 + deltaRho_rem_G)
                   /(1.0 - myEWSM.getDeltaRho(EWSM::orders_EW_size)
                            *(1.0 - myEWSM.getDeltaRbar_rem()))
                   + deltaRho_rem[EWSM::EW2];            
            break;        
        case OMSII:
            rhoZ = 1.0 + myEWSM.getDeltaRho(EWSM::orders_EW_size)
                   + pow(myEWSM.getDeltaRho(EWSM::EW1), 2.0) 
                   - myEWSM.getDeltaRho(EWSM::EW1)*myEWSM.getDeltaRbar_rem()
                   + deltaRho_rem_G*(1.0 + myEWSM.getDeltaRho(EWSM::EW1))
                   + deltaRho_rem[EWSM::EW2];  
            break;
        default:
            throw "Error in EW::resumRhoZ()";
            break;
    }

    return rhoZ;
}

double EW::resumKappaZ(const schemes_EW schemeKappaZ, 
                       const double deltaKappa_rem[EWSM::orders_EW_size]) {
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double cW2 = 1.0 - sW2;
    double cW2_to_sW2 = cW2/sW2;
    double kappaZ;
    double deltaKappa_rem_G = deltaKappa_rem[EWSM::EW1] 
                              + deltaKappa_rem[EWSM::EW1QCD1] 
                              + deltaKappa_rem[EWSM::EW1QCD2];
    if (deltaKappa_rem[EWSM::EW2QCD1]!=0.0) throw "Error in EW::resumKappaZ()";    
    if (deltaKappa_rem[EWSM::EW3]!=0.0) throw "Error in EW::resumKappaZ()";     
    
    /* Real parts */
    switch (schemeKappaZ) {
        case NORESUM: 
            kappaZ = 1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size)
                     + deltaKappa_rem[EWSM::orders_EW_size];
            break;
        case OMSI:
            kappaZ = (1.0 + deltaKappa_rem_G 
                      + deltaKappa_rem[EWSM::EW2])
                     *(1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size)
                             *(1.0 - myEWSM.getDeltaRbar_rem()));
            break;
        case INTERMEDIATE:
            kappaZ = (1.0 + deltaKappa_rem_G)
                     *(1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size)
                             *(1.0 - myEWSM.getDeltaRbar_rem()))
                     + deltaKappa_rem[EWSM::EW2];
            break;        
        case OMSII:
            kappaZ = 1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::orders_EW_size)
                     - cW2_to_sW2*myEWSM.getDeltaRho(EWSM::EW1)
                       *myEWSM.getDeltaRbar_rem()
                     + deltaKappa_rem_G
                       *(1.0 + cW2_to_sW2*myEWSM.getDeltaRho(EWSM::EW1))
                     + deltaKappa_rem[EWSM::EW2];
            break;
        case APPROXIMATEFORMULA:
            /* The real parts are given by the approximate formulae. 
             * See ComputeKappaZ() */
            kappaZ = 0.0; // dummy
            break;
        default:
            throw "Error in EW::resumKappaZ()";
            break;
    }

    return kappaZ;
}

void EW::SetZFitterFlags(const schemes_EW schemeMw, 
                         const schemes_EW schemeRhoZ,
                         const schemes_EW schemeKappaZ,
                         const bool flag_order[EWSM::orders_EW_size]) {
    // DAL5H is supplied by the user as input. 
    myZFitter.flag("ALEM", 2); 

    if (schemeMw==NORESUM 
        && schemeRhoZ==NORESUM && schemeKappaZ==NORESUM) {
            myZFitter.flag("AMT4", 0); // Does this option work correctly? 
    } else if (schemeMw==OMSI   
               && schemeRhoZ==OMSI && schemeKappaZ==OMSI) {    
            myZFitter.flag("AMT4", 4);
            myZFitter.flag("IFACR", 0);    
            myZFitter.flag("IFACT", 0);
    } else if (schemeMw==INTERMEDIATE
               && schemeRhoZ==INTERMEDIATE && schemeKappaZ==INTERMEDIATE) {    
            myZFitter.flag("AMT4", 4);
            myZFitter.flag("IFACR", 1);    
            myZFitter.flag("IFACT", 1);            
    } else if (schemeMw==OMSII  
               && schemeRhoZ==OMSII && schemeKappaZ==OMSII) {    
            myZFitter.flag("AMT4", 4);
            myZFitter.flag("IFACR", 2);    
            myZFitter.flag("IFACT", 2);
     } else if (schemeMw==APPROXIMATEFORMULA
               && schemeKappaZ==APPROXIMATEFORMULA) {
            myZFitter.flag("AMT4", 6);
            if (schemeRhoZ==OMSI) myZFitter.flag("IFACT", 0);
            if (schemeRhoZ==INTERMEDIATE) myZFitter.flag("IFACT", 1);
            if (schemeRhoZ==OMSII) myZFitter.flag("IFACT", 2);    
    } else {
        throw "Write codes in EW::ComputeZFitter()";
    }

    for (int i=0; i<EWSM::orders_EW_size; i++) {
        if (flag_order[i]!=true) 
        throw "Invalid flag_order[] in EW::SetZFitterFlags()";
    }    
}