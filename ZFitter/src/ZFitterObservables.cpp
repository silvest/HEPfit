/* 
 * File:   ZFitterObservables.cpp
 * Author: mishima
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ZFitterObservables.h"


ZFitterObservables::ZFitterObservables(const ZFitter& ZFitter_i) : ZFitter(ZFitter_i) {
}

//ZFitterObservables::ZFitterObservables(const ZFitterObservables& orig) {
//}

ZFitterObservables::~ZFitterObservables() {
}


////////////////////////////////////////////////////////////////////////

double ZFitterObservables::getCommonWIDTHS(const int INDF) {
    if (INDF < 12) {
        return zupars_.WIDTHS[INDF]*0.001;
    } else {
        std::cout << "WIDTHS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterObservables::getCommonPARTZ(const int INDF) {
    if (INDF < 12) {
        return partzw_.PARTZ[INDF]*0.001;
    } else {
        std::cout << "PARTZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterObservables::getCommonPARTW(const int i) {
    if (i < 3) {
        return partzw_.PARTW[i]*0.001;
    } else {
        std::cout << "PARTW[i < 3]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////

void ZFitterObservables::setCuts(const int INDF, const int ICUT, const double ACOL,
                                 const double EMIN, const double S_PR, const double ANG0,
                                 const double ANG1, const double SIPP) {
    zucuts_(&INDF, &ICUT, &ACOL, &EMIN, &S_PR, &ANG0, &ANG1, &SIPP);
}

void ZFitterObservables::setAllCuts(const int ICUT[12], const double ACOL[12],
                                    const double EMIN[12], const double S_PR[12],
                                    const double ANG0[12], const double ANG1[12],
                                    const double SPP[12], const int flagPrint) {
    /* set cuts */
    int indexFermion;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        setCuts(indexFermion, ICUT[indexFermion], ACOL[indexFermion],
                EMIN[indexFermion], S_PR[indexFermion], ANG0[indexFermion],
                ANG1[indexFermion], SPP[indexFermion]);
    }

    /* print cuts if flagPrint=1 */
    if (flagPrint==1) {
        std::cout << "  Channel ICUT  ACOL   EMIN    S_PR  ANG0  ANG1   SPP"
                  << std::endl;
        for (indexFermion=0; indexFermion<12; indexFermion++) {
            std::cout << std::setw(9) << convertINDF(indexFermion)
                      << std::setw(4) << ICUT[indexFermion]
                      << std::setw(7) << ACOL[indexFermion]
                      << std::setw(7) << EMIN[indexFermion]
                      << std::setw(9) << S_PR[indexFermion]
                      << std::setw(5) << ANG0[indexFermion]
                      << std::setw(5) << ANG1[indexFermion]
                      << std::setw(9) << SPP[indexFermion] << std::endl;
        }
    }
}

void ZFitterObservables::CutInfo() { 
    const int mode = 1;
    zuinfo_(&mode); 
}


////////////////////////////////////////////////////////////////////////

void ZFitterObservables::calcXS_AFB(const int INDF, const double SQRS,
                                    double *XS, double *AFB) {
    zvthsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            XS, AFB);
}

void ZFitterObservables::calcDXS(const int INDF, const double SQRS, 
                                 const double CSA, double *DXS) {
    zvatsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &CSA,
            DXS);
}

void ZFitterObservables::calcTauPol(const double SQRS, double *TAUPOL, double *TAUAFB) {
    zvtpsm_(&SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            TAUPOL, TAUAFB);
}

void ZFitterObservables::calcALR(const int INDF, const double SQRS, const double POL,
                                 double *XSPL, double *XSMI) {
    zvlrsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &POL, 
            XSPL, XSMI);
}

void ZFitterObservables::calcAPV(double *C1U, double *C1D, double *C2U, double *C2D) {
                                 double SIN2TW = getCommonSIN2TW();
    zu_apv_(&ZMASS, &TMASS, &HMASS, &SIN2TW, &UMASS, &DMASS,
            C1U, C1D, C2U, C2D);
}

void ZFitterObservables::calcXS(const int INDF, const double SQRS, const double GAMZ0,
                                const double GAMEE, const double GAMFF, double *XS) {
    zuxsec_(&INDF, &SQRS, &ZMASS, &GAMZ0, &GAMEE, &GAMFF, XS);
}

void ZFitterObservables::calcXS_AFB_2(const int INDF, const double SQRS, const double GAMZ0,
                                      const int MODE, const double GVE, const double XE,
                                      const double GVF, const double XF,
                                      double *XS, double *AFB) {
    zuxsa_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, XS, AFB);
}

void ZFitterObservables::calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                                      const int MODE, const double GV2, const double X2,
                                      double *XS, double *AFB) {
    zuxsa2_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GV2, &X2, XS, AFB);
}

void ZFitterObservables::calcXS_AFB_4(const int INDF, const double SQRS, const double GAMZ0,
                                      const double PFOUR, const double PVAE2,
                                      const double PVAF2, double *XS, double *AFB) {
    zuxafb_(&INDF, &SQRS, &ZMASS, &GAMZ0, &PFOUR, &PVAE2, &PVAF2, XS, AFB);
}

void ZFitterObservables::calcTauPol_2(const double SQRS, const double GAMZ0, const int MODE,
                                      const double GVE, const double XE, const double GVF,
                                      const double XF, double *TAUPOL, double *TAUAFB) {
    zutau_(&SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, TAUPOL, TAUAFB);
}

void ZFitterObservables::test(const int IMISC) { zftest_(&IMISC); }


////////////////////////////////////////////////////////////////////////

double ZFitterObservables::Gamma_W() {
    return ( getCommonPARTW(2) );
}

double ZFitterObservables::sw2() {
    return ( getCommonSIN2TW() );
}

double ZFitterObservables::s2teff_f(const int INDF) {
    return ( getCommonARSEFZ(INDF) );
}

double ZFitterObservables::Gamma_f(const int INDF) {
    return ( getCommonWIDTHS(INDF) );
}

double ZFitterObservables::Gamma_inv() {
    return ( 3.0*getCommonWIDTHS(0) );
}

double ZFitterObservables::Gamma_had() {
    return ( getCommonWIDTHS(10) );
}

double ZFitterObservables::Gamma_Z() {
    return ( getCommonWIDTHS(11) );
}


////////////////////////////////////////////////////////////////////////

std::string ZFitterObservables::convertINDF(const int INDF) {
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

void ZFitterObservables::printConstants() {
    std::cout << "------------ Constants ------------" << std::endl;
    std::cout << "  Channel    charges       masses" << std::endl;
    int indexFermion;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        std::cout << std::setw(9) << convertINDF(indexFermion)
                  << std::setw(11) << getCommonALLCH(indexFermion)
                  << std::setw(13) << getCommonALLMS(indexFermion) << std::endl;
    }
    std::cout << std::endl;
}

void ZFitterObservables::printInputs() {
    std::cout << "----- Input parameters -----" << std::endl;
    std::cout << "  ZMASS = " << ZMASS << std::endl;
    std::cout << "  TMASS = " << TMASS << std::endl;
    std::cout << "  HMASS = " << HMASS << std::endl;
    std::cout << "  ALFAS = " << ALFAS << std::endl;
    std::cout << "  DAL5H = " << DAL5H << std::endl;
    std::cout << "  V_TB  = " << V_TB << std::endl;
    std::cout << "  UMASS = " << UMASS << std::endl;
    std::cout << "  DMASS = " << DMASS << std::endl;
    std::cout << std::endl;
}

void ZFitterObservables::printIntermediateResults() {
    std::cout << "----- Intermediate Results -----" << std::endl;

    std::cout << "Mw = " << Mw << "  sin^2(theta_W) = "
              << getCommonSIN2TW() << std::endl << std::endl;

    std::cout << "Z decays:" << std::endl;
    std::cout << "  Channel Gamma (rho_Z^f)' Re[rho_Z^f] Im[rho_Z^f] "
              << "Re[g_Z^f] Im[g_Z^f] Re[k_Z^f] Im[k_Z^f] "
              << "sin^2(th_eff^f)[ARSEFZ] sin^2(th_eff^f)[S2TEFF]"
              << std::endl;
    int indexFermion;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        std::cout << std::setw(9) << convertINDF(indexFermion)
                  //<< std::setw(11) << getCommonPARTZ(indexFermion) // TEST
                  << std::setw(11) << getCommonWIDTHS(indexFermion)
                  << std::setw(9)  << getCommonARROFZ(indexFermion)
                  << std::setw(9)  << getCommonAROTFZ(indexFermion)
                  << std::setw(13) << getCommonAIROFZ(indexFermion)
                  << std::setw(11) << getCommonARVEFZ(indexFermion)
                  << std::setw(12) << getCommonAIVEFZ(indexFermion)
                  << std::setw(9)  << getCommonARKAFZ(indexFermion)
                  << std::setw(11) << getCommonAIKAFZ(indexFermion)
                  << std::setw(10) << getCommonARSEFZ(indexFermion)
                  //<< std::setw(10) << getCommonS2TEFF(indexFermion) // TEST
                  << std::endl;
    }
    for (indexFermion=10; indexFermion<12; indexFermion++) {
        std::cout << std::setw(9) << convertINDF(indexFermion)
                  << std::setw(10) << getCommonWIDTHS(indexFermion)
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "W decays:" << std::endl;
    std::cout << "  Gamma(W->leptons) = " << getCommonPARTW(0) << std::endl;
    std::cout << "  Gamma(W->quarks)  = " << getCommonPARTW(1) << std::endl;
    std::cout << "  Gamma(W->total)   = " << getCommonPARTW(2) << std::endl;
    std::cout << std::endl;

    /* test for sin^2(theta_eff^f) */
    std::cout << "  Channel sin^2(th_eff^f)[ARSEFZ] sin^2(th_eff^f)[S2TEFF] "
              << "Re[k_Z^f]*sin^2(theta_W) " << std::endl;
    for (indexFermion = 0; indexFermion < 10; indexFermion++) {
        std::cout << std::setw(9) << convertINDF(indexFermion)
                  << std::setw(10) << getCommonARSEFZ(indexFermion)
                  << std::setw(10) << getCommonS2TEFF(indexFermion) 
                  << std::setw(10)
                  << getCommonARKAFZ(indexFermion)*getCommonSIN2TW()
                  << std::endl;
    }
    std::cout << std::endl;

}



