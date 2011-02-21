/* 
 * File:   ZFitter.cpp
 * Author: mishima
 * 
 * Created on February 15, 2011, 11:15 PM
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ZFitter.h"


ZFitter::ZFitter(double ZMASS_i, double TMASS_i, double HMASS_i, double ALFAS_i,
                 double DAL5H_i, double V_TB_i, double UMASS_i, double DMASS_i) {
    ZMASS = ZMASS_i;
    TMASS = TMASS_i;
    HMASS = HMASS_i;
    ALFAS = ALFAS_i;
    DAL5H = DAL5H_i;
    V_TB  = V_TB_i;
    UMASS = UMASS_i;
    DMASS = DMASS_i;

    init(0);
}

ZFitter::ZFitter(const ZFitter& orig) {
    ZMASS = orig.getZMASS();
    TMASS = orig.getTMASS();
    HMASS = orig.getHMASS();
    ALFAS = orig.getALFAS();
    DAL5H = orig.getDAL5H();
    V_TB  = orig.getV_TB();
    UMASS = orig.getUMASS();
    DMASS = orig.getDMASS();

    init(0);
}

ZFitter::~ZFitter() {
}


///////////////////////////////////////////////////////////////////////////

double ZFitter::getCommonALPHST() { return zupars_.ALPHST; }

double ZFitter::getCommonSIN2TW() { return zupars_.SIN2TW; }

double ZFitter::getCommonS2TEFF(const int INDF) {
    if (INDF < 12) {
        return zupars_.S2TEFF[INDF];
    } else {
        std::cout << "S2TEFF[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonWIDTHS(const int INDF) {
    if (INDF < 12) {
        return zupars_.WIDTHS[INDF]*0.001;
    } else {
        std::cout << "WIDTHS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonALLCH(const int INDF) {
    if (INDF < 12) {
        return zfchms_.ALLCH[INDF];
    } else {
        std::cout << "ALLCH[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonALLMS(const int INDF) {
    if (INDF < 12) {
        return zfchms_.ALLMS[INDF];
    } else {
        std::cout << "ALLMS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonPARTZ(const int INDF) {
    if (INDF < 12) {
        return partzw_.PARTZ[INDF]*0.001;
    } else {
        std::cout << "PARTZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonPARTW(const int i) {
    if (i < 3) {
        return partzw_.PARTW[i]*0.001;
    } else {
        std::cout << "PARTW[i < 3]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


///////////////////////////////////////////////////////////////////////////

double ZFitter::getCommonARROFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARROFZ[INDF];
    } else {
        std::cout << "ARROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARKAFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARKAFZ[INDF];
    } else {
        std::cout << "ARKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARVEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARVEFZ[INDF];
    } else {
        std::cout << "ARVEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonARSEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.ARSEFZ[INDF];
    } else {
        std::cout << "ARSEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAROTFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AROTFZ[INDF];
    } else {
        std::cout << "AROTFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIROFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIROFZ[INDF];
    } else {
        std::cout << "AIROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIKAFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIKAFZ[INDF];
    } else {
        std::cout << "AIKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonAIVEFZ(const int INDF) {
    if (INDF < 11) {
        return cdzrkz_.AIVEFZ[INDF];
    } else {
        std::cout << "AIVEFZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


///////////////////////////////////////////////////////////////////////////

void ZFitter::test(const int IMISC) { zftest_(&IMISC); }

void ZFitter::init(const int IPRINT) { zuinit_(&IPRINT); }

void ZFitter::info(const int MODE) { zuinfo_(&MODE); }

void ZFitter::flag(const std::string CHFLAG, const int IVALUE) {
    zuflag_(CHFLAG.c_str(), &IVALUE, CHFLAG.size());
}

void ZFitter::calcCommonBlocks() {
    zvweak_(&ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS);
}

void ZFitter::setCuts(const int INDF, const int ICUT, const double ACOL,
                      const double EMIN, const double S_PR, const double ANG0,
                      const double ANG1, const double SIPP) {
    zucuts_(&INDF, &ICUT, &ACOL, &EMIN, &S_PR, &ANG0, &ANG1, &SIPP);
}

void ZFitter::calcXS_AFB(const int INDF, const double SQRS,
                         double *XS, double *AFB) {
    zvthsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            XS, AFB);
}

void ZFitter::calcDXS(const int INDF, const double SQRS, const double CSA,
                      double *DXS) {
    zvatsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &CSA,
            DXS);
}

void ZFitter::calcTauPol(const double SQRS, double *TAUPOL, double *TAUAFB) {
    zvtpsm_(&SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            TAUPOL, TAUAFB);
}

void ZFitter::calcALR(const int INDF, const double SQRS, const double POL,
                      double *XSPL, double *XSMI) {
    zvlrsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &POL, 
            XSPL, XSMI);
}

void ZFitter::calcAPV(double *C1U, double *C1D, double *C2U, double *C2D) {
    double SIN2TW = getCommonSIN2TW();
    zu_apv_(&ZMASS, &TMASS, &HMASS, &SIN2TW, &UMASS, &DMASS,
            C1U, C1D, C2U, C2D);
}

void ZFitter::calcXS(const int INDF, const double SQRS, const double GAMZ0,
                     const double GAMEE, const double GAMFF, double *XS) {
    zuxsec_(&INDF, &SQRS, &ZMASS, &GAMZ0, &GAMEE, &GAMFF, XS);
}

void ZFitter::calcXS_AFB_2(const int INDF, const double SQRS, const double GAMZ0,
                           const int MODE, const double GVE, const double XE,
                           const double GVF, const double XF,
                           double *XS, double *AFB) {
    zuxsa_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, XS, AFB);
}

void ZFitter::calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                           const int MODE, const double GV2, const double X2,
                           double *XS, double *AFB) {
    zuxsa2_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GV2, &X2, XS, AFB);
}

void ZFitter::calcXS_AFB_4(const int INDF, const double SQRS, const double GAMZ0,
                           const double PFOUR, const double PVAE2,
                           const double PVAF2, double *XS, double *AFB) {
    zuxafb_(&INDF, &SQRS, &ZMASS, &GAMZ0, &PFOUR, &PVAE2, &PVAF2, XS, AFB);
}

void ZFitter::calcTauPol_2(const double SQRS, const double GAMZ0, const int MODE,
                           const double GVE, const double XE, const double GVF,
                           const double XF, double *TAUPOL, double *TAUAFB) {
    zutau_(&SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, TAUPOL, TAUAFB);
}


///////////////////////////////////////////////////////////////////////////

void ZFitter::setAllFlags(const int flags[46], const int flagPrint) {

    std::string flagNames[46]
            = {"AFBC", "SCAL", "SCRE", "AMT4", "BORN",
               "BOXD", "CONV", "FINR", "FOT2", "GAMS",
               "DIAG", "INTF", "BARB", "PART", "POWR",
               "PRNT", "ALEM", "QCDC", "VPOL", "WEAK",
               "FTJR", "EXPR", "EXPF", "HIGS", "AFMT",
               "CZAK", "PREC", "HIG2", "ALE2", "GFER",
               "ISPP", "FSRS", "MISC", "MISD", "IPFC",
               "IPSC", "IPTO", "FBHO", "FSPP", "FUNA",
               "ASCR", "SFSR", "ENUE", "TUPV", "DMWW",
               "DSWW" };

    /* set flags */
    int i;
    for (i=0; i<46; i++) {
        flag(flagNames[i].c_str(), flags[i]);
    }

    /* print flags if flagPrint=1 */
    if (flagPrint==1) {
        for (i=0; i<46; i++) {
            std::cout << "  " << flagNames[i].c_str() << ": " << flags[i];
            if ((i+1)%5==0) std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
}

void ZFitter::setAllCuts(const int ICUT[12], const double ACOL[12],
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

std::string ZFitter::convertINDF(const int INDF) {
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

void ZFitter::printConstants() {
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

void ZFitter::printInputs() {
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

void ZFitter::printResults() {
    std::cout << "----- Results -----" << std::endl;

    double mW = ZMASS*sqrt(1.0-getCommonSIN2TW());
    std::cout << "mW = " << mW << "  sin^2(theta_W) = "
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

void ZFitter::calcPO(double *mW, double *Gamma_W, double *sw2,
                     double *Gamma_inv, double *Gamma_had, double *Gamma_Z,
                     double *sigma0_e, double *sigma0_mu, double *sigma0_tau,
                     double *sigma0_had,
                     double *R0_e, double *R0_mu, double *R0_tau,
                     double *R0_b, double *R0_c,
                     double *A_e, double *A_mu, double *A_tau,
                     double *A_b, double *A_c,
                     double *AFB0_e, double *AFB0_mu, double *AFB0_tau,
                     double *AFB0_b, double *AFB0_c,
                     double *s2teff_e, double *s2teff_mu, double *s2teff_tau,
                     double *s2teff_b, double *s2teff_c) {

    *mW = ZMASS*sqrt(1.0-getCommonSIN2TW());
    *Gamma_W =getCommonPARTW(2);
    *sw2 = getCommonSIN2TW();

    double Gamma_e = getCommonWIDTHS(1);
    double Gamma_mu = getCommonWIDTHS(2);
    double Gamma_tau = getCommonWIDTHS(3);
    double Gamma_hadron = getCommonWIDTHS(10);
    double Gamma_total = getCommonWIDTHS(11);

    *Gamma_inv = 3.0*getCommonWIDTHS(0);
    *Gamma_had = Gamma_hadron;
    *Gamma_Z = Gamma_total;

    *sigma0_e = 12.0*M_PI*Gamma_e*Gamma_e
                /ZMASS/ZMASS/Gamma_total/Gamma_total;
    *sigma0_mu = 12.0*M_PI*Gamma_e*Gamma_mu
                 /ZMASS/ZMASS/Gamma_total/Gamma_total;
    *sigma0_tau = 12.0*M_PI*Gamma_e*Gamma_tau
                  /ZMASS/ZMASS/Gamma_total/Gamma_total;
    *sigma0_had = 12.0*M_PI*Gamma_e*Gamma_hadron
                  /ZMASS/ZMASS/Gamma_total/Gamma_total;

    *R0_e = Gamma_hadron/Gamma_e;
    *R0_mu = Gamma_hadron/Gamma_mu;
    *R0_tau = Gamma_hadron/Gamma_tau;
    *R0_b = getCommonWIDTHS(9)/Gamma_hadron;
    *R0_c = getCommonWIDTHS(6) / Gamma_hadron;

    double Re_gZ_e = getCommonARVEFZ(1);
    double Re_gZ_mu = getCommonARVEFZ(2);
    double Re_gZ_tau = getCommonARVEFZ(3);
    double Re_gZ_b = getCommonARVEFZ(9);
    double Re_gZ_c = getCommonARVEFZ(6);
    *A_e = 2.0*Re_gZ_e/(Re_gZ_e*Re_gZ_e + 1.0);
    *A_mu = 2.0*Re_gZ_mu/(Re_gZ_mu*Re_gZ_mu + 1.0);
    *A_tau = 2.0*Re_gZ_tau/(Re_gZ_tau*Re_gZ_tau + 1.0);
    *A_b = 2.0*Re_gZ_b/(Re_gZ_b*Re_gZ_b + 1.0);
    *A_c = 2.0*Re_gZ_c/(Re_gZ_c*Re_gZ_c + 1.0);

    *AFB0_e = 3.0/4.0* (*A_e) * (*A_e);
    *AFB0_mu = 3.0/4.0* (*A_e) * (*A_mu);
    *AFB0_tau = 3.0/4.0* (*A_e) * (*A_tau);
    *AFB0_b = 3.0/4.0* (*A_e) * (*A_b);
    *AFB0_c = 3.0/4.0* (*A_e) * (*A_c);

    *s2teff_e = getCommonARSEFZ(1);
    *s2teff_mu = getCommonARSEFZ(2);
    *s2teff_tau = getCommonARSEFZ(3);
    *s2teff_b = getCommonARSEFZ(9);
    *s2teff_c = getCommonARSEFZ(6);
}

void ZFitter::calcSTU(double *S, double *T, double *U) {

    

}

