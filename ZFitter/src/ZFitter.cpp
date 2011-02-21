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

void ZFitter::printIntermediateResults() {
    std::cout << "----- Intermediate Results -----" << std::endl;

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


///////////////////////////////////////////////////////////////////////////

int ZFitter::flavour_st_to_int(const std::string flavour) {
    if (flavour == "e") {
        return 1;
    } else if (flavour == "mu") {
        return 2;
    } else if (flavour == "tau") {
        return 3;
    } else if (flavour == "b") {
        return 9;
    } else if (flavour == "c") {
        return 6;
    } else if (flavour == "s") {
        return 7;
    } else {
        std::cout << "flavour = e, mu, tau, b, c, s" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::mW() {
    return ( ZMASS*sqrt(1.0 - getCommonSIN2TW()) );
}

double ZFitter::Gamma_W() {
    return ( getCommonPARTW(2) );
}

double ZFitter::sw2() {
    return ( getCommonSIN2TW() );
}

double ZFitter::Gamma_f(const std::string flavour) {
    return ( getCommonWIDTHS(flavour_st_to_int(flavour)) );
}

double ZFitter::Gamma_inv() {
    return ( 3.0*getCommonWIDTHS(0) );
}

double ZFitter::Gamma_had() {
    return ( getCommonWIDTHS(10) );
}

double ZFitter::Gamma_Z() {
    return ( getCommonWIDTHS(11) );
}

double ZFitter::sigma0_l(const std::string flavour_l) {
    if (flavour_l!="e" && flavour_l!="mu" && flavour_l!="tau") {
        std::cout << "flavour_l = e, mu, tau" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( 12.0*M_PI*Gamma_f("e")*Gamma_f(flavour_l)
            /ZMASS/ZMASS/Gamma_Z()/Gamma_Z() );
}

double ZFitter::sigma0_had() {
    return ( 12.0*M_PI*Gamma_f("e")*Gamma_had()
             /ZMASS/ZMASS/Gamma_Z()/Gamma_Z() );
}

double ZFitter::R0_l(const std::string flavour_l) {
    if (flavour_l!="e" && flavour_l!="mu" && flavour_l!="tau") {
        std::cout << "flavour_l = e, mu, tau" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( Gamma_had()/Gamma_f(flavour_l) );
}

double ZFitter::R0_q(const std::string flavour_q) {
    if (flavour_q!="b" && flavour_q!="c" && flavour_q!="s") {
        std::cout << "flavour_q = b, c, s" << std::endl;
        exit(EXIT_FAILURE);
    }
    return ( Gamma_f(flavour_q)/Gamma_had() );
}

double ZFitter::A_f(const std::string flavour) {
    double Re_gZ_f = getCommonARVEFZ(flavour_st_to_int(flavour));
    return ( 2.0*Re_gZ_f/(Re_gZ_f*Re_gZ_f + 1.0) );
}

double ZFitter::AFB0_f(const std::string flavour) {
    return ( 3.0/4.0*A_f("e")*A_f(flavour) );
}

double ZFitter::s2teff_f(const std::string flavour) {
    return ( getCommonARSEFZ(flavour_st_to_int(flavour)) );
}

double ZFitter::obliqueEpsilon1() {
    double Delta_rho = 2.0*(sqrt(getCommonAROTFZ(1)) - 1.0);

    return Delta_rho;
}

double ZFitter::obliqueEpsilon2() {
    double alpha_mZ = calqed_.ALQEDZ;
    double Delta_rW = 1.0 - M_PI*alpha_mZ/sqrt(2.0)/GF()
                            /( 1.0 - mW()*mW()/ZMASS/ZMASS )/mW()/mW();
    double Delta_rho = 2.0*(sqrt(getCommonAROTFZ(1)) - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*alpha_mZ/sqrt(2.0)/GF()/ZMASS/ZMASS);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f("e")/s02 - 1.0;

    //std::cout << "  Delta_rW = " << Delta_rW << std::endl; // TEST
    //std::cout << "  Delta_rho = " << Delta_rho << std::endl; // TEST
    //std::cout << "  Delta_k = " << Delta_k << std::endl; // TEST
    //std::cout << "  s02 = " << s02 << std::endl << std::endl; // TEST

    return ( c02*Delta_rho + s02*Delta_rW/(c02-s02) - 2.0*s02*Delta_k );
}

double ZFitter::obliqueEpsilon3() {
    double alpha_mZ = calqed_.ALQEDZ; 
    double Delta_rho = 2.0*(sqrt(getCommonAROTFZ(1)) - 1.0);
    double s02 = 0.5 - sqrt(0.25 - M_PI*alpha_mZ/sqrt(2.0)/GF()/ZMASS/ZMASS);
    double c02 = 1.0 - s02;
    double Delta_k = s2teff_f("e")/s02 - 1.0;

    return ( c02*Delta_rho + (c02-s02)*Delta_k );
}

double ZFitter::obliqueS() {
    double alpha_mZ = calqed_.ALQEDZ;
    double s02 = 0.5 - sqrt(0.25 - M_PI*alpha_mZ/sqrt(2.0)/GF()/ZMASS/ZMASS);
    
    return ( obliqueEpsilon3()/alpha()*4.0*s02 );
}

double ZFitter::obliqueT() {
    return ( obliqueEpsilon1()/alpha() );
}

double ZFitter::obliqueU() {
    double alpha_mZ = calqed_.ALQEDZ;
    double s02 = 0.5 - sqrt(0.25 - M_PI*alpha_mZ/sqrt(2.0)/GF()/ZMASS/ZMASS);
    
    return ( - obliqueEpsilon2()/alpha()*4.0*s02 );
}

void ZFitter::printPO() {

    // GeV^{-2} --> nb
    const double GeVminus2_to_nb = pow(10.0, -6.0)
                                   / pow(10.0, -28.0)
                                   / pow(299792458.0, -2.0)
                                   / pow(6.58211899 * pow(10.0, -22.0), -2.0)
                                   * pow(10.0, 9.0);

    std::cout << std::setw(15) << "m_W [GeV]" << std::setw(13)
              << mW() << std::endl
              << std::setw(15) << "Gamma_W [GeV]" << std::setw(13)
              << Gamma_W() << std::endl
              << std::setw(15) << "sin^2(th_W)" << std::setw(13)
              << sw2() << std::endl
              << std::setw(15) << "Gamma_inv [GeV]" << std::setw(13)
              << Gamma_inv() << std::endl
              << std::setw(15) << "Gamma_had [GeV]" << std::setw(13)
              << Gamma_had() << std::endl
              << std::setw(15) << "Gamma_Z [GeV]" << std::setw(13)
              << Gamma_Z() << std::endl
              << std::setw(15) << "sigma0_e [nb]" << std::setw(13)
              << sigma0_l("e")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_mu [nb]" << std::setw(13)
              << sigma0_l("mu")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_tau [nb]" << std::setw(13)
              << sigma0_l("tau")*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_had [nb]" << std::setw(13)
              << sigma0_had()*GeVminus2_to_nb << std::endl
              << std::setw(15) << "R0_e" << std::setw(13)
              << R0_l("e") << std::endl
              << std::setw(15) << "R0_mu" << std::setw(13)
              << R0_l("mu") << std::endl
              << std::setw(15) << "R0_tau" << std::setw(13)
              << R0_l("tau") << std::endl
              << std::setw(15) << "R0_b" << std::setw(13)
              << R0_q("b") << std::endl
              << std::setw(15) << "R0_c" << std::setw(13)
              << R0_q("c") << std::endl
              << std::setw(15) << "R0_s" << std::setw(13)
              << R0_q("s") << std::endl
              << std::setw(15) << "A_e" << std::setw(13)
              << A_f("e") << std::endl
              << std::setw(15) << "A_mu" << std::setw(13)
              << A_f("mu") << std::endl
              << std::setw(15) << "A_tau" << std::setw(13)
              << A_f("tau") << std::endl
              << std::setw(15) << "A_b" << std::setw(13)
              << A_f("b") << std::endl
              << std::setw(15) << "A_c" << std::setw(13)
              << A_f("c") << std::endl
              << std::setw(15) << "A_s" << std::setw(13)
              << A_f("s") << std::endl
              << std::setw(15) << "AFB0_e" << std::setw(13)
              << AFB0_f("e") << std::endl
              << std::setw(15) << "AFB0_mu" << std::setw(13)
              << AFB0_f("mu") << std::endl
              << std::setw(15) << "AFB0_tau" << std::setw(13)
              << AFB0_f("tau") << std::endl
              << std::setw(15) << "AFB0_b" << std::setw(13)
              << AFB0_f("b") << std::endl
              << std::setw(15) << "AFB0_c" << std::setw(13)
              << AFB0_f("c") << std::endl
              << std::setw(15) << "AFB0_s" << std::setw(13)
              << AFB0_f("s") << std::endl
              << std::setw(15) << "sin^2(teff_e)" << std::setw(13)
              << s2teff_f("e") << std::endl
              << std::setw(15) << "sin^2(teff_mu)" << std::setw(13)
              << s2teff_f("mu") << std::endl
              << std::setw(15) << "sin^2(teff_tau)" << std::setw(13)
              << s2teff_f("tau") << std::endl
              << std::setw(15) << "sin^2(teff_b)" << std::setw(13)
              << s2teff_f("b") << std::endl
              << std::setw(15) << "sin^2(teff_c)" << std::setw(13)
              << s2teff_f("c") << std::endl
              << std::setw(15) << "sin^2(teff_s)" << std::setw(13)
              << s2teff_f("s") << std::endl
              << std::endl;

    std::cout << "  epsilon1 = " << obliqueEpsilon1() << std::endl
              << "  epsilon2 = " << obliqueEpsilon2() << std::endl
              << "  epsilon3 = " << obliqueEpsilon3() << std::endl
              << "         S = " << obliqueS() << std::endl
              << "         T = " << obliqueT() << std::endl
              << "         U = " << obliqueU() << std::endl
              << std::endl;
}


///////////////////////////////////////////////////////////////////////////

double ZFitter::GF() const {
    //return (1.166388e-5); // for GFER=0
    //return (1.16639e-5);  // for GFER=1
    return (1.16637e-5);  // for GFER=2
}

double ZFitter::alpha() const {
    return (1.0/137.0359895);
}


