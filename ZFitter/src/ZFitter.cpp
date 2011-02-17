/* 
 * File:   ZFitter.cpp
 * Author: mishima
 * 
 * Created on February 15, 2011, 11:15 PM
 */

#include <stdlib.h>
#include <iostream>
#include "ZFitter.h"


ZFitter::ZFitter(double ZMASS_i, double TMASS_i, double HMASS_i, double ALFAS_i,
                 double DAL5H_i, double ALFA_i, double V_TB_i, double UMASS_i,
                 double DMASS_i) {
    ZMASS = ZMASS_i;
    TMASS = TMASS_i;
    HMASS = HMASS_i;
    ALFAS = ALFAS_i;
    DAL5H = DAL5H_i;
    ALFA = ALFA_i;
    V_TB = V_TB_i;
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
    ALFA = orig.getALFA();
    V_TB = orig.getV_TB();
    UMASS = orig.getUMASS();
    DMASS = orig.getDMASS();

    init(0);
}

ZFitter::~ZFitter() {
}


///////////////////////////////////////////////////////////////////////////

double ZFitter::getCommonAlphas() { return zupars_.ALPHST; }

double ZFitter::getCommonSin2thetaW() { return zupars_.SIN2TW; }

double ZFitter::getCommonS2TEFF(int INDF) {
    if (INDF < 12) {
        return zupars_.S2TEFF[INDF];
    } else {
        std::cout << "S2TEFF[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitter::getCommonWidths(int INDF) {
    if (INDF < 12) {
        return zupars_.WIDTHS[INDF];
    } else {
        std::cout << "WIDTHS[INDF < 12]" << std::endl;
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

void ZFitter::calcXS_AFB(const int INDF, const double SQRS) {
    zvthsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            &XS[INDF], &AFB[INDF]);
}

void ZFitter::calcDXS(const int INDF, const double SQRS, const double CSA) {
    zvatsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &CSA,
            &DXS[INDF]);
}

void ZFitter::calcTauPol(const double SQRS) {
    zvtpsm_(&SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS,
            &TAUPOL, &TAUAFB);
}

void ZFitter::calcALR(const int INDF, const double SQRS, const double POL) {
    zvlrsm_(&INDF, &SQRS, &ZMASS, &TMASS, &HMASS, &DAL5H, &V_TB, &ALFAS, &POL, 
            &XSPL[INDF], &XSMI[INDF]);
}

void ZFitter::calcAPV() {
    double SIN2TW = getCommonSin2thetaW();
    zu_apv_(&ZMASS, &TMASS, &HMASS, &SIN2TW, &UMASS, &DMASS,
            &C1U, &C1D, &C2U, &C2D);
}

void ZFitter::calcXS(const int INDF, const double SQRS, const double GAMZ0,
                     const double GAMEE, const double GAMFF) {
    zuxsec_(&INDF, &SQRS, &ZMASS, &GAMZ0, &GAMEE, &GAMFF, &XS[INDF]);
}

void ZFitter::calcXS_AFB(const int INDF, const double SQRS, const double GAMZ0,
                         const int MODE, const double GVE, const double XE,
                         const double GVF, const double XF) {
    zuxsa_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, 
            &XS[INDF], &AFB[INDF]);
}

void ZFitter::calcXS_AFB(const int INDF, const double SQRS, const double GAMZ0,
                         const int MODE, const double GV2, const double X2) {
    zuxsa2_(&INDF, &SQRS, &ZMASS, &GAMZ0, &MODE, &GV2, &X2, 
            &XS[INDF], &AFB[INDF]);
}

void ZFitter::calcXS_AFB(const int INDF, const double SQRS, const double GAMZ0,
                         const double PFOUR, const double PVAE2,
                         const double PVAF2) {
    zuxafb_(&INDF, &SQRS, &ZMASS, &GAMZ0, &PFOUR, &PVAE2, &PVAF2, 
            &XS[INDF], &AFB[INDF]);
}

void ZFitter::calcTauPol(const double SQRS, const double GAMZ0, const int MODE,
                         const double GVE, const double XE, const double GVF,
                         const double XF) {
    zutau_(&SQRS, &ZMASS, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, &TAUPOL, &TAUAFB);
}


