/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "ZFitter-f2c.h"
#include "ZFitterWrapper.h"
#include <stdexcept>


ZFitterWrapper::ZFitterWrapper(const StandardModel& SM_i)
: ThObsType(SM_i), SM(SM_i)
{
    MzCache = 0.0;
    MtCache = 0.0;
    MhCache = 0.0; 
    AlsMzCache = 0.0; 
    DAle5MzCache = 0.0;

    // constants
    V_TB = 1.0;
    UMASS = 0.1;
    DMASS = 0.1;

    init(0);

    // Flags [default: AMT4=4, ISPP=2, ALEM = 3, FINR = 1]
    int  AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
         BOXD = 1, CONV = 1, FINR = 0, FOT2 = 3, GAMS = 1,
         DIAG = 1, INTF = 0, BARB = 2, PART = 0, POWR = 1,
         PRNT = 0, ALEM = 2, QCDC = 3, VPOL = 1, WEAK = 1,
         FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
         CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
         ISPP = 2, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
         IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
         ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
         DSWW = 0;

    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    setAllFlags(flags, 0);
}


////////////////////////////////////////////////////////////////////////

void ZFitterWrapper::setFlag(const std::string CHFLAG, const int IVALUE) const
{
    zuflag_(CHFLAG.c_str(), &IVALUE, CHFLAG.size());
    IsFlagChanged = true;
}

void ZFitterWrapper::setAllFlags(const int flags[46], const int flagPrint) const
{

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
        setFlag(flagNames[i].c_str(), flags[i]);
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

void ZFitterWrapper::FlagInfo() const
{ 
    const int mode = 0;
    zuinfo_(&mode); 
}


////////////////////////////////////////////////////////////////////////

void ZFitterWrapper::setCuts(const int INDF, const int ICUT, const double ACOL,
                             const double EMIN, const double S_PR, const double ANG0,
                             const double ANG1, const double SIPP) const
{
    zucuts_(&INDF, &ICUT, &ACOL, &EMIN, &S_PR, &ANG0, &ANG1, &SIPP);
}

void ZFitterWrapper::setAllCuts(const int ICUT[12], const double ACOL[12],
                                const double EMIN[12], const double S_PR[12],
                                const double ANG0[12], const double ANG1[12],
                                const double SPP[12], const int flagPrint) const
{
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

void ZFitterWrapper::setSprimeCut(const int INDF, const double s) const
{
    // Cuts
    int ICUT;
    double ACOL, EMIN, S_PR, ANG0, ANG1, SPP;
    ICUT = -1;
    //ICUT = 1;
    ACOL = 0.0;// for ICUT=0,2,3
    EMIN = 0.0;// for ICUT=0,2,3
    S_PR = 0.85*0.85*s;// for ICUT=-1,1
    ANG0 = 0.0;
    ANG1 = 180.0;
    SPP = 0.0;// for FSPP=1,2
    setCuts(INDF, ICUT, ACOL, EMIN, S_PR, ANG0, ANG1, SPP);
}

void ZFitterWrapper::CutInfo() const
{ 
    const int mode = 1;
    zuinfo_(&mode); 
}


////////////////////////////////////////////////////////////////////////

complex ZFitterWrapper::rhoZ_f(const int INDF) const
{
    if (INDF<0 || INDF>11)
        throw std::runtime_error("Error in ZFitterWrapper::rhoZ_f");
    calcCommonBlocks();
    return complex(getCommonAROTFZ(INDF), getCommonAIROFZ(INDF), false);
}

complex ZFitterWrapper::kappaZ_f(const int INDF) const
{
    if (INDF<0 || INDF>11)
        throw std::runtime_error("Error in ZFitterWrapper::kappaZ_f");
    calcCommonBlocks();    
    return complex(getCommonARKAFZ(INDF), getCommonAIKAFZ(INDF), false);
}

complex ZFitterWrapper::gZ_f(const int INDF) const
{
    if (INDF<0 || INDF>11)
        throw std::runtime_error("Error in ZFitterWrapper::gZ_f");
    calcCommonBlocks();
    return complex(getCommonARVEFZ(INDF), getCommonAIVEFZ(INDF), false);
}


////////////////////////////////////////////////////////////////////////

double ZFitterWrapper::Af(const int INDF) const
{
    if (INDF<0 || INDF>9)
        throw std::runtime_error("Error in ZFitterWrapper::Af");
    double Qf = getCommonALLCH(INDF);
    double Re_gVf_over_gAf = 1.0 - 4.0*fabs(Qf)*kappaZ_f(INDF).real()*sw2();
    return ( 2.0*Re_gVf_over_gAf/(1.0 + Re_gVf_over_gAf*Re_gVf_over_gAf) );
}


////////////////////////////////////////////////////////////////////////
    
double ZFitterWrapper::alphaMZ() const
{
    calcCommonBlocks();    
    return ( calqed_.ALQEDZ );
}

double ZFitterWrapper::Mw() const
{
    calcCommonBlocks();    
    return ( MzCache*sqrt(1.0 - getCommonSIN2TW()) );
}

double ZFitterWrapper::Gamma_W() const
{
    calcCommonBlocks();    
    return ( getCommonPARTW(2) );
}

double ZFitterWrapper::sw2() const
{
    calcCommonBlocks();    
    return ( getCommonSIN2TW() );
}

double ZFitterWrapper::s2teff_f(const int INDF) const
{
    calcCommonBlocks();    
    return ( getCommonARSEFZ(INDF) );
}

double ZFitterWrapper::Gamma_f(const int INDF) const
{
    calcCommonBlocks();    
    return ( getCommonWIDTHS(INDF) );
}

double ZFitterWrapper::Gamma_inv() const
{
    calcCommonBlocks();    
    return ( 3.0*getCommonWIDTHS(0) );
}

double ZFitterWrapper::Gamma_had() const
{
    calcCommonBlocks();    
    return ( getCommonWIDTHS(10) );
}

double ZFitterWrapper::Gamma_Z() const
{
    calcCommonBlocks();    
    return ( getCommonWIDTHS(11) );
}

////////////////////////////////////////////////////////////////////////

void ZFitterWrapper::calcXS_AFB(const int INDF, const double SQRS,
                                double *XS, double *AFB) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zvthsm_(&INDF, &SQRS, &MzCache, &MtCache, &MhCache, &DAle5MzCache, &V_TB, 
            &AlsMzCache, XS, AFB);
}

void ZFitterWrapper::calcDXS(const int INDF, const double SQRS,
                             const double CSA, double *DXS) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zvatsm_(&INDF, &SQRS, &MzCache, &MtCache, &MhCache, &DAle5MzCache, &V_TB, 
            &AlsMzCache, &CSA, DXS);
}

void ZFitterWrapper::calcTauPol(const double SQRS, double *TAUPOL, double *TAUAFB) const
{
    calcCommonBlocks();    
    setSprimeCut(3, SQRS*SQRS);
    zvtpsm_(&SQRS, &MzCache, &MtCache, &MhCache, &DAle5MzCache, &V_TB, 
            &AlsMzCache, TAUPOL, TAUAFB);
}

void ZFitterWrapper::calcALR(const int INDF, const double SQRS, const double POL,
                             double *XSPL, double *XSMI) const {
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zvlrsm_(&INDF, &SQRS, &MzCache, &MtCache, &MhCache, &DAle5MzCache, &V_TB, 
            &AlsMzCache, &POL, XSPL, XSMI);
}

void ZFitterWrapper::calcAPV(double *C1U, double *C1D, double *C2U, double *C2D) const
{
    calcCommonBlocks();    
    double SIN2TW = getCommonSIN2TW();
    zu_apv_(&MzCache, &MtCache, &MhCache, &SIN2TW, &UMASS, &DMASS,
            C1U, C1D, C2U, C2D);
}

void ZFitterWrapper::calcXS(const int INDF, const double SQRS, const double GAMZ0,
                            const double GAMEE, const double GAMFF, double *XS) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zuxsec_(&INDF, &SQRS, &MzCache, &GAMZ0, &GAMEE, &GAMFF, XS);
}

void ZFitterWrapper::calcXS_AFB_2(const int INDF, const double SQRS, const double GAMZ0,
                                  const int MODE, const double GVE, const double XE,
                                  const double GVF, const double XF,
                                  double *XS, double *AFB) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zuxsa_(&INDF, &SQRS, &MzCache, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, XS, AFB);
}

void ZFitterWrapper::calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                                  const int MODE, const double GV2, const double X2,
                                  double *XS, double *AFB) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zuxsa2_(&INDF, &SQRS, &MzCache, &GAMZ0, &MODE, &GV2, &X2, XS, AFB);
}

void ZFitterWrapper::calcXS_AFB_4(const int INDF, const double SQRS, const double GAMZ0,
                                  const double PFOUR, const double PVAE2,
                                  const double PVAF2, double *XS, double *AFB) const
{
    calcCommonBlocks();    
    setSprimeCut(INDF, SQRS*SQRS);
    zuxafb_(&INDF, &SQRS, &MzCache, &GAMZ0, &PFOUR, &PVAE2, &PVAF2, XS, AFB);
}

void ZFitterWrapper::calcTauPol_2(const double SQRS, const double GAMZ0, const int MODE,
                                  const double GVE, const double XE, const double GVF,
                                  const double XF, double *TAUPOL, double *TAUAFB) const
{
    calcCommonBlocks();    
    setSprimeCut(3, SQRS*SQRS);
    zutau_(&SQRS, &MzCache, &GAMZ0, &MODE, &GVE, &XE, &GVF, &XF, TAUPOL, TAUAFB);
}


////////////////////////////////////////////////////////////////////////

void ZFitterWrapper::printConstants() const
{
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

void ZFitterWrapper::printInputs() const
{
    std::cout << "----- Input parameters -----" << std::endl;
    std::cout << "  ZMASS = " << SM.getMz() << std::endl;
    std::cout << "  TMASS = " << SM.getMtpole() << std::endl;
    std::cout << "  HMASS = " << SM.getMHl() << std::endl;
    std::cout << "  ALFAS = " << SM.getAlsMz() << std::endl;
    std::cout << "  DAL5H = " << SM.getDAle5Mz() << std::endl;
    std::cout << "  V_TB  = " << V_TB << std::endl;
    std::cout << "  UMASS = " << UMASS << std::endl;
    std::cout << "  DMASS = " << DMASS << std::endl;
    std::cout << std::endl;
}

void ZFitterWrapper::printIntermediateResults() const
{
    std::cout << "----- Intermediate Results -----" << std::endl;

    calcCommonBlocks();    
    
    std::cout << "Mw = " << Mw() << "  sin^2(theta_W) = "
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

std::string ZFitterWrapper::convertINDF(const int INDF) const
{
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

void ZFitterWrapper::test(const int IMISC) const
{ 
    zftest_(&IMISC);
}


////////////////////////////////////////////////////////////////////////

void ZFitterWrapper::init(const int IPRINT) const
{ 
    zuinit_(&IPRINT);
}

void ZFitterWrapper::calcCommonBlocks() const
{
    double MzTMP = SM.getMz();
    double MtTMP = SM.getMtpole();
    double MhTMP = SM.getMHl();
    double AlsMzTMP = SM.getAlsMz();
    double DAle5MzTMP = SM.getDAle5Mz();
    
    // TEST
    //std::cout << MzTMP << " " << MtTMP << " " << MhTMP << " " 
    //          << AlsMzTMP << " " << DAle5MzTMP << std::endl;
    //std::cout << MzCache << " " << MtCache << " " << MhCache << " " 
    //          << AlsMzCache << " " << DAle5MzCache << std::endl;
    
    if (MzTMP != MzCache || MtTMP != MtCache || MhTMP != MhCache || 
        AlsMzTMP != AlsMzCache || DAle5MzTMP != DAle5MzCache || 
        IsFlagChanged == true) { 
        MzCache = MzTMP;
        MtCache = MtTMP;
        MhCache = MhTMP;
        AlsMzCache = AlsMzTMP;
        DAle5MzCache = DAle5MzTMP;        
        //std::cout << "ZFitterWrapper::calcCommonBlocks() is working" << std::endl;
        zvweak_(&MzCache, &MtCache, &MhCache, &DAle5MzCache, &V_TB, &AlsMzCache); 
        IsFlagChanged = false;
    }
}


////////////////////////////////////////////////////////////////////////

double ZFitterWrapper::getCommonALPHST() const
{
    return zupars_.ALPHST;
}

double ZFitterWrapper::getCommonSIN2TW() const
{
    return zupars_.SIN2TW;
}

double ZFitterWrapper::getCommonS2TEFF(const int INDF) const
{
    if (INDF < 12) {
        return zupars_.S2TEFF[INDF];
    } else {
        std::cout << "S2TEFF[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonALLCH(const int INDF) const
{
    if (INDF < 12) {
        return zfchms_.ALLCH[INDF];
    } else {
        std::cout << "ALLCH[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonALLMS(const int INDF) const
{
    if (INDF < 12) {
        return zfchms_.ALLMS[INDF];
    } else {
        std::cout << "ALLMS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonWIDTHS(const int INDF) const
{
    if (INDF < 12) {
        return zupars_.WIDTHS[INDF]*0.001;
    } else {
        std::cout << "WIDTHS[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonPARTZ(const int INDF) const
{
    if (INDF < 12) {
        return partzw_.PARTZ[INDF]*0.001;
    } else {
        std::cout << "PARTZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonPARTW(const int i) const
{
    if (i < 3) {
        return partzw_.PARTW[i]*0.001;
    } else {
        std::cout << "PARTW[i < 3]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////

double ZFitterWrapper::getCommonARROFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.ARROFZ[INDF];
    } else {
        std::cout << "ARROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonARKAFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.ARKAFZ[INDF];
    } else {
        std::cout << "ARKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonARVEFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.ARVEFZ[INDF];
    } else {
        std::cout << "ARVEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonARSEFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.ARSEFZ[INDF];
    } else {
        std::cout << "ARSEFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonAROTFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.AROTFZ[INDF];
    } else {
        std::cout << "AROTFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonAIROFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.AIROFZ[INDF];
    } else {
        std::cout << "AIROFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonAIKAFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.AIKAFZ[INDF];
    } else {
        std::cout << "AIKAFZ[INDF < 11]" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double ZFitterWrapper::getCommonAIVEFZ(const int INDF) const
{
    if (INDF < 11) {
        return cdzrkz_.AIVEFZ[INDF];
    } else {
        std::cout << "AIVEFZ[INDF < 12]" << std::endl;
        exit(EXIT_FAILURE);
    }
}


////////////////////////////////////////////////////////////////////////

