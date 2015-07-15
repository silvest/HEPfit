/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ZFITTER_F2C_H
#define	ZFITTER_F2C_H


/* complex double structure defined for ewform_ */
//struct dcomplex { double real, imag; };


/* ZFITTER common blocks (global variables in ZFITTER)
 *   see Appendix A.1 in hep-ph/9412201 for detail */
extern "C" {

    /* parameters */
    extern struct {
        double QDF;
        double QCDCOR[15];
        double ALPHST;
        double SIN2TW;
        double S2TEFF[12]; // do not correspond to the usual sin^2(th_eff^f)
        double WIDTHS[12]; // in MeV
    } zupars_;

    /* electroweak form factors */
    //extern struct {
    //    dcomplex XALLCH[4][5]; // form factors
    //    dcomplex XFOTF;        // 1+\Delta\alpha(s)
    //} ewform_;

    /* the charges and masses of fermions 
     * note: ALLCH[10,11] and ALLMS[10,11] are not used */
    extern struct {
        double ALLCH[12]; // the fermion charges
        double ALLMS[12]; // the fermion masses in GeV
    } zfchms_;

    /* widths */
    extern struct {
        double PARTZ[12]; // in MeV
        double PARTW[3];  // in MeV
    } partzw_;
            
}


/* DIZET common blocks (global variables used in DIZET)
 *   see Sec. 2.5.1 in hep-ph/9908433 and Sec. 6.1 in hep-ph/0507146 for detail */
extern "C" {

    /* effective couplings of Z decays
     * note: The LEP EWWG used the real parameter \rho_f=AROTFZ/. */
    extern struct {
        double ARROFZ[11]; // (\rho_Z^f)'
        double ARKAFZ[11]; // Re[\kappa_Z^f]
        double ARVEFZ[11]; // Re[g_Z^f]
        double ARSEFZ[11]; // effective weak mixing angles \sin^2\theta_{eff}^f
        double AROTFZ[11]; // Re[\rho_Z^f]
        double AIROFZ[11]; // Im[\rho_Z^f]
        double AIKAFZ[11]; // Im[\kappa_Z^f]
        double AIVEFZ[11]; // Im[g_Z^f]
    } cdzrkz_;

    //extern struct {
    //    double PARTZA[11];
    //    double PARTZI[11];
    //    double RENFAC[11]; // R_f
    //    double SRENFC[11]; // \sqrt{R_f}
    //} cdzaux_;

    extern struct {
         double ALQEDZ; // alpha(Mz)
         double ALQEDS;
    } calqed_;

}


/* Subroutines in zfusr6_43.f */
extern "C" {

    void zftest_(const int *IMISC);
    void zuinit_(const int *IPRINT);
    void zuinfo_(const int *MODE);
    void zuflag_(const char *CHFLAG, const int *IVALUE, long int CHFLAG_SIZE);

    void zuweak_(const double *ZMASS, const double *TMASS, const double *HMASS,
                 const double *DAL5H, const double *ALFAS);
    void zvweak_(const double *ZMASS, const double *TMASS, const double *HMASS,
                 const double *DAL5H, const double *V_TB, const double *ALFAS);

    void zucuts_(const int *INDF, const int *ICUT, const double *ACOL,
                 const double *EMIN, const double *S_PR, const double *ANG0,
                 const double *ANG1, const double *SIPP);

    void zuthsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *ALFAS, double *XS, double *AFB);
    void zvthsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *V_TB, const double *ALFAS, double *XS,
                 double *AFB);

    void zuatsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *ALFAS, const double *CSA, double *DXS);
    void zvatsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *V_TB, const double *ALFAS, const double *CSA,
                 double *DXS);

    void zutpsm_(const double *SQRS, const double *ZMASS, const double *TMASS,
                 const double *HMASS, const double *DAL5H, const double *ALFAS,
                 double *TAUPOL, double *TAUAFB);
    void zvtpsm_(const double *SQRS, const double *ZMASS, const double *TMASS,
                 const double *HMASS, const double *DAL5H, const double *V_TB,
                 const double *ALFAS, double *TAUPOL, double *TAUAFB);

    void zulrsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *ALFAS, const double *POL, double *XSPL,
                 double *XSMI);
    void zvlrsm_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *TMASS, const double *HMASS, const double *DAL5H,
                 const double *V_TB, const double *ALFAS, const double *POL,
                 double *XSPL, double *XSMI);

    void zu_apv_(const double *ZMASS,const double *TMASS, const double *HMASS,
                 const double *SIN2TW, const double *UMASS, const double *DMASS,
                 double *C1U, double *C1D, double *C2U, double *C2D);

    /*
    // CPU saving interface to DIZET, called by ZUTHSM, ZUATSM and ZUTPSM
    void zdizet_(const double *ZMASS, const double *TMASS, const double *HMASS,
                 const double *DAL5H, const double *V_TB, const double *ALFAS,
                 double *ALQED, double *ALFAT);

    // CPU saving interface to QCDCOF, called by ZUATSM
    void zqcdco_(const double *SQRS, const double *TMASS, const double *SINEFF,
                 const double *ALQED, const double *ALFAST, const double *ALFATT,
                 const double *ALPHXI, const double *QCDCCR);
     */

    void zuxsec_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *GAMZ0, const double *GAMEE, const double *GAMFF,
                 double *XS);

    void zuxsa_(const int *INDF, const double *SQRS, const double *ZMASS,
                const double *GAMZ0, const int *MODE, const double *GVE,
                const double *XE, const double *GVF, const double *XF,
                double *XS, double *AFB);
    void zuxsa2_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *GAMZ0, const int *MODE, const double *GV2,
                 const double *X2, double *XS, double *AFB);
    void zuxafb_(const int *INDF, const double *SQRS, const double *ZMASS,
                 const double *GAMZ0, const double *PFOUR, const double *PVAE2,
                 const double *PVAF2, double *XS, double *AFB);

    void zutau_(const double *SQRS, const double *ZMASS, const double *GAMZ0,
                const int *MODE, const double *GVE, const double *XE,
                const double *GVF, const double *XF, double *TAUPOL,
                double *TAUAFB);
    void zualr_(const double *SQRS, const double *ZMASS, const double *GAMZ0,
                const int *MODE, const double *GVE, const double *XE,
                const double *GVF, const double *XF, double *TAUPOL,
                double *TAUAFB);

}


#endif	/* ZFITTER_F2C_H */

