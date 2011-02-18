/* 
 * File:   ZFitter.h
 * Author: mishima
 *
 * Created on February 15, 2011, 11:15 PM
 */

#ifndef ZFITTER_H
#define	ZFITTER_H

#include <iostream>
#include <string>


/*  Documents of ZFITTER:
 *     hep-ph/9412201
 *     hep-ph/9908433 for v6.21
 *     hep-ph/0507146 for v6.42  */


/* complex double structure defined for ewform_ */
struct dcomplex { double real, imag; };


/* ZFITTER common blocks (glocal variables in ZFITTER)
 *   see Appendix A.1 in hep-ph/9412201 for detail */
extern "C" {

    /* parameters */
    extern struct {
        double QDF;
        double QCDCOR[15];
        double ALPHST;
        double SIN2TW;
        double S2TEFF[12];
        double WIDTHS[12];
    } zupars_;

    /* electroweak form factors */
    //extern struct {
    //    dcomplex XALLCH[4][5]; // form factors
    //    dcomplex XFOTF;        // 1+\Delta\alpha(s)
    //} ewform_;

    /* the charges and masses of fermions 
     * note: ALLCH[10,11] and ALLMS[10,11] are not used */
    //extern struct {
    //    double ALLCH[12]; // the fermion charges
    //    double ALLMS[12]; // the fermion masses
    //} zfchms_;

    //extern struct {
    //    double PARTZ[12]; // \Gamma_f : partial Z-decay widths
    //    double PARTW[3];
    //} partzw_;
            
}


/* DIZET common blocks (glocal variables used in DIZET)
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


/**
 * @class ZFitter
 * @brief ZFITTER C++ Interface Class
 */
class ZFitter {
public:

    /**
     * @brief ZFitter constructor
     * @param[in] ZMASS_i the pole mass of Z in GeV
     * @param[in] TMASS_i the pole mass of the top quark in GeV [10-400]
     * @param[in] HMASS_i the pole mass of the Higgs in GeV [10-1000]
     * @param[in] ALFAS_i the strong coupling constant at mZ
     * @param[in] DAL5H_i the 5-quark flavour hadronic vacuum polarization at mZ
     * @param[in] ALFA_i the fine-structure constant
     * @param[in] V_TB_i the CKM element V_tb
     * @param[in] UMASS_i the constituent mass of the up quark (for calcAPV())
     * @param[in] DMASS_i the constituent mass of the down quark (for calcAPV())
     */
    ZFitter(double ZMASS_i, double TMASS_i, double HMASS_i, double ALFAS_i,
            double DAL5H_i, double ALFA_i, double V_TB_i, double UMASS_i,
            double DMASS_i);

    /**
     * @brief ZFitter copy constructor
     * @param[in] orig reference to a ZFitter object
     */
    ZFitter(const ZFitter& orig);

    /**
     * @brief ZFitter destructor
     */
    virtual ~ZFitter();


    ///////////////////////////////////////////////////////////////////////////

    
    /**
     * @return the pole mass of Z in GeV
     */
    double getZMASS() const { return ZMASS; }

    /**
     * @brief set the pole mass of Z in GeV
     * @param[in] ZMASS the Z pole mass of Z in GeV
     */
    void setZMASS(double ZMASS) { this->ZMASS = ZMASS; }

    /**
     * @return the pole mass of the top quark in GeV [10-400]
     */
    double getTMASS() const { return TMASS; }

    /**
     * @brief set the pole mass of the top quark in GeV [10-400]
     * @param[in] TMASS the pole mass of the top quark in GeV [10-400]
     */
    void setTMASS(double TMASS) { this->TMASS = TMASS; }

    /**
     * @return the pole mass of the Higgs in GeV [10-1000]
     */
    double getHMASS() const { return HMASS; }

    /**
     * @brief set the pole mass of the Higgs in GeV [10-1000]
     * @param[in] HMASS the pole mass of the Higgs in GeV [10-1000]
     */
    void setHMASS(double HMASS) { this->HMASS = HMASS; }

    /**
     * @return the strong coupling constant at mZ
     */
    double getALFAS() const { return ALFAS; }

    /**
     * @brief set the strong coupling constant at mZ
     * @param[in] ALFAS the strong coupling constant at mZ
     */
    void setALFAS(double ALFAS) { this->ALFAS = ALFAS; }

    /**
     * @return the 5-quark flavour hadronic vacuum polarization at mZ
     */
    double getDAL5H() const { return DAL5H; }

    /**
     * @brief set the 5-quark flavour hadronic vacuum polarization at mZ
     * @param[in] DAL5H the 5-quark flavour hadronic vacuum polarization at mZ
     */
    void setDAL5H(double DAL5H) { this->DAL5H = DAL5H; }

    /**
     * @return the fine-structure constant
     */
    double getALFA() const { return ALFA; }

    /**
     * @brief set the fine-structure constant
     * @param[in] ALFA the fine-structure constant
     */
    void setALFA(double ALFA) { this->ALFA = ALFA; }

    /**
     * @return the CKM element V_tb
     */
    double getV_TB() const { return V_TB; }

    /**
     * @brief set the CKM element V_tb
     * @param[in] V_TB the CKM element V_tb
     */
    void setV_TB(double V_TB) { this->V_TB = V_TB; }

    /**
     * @return the constituent mass of the up quark
     */
    double getUMASS() const { return UMASS; }

    /**
     * @brief set the constituent mass of the up quark
     * @param[in] UMASS the constituent mass of the up quark
     */
    void setUMASS(double UMASS) { this->UMASS = UMASS; }

    /**
     * @return the constituent mass of the down quark
     */
    double getDMASS() const { return DMASS; }

    /**
     * @brief set the constituent mass of the down quark
     * @param[in] DMASS the constituent mass of the down quark
     */
    void setDMASS(double DMASS) { this->DMASS = DMASS; }


    ///////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of ZFITTER


    /**
     * @return alpha_s(mZ)
     */
    double getCommonAlphas();

    /**
     * @return sin^2theta_w = 1 - mW^2/mZ^2
     */
    double getCommonSin2thetaW();

    /**
     * @param[in] INDF fermion index
     * @return effective weak mixing angles @f$\sin^2\theta_{\mathrm{eff}}^f@f$
     */
    double getCommonS2TEFF(int INDF);

    /**
     * @param[in] INDF fermion index
     * @return partial and total decay widths of Z in MeV
     */
    double getCommonWidths(int INDF);


    ///////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of DIZET


    /**
     * @param INDF[in] INDF fermion index
     * @return (\rho_Z^f)'
     */
    double getCommonARROFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[\kappa_Z^f]
     */
    double getCommonARKAFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[g_Z^f]
     */
    double getCommonARVEFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return effective weak mixing angles \sin^2\theta_{eff}^f
     */
    double getCommonARSEFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[\rho_Z^f]
     */
    double getCommonAROTFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[\rho_Z^f]
     */
    double getCommonAIROFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[\kappa_Z^f]
     */
    double getCommonAIKAFZ(int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[g_Z^f]
     */
    double getCommonAIVEFZ(int INDF);


    ///////////////////////////////////////////////////////////////////////////

    
    /**
     * @brief Test program calculates cross sections and asymmetries as functions of @f$\sqrt{s}@f$
     * @param[in] IMISC 0 for non-scaled @f$\rho_Z^f@f$'s and 1 for scaled @f$\rho_Z^f@f$'s
     */
    void test(const int IMISC);

    /**
     * @brief initializes flag and cut variables with their default values. called by the constructors
     * @attention must be called before any other ZFITTER routine.
     * @param[in] IPRINT 0 for no output
     */
    void init(const int IPRINT);

    /**
     * @brief prints flag or cut info
     * @param[in] MODE 0(1) for printing flag(cut) info
     */
    void info(const int MODE);

    /**
     * @brief modifies the default values of flags which control various options
     * @attention must be called before ZUWEAK
     * @param[in] CHFLAG the string identifer of a ZFITTER flag
     * @param[in] IVALUE the value to assign to the flag
     */
    void flag(const std::string CHFLAG, const int IVALUE);

    /**
     * @brief calculates a number of important EW parameters which are stored in common blocks for later use
     */
    void calcCommonBlocks();

    /**
     * @bfief defines kinematic and geometric cuts for each fermion channel
     * @param[in] INDF fermion index
     * @param[in] ICUT flag for cuts, -1, 0, 1, 2 or 3
     * @param[in] ACOL acolinearity cut
     * @param[in] EMIN minimum fermion energy
     * @param[in] S_PR sprime, related to maximum photon energy
     * @param[in] ANG0 minimum polar angle theta (deg)
     * @param[in] ANG1 maximum polar angle theta (deg)
     * @param[in] SIPP s cut for the secondary pairs, ISPP+FSPP
     */
    void setCuts(const int INDF, const int ICUT, const double ACOL,
                 const double EMIN, const double S_PR, const double ANG0,
                 const double ANG1, const double SIPP);

    /**
     * @brief calculates a cross section and a forward-backward asymmetry
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s) 
     */
    void calcXS_AFB(const int INDF, const double SQRS);

    /**
     * @brief calculates a differential cross section w.r.t cos(theta)
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] CSA cosine of scattering angle
     */
    void calcDXS(const int INDF, const double SQRS, const double CSA);

    /**
     * @brief calculates the tau polarization and tau polarization asymmetry
     * @param[in] SQRS sqrt(s)
     */
    void calcTauPol(const double SQRS);

    /**
     * @brief calculates a left-right polarization asymmetry
     * @param[in] INDF fermion index 
     * @param[in] SQRS sqrt(s)
     * @param[in] POL degree of the longitudinal e^- polarization
     */
    void calcALR(const int INDF, const double SQRS, const double POL);


    /**
     * @brief calculates C_1u, C_1d, C_2u and C_2d with an exp. input from atomic parity violation measurements in heavy atoms
     */
    void calcAPV();

    /** 
     * @brief calculates a cross section as a function of s, mZ, Gamma_Z, Gamma_e, Gamma_f
     * @param[in] INDF INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] GAMEE the e^+e^- partial width of Z
     * @param[in] GAMFF the f\bar{f} partial width of Z
     */
    void calcXS(const int INDF, const double SQRS, const double GAMZ0,
                const double GAMEE, const double GAMFF);

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, mZ, Gamma_Z and the weak couplings
     * @param[in] INDF fermion index [1-9, 11]
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of XE/XF (electron/fermion), 0 for axial vector couplings, 1 for rho
     * @param[in] GVE effective vector coupling (electron)
     * @param[in] XE axial vector coupling or rho (electron)
     * @param[in] GVF effective vector coupling (fermion)
     * @param[in] XF axial vector coupling or rho (fermion)
     */
    void calcXS_AFB_2(const int INDF, const double SQRS, const double GAMZ0,
                      const int MODE, const double GVE, const double XE,
                      const double GVF, const double XF);

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, mZ, Gamma_Z and the weak couplings, assuming lepton universality
     * @attention The couplings are squared.
     * @param[in] INDF fermion index [1-3, 11]
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of X2, 0 for axial vector coupling squared, 1 for rho squared
     * @param[in] GV2 effective vector coupling squared
     * @param[in] X2 eff. axial vector coupling or rho, squared
     */
    void calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                      const int MODE, const double GV2, const double X2);

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, mZ, Gamma_Z and the weak couplings
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] PFOUR VE*AE*VF*AF
     * @param[in] PVAE2 VE**2+AE**2
     * @param[in] PVAF2 VF**2+AF**2
     */
    void calcXS_AFB_4(const int INDF, const double SQRS, const double GAMZ0,
                      const double PFOUR, const double PVAE2, const double PVAF2);

    /**
     * @brief calculates the tau polarization and tau polarization asymmetry as functions of s, mZ, Gamma_Z and the weak couplings
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of XE/XF (electron/fermion), 0 for axial vector couplings, 1 for rho
     * @param[in] GVE effective vector coupling (electron)
     * @param[in] XE axial vector coupling or rho (electron)
     * @param[in] GVF effective vector coupling (fermion)
     * @param[in] XF axial vector coupling or rho (fermion)
     */
    void calcTauPol_2(const double SQRS, const double GAMZ0, const int MODE,
                      const double GVE, const double XE, const double GVF,
                      const double XF);


    ///////////////////////////////////////////////////////////////////////////

    /* Outputs */
    double XS[12]; // cross section in nb
    double DXS[12]; // differential cross section
    double AFB[12]; // forward-backward asymmetry
    double TAUPOL; // tau polarization
    double TAUAFB; // tau polarization forward-backward asymmetry
    double XSPL[12]; // cross section for a given +POL
    double XSMI[12]; // cross section for a given -POL
    double C1U, C1D, C2U, C2D; // couplings in the electron-quark parity-violating Hamiltonian


    ///////////////////////////////////////////////////////////////////////////

private:

    /* Inputs */
    double ZMASS; // the pole mass of Z in GeV
    double TMASS; // the pole mass of the top quark in GeV [10-400]
    double HMASS; // the pole mass of the Higgs in GeV [10-1000]
    double ALFAS; // the strong coupling constant at mZ
    double DAL5H; // the 5-quark flavour hadronic vacuum polarization at mZ
    double ALFA;  // the fine-structure constant
    double V_TB;  // the CKM element V_tb
    double UMASS; // the constituent mass of the up quark
    double DMASS; // the constituent mass of the down quark

};

#endif	/* ZFITTER_H */

