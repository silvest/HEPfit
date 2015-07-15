/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

/*
 *  Documents of ZFITTER:
 *    hep-ph/9412201
 *    hep-ph/9908433 for v6.21
 *    hep-ph/0507146 for v6.42
 *
 *  Constants defined in the ZFITTER Fortran codes:
 *   - G_F^mu and alpha(0) are set in CONST1() in zfbib6_40.f and EWINIT()
 *     in zfbib6_40.f (and other places?). 
 *   - All lepton and quark masses (u,d,s,c,b) are set in CONST1().
 *   - Since it is hard to access some of these constants from the current
 *     class, the constants G_F^mu and alpha(0), which are necessary for
 *     the computations of the oblique parameters, are also set here in
 *     GF() and alpha().
 *
 *  Usage:
 *   1. call a constructor
 *        ZFitter()
 *   2. set flags (if necessary)
 *        setAllFlags() or setFlag()
 *   3. set cuts (if necessary)
 *   4. compute observables
 *
 *  Important flags in ZFITTER:
 *       AMT4 = 6 : with complete two-loop corrections to Mw and fermionic
 *                  two-loop corrections to sin^2(theta_eff^lept)
 *       ALEM = 2 : DAL5H is supplied by the user as input. 
 *       FOT2 = 3 : controls Initial-state radiation
 *       FINR = 0 : controls Final-state radiation 
 *       INTF = 0 : ISR/FSR interference is ignored. 
 *                  (INTF=1: with O(alpha), INTF=2: with the exponentiation)
 *       ISPP = 2 : with ISR pairs (ISPP=0: without the pairs)
 *       ICOL = -1: controls cuts
 *
 *  Fitted parameters in ZFITTER:
 *    ZMASS: the pole mass of Z in GeV
 *    TMASS: the pole mass of the top quark in GeV [10-400]
 *    HMASS: the pole mass of the Higgs in GeV [10-1000]
 *    ALFAS: the strong coupling constant at Mz
 *    DAL5H: the 5-quark flavor hadronic vacuum polarization at Mz
 *    UMASS: the constituent mass of the up quark (fixed)
 *    DMASS: the constituent mass of the down quark (fixed)
 *    V_TB: the CKM matrix elements V_{tb} (fixed)
 *
 *  INDF = Fermion index:
 *    INDF==0 neutrinos
 *    INDF==1 electron
 *    INDF==2 mu
 *    INDF==3 tau
 *    INDF==4 u
 *    INDF==5 d
 *    INDF==6 c
 *    INDF==7 s
 *    INDF==8 t
 *    INDF==9 b
 *    INDF==10 hadronic
 *    INDF==11 total
 * 
 */

#ifndef ZFITTERWRAPPER_H
#define	ZFITTERWRAPPER_H

#include <iostream>
#include <cstring>
#include <gslpp.h>


#include <StandardModel.h>

using namespace gslpp;


/**
 * @class ZFitterWrapper
 * @brief ZFITTER C++ Interface Class
 */
class ZFitterWrapper : public ThObsType {
public:

    /**
     * @brief ZFitterWrapper constructor
     * @param[in] mySM reference to a StandardModel object
     */
    ZFitterWrapper(const StandardModel& mySM);


    ////////////////////////////////////////////////////////////////////////
    // Flags
    
    /**
     * @brief modifies the default values of flags which control various options
     * @param[in] CHFLAG the string identifer of a ZFITTER flag
     * @param[in] IVALUE the value to assign to the flag
     *
     * @attention must be called before ZUWEAK
     */
    void setFlag(const std::string CHFLAG, const int IVALUE) const;

    /**
     * @brief sets flags (see Appendix B.2 in hep-ph/0507146)
     * @param[in] flags the array of the flags
     * @param[in] flagPrint prints flags if flagPrint=1
     */
    void setAllFlags(const int flags[46], const int flagPrint) const;    
    
    /**
     * @brief prints flag info with a ZFITTER subroutine
     */
    void FlagInfo() const;

    
    ////////////////////////////////////////////////////////////////////////    
    // Cuts
    
    /**
     * @bfief defines kinematic and geometric cuts for each fermion channel
     * @param[in] INDF fermion index
     * @param[in] ICUT flag for cuts, -1, 0, 1, 2 or 3
     * @param[in] ACOL acolinearity cut
     * @param[in] EMIN minimum fermion energy
     * @param[in] S_PR s'
     * @param[in] ANG0 minimum polar angle theta (deg)
     * @param[in] ANG1 maximum polar angle theta (deg)
     * @param[in] SIPP s cut for the secondary pairs, ISPP+FSPP
     */
    void setCuts(const int INDF, const int ICUT, const double ACOL,
                 const double EMIN, const double S_PR, const double ANG0,
                 const double ANG1, const double SIPP) const;

    /**
     * @brief sets cuts
     * @param[in] ICUT controls thekinds of cuts to be used
     * @param[in] ACOL acolinearity cut
     * @param[in] EMIN minimum fermion energy
     * @param[in] S_PR s'
     * @param[in] ANG0 minimum polar angle theta (deg)
     * @param[in] ANG1 maximum polar angle theta (deg)
     * @param[in] SPP s cut for for the secondary pairs, ISPP+FSPP
     * @param[in] flagPrint prints flags if flagPrint=1
     */
    void setAllCuts(const int ICUT[12], const double ACOL[12],
                    const double EMIN[12], const double S_PR[12],
                    const double ANG0[12], const double ANG1[12],
                    const double SPP[12], const int flagPrint) const;

    /**
     * @brief sets the s' cut
     * @param[in] INDF fermion index
     * @param[in] s invariant-mass-squared for e+ e- in the initial state
     */
    void setSprimeCut(const int INDF, const double s) const;
    
    /**
     * @brief prints cut info with a ZFITTER subroutine
     */
    void CutInfo() const;
    
    
    ////////////////////////////////////////////////////////////////////////
    // Effective couplings
    
    /**
     * @brief effective coupling rho_Z^f
     * @param[in] INDF fermion index
     * @return rho_Z^f
     */
    complex rhoZ_f(const int INDF) const;
    
    /**
     * @brief effective coupling kappa_Z^f
     * @param[in] INDF fermion index
     * @return kappa_Z^f
     */
    complex kappaZ_f(const int INDF) const;

    /**
     * @brief effective coupling g_Z^f
     * @param[in] INDF fermion index
     * @return g_Z^f
     */
    complex gZ_f(const int INDF) const;

    
    ////////////////////////////////////////////////////////////////////////
    // LR-asymmetry parameter
    
    /**
     * @param[in] INDF fermion index
     * @return LR-asymmetry parameter A_f
     */
    double Af(const int INDF) const;
    
    
    ////////////////////////////////////////////////////////////////////////
    // Z-pole observables
    
    /**
     * @return the electromagnetic coupling at Mz
     */
    double alphaMZ() const;
    
    /**
     * @return the W boson mass 
     */
    double Mw() const;
    
    /**
     * @return the total width of the W boson in GeV
     */
    double Gamma_W() const;

    /**
     * @return the weak mixing angle
     */
    double sw2() const;

    /**
     * @param[in] INDF fermion index
     * @return the effective weak mixing angle
     */
    double s2teff_f(const int INDF) const;
    
    /**
     * @param[in] INDF fermion index
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_f(const int INDF) const;

    /**
     * @return the invisible width of the Z boson in GeV
     */
    double Gamma_inv() const;

    /**
     * @return the hadronic width of the Z boson in GeV
     */
    double Gamma_had() const;

    /**
     * @return the total width of the Z boson in GeV
     */
    double Gamma_Z() const;
    
    
    ////////////////////////////////////////////////////////////////////////  
    
    /**
     * @brief calculates a cross section and a forward-backward asymmetry
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     */
    void calcXS_AFB(const int INDF, const double SQRS, double *XS, double *AFB) const;

    /**
     * @brief calculates a differential cross section w.r.t cos(theta)
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] CSA cosine of scattering angle
     * @param[out] DXS differential cross section
     */
    void calcDXS(const int INDF, const double SQRS, const double CSA, double *DXS) const;

    /**
     * @brief calculates the tau polarization and tau polarization asymmetry
     * @param[in] SQRS sqrt(s)
     * @param[out] TAUPOL tau polarization
     * @param[out] TAUAFB tau polarization forward-backward asymmetry
     */
    void calcTauPol(const double SQRS, double *TAUPOL, double *TAUAFB) const;

    /**
     * @brief calculates a left-right polarization asymmetry
     * @param[in] INDF fermion index 
     * @param[in] SQRS sqrt(s)
     * @param[in] POL degree of the longitudinal e^- polarization
     * @param[out] XSPL cross section for a given +POL
     * @param[out] XSMI cross section for a given -POL
     */
    void calcALR(const int INDF, const double SQRS, const double POL,
                 double *XSPL, double *XSMI) const;

    /**
     * @brief calculates C_1u, C_1d, C_2u and C_2d with an exp. input from atomic parity violation measurements in heavy atoms
     * @param[out] C1U coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C1D coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C2U coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C2D coupling in the electron-quark parity-violating Hamiltonian
     */
    void calcAPV(double *C1U, double *C1D, double *C2U, double *C2D) const;

    /** 
     * @brief calculates a cross section as a function of s, Mz, Gamma_Z, Gamma_e, Gamma_f
     * @param[in] INDF INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] GAMEE the e^+e^- partial width of Z
     * @param[in] GAMFF the f\bar{f} partial width of Z
     * @param[out] XS cross section in nb
     */
    void calcXS(const int INDF, const double SQRS, const double GAMZ0,
                const double GAMEE, const double GAMFF, double *XS) const;

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, Mz, Gamma_Z and the weak couplings
     * @param[in] INDF fermion index [1-9, 11]
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of XE/XF (electron/fermion), 0 for axial vector couplings, 1 for rho
     * @param[in] GVE effective vector coupling (electron)
     * @param[in] XE axial vector coupling or rho (electron)
     * @param[in] GVF effective vector coupling (fermion)
     * @param[in] XF axial vector coupling or rho (fermion)
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     */
    void calcXS_AFB_2(const int INDF, const double SQRS, const double GAMZ0,
                      const int MODE, const double GVE, const double XE,
                      const double GVF, const double XF, double *XS, double *AFB) const;

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, Mz, Gamma_Z and the weak couplings, assuming lepton universality
     * @param[in] INDF fermion index [1-3, 11]
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of X2, 0 for axial vector coupling squared, 1 for rho squared
     * @param[in] GV2 effective vector coupling squared
     * @param[in] X2 eff. axial vector coupling or rho, squared
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     *
     * @attention The couplings are squared.
     */
    void calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                      const int MODE, const double GV2, const double X2,
                      double *XS, double *AFB) const;

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, Mz, Gamma_Z and the weak couplings
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] PFOUR VE*AE*VF*AF
     * @param[in] PVAE2 VE**2+AE**2
     * @param[in] PVAF2 VF**2+AF**2
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     */
    void calcXS_AFB_4(const int INDF, const double SQRS, const double GAMZ0,
                      const double PFOUR, const double PVAE2, const double PVAF2,
                      double *XS, double *AFB) const;

    /**
     * @brief calculates the tau polarization and tau polarization asymmetry as functions of s, Mz, Gamma_Z and the weak couplings
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of XE/XF (electron/fermion), 0 for axial vector couplings, 1 for rho
     * @param[in] GVE effective vector coupling (electron)
     * @param[in] XE axial vector coupling or rho (electron)
     * @param[in] GVF effective vector coupling (fermion)
     * @param[in] XF axial vector coupling or rho (fermion)
     * @param[out] TAUPOL tau polarization
     * @param[out] TAUAFB tau polarization forward-backward asymmetry
     */
    void calcTauPol_2(const double SQRS, const double GAMZ0, const int MODE,
                      const double GVE, const double XE, const double GVF,
                      const double XF, double *TAUPOL, double *TAUAFB) const;

    
    ////////////////////////////////////////////////////////////////////////       

    /**
     * @brief prints constants defined in ZFITTER
     */
    void printConstants() const;

    /**
     * @brief prints input parameters
     */
    void printInputs() const;

    /**
     * @brief prints intermediate results
     */
    void printIntermediateResults() const;

    /**
     * @param[in] INDF fermion index
     * @return string of the channel corresponding to INDF
     */
    std::string convertINDF(const int INDF) const;
    
    /**
     * @brief Test program calculates cross sections and asymmetries as functions of @f$\sqrt{s}@f$
     * @param[in] IMISC 0 for non-scaled @f$\rho_Z^f@f$'s and 1 for scaled @f$\rho_Z^f@f$'s
     */
    void test(const int IMISC) const;

    
    ////////////////////////////////////////////////////////////////////////       
    
private:
    const StandardModel& SM;
    
    mutable bool IsFlagChanged;

    // Caches
    mutable double MzCache, MtCache, MhCache, AlsMzCache, DAle5MzCache;

    // Constants (set in the constructor)
    double V_TB;
    double UMASS; // the constituent mass of the up quark
    double DMASS; // the constituent mass of the down quark 

    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief initializes flag and cut variables with their default values. called by the constructors
     * @param[in] IPRINT 0 for no output
     *
     * @attention This method must be called before any other ZFITTER routine.
     */
    void init(const int IPRINT) const;

    /**
     * @brief calculates a number of important EW parameters which are stored in common blocks for later use
     */
    void calcCommonBlocks() const;

        
    ////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of ZFITTER

    /**
     * @return alpha_s(Mz)
     */
    double getCommonALPHST() const;

    /**
     * @return sin^2theta_w = 1 - Mw^2/Mz^2
     */
    double getCommonSIN2TW() const;

    /**
     * @param[in] INDF fermion index
     * @return an auxiliary quantity, which do not coincide with the effective weak mixing angle
     */
    double getCommonS2TEFF(const int INDF) const;

    /**
     * @param INDF fermion index
     * @return fermion charge
     */
    double getCommonALLCH(const int INDF) const;

    /**
     * @param INDF fermion index
     * @return fermion mass in GeV
     */
    double getCommonALLMS(const int INDF) const;
    
    /**
     * @param[in] INDF fermion index
     * @return partial and total decay widths of Z in GeV
     */
    double getCommonWIDTHS(const int INDF) const;
    
    /**
     * @param[in] INDF fermion index
     * @return partial or total Z widths in GeV
     */
    double getCommonPARTZ(const int INDF) const;
   
    /**
     * @param[in] i 0:leptons, 1:quarks, 3:total
     * @return partial or total W widths in GeV
     */
    double getCommonPARTW(const int i) const;   
    

    ////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of DIZET

    /**
     * @param INDF[in] INDF fermion index
     * @return (rho_Z^f)'
     */
    double getCommonARROFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[kappa_Z^f]
     */
    double getCommonARKAFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[g_Z^f]
     */
    double getCommonARVEFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return effective weak mixing angles \sin^2\theta_{eff}^f
     */
    double getCommonARSEFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[rho_Z^f]
     */
    double getCommonAROTFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[rho_Z^f]
     */
    double getCommonAIROFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[kappa_Z^f]
     */
    double getCommonAIKAFZ(const int INDF) const;

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[g_Z^f]
     */
    double getCommonAIVEFZ(const int INDF) const;    
    
    
};
    
#endif	/* ZFITTERWRAPPER_H */
    
