/* 
 * File:   ZFitterObservables.h
 * Author: mishima
 */

#ifndef ZFITTEROBSERVABLES_H
#define	ZFITTEROBSERVABLES_H

#include "ZFitter.h"


class ZFitterObservables : public ZFitter {
public:

    ZFitterObservables(const ZFitter& ZFitter_i);

    //ZFitterObservables(const ZFitterObservables& orig);

    virtual ~ZFitterObservables();
    
    
    ////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of ZFITTER
 
    /**
     * @param[in] INDF fermion index
     * @return partial and total decay widths of Z in GeV
     */
    double getCommonWIDTHS(const int INDF);
    
    /**
     * @param[in] INDF fermion index
     * @return partial or total Z widths in GeV
     */
    double getCommonPARTZ(const int INDF);
   
    /**
     * @param[in] i 0:leptons, 1:quarks, 3:total
     * @return partial or total W widths in GeV
     */
    double getCommonPARTW(const int i);    
    
    
    ////////////////////////////////////////////////////////////////////////    
    
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
     * @brief sets cuts
     * @param[in] ICUT controls thekinds of cuts to be used
     * @param[in] ACOL acolinearity cut
     * @param[in] EMIN minimum fermion energy
     * @param[in] S_PR s' (variable related to maximum photon energy)
     * @param[in] ANG0 minimum polar angle theta (deg)
     * @param[in] ANG1 maximum polar angle theta (deg)
     * @param[in] SPP s cut for for the secondary pairs, ISPP+FSPP
     * @param[in] flagPrint prints flags if flagPrint=1
     */
    void setAllCuts(const int ICUT[12], const double ACOL[12],
                    const double EMIN[12], const double S_PR[12],
                    const double ANG0[12], const double ANG1[12],
                    const double SPP[12], const int flagPrint);
    
    /**
     * @brief prints cut info using a ZFITTER subroutine
     */
    void CutInfo();

    
    ////////////////////////////////////////////////////////////////////////    
    
    /**
     * @brief calculates a cross section and a forward-backward asymmetry
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     */
    void calcXS_AFB(const int INDF, const double SQRS, double *XS, double *AFB);

    /**
     * @brief calculates a differential cross section w.r.t cos(theta)
     * @param[in] INDF fermion index
     * @param[in] SQRS sqrt(s)
     * @param[in] CSA cosine of scattering angle
     * @param[out] DXS differential cross section
     */
    void calcDXS(const int INDF, const double SQRS, const double CSA, double *DXS);

    /**
     * @brief calculates the tau polarization and tau polarization asymmetry
     * @param[in] SQRS sqrt(s)
     * @param[out] TAUPOL tau polarization
     * @param[out] TAUAFB tau polarization forward-backward asymmetry
     */
    void calcTauPol(const double SQRS, double *TAUPOL, double *TAUAFB);

    /**
     * @brief calculates a left-right polarization asymmetry
     * @param[in] INDF fermion index 
     * @param[in] SQRS sqrt(s)
     * @param[in] POL degree of the longitudinal e^- polarization
     * @param[out] XSPL cross section for a given +POL
     * @param[out] XSMI cross section for a given -POL
     */
    void calcALR(const int INDF, const double SQRS, const double POL,
                 double *XSPL, double *XSMI);

    /**
     * @brief calculates C_1u, C_1d, C_2u and C_2d with an exp. input from atomic parity violation measurements in heavy atoms
     * @param[out] C1U coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C1D coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C2U coupling in the electron-quark parity-violating Hamiltonian
     * @param[out] C2D coupling in the electron-quark parity-violating Hamiltonian
     */
    void calcAPV(double *C1U, double *C1D, double *C2U, double *C2D);

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
                const double GAMEE, const double GAMFF, double *XS);

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
                      const double GVF, const double XF, double *XS, double *AFB);

    /**
     * @brief calculates a cross section and a forward-backward asymmetry as functions of s, Mz, Gamma_Z and the weak couplings, assuming lepton universality
     * @attention The couplings are squared.
     * @param[in] INDF fermion index [1-3, 11]
     * @param[in] SQRS sqrt(s)
     * @param[in] GAMZ0 the total width of Z
     * @param[in] MODE meaning of X2, 0 for axial vector coupling squared, 1 for rho squared
     * @param[in] GV2 effective vector coupling squared
     * @param[in] X2 eff. axial vector coupling or rho, squared
     * @param[out] XS cross section in nb
     * @param[out] AFB forward-backward asymmetry
     */
    void calcXS_AFB_3(const int INDF, const double SQRS, const double GAMZ0,
                      const int MODE, const double GV2, const double X2,
                      double *XS, double *AFB);

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
                      double *XS, double *AFB);

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
                      const double XF, double *TAUPOL, double *TAUAFB);
    
    /**
     * @brief Test program calculates cross sections and asymmetries as functions of @f$\sqrt{s}@f$
     * @param[in] IMISC 0 for non-scaled @f$\rho_Z^f@f$'s and 1 for scaled @f$\rho_Z^f@f$'s
     */
    void test(const int IMISC);

    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the total width of the W boson in GeV
     */
    double Gamma_W();

    /**
     * @return the weak mixing angle
     */
    double sw2();

    /**
     * @param[in] INDF fermion index
     * @return the effective weak mixing angle
     */
    double s2teff_f(const int INDF);
    
    /**
     * @param[in] INDF fermion index
     * @return the partial decay width of the Z boson in GeV
     */
    double Gamma_f(const int INDF);

    /**
     * @return the invisible width of the Z boson in GeV
     */
    double Gamma_inv();

    /**
     * @return the hadronic width of the Z boson in GeV
     */
    double Gamma_had();

    /**
     * @return the total width of the Z boson in GeV
     */
    double Gamma_Z();

 
    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] INDF fermion index
     * @return string of the channel corresponding to INDF
     */
    std::string convertINDF(const int INDF);

    /**
     * @brief prints constants defined in ZFITTER
     */
    void printConstants();

    /**
     * @brief prints input parameters
     */
    void printInputs();

    /**
     * @brief prints intermediate results
     */
    void printIntermediateResults();

    
    ////////////////////////////////////////////////////////////////////////    
    
private:

};

#endif	/* ZFITTEROBSERVABLES_H */

