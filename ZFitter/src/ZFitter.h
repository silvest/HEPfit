/* 
 * File:   ZFitter.h
 * Author: mishima
 */

/*
 *  Documents of ZFITTER:
 *    hep-ph/9412201
 *    hep-ph/9908433 for v6.21
 *    hep-ph/0507146 for v6.42
 *
 *  Oblique parameters:
 *    NPB369, 3 (1992)
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
 *  How to use (see ZFitterTest-2.cpp for detail):
 *   1. call a constructor
 *        ZFitter(const StandardModel& mySM)
 *   2. set flags, if necessary to change from their default values
 *        setAllFlags(const int flags[46], const int flagPrint)
 *          or
 *        flag(const std::string CHFLAG, const int IVALUE)
 *   3. compute the effective couplings, the partial widths, the effective
 *      weak mixing angle, etc.
 *        calcCommonBlocks()
 *   4. calculale and get EW precision observables
 *        ---> use ZFitterObservables class or another method
 *
 *  ZFITTER flags:
 *   - flags which should be changed from their defaults
 *       AMT4 = 6 : with complete two-loop corrections to Mw and fermionic
 *                  two-loop corrections to sin^2(theta_eff^lept)
 *       ALEM = 2 : DAL5H is supplied by the user as input. 
 *   - flags used with their default value:
 *       AFBC, SCAL, SCRE, BORN, BOXD, CONV, FINR, FOT2, GAMS, DIAG, INTF, 
 *       BARB, PART, POWR, PRNT, QCDC, VPOL, WEAK, FTJR, EXPR, EXPF, HIGS,
 *       AFMT, CZAK, PREC, HIG2, ALE2, GFER, ISPP, FSRS, MISC, MISD, IPFC,
 *       IPSC, IPTO, FBHO, FSPP, FUNA, ASCR, SFSR, ENUE, TUPV, DMWW, DSWW
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

#ifndef ZFITTER_H
#define	ZFITTER_H

#include <iostream>
#include <cstring>
#include <gslpp.h>
#include <StandardModel.h>
#include "ZFitter-f2c.h"

using namespace gslpp;


/**
 * @class ZFitter
 * @brief ZFITTER C++ Interface Class
 */
class ZFitter {
public:

    /**
     * @brief ZFitter constructor
     * @param[in] mySM reference to a StandardModel object
     */
    ZFitter(const StandardModel& mySM);
    
    /**
     * @brief ZFitter copy constructor
     * @param[in] orig reference to a ZFitter object
     */
    //ZFitter(const ZFitter& orig);

    /**
     * @brief ZFitter destructor
     */
    virtual ~ZFitter();


    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @return the pole mass of Z in GeV
     */
    double getZMASS() const { return ZMASS; }

    /**
     * @return the pole mass of the top quark in GeV [10-400]
     */
    double getTMASS() const { return TMASS; }

    /**
     * @return the pole mass of the Higgs in GeV [10-1000]
     */
    double getHMASS() const { return HMASS; }

    /**
     * @return the strong coupling constant at Mz
     */
    double getALFAS() const { return ALFAS; }

    /**
     * @return the 5-quark flavour hadronic vacuum polarization at Mz
     */
    double getDAL5H() const { return DAL5H; }

    /**
     * @return the CKM element V_tb
     */
    double getV_TB() const { return V_TB; }

    /**
     * @return the constituent mass of the up quark
     */
    double getUMASS() const { return UMASS; }

    /**
     * @return the constituent mass of the down quark
     */
    double getDMASS() const { return DMASS; }

    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double getAlphaMZ() const {
        return alphaMZ;
    }

    /**
     * @brief the W bosn mass
     * @return M_W
     */
    double getMw() const {
        return Mw;
    }    
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l name of a lepton 
     * @return rho_Z^l for lepton "l" 
     */
    complex getRhoZ_l(const StandardModel::lepton l) const {
        return rhoZ_l[l];
    }

    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q name of a quark
     * @return rho_Z^q for quark "q" 
     */
    complex getRhoZ_q(const StandardModel::quark q) const {
        return rhoZ_q[q];
    }
    
    /**
     * @brief effective coupling kappa_Z^l
     * @param[in] l name of a lepton 
     * @return kappa_Z^l for lepton "l" 
     */
    complex getKappaZ_l(const StandardModel::lepton l) const {
        return kappaZ_l[l];
    }
    
    /**
     * @brief effective coupling kappa_Z^q
     * @param[in] q name of a quark
     * @return kappa_Z^q for quark "q" 
     */
    complex getKappaZ_q(const StandardModel::quark q) const {
        return kappaZ_q[q];
    }


    ////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of ZFITTER

    /**
     * @return alpha_s(Mz)
     */
    double getCommonALPHST();

    /**
     * @return sin^2theta_w = 1 - Mw^2/Mz^2
     */
    double getCommonSIN2TW();

    /**
     * @param[in] INDF fermion index
     * @return an auxiliary quantity, which do not coincide with the effective weak mixing angle
     */
    double getCommonS2TEFF(const int INDF);

    /**
     * @param INDF fermion index
     * @return fermion charge
     */
    double getCommonALLCH(const int INDF);

    /**
     * @param INDF fermion index
     * @return fermion mass in GeV
     */
    double getCommonALLMS(const int INDF);


    ////////////////////////////////////////////////////////////////////////
    // get variables from Common blocks of DIZET

    /**
     * @param INDF[in] INDF fermion index
     * @return (rho_Z^f)'
     */
    double getCommonARROFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[kappa_Z^f]
     */
    double getCommonARKAFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[g_Z^f]
     */
    double getCommonARVEFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return effective weak mixing angles \sin^2\theta_{eff}^f
     */
    double getCommonARSEFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Re[rho_Z^f]
     */
    double getCommonAROTFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[rho_Z^f]
     */
    double getCommonAIROFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[kappa_Z^f]
     */
    double getCommonAIKAFZ(const int INDF);

    /**
     * @param INDF[in] INDF fermion index
     * @return Im[g_Z^f]
     */
    double getCommonAIVEFZ(const int INDF);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief initializes flag and cut variables with their default values. called by the constructors
     * @attention must be called before any other ZFITTER routine.
     * @param[in] IPRINT 0 for no output
     */
    void init(const int IPRINT);

    /**
     * @brief modifies the default values of flags which control various options
     * @attention must be called before ZUWEAK
     * @param[in] CHFLAG the string identifer of a ZFITTER flag
     * @param[in] IVALUE the value to assign to the flag
     */
    void flag(const std::string CHFLAG, const int IVALUE);

    /**
     * @brief sets flags (see Appendix B.2 in hep-ph/0507146)
     * @param[in] flags the array of the flags
     * @param[in] flagPrint prints flags if flagPrint=1
     */
    void setAllFlags(const int flags[46], const int flagPrint);    
    
    /**
     * @brief prints flag info using a ZFITTER subroutine
     */
    void FlagInfo();
    
    /**
     * @brief calculates a number of important EW parameters which are stored in common blocks for later use
     */
    void calcCommonBlocks();

    
    ////////////////////////////////////////////////////////////////////////

protected:
    const StandardModel& SM;
    
    /* Inputs */
    double ZMASS; // the pole mass of Z in GeV
    double TMASS; // the pole mass of the top quark in GeV [10-400]
    double HMASS; // the pole mass of the Higgs in GeV [10-1000]
    double ALFAS; // the strong coupling constant at Mz
    double DAL5H; // the 5-quark flavour hadronic vacuum polarization at Mz
    double V_TB;  // the CKM element V_tb
    double UMASS; // the constituent mass of the up quark
    double DMASS; // the constituent mass of the down quark

    /* Outputs */
    double alphaMZ;
    double Mw;
    complex rhoZ_l[6], rhoZ_q[6];
    complex kappaZ_l[6], kappaZ_q[6];   
    
};


#endif	/* ZFITTER_H */
    