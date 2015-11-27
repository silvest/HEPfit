/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SUSYMATCHING_H
#define	SUSYMATCHING_H

#include <gslpp.h>
#include <complex>
#include <PVfunctions.h>
#include "StandardModelMatching.h"

#define LEPS 1.e-5     // tolerance in the limit of S(x,y) to S(x) 
#define SUSYLEPS 5.e-3 // tolerance in the limits of D0[x,y,z,t] and D2[x,y,z,t]
#define SUSYLEPS2 1.e-10 // tolerance in the limits of B0[x,y]
#define SUSYLEPS3 5.e-3 // tolerance in the limits of D0N

class SUSY;

/**
 * @class SUSYMatching
 * @ingroup SUSY
 * @brief A class for the matching in the MSSM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class SUSYMatching : public StandardModelMatching {
public:
    SUSYMatching(const SUSY & SUSY_i);
   
    /**
     *
     * @brief Wilson coefficient for the process \f$ \ell_j \to \ell_i \gamma \f$.
     * @param[in] li_lj determines the process, e.g., 1 = \f$ \mu \to e \gamma \f$, 2 = \f$ \tau \to \mu \gamma \f$, 3 = \f$ \tau \to e \gamma \f$
     * @return returns the vector of SUSY Wilson coefficients
     */
    virtual std::vector<WilsonCoefficient>& CMDLij(int li_lj);

    /**
     *
     * @brief Wilson coefficient for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] li_lj determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of SUSY Wilson coefficients
     */
    virtual std::vector<WilsonCoefficient>& CMDLi3j(int li_lj);

    /**
     *
     * @brief Wilson coefficient for the process \f$ \mu \to e \f$ conversion in Nuclei.
     * @return returns the vector of SUSY Wilson coefficients
     */
    virtual std::vector<WilsonCoefficient>& CMmueconv();

    /**
     *
     * @brief Wilson coefficient for \f$ (g-2)_{\mu} \f$ at one-loop.
     * @return returns the vector of SUSY Wilson coefficients
     */
    virtual std::vector<WilsonCoefficient>& CMgminus2mu();

//   /** Calculates gamma penguin amplitudes for m->(3)e (1), t->(3)m (2) and t->(3)e (3)**/
    /**
     * 
     * @brief Calculates gamma penguin amplitudes for the process \f$ \ell_j \to \ell_i \gamma \f$ from \cite Hisano:1995cp \cite Arganda:2005ji.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to e \gamma \f$, 2 = \f$ \tau \to \mu \gamma \f$, 3 = \f$ \tau \to e \gamma \f$
     * @return returns the vector of gamma penguin amplitude
     */
    virtual gslpp::vector<gslpp::complex> AFunctions(int n);
    
    //   /** Calculates Z penguin amplitudes for m->3e (1), t->3m (2) and t->3e (3)**/
    /**
     * 
     * @brief Calculates Z penguin amplitudes for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$ from \cite Hisano:1995cp \cite Arganda:2005ji.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of Z penguin amplitude
     */
    virtual gslpp::vector<gslpp::complex> FFunctions(int n);


//   /** Calculates box amplitudes for m->3e (1), t->3m (2) and t->3e (3) **/
    /**
     * @brief Calculates box diagram amplitudes for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$ from from \cite Hisano:1995cp \cite Arganda:2005ji.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of gamma diagram amplitude
     */
    virtual gslpp::vector<gslpp::complex> BFunctions(int n);

//   /** Calculates Higgs penguin amplitudes for m->3e (1), t->3m (2) and t->3e (3)**/
    /**
     * @brief Calculates Higgs penguin amplitudes for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$ from \cite Arganda:2005ji.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of Higgs penguin amplitude
     */
    virtual gslpp::vector<gslpp::complex> BHFunctions(int n);

//   /** Calculates box amplitudes for m->e conversion**/
    /**
     * @brief Calculates box diagram amplitudes for the process \f$ \mu \to e \f$ conversion in Nuclei from from \cite Hisano:1995cp.
     * @return returns the vector of box diagram amplitude
     */    
    virtual gslpp::vector<gslpp::complex> DFunctions();

//   /** Calculates the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at one-loop from \cite Hisano:1995cp.
     * @return returns the vector of \f$ (g-2)_{\mu} \f$ amplitude
     */    
    virtual gslpp::vector<gslpp::complex> gminus2mu();

   /** Calculates C7 and C7' for m->e (1), t->m (2) and t->e (3)**/
    /**
     * 
     * @brief Calculates \f$ C_{7} \f$ and \f$ C^{\prime}_{7} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \gamma \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to e \gamma \f$, 2 = \f$ \tau \to \mu \gamma \f$, 3 = \f$ \tau \to e \gamma \f$
     * @return returns the vector of \f$ C_{7} \f$ and \f$ C^{\prime}_{7} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> C7_Lepton(int n);

   /** Calculates C9 and C9' for m->eee (1), t->mmm (2), t->eee (3) and t->muee (4)**/
    /**
     * 
     * @brief Calculates \f$ C_{9} \f$ and \f$ C^{\prime}_{9} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of \f$ C_{9} \f$ and \f$ C^{\prime}_{9} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> C9_Lepton(int n);

   /** Calculates C10 and C10' for m->eee (1), t->mmm (2), t->eee (3) and t->muee (4)**/
    /**
     * 
     * @brief Calculates \f$ C_{10} \f$ and \f$ C^{\prime}_{10} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of \f$ C_{10} \f$ and \f$ C^{\prime}_{10} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> C10_Lepton(int n);

   /** Calculates CS and CS' for m->eee (1), t->mmm (2), t->eee (3) and t->muee (4)**/
    /**
     * 
     * @brief Calculates \f$ C_{S} \f$ and \f$ C^{\prime}_{S} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of \f$ C_{S} \f$ and \f$ C^{\prime}_{S} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> CS_Lepton(int n);

   /** Calculates CP and CP' for m->eee (1), t->mmm (2), t->eee (3) and t->muee (4)**/
    /**
     * 
     * @brief Calculates \f$ C_{P} \f$ and \f$ C^{\prime}_{P} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of \f$ C_{P} \f$ and \f$ C^{\prime}_{P} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> CP_Lepton(int n);

   /** Calculates CT and CT5 for m->eee (1), t->mmm (2), t->eee (3) and t->muee (4)**/
    /**
     * 
     * @brief Calculates \f$ C_{T} \f$ and \f$ C_{T5} \f$ Wilson coefficients for the process \f$ \ell_j \to \ell_i \ell_i \ell_i \f$.
     * @param[in] n determines the process, e.g., 1 = \f$ \mu \to eee \f$, 2 = \f$ \tau \to \mu \mu \mu \f$, 3 = \f$ \tau \to eee \f$
     * @return returns the vector of \f$ C_{T} \f$ and \f$ C_{T5} \f$ Wilson coefficients
     */
    virtual gslpp::vector<gslpp::complex> CT_Lepton(int n);
    
    /**
     *
     * @brief Updates to new SUSY parameter sets.
     * @return
     */
    
    void updateSUSYParameters();
    
private:
    const SUSY & mySUSY;
    const PVfunctions PV;

    WilsonCoefficient mcDLij,mcDLi3j,mcmueconv,mcgminus2mu;

    /**
     * @brief Sfermion mass-eigenvalue squared.
     */
    gslpp::vector<double> mym_su_sq,mym_sd_sq,mym_se_sq,mym_sn_sq;
    
    /**
     * @brief Sfermion mixing matrices.
     */
    gslpp::matrix<gslpp::complex> myRu, myRd, myRl, myRn;
    
    /**
     * @brief Chargino mass-eigenvalue
     */
    gslpp::vector<double> MChi;
    
    /**
     * @brief Chargino mixing matrices
     */
    gslpp::matrix<gslpp::complex> myV, myU;

    /**
     * @brief Neutralino mass-eigenvalue
     */
    gslpp::vector<double> MChi0, MNeig;
    /**
     * @brief Neutralino mixing matrix
     */
    gslpp::matrix<gslpp::complex> myN, ON;
    
    /**
     *@brief Functions needed to calculate various LFV observables
     */
    gslpp::matrix<double> Lepty,Leptz,Leptfa1,Leptfa2,Leptf1,Leptf2,Leptf3,Leptf4;
    std::complex<double> Leptfzn[6][4][4],Leptgzn[6][4][4],Leptfzc[3][2][2],Leptgzc[3][2][2];

    /**
     *@brief Chargino and Neutralino couplings to sfermions
     */    
    gslpp::matrix<gslpp::complex> CRlE, CRlMU, CRlTAU, CRqUP, CRqDOWN, CLlE, CLlMU, CLlTAU, CLqUP, CLqDOWN, NRlE, NRlMU, NRlTAU, NRqUP, NRqDOWN, NLlE, NLlMU, NLlTAU, NLqUP, NLqDOWN;
    
    /**
     *@brief Amplitudes of Chargino and Neutralino contribution to various LFV observables
     */    
    gslpp::matrix<gslpp::complex> AmpA1LN, AmpA1RN, AmpA1LC, AmpA1RC, AmpTauA1LN, AmpTauA1RN, AmpTauA1LC, AmpTauA1RC, AmpTEA1LN, AmpTEA1RN, AmpTEA1LC, AmpTEA1RC, AmpALN, AmpARN, AmpALC, AmpARC, AmpTauALN, AmpTauARN, AmpTauALC, AmpTauARC, AmpTEALN, AmpTEARN, AmpTEALC, AmpTEARC;

    /**
     *@brief Slepton tri-linear coupling matrix
     */
    gslpp::matrix<gslpp::complex> TEhat;

    double tanb;
    double sinb;
    double cosb;
    double v;
    double v1;
    double v2;
    double gW;

    /** Re-diagonalisation of the Neutralino mass matrix **/
    void NeutralinoRemixing();

    /** Kronecker delta **/
    /**
     * @brief Kronecker delta
     * @param[in] a is the first index of the Kronecker delta function
     * @param[in] b is the second index of the Kronecker delta function
     */
    int delta_ab(int a, int b);
};

#endif	/* SUSYMATCHING_H */
