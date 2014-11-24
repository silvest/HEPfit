/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLL_H
#define	MPLL_H

#include <math.h>
#include "Flavour.h"
#include "MVll.h"
#include <StandardModel.h>
#include <ThObservable.h>
#include <gsl/gsl_integration.h>
#include <assert.h>


#define CUTOFF 10    //cutoff between LCSR and lattice values for Form Factors, in GeV^2

/**
 * @class MPll
 * @ingroup flavour
 * @brief A class for the decay B -> K^*ll. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class MPll : public ThObservable {
public:
    MPll(const StandardModel& SM_i);
    virtual ~MPll();
    void updateParameters();
    //void checkCache( double qmin, double qmax);
    virtual double computeThValue()=0;
    
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Me;            //electron mass
    double Mmu;           //muon mass
    double MB;            //B meson mass
    double MK;            //K star meson mass
    double Mb;            //b quark mass
    double mu_b;          //b mass scale
    double width_Bd;      //B meson width
    double Ms;            //s quark mass
    double MW;            //W boson mass
    complex lambda_t;     //Vckm factor
    double b;             //BF of the decay K^* -> K pi
    complex h_0;          //parameter that contains the contribution from the hadronic hamiltonian
    complex h_0_1;          //parameter that contains the contribution from the hadronic hamiltonian
    double q2;            //q^2 of the decay
    
    /*LCSR fit parameters*/
    double r_1_fplus, r_2_fplus, m_fit2_fplus;
    double r_1_fT, r_2_fT, m_fit2_fT;
    double r_2_f0, m_fit2_f0;
    

    vector<complex> ** allcoeff;
    vector<complex> ** allcoeffprime;
    
    complex C_7;
    complex C_9;
    complex C_10;
    complex C_S;
    complex C_P;
    
    complex C_7p;
    complex C_9p;
    complex C_10p;
    complex C_Sp;
    complex C_Pp;
    
    
    /**
    * @brief \f$ f_+ \f$
    * @param[in] q2 q^2 of the decay
    * @return return the form factor f_+
    */
    double f_plus(double q2);
    
    
    /**
    * @brief \f$ f_T \f$
    * @param[in] q2 q^2 of the decay
    * @return return the form factor f_T
    */
    double f_T(double q2);
    
    
    /**
    * @brief \f$ f_0 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the form factor f_0
    */
    double f_0(double q2);
    
    /**
    * @brief \f$ V_L \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_L
    */
    gslpp::complex V_L(double q2);

    
    /**
    * @brief \f$ V_R \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_R
    */
    gslpp::complex V_R(double q2);


    /**
    * @brief \f$ T_L \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_L
    */
    gslpp::complex T_L(double q2);


    /**
    * @brief \f$ T_R \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_R
    */
    gslpp::complex T_R(double q2);


    /**
    * 
    * @brief \f$ S_L \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_L
    */
    double S_L(double q2);


    /**
    * 
    * @brief \f$ S_R \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_R
    */
    double S_R(double q2);
    
    
    /**
    * @brief \f$ N \f$ 
    * @return return the helicity amplitude normalization factor N
    */
    
    complex N();
    
    
    /**
    * @brief \f$ H_V(\lambda) \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_V(lambda)
    */
    gslpp::complex H_V(double q2, int bar);


    /**
    * @brief \f$ H_A(\lambda) \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_A(lambda)
    */
    gslpp::complex H_A(double q2, int bar);


    /**
    * @brief \f$ H_S \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_S
    */
    gslpp::complex H_S(double q2, int bar);


    /**
    * @brief \f$ H_P \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @param[in] Mlep mass of the lepton
    * @return return the helicity amplitude H_P
    */
    gslpp::complex H_P(double q2, int bar, double Mlep);
    
    
    /**
    * @brief \f$ k^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the 3-momentum of the recoiling meson in the B rest frame
    */
    double k2 (double q2);
    
    
    /**
    * @brief \f$ beta \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] Mlep mass of the lepton
    * @return return the factor beta used in the angular coefficients I_i
    */
    double beta (double q2, double Mlep);
    
    
    /**
    * @brief \f$ lambda \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor lambda used in the angular coefficients I_i
    */
    double lambda(double q2);

    
    /**
    * @brief \f$ F \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] Mlep mass of the lepton
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2, double Mlep);
    
    
    /**
    * i values:
    * 0 = 1c
    * 2 = 2c
    * 8 = 6c
    */
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @param[in] Mlep mass of the lepton
    * @return return the angular coefficient I_i
    */
    double  I(int i, double q2, int bar, double Mlep);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @param[in] Mlep mass of the lepton
    * @return return the CP average Sigma_i
    */
    double Sigma(int i, double q2, double Mlep);
    
    
    /**
    * @brief \f$ Delta_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @param[in] Mlep mass of the lepton
    * @return return the CP asymmetry Delta_i
    */
    double Delta(int i, double q2, double Mlep);
    
    /**
    * @brief \f$ Sigma_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1c for electron channel
    */
    double getSigma0e(double q2){
        return Sigma(0, q2, Me);
    };
    
    /**
    * @brief \f$ Sigma_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1c for muon channel
    */
    double getSigma0mu(double q2){
        return Sigma(0, q2, Mmu);
    };
    
    /**
    * @brief \f$ Sigma_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2c for electron channel
    */
    double getSigma2e(double q2){
        return Sigma(2, q2, Me);
    };
    
    /**
    * @brief \f$ Sigma_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2c for muon channel
    */
    double getSigma2mu(double q2){
        return Sigma(2, q2, Mmu);
    };
    
    /**
    * @brief \f$ Delta_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1c for electron channel
    */
    double getDelta0e(double q2){
        return Delta(0, q2, Me);
    };
    
    /**
    * @brief \f$ Delta_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1c for muon channel
    */
    double getDelta0mu(double q2){
        return Delta(0, q2, Mmu);
    };
    
    /**
    * @brief \f$ Delta_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2c for electron channel
    */    
    double getDelta2e(double q2){
        return Delta(2, q2, Me);
    };
    
    /**
    * @brief \f$ Delta_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2c for muon channel
    */    
    double getDelta2mu(double q2){
        return Delta(2, q2, Mmu);
    };
    

private:
    const StandardModel& mySM;
    MVll myMVll;
};



/**
 * @class Branching Fraction for electron channel
 * @ingroup flavour
 * @brief A class for the clean observable BR_e. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MPll_e : public MPll{
public:
    
    /**
    * @brief \f$ BR_e \f$ 
    */
    BR_MPll_e(const StandardModel& SM_i);
    
    /**
    * @return return the clean observable BR_e
    */
    double computeBR_MPll_e(double qmin, double qmax);
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2;
    double avaSigma0, errSigma0, avaSigma2, errSigma2;
};


/**
 * @class Branching Fraction for muon channel
 * @ingroup flavour
 * @brief A class for the clean observable BR_mu. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BR_MPll_mu : public BR_MPll_e{
public:
    
    /**
    * @brief \f$ BR_mu \f$ 
    */
    BR_MPll_mu(const StandardModel& SM_i);
    
    /**
    * @return return the clean observable BR_mu
    */
    double computeBR_MPll_mu(double qmin, double qmax);
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2;
    double avaSigma0, errSigma0, avaSigma2, errSigma2;
};


/**
 * @class ratio between BR for electron and muon channels
 * @ingroup flavour
 * @brief A class for the Branching Fraction ratio. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class R_MPll : public BR_MPll_mu{
public:
    
    /**
    * @brief \f$ BR_e/BR_mu \f$ 
    */
    R_MPll(const StandardModel& SM_i);
    
    /**
    * @return the ratio between branching fractions of \f$ B\to K e^+ e^- \f$ and \f$ B\to K \mu^+ \mu^- \f$
    */
    double computeThValue ();
};


/**
 * @class ACP for electron channel
 * @ingroup flavour
 * @brief A class for the clean observable ACP_e. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP_e : public BR_MPll_e{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP_e(const StandardModel& SM_i);
    
    /**
    * @return return the clean observable ACP_e
    */
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2;
    double avaDelta0, errDelta0, avaDelta2, errDelta2;
};


/**
 * @class ACP for muon channel
 * @ingroup flavour
 * @brief A class for the clean observable ACP_mu. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP_mu : public BR_MPll_mu{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP_mu(const StandardModel& SM_i);
    
    /**
    * @return return the clean observable ACP_mu
    */
    double computeThValue ();
    
protected:
    
    
private:
    gsl_function F1, F2;
    double avaDelta0, errDelta0, avaDelta2, errDelta2;
};

#endif	/* MPLL_H */

    