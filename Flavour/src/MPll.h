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
    void checkCache();
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
    
    /**
    * @brief \f$ <Sigma_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Sigma_i from q_min to q_max
    */
    double integrateSigma_e(int i, double q_min, double q_max);
    
    /**
    * @brief \f$ <Delta_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Delta_i from q_min to q_max
    */
    double integrateDelta_e(int i, double q_min, double q_max);
    
    /**
    * @brief \f$ <Sigma_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Sigma_i from q_min to q_max
    */
    double integrateSigma_mu(int i, double q_min, double q_max);
    
    /**
    * @brief \f$ <Delta_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Delta_i from q_min to q_max
    */
    double integrateDelta_mu(int i, double q_min, double q_max);
    

private:
    const StandardModel& mySM;
    int iter;
    
    unsigned int fplus_updated;
    gslpp::vector<double> fplus_cache;
    
    unsigned int fT_updated;
    gslpp::vector<double> fT_cache;
    
    unsigned int f0_updated;
    double f0_cache;
    
    unsigned int k2_updated;
    gslpp::vector<double> k2_cache;
    
    unsigned int beta_e_updated;
    unsigned int beta_mu_updated;
    gslpp::vector<double> beta_cache;
    
    unsigned int lambda_updated;
    double lambda_cache;
    
    unsigned int F_e_updated;
    unsigned int F_mu_updated;
    
    unsigned int VL_updated;
    
    unsigned int VR_updated;
    
    unsigned int TL_updated;
    
    unsigned int TR_updated;
    
    unsigned int SL_updated;
    gslpp::vector<double> SL_cache;
    
    unsigned int SR_updated;
    
    unsigned int N_updated;
    gslpp::vector<double> N_cache;
    gslpp::complex Nc_cache;
    
    unsigned int C_7_updated;
    gslpp::complex C_7_cache;

    unsigned int C_9_updated;
    gslpp::complex C_9_cache;
    
    unsigned int C_10_updated;
    gslpp::complex C_10_cache;
    
    unsigned int C_7p_updated;
    gslpp::complex C_7p_cache;
    
    unsigned int C_9p_updated;
    gslpp::complex C_9p_cache;
    
    unsigned int C_10p_updated;
    gslpp::complex C_10p_cache;
    
    unsigned int C_S_updated;
    gslpp::complex C_S_cache;
    
    unsigned int C_P_updated;
    gslpp::complex C_P_cache;
    
    unsigned int C_Sp_updated;
    gslpp::complex C_Sp_cache;
    
    unsigned int C_Pp_updated;
    gslpp::complex C_Pp_cache;
    
    unsigned int H_V0updated;
    gslpp::vector<double> H_V0cache;
    gslpp::complex H_V0Ccache[2];
    
    unsigned int H_A0updated;
    
    unsigned int H_Supdated;
    gslpp::vector<double> H_Scache;
    
    unsigned int H_Pe_updated;
    gslpp::vector<double> H_Pe_cache;
    
    unsigned int H_Pmu_updated;
    gslpp::vector<double> H_Pmu_cache;
    
    unsigned int I0e_updated;
    unsigned int I2e_updated;
    unsigned int I8e_updated;
    
    
    unsigned int I0mu_updated;
    unsigned int I2mu_updated;
    unsigned int I8mu_updated;
    
    std::map<std::pair<double, double>, unsigned int > sigma0eCached;
    std::map<std::pair<double, double>, unsigned int > sigma2eCached;
    
    std::map<std::pair<double, double>, unsigned int > delta0eCached;
    std::map<std::pair<double, double>, unsigned int > delta2eCached;
    
    std::map<std::pair<double, double>, unsigned int > sigma0muCached;
    std::map<std::pair<double, double>, unsigned int > sigma2muCached;
    
    std::map<std::pair<double, double>, unsigned int > delta0muCached;
    std::map<std::pair<double, double>, unsigned int > delta2muCached;
    
    double avaSigma0e;
    double avaSigma2e;
    
    double errSigma0e;
    double errSigma2e;
    
    double avaDelta0e;
    double avaDelta2e;
    
    double errDelta0e;
    double errDelta2e;
    
    double avaSigma0mu;
    double avaSigma2mu;
    
    double errSigma0mu;
    double errSigma2mu;
    
    double avaDelta0mu;
    double avaDelta2mu;
    
    double errDelta0mu;
    double errDelta2mu;
    
    gsl_function FS0e;
    gsl_function FS2e;
    
    gsl_function FD0e;
    gsl_function FD2e;
    
    gsl_function FS0mu;
    gsl_function FS2mu;
    
    gsl_function FD0mu;
    gsl_function FD2mu;
    
    gsl_integration_workspace * w_sigma0e;
    gsl_integration_workspace * w_sigma2e;
    
    gsl_integration_workspace * w_delta0e;
    gsl_integration_workspace * w_delta2e;
    
    gsl_integration_workspace * w_sigma0mu;
    gsl_integration_workspace * w_sigma2mu;
    
    gsl_integration_workspace * w_delta0mu;
    gsl_integration_workspace * w_delta2mu;
    
    std::map<std::pair<double, double>, double > cacheSigma0e;
    std::map<std::pair<double, double>, double > cacheSigma2e;
    
    std::map<std::pair<double, double>, double > cacheDelta0e;
    std::map<std::pair<double, double>, double > cacheDelta2e;
    
    std::map<std::pair<double, double>, double > cacheSigma0mu;
    std::map<std::pair<double, double>, double > cacheSigma2mu;
    
    std::map<std::pair<double, double>, double > cacheDelta0mu;
    std::map<std::pair<double, double>, double > cacheDelta2mu;
    
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

    