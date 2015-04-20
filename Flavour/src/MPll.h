/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLL_H
#define	MPLL_H

#include <math.h>
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
class MPll{
public:
    MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    virtual ~MPll();
    void updateParameters();
    void checkCache();
    
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Mlep;          //muon mass
    double MM;            //initial meson mass
    double MP;            //final pseudoscalar meson mass
    double Mb;            //b quark mass
    double Mc;            //c quark mass
    double mu_b;          //b mass scale
    double width;         //initial meson width
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
    
    gslpp::complex C_1;
    gslpp::complex C_2;
    gslpp::complex C_3;
    gslpp::complex C_4;
    gslpp::complex C_5;
    gslpp::complex C_6;
    gslpp::complex C_7;
    gslpp::complex C_9;
    gslpp::complex C_10;
    gslpp::complex C_S;
    gslpp::complex C_P;
    
    gslpp::complex C_7p;
    gslpp::complex C_9p;
    gslpp::complex C_10p;
    gslpp::complex C_Sp;
    gslpp::complex C_Pp;
    
    
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
    
    gslpp::complex N();
    
    
    /**
    * 
    * @brief \f$ h(q^2,m) \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] m mass
    * @return return the h(q^2,m) function involved into C_9^eff
    */
    gslpp::complex H(double q2, double m);
    
    
    /**
    * 
    * @brief \f$ Y(q^2) \f$
    * @param[in] q2 q^2 of the decay
    * @return return the Y(q^2) function involved into C_9^eff
    */
    gslpp::complex Y(double q2);
    
    
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
    * @return return the helicity amplitude H_P
    */
    gslpp::complex H_P(double q2, int bar);
    
    
    /**
    * @brief \f$ k^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the 3-momentum of the recoiling meson in the B rest frame
    */
    double k2 (double q2);
    
    
    /**
    * @brief \f$ beta \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor beta used in the angular coefficients I_i
    */
    double beta (double q2);
    
    
    /**
    * @brief \f$ lambda \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor lambda used in the angular coefficients I_i
    */
    double lambda(double q2);

    
    /**
    * @brief \f$ F \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2);
    
    
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
    * @return return the angular coefficient I_i
    */
    double  I(int i, double q2, int bar);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_i
    */
    double Sigma(int i, double q2);
    
    
    /**
    * @brief \f$ Delta_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_i
    */
    double Delta(int i, double q2);
    
    /**
    * @brief \f$ Sigma_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1c
    */
    double getSigma0(double q2){
        return Sigma(0, q2);
    };
    
    /**
    * @brief \f$ Sigma_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2c
    */
    double getSigma2(double q2){
        return Sigma(2, q2);
    };
    
    /**
    * @brief \f$ Delta_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1c
    */
    double getDelta0(double q2){
        return Delta(0, q2);
    };
    
    /**
    * @brief \f$ Delta_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2c
    */    
    double getDelta2(double q2){
        return Delta(2, q2);
    };
    
    /**
    * @brief \f$ <Sigma_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Sigma_i from q_min to q_max
    */
    double integrateSigma(int i, double q_min, double q_max);
    
    /**
    * @brief \f$ <Delta_{i}> \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return return the CP average integral of Delta_i from q_min to q_max
    */
    double integrateDelta(int i, double q_min, double q_max);
    

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson pseudoscalar;
    
    unsigned int fplus_updated;
    gslpp::vector<double> fplus_cache;
    
    unsigned int fT_updated;
    gslpp::vector<double> fT_cache;
    
    unsigned int f0_updated;
    double f0_cache;
    
    unsigned int k2_updated;
    gslpp::vector<double> k2_cache;
    
    unsigned int beta_updated;
    double beta_cache;
    
    unsigned int lambda_updated;
    double lambda_cache;
    
    unsigned int F_updated;
    
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
    
    unsigned int C_1_updated;
    gslpp::complex C_1_cache;

    unsigned int C_2_updated;
    gslpp::complex C_2_cache;
    
    unsigned int C_3_updated;
    gslpp::complex C_3_cache;
    
    unsigned int C_4_updated;
    gslpp::complex C_4_cache;
    
    unsigned int C_5_updated;
    gslpp::complex C_5_cache;
    
    unsigned int C_6_updated;
    gslpp::complex C_6_cache;
    
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
    
    unsigned int Yupdated;
    gslpp::vector<double> Ycache;
    
    unsigned int H_V0updated;
    gslpp::vector<double> H_V0cache;
    gslpp::complex H_V0Ccache[2];
    
    unsigned int H_A0updated;
    
    unsigned int H_Supdated;
    gslpp::vector<double> H_Scache;
    
    unsigned int H_P_updated;
    gslpp::vector<double> H_P_cache;
    
    unsigned int I0_updated;
    unsigned int I2_updated;
    unsigned int I8_updated;
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;
    std::map<std::pair<double, double>, unsigned int > delta2Cached;
    
    double avaSigma0;
    double avaSigma2;
    
    double errSigma0;
    double errSigma2;
    
    double avaDelta0;
    double avaDelta2;
    
    double errDelta0;
    double errDelta2;
    
    gsl_function FS0;
    gsl_function FS2;
    
    gsl_function FD0;
    gsl_function FD2;
    
    gsl_integration_workspace * w_sigma0;
    gsl_integration_workspace * w_sigma2;
    
    gsl_integration_workspace * w_delta0;
    gsl_integration_workspace * w_delta2;
    
    std::map<std::pair<double, double>, double > cacheSigma0;
    std::map<std::pair<double, double>, double > cacheSigma2;
    
    std::map<std::pair<double, double>, double > cacheDelta0;
    std::map<std::pair<double, double>, double > cacheDelta2;
    
};

#endif	/* MPLL_H */

    