/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVLL_H
#define	MVLL_H

#include <math.h>
#include <StandardModel.h>
#include <ThObservable.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>


#define CUTOFF 10    //cutoff between LCSR and lattice values for Form Factors, in GeV^2

/*******************************************************************************
 * GSL Function Conversion BEGIN                                                  *
 * ****************************************************************************/

template<class F>
static double gslFunctionAdapter( double x, void* p)
{
    // Here I do recover the "right" pointer, safer to use static_cast
    // than reinterpret_cast.
    F* function = static_cast<F*>( p );
    return (*function)( x );
}

template<class F>
gsl_function convertToGslFunction( const F& f )
{
    gsl_function gslFunction;
    
    const void* p = &f;
    assert (p != 0);
    
    gslFunction.function = &gslFunctionAdapter<F>;
    // Just to eliminate the const.
    gslFunction.params = const_cast<void*>( p );
    
    return gslFunction;
}

/*******************************************************************************
 * GSL Function conversion END                                                     *
 * ****************************************************************************/

/**
 * @class MVll
 * @ingroup flavour
 * @brief A class for the decay B -> K^*ll. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class MVll {
public:
    MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    virtual ~MVll();
    void updateParameters();
    void checkCache();
    //virtual double computeThValue()=0;
    
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Mlep;            //muon mass
    double MB;            //B meson mass
    double MKstar;        //K star meson mass
    double Mb;            //b quark mass
    double mu_b;          //b mass scale
    double width_Bd;      //B meosn width
    double Ms;            //s quark mass
    double MW;            //W boson mass
    complex lambda_t;     //Vckm factor
    double b;             //BF of the decay K^* -> K pi
    complex h[3];         //parameter that contains the contribution from the hadronic hamiltonian  
    double q2;            //q^2 of the decay
    
    /*lattice fit parameters*/
    double a_0V, a_1V, dmV;
    double a_0A0, a_1A0, dmA0;
    double a_0A1, a_1A1, dmA1;
    double a_0A12, a_1A12, dmA12;
    double a_0T1, a_1T1, dmT1;
    double a_0T2, a_1T2, dmT2;
    double a_0T23, a_1T23, dmT23;
    
    /*LCSR fit parameters*/
    double r_1V, r_2V, m_RV, m_fit2V;
    double r_1A0, r_2A0, m_RA0, m_fit2A0;
    double r_2A1, m_fit2A1;
    double r_1A2, r_2A2, m_fit2A2;
    double r_1T1, r_2T1, m_RT1, m_fit2T1;
    double r_2T2, m_fit2T2;
    double r_1T3t, r_2T3t, m_fit2T3t;

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
    * @brief \f$ LCSR_fit1 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_R2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the first fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit1(double q2, double r_1, double r_2, double m_R2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit2 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the second fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit2(double q2, double r_1, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit3 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the third fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit3(double q2, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ z \f$
    * @param[in] q2 q^2 of the decay
    * @return return the lattice parameter z from arXiv:1310.3722v3
    */
    double z(double q2);
    
    
    /**
    * @brief \f$ lat_fit(q^2) \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] c_01 fit parameter
    * @param[in] c_01s fit parameter
    * @param[in] q2 q^2 of the decay
    * @return return the lattice pole factor P(q^2,dm) from arXiv:1310.3722v3
    */
    double lat_fit(double q2, double a_0, double a_1, double dm);
    
    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double V(double q2);

    
    /**
    * @brief \f$ A_0 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_0(q^2)
    */
    double A_0(double q2);

    
    /**
    * @brief \f$ A_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_1(q^2)
    */
    double A_1(double q2);

    
    /**
    * @brief \f$ A_2 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_2(q^2)
    */
    double A_2(double q2);

    
    /**
    * @brief \f$ T_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_1(q^2)
    */
    double T_1(double q2);

    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double T_2(double q2);

    
    /**
    * @brief \f$ T_3tilde \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3tilde(q^2)
    */
    double T_3tilde(double q2);

    
    /**
    * @brief \f$ T_3 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3(q^2)
    */
    double T_3(double q2);
    
    
    /**
    * polarization   i
    *      0         0
    *      +         1
    *      -         2
    */

    
    /**
    * @brief \f$ V_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_L(lambda)
    */
    double V_L(int i, double q2);

    
    /**
    * @brief \f$ V_R \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_R(lambda)
    */
    double V_R(int i, double q2);


    /**
    * @brief \f$ T_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_L(lambda)
    */
    double T_L(int i, double q2);


    /**
    * @brief \f$ T_R \f$ 
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_R(lambda)
    */
    double T_R(int i, double q2);


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
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_V(lambda)
    */
    gslpp::complex H_V(int i, double q2, int bar);


    /**
    * @brief \f$ H_A(\lambda) \f$ 
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return return the helicity amplitude H_A(lambda)
    */
    gslpp::complex H_A(int i, double q2, int bar);


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
    * @param[in] b_i BF of the decay K* -> K pi
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2, double b_i);
    
    
    /**
    * i values:
    * 0 = 1c
    * 1 = 1s
    * 2 = 2c
    * 3 = 2s
    * 4 = 3
    * 5 = 4
    * 6 = 5
    * 7 = 6s
    * 8 = 6c
    * 9 = 7
    * 10 = 8
    * 11 = 9
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
    
    
    /**
    * @brief \f$ Delta_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_i
    */
    double Delta(int i, double q2);
    
    /**
    * @brief \f$ Sigma_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1c
    */
    double getSigma0(double q2){
        return Sigma(0, q2);
    };
    
    /**
    * @brief \f$ Sigma_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_1s
    */
    double getSigma1(double q2){
        return Sigma(1, q2);
    };
    
    /**
    * @brief \f$ Sigma_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2c
    */
    double getSigma2(double q2){
        return Sigma(2, q2);
    };
    
    /**
    * @brief \f$ Sigma_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_2s
    */
    double getSigma3(double q2){
        return Sigma(3, q2);
    };
    
    /**
    * @brief \f$ Sigma_{3} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_3
    */    
    double getSigma4(double q2){
        return Sigma(4, q2);
    };
    
    /**
    * @brief \f$ Sigma_{4} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_4
    */
    double getSigma5(double q2){
        return Sigma(5, q2);
    };
    
    /**
    * @brief \f$ Sigma_{5} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_5
    */
    double getSigma6(double q2){
        return Sigma(6, q2);
    };
    
    /**
    * @brief \f$ Sigma_{6s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_6s
    */
    double getSigma7(double q2){
        return Sigma(7, q2);
    };
    
    /**
    * @brief \f$ Sigma_{7} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_7
    */
    double getSigma9(double q2){
        return Sigma(9, q2);
    };
    
    /**
    * @brief \f$ Sigma_{9} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_9
    */
    double getSigma11(double q2){
        return Sigma(11, q2);
    };
    
    /**
    * @brief \f$ Delta_{1s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1s
    */
    double getDelta0(double q2){
        return Delta(0, q2);
    };
    
    /**
    * @brief \f$ Delta_{1c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_1c
    */
    double getDelta1(double q2){
        return Delta(1, q2);
    };
    
    /**
    * @brief \f$ Delta_{2s} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2s
    */    
    double getDelta2(double q2){
        return Delta(2, q2);
    };
    
    /**
    * @brief \f$ Delta_{2c} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_2c
    */
    double getDelta3(double q2){
        return Delta(3, q2);
    };
    
    /**
    * @brief \f$ Delta_{9} \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_9
    */
    double getDelta11(double q2){
        return Delta(11, q2);
    };
    
    /**
    * @brief \f$ |H_V(0)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(0)
    */
    double getHV0_abs2(double q2){
        return H_V(0, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_V(1)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(1)
    */
    double getHV1_abs2(double q2){
        return H_V(1, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_V(2)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_V(2)
    */
    double getHV2_abs2(double q2){
        return H_V(2, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(0)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(0)
    */
    double getHA0_abs2(double q2){
        return H_A(0, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(1)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(1)
    */
    double getHA1_abs2(double q2){
        return H_A(1, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_A(2)|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_A(2)
    */
    double getHA2_abs2(double q2){
        return H_A(2, q2, 0).abs2();
    };
    
    /**
    * @brief \f$ |H_S|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_S
    */
    double getHS_abs2(double q2){
        return H_S(q2, 0).abs2();
    }
    
    /**
    * @brief \f$ |H_P|^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the absolute value of the helicity amplitude H_P
    */
    double getHP_abs2(double q2){
        return H_P(q2, 0).abs2();
    }
    
    
    double getFactor(double q2){
        return q2/(2*Mlep*Mlep);
    }

private:
    const StandardModel& mySM;
    StandardModel::lepton lep;
    StandardModel::meson meson;
    StandardModel::meson vectorM;
    int iter;
    
    std::map<std::pair<double, double>, double > cacheSigma0;
    std::map<std::pair<double, double>, double > cacheSigma1;
    std::map<std::pair<double, double>, double > cacheSigma2;
    std::map<std::pair<double, double>, double > cacheSigma3;
    std::map<std::pair<double, double>, double > cacheSigma4;
    std::map<std::pair<double, double>, double > cacheSigma5;
    std::map<std::pair<double, double>, double > cacheSigma6;
    std::map<std::pair<double, double>, double > cacheSigma7;
    std::map<std::pair<double, double>, double > cacheSigma9;
    std::map<std::pair<double, double>, double > cacheSigma11;
    
    std::map<std::pair<double, double>, double > cacheDelta0;
    std::map<std::pair<double, double>, double > cacheDelta1;
    std::map<std::pair<double, double>, double > cacheDelta2;
    std::map<std::pair<double, double>, double > cacheDelta3;
    std::map<std::pair<double, double>, double > cacheDelta11;
    
    double avaSigma0;
    double avaSigma1;
    double avaSigma2;
    double avaSigma3;
    double avaSigma4;
    double avaSigma5;
    double avaSigma6;
    double avaSigma7;
    double avaSigma9;
    double avaSigma11;
    
    double errSigma0;
    double errSigma1;
    double errSigma2;
    double errSigma3;
    double errSigma4;
    double errSigma5;
    double errSigma6;
    double errSigma7;
    double errSigma9;
    double errSigma11;
    
    double avaDelta0;
    double avaDelta1;
    double avaDelta2;
    double avaDelta3;
    double avaDelta11;
    
    double errDelta0;
    double errDelta1;
    double errDelta2;
    double errDelta3;
    double errDelta11;
    
    
    gsl_function FS0;
    gsl_function FS1;
    gsl_function FS2;
    gsl_function FS3;
    gsl_function FS4;
    gsl_function FS5;
    gsl_function FS6;
    gsl_function FS7;
    gsl_function FS9;
    gsl_function FS11;
    
    gsl_function FD0;
    gsl_function FD1;
    gsl_function FD2;
    gsl_function FD3;
    gsl_function FD11;
    
    gsl_integration_workspace * w_sigma0;
    gsl_integration_workspace * w_sigma1;
    gsl_integration_workspace * w_sigma2;
    gsl_integration_workspace * w_sigma3;
    gsl_integration_workspace * w_sigma4;
    gsl_integration_workspace * w_sigma5;
    gsl_integration_workspace * w_sigma6;
    gsl_integration_workspace * w_sigma7;
    gsl_integration_workspace * w_sigma9;
    gsl_integration_workspace * w_sigma11;
    
    gsl_integration_workspace * w_delta0;
    gsl_integration_workspace * w_delta1;
    gsl_integration_workspace * w_delta2;
    gsl_integration_workspace * w_delta3;
    gsl_integration_workspace * w_delta11;
    
    
    unsigned int N_updated;
    vector<double> N_cache;
    complex Nc_cache;
    
    unsigned int V_updated;
    vector<double> V_cache;
    
    unsigned int A0_updated;
    vector<double> A0_cache;
    
    unsigned int A1_updated;
    vector<double> A1_cache;
    
    unsigned int A2_updated;
    vector<double> A2_cache;
    
    unsigned int T1_updated;
    vector<double> T1_cache;
    
    unsigned int T2_updated;
    vector<double> T2_cache;
    
    unsigned int T3t_updated;
    vector<double> T3t_cache;
    
    unsigned int T3_updated;
    
    unsigned int k2_updated;
    vector<double> k2_cache;
    
    unsigned int z_updated;
    
    unsigned int lambda_updated;
    double lambda_cache;
    
    unsigned int beta_updated;
    double beta_cache;
    
    unsigned int F_updated;
    
    unsigned int VL1_updated;
    unsigned int VL2_updated;
    
    unsigned int TL1_updated;
    unsigned int TL2_updated;

    unsigned int VR1_updated;
    unsigned int VR2_updated;
    
    unsigned int TR1_updated;
    unsigned int TR2_updated;
    
    unsigned int VL0_updated;
    vector<double> VL0_cache;
    
    unsigned int TL0_updated;
    vector<double> TL0_cache;
    
    unsigned int VR0_updated;
    
    unsigned int TR0_updated;
    
    unsigned int SL_updated;
    vector<double> SL_cache;
    
    unsigned int SR_updated;
    
    unsigned int C_7_updated;
    complex C_7_cache;

    unsigned int C_9_updated;
    complex C_9_cache;
    
    unsigned int C_10_updated;
    complex C_10_cache;
    
    unsigned int C_7p_updated;
    complex C_7p_cache;
    
    unsigned int C_9p_updated;
    complex C_9p_cache;
    
    unsigned int C_10p_updated;
    complex C_10p_cache;
    
    unsigned int C_S_updated;
    complex C_S_cache;
    
    unsigned int C_P_updated;
    complex C_P_cache;
    
    unsigned int C_Sp_updated;
    complex C_Sp_cache;
    
    unsigned int C_Pp_updated;
    complex C_Pp_cache;
    
    unsigned int H_V0updated;
    vector<double> H_V0cache;
    complex H_V0Ccache;
    
    unsigned int H_V1updated;
    vector<double> H_V1cache;
    complex H_V1Ccache;
    
    unsigned int H_V2updated;
    vector<double> H_V2cache;
    complex H_V2Ccache;
    
    unsigned int H_A0updated;
    unsigned int H_A1updated;
    unsigned int H_A2updated;
    
    unsigned int H_Supdated;
    vector<double> H_Scache;
    
    unsigned int H_Pupdated;
    vector<double> H_Pcache;
    
    unsigned int I0_updated;
    unsigned int I1_updated;
    unsigned int I2_updated;
    unsigned int I3_updated;
    unsigned int I4_updated;
    unsigned int I5_updated;
    unsigned int I6_updated;
    unsigned int I7_updated;
    unsigned int I8_updated;
    unsigned int I9_updated;
    unsigned int I10_updated;
    unsigned int I11_updated;
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;
    std::map<std::pair<double, double>, unsigned int > sigma1Cached;
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;
    std::map<std::pair<double, double>, unsigned int > sigma3Cached;
    std::map<std::pair<double, double>, unsigned int > sigma4Cached;
    std::map<std::pair<double, double>, unsigned int > sigma5Cached;
    std::map<std::pair<double, double>, unsigned int > sigma6Cached;
    std::map<std::pair<double, double>, unsigned int > sigma7Cached;
    std::map<std::pair<double, double>, unsigned int > sigma9Cached;
    std::map<std::pair<double, double>, unsigned int > sigma11Cached;
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;
    std::map<std::pair<double, double>, unsigned int > delta1Cached;
    std::map<std::pair<double, double>, unsigned int > delta2Cached;
    std::map<std::pair<double, double>, unsigned int > delta3Cached;
    std::map<std::pair<double, double>, unsigned int > delta11Cached;
    
    
    
};

#endif	/* MVLL_H */

