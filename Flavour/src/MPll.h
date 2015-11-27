/* 
 * Copyright (C) 2014 HEPfit Collaboration
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
 * @ingroup Flavour
 * @brief A class for the @f$M \to P l^+ l^-@f$ decay.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$M \to P l^+ l^-@f$ decay. After the
 * parameters are updated in updateParameters() and the cache is checked in 
 * checkCache(), the form factor are build in the transverse basis in the functions
 * f_plus(), f_0() and f_T() @cite Ball:2004ye. Following @cite Jager:2012uw, 
 * they are consequentely translated in the helicity basis through the
 * functions V_L(), V_R(), T_L(), T_R(), S_L() and S_R(). Form factors and parameters
 * are combined together in the functions H_V(), H_A(), H_S() and H_P() in order
 * to build the helicity amplitudes, which are consequentely combined to create
 * the angular coefficients in the function I(). Those coefficients are used to
 * create the CP averaged coefficients in the function Sigma() ad the CP asymmetric
 * coefficients in the function Delta(). Form factors, CP averaged and asymmetric
 * coefficients and hadronic contributions are integrated in the functions 
 * integrateSigma() and integrateDelta() in order to be further used to build the observables.
 */
class MPll{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~MPll();
    
    /**
     * @brief The update parameter method for MPll.
     */
    void updateParameters();
    
    /**
     * @brief The caching method for MPll.
     */
    void checkCache();
    
    double GF;            /**<Fermi constant */
    double ale;           /**<alpha electromagnetic */
    double Mlep;          /**<muon mass */
    double MM;            /**<initial meson mass */
    double MP;            /**<final pseudoscalar meson mass */
    double Mb;            /**<b quark mass */
    double mu_b;          /**<b mass scale */
    double Mc;            /**<c quark mass */
    double Ms;            /**<s quark mass */
    double width;         /**<initial meson width */
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor */
    double b;             /**<BF of the decay V -> final states */
    gslpp::complex h_0;          /**<parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_0_1;        /**<parameter that contains the contribution from the hadronic hamiltonian */
    double q2;            /**<\f$q^2\f$ of the decay */
    
    /*LCSR fit parameters*/
    double r_1_fplus;/**<LCSR fit parameter */
    double r_2_fplus;/**<LCSR fit parameter */
    double m_fit2_fplus;/**<LCSR fit parameter */
    double r_1_fT;/**<LCSR fit parameter */
    double r_2_fT;/**<LCSR fit parameter */
    double m_fit2_fT;/**<LCSR fit parameter */
    double r_2_f0;/**<LCSR fit parameter */
    double m_fit2_f0;/**<LCSR fit parameter */
    

    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients */
    
    gslpp::complex C_1;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_3;/**<Wilson coeffients @f$C_3@f$*/
    gslpp::complex C_4;/**<Wilson coeffients @f$C_4@f$*/
    gslpp::complex C_5;/**<Wilson coeffients @f$C_5@f$*/
    gslpp::complex C_6;/**<Wilson coeffients @f$C_6@f$*/
    gslpp::complex C_7;/**<Wilson coeffients @f$C_7@f$*/
    gslpp::complex C_9;/**<Wilson coeffients @f$C_9@f$*/
    gslpp::complex C_10;/**<Wilson coeffients @f$C_{10}@f$*/
    gslpp::complex C_S;/**<Wilson coeffients @f$C_S@f$*/
    gslpp::complex C_P;/**<Wilson coeffients @f$C_P@f$*/
    
    gslpp::complex C_7p;/**<Wilson coeffients @f$C_7'@f$*/
    gslpp::complex C_9p;/**<Wilson coeffients @f$C_9'@f$*/
    gslpp::complex C_10p;/**<Wilson coeffients @f$C_{10}'@f$*/
    gslpp::complex C_Sp;/**<Wilson coeffients @f$C_S'@f$*/
    gslpp::complex C_Pp;/**<Wilson coeffients @f$C_P'@f$*/
    
    /**
    * @brief The second fit function from arXiv:hep-ph/0412079v1,\f$ f_2^{LCSR} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return \f$ f_2^{LCSR} \f$
    */
    double LCSR_fit1(double q2, double r_1, double r_2, double m_fit2);
    
    
    /**
    * @brief The third fit function from arXiv:hep-ph/0412079v1, \f$ f_3^{LCSR} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return \f$ f_3^{LCSR} \f$
    */
    double LCSR_fit2(double q2, double r_2, double m_fit2);
    
    
    /**
    * @brief The form factor \f$ f_+ \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_+ \f$
    */
    double f_plus(double q2);
    
    
    /**
    * @brief The form factor \f$ f_T \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_T \f$
    */
    double f_T(double q2);
    
    
    /**
    * @brief The form factor \f$ f_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ f_0 \f$
    */
    double f_0(double q2);
    
    /**
    * @brief The helicity form factor \f$ V_L^0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_L^{\lambda} \f$
    */
    gslpp::complex V_L(double q2);

    
    /**
    * @brief The helicity form factor \f$ V_R^0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_R^{\lambda} \f$
    */
    gslpp::complex V_R(double q2);


    /**
    * @brief The helicity form factor \f$ T_L^0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_L^{\lambda} \f$
    */
    gslpp::complex T_L(double q2);


    /**
    * @brief The helicity form factor \f$ T_R^0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_R^{\lambda} \f$
    */
    gslpp::complex T_R(double q2);


    /**
    * @brief The helicity form factor \f$ S_L \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S_L \f$
    */
    double S_L(double q2);


    /**
    * @brief The helicity form factor \f$ S_R \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S_R \f$
    */
    double S_R(double q2);


    /**
    * @brief The helicity amplitudes normalization factor \f$ N \f$ .
    * @return \f$ N \f$
    */
    gslpp::complex N();
    
    
    /**
    * @brief The \f$ h(q^2,m) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] m mass
    * @return \f$ h(q^2,m) \f$
    */
    gslpp::complex H(double q2, double m);
    
    
    /**
    * @brief The \f$ Y(q^2) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ Y(q^2) \f$
    */
    gslpp::complex Y(double q2);
    
    
    /**
    * @brief The helicity amplitude \f$ H_V^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_V^{\lambda} \f$
    */
    gslpp::complex H_V(double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_A^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_A^{\lambda} \f$
    */
    gslpp::complex H_A(double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_S^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_S^{\lambda} \f$
    */
    gslpp::complex H_S(double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_P^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_P^{\lambda} \f$
    */
    gslpp::complex H_P(double q2, int bar);
    
    
    /**
    * @brief The square of the 3-momentum of the recoiling meson in the M rest frame, \f$ k^2 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ k^2 \f$ 
    */
    double k2 (double q2);
    
    
    /**
    * @brief The factor \f$ \beta \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta \f$
    */
    double beta (double q2);
    
    
    /**
    * @brief The factor \f$ \lambda \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \lambda \f$
    */
    double lambda(double q2);

    
    /**
    * @brief The factor \f$ F \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F \f$
    */
    double F(double q2);
    
    
    /**
    * @brief The angular coefficient \f$ I_{i} \f$ .
    * @param[in] i index of the angular coefficient: 0 for 1c, 2 for 2c, 8 for 6c
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ I_{i} \f$
    */
    double  I(int i, double q2, int bar);
    
    
    /**
    * @brief The CP average \f$ \Sigma_{i} \f$ .
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{i} \f$
    */
    double Sigma(int i, double q2);
    
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{i} \f$ .
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{i} \f$
    */
    double Delta(int i, double q2);
    
    /**
    * @brief The CP average \f$ \Sigma_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{1s} \f$
    */
    double getSigma0(double q2)
    {
        return Sigma(0, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2s} \f$
    */
    double getSigma2(double q2)
    {
        return Sigma(2, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1s} \f$
    */
    double getDelta0(double q2)
    {
        return Delta(0, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2s} \f$
    */
    double getDelta2(double q2)
    {
        return Delta(2, q2);
    };
    
    /**
    * @brief The integral of \f$ \Sigma_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Sigma_{i}> \f$ 
    */
    double integrateSigma(int i, double q_min, double q_max);
    
    /**
    * @brief The integral of \f$ \Delta_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Delta_{i}> \f$ 
    */
    double integrateDelta(int i, double q_min, double q_max);
    

private:
    const StandardModel& mySM;/**< Model type */
    StandardModel::lepton lep;/**< Final leptons type */
    StandardModel::meson meson;/**< Initial meson type */
    StandardModel::meson pseudoscalar;/**< Final pseudoscalar meson type */
    
    unsigned int fplus_updated;/**< Cache variable */
    gslpp::vector<double> fplus_cache;/**< Cache variable */
    
    unsigned int fT_updated;/**< Cache variable */
    gslpp::vector<double> fT_cache;/**< Cache variable */
    
    unsigned int f0_updated;/**< Cache variable */
    double f0_cache;/**< Cache variable */
    
    unsigned int k2_updated;/**< Cache variable */
    gslpp::vector<double> k2_cache;/**< Cache variable */
    
    unsigned int beta_updated;/**< Cache variable */
    double beta_cache;/**< Cache variable */
    
    unsigned int lambda_updated;/**< Cache variable */
    double lambda_cache;/**< Cache variable */
    
    unsigned int F_updated;/**< Cache variable */
    
    unsigned int VL_updated;/**< Cache variable */
    
    unsigned int VR_updated;/**< Cache variable */
    
    unsigned int TL_updated;/**< Cache variable */
    
    unsigned int TR_updated;/**< Cache variable */
    
    unsigned int SL_updated;/**< Cache variable */
    gslpp::vector<double> SL_cache;/**< Cache variable */
    
    unsigned int SR_updated;/**< Cache variable */
    
    unsigned int N_updated;/**< Cache variable */
    gslpp::vector<double> N_cache;/**< Cache variable */
    gslpp::complex Nc_cache;/**< Cache variable */
    
    unsigned int C_1_updated;/**< Cache variable */
    gslpp::complex C_1_cache;/**< Cache variable */

    unsigned int C_2_updated;/**< Cache variable */
    gslpp::complex C_2_cache;/**< Cache variable */
    
    unsigned int C_3_updated;/**< Cache variable */
    gslpp::complex C_3_cache;/**< Cache variable */
    
    unsigned int C_4_updated;/**< Cache variable */
    gslpp::complex C_4_cache;/**< Cache variable */
    
    unsigned int C_5_updated;/**< Cache variable */
    gslpp::complex C_5_cache;/**< Cache variable */
    
    unsigned int C_6_updated;/**< Cache variable */
    gslpp::complex C_6_cache;/**< Cache variable */
    
    unsigned int C_7_updated;/**< Cache variable */
    gslpp::complex C_7_cache;/**< Cache variable */

    unsigned int C_9_updated;/**< Cache variable */
    gslpp::complex C_9_cache;/**< Cache variable */
    
    unsigned int C_10_updated;/**< Cache variable */
    gslpp::complex C_10_cache;/**< Cache variable */
    
    unsigned int C_7p_updated;/**< Cache variable */
    gslpp::complex C_7p_cache;/**< Cache variable */
    
    unsigned int C_9p_updated;/**< Cache variable */
    gslpp::complex C_9p_cache;/**< Cache variable */
    
    unsigned int C_10p_updated;/**< Cache variable */
    gslpp::complex C_10p_cache;/**< Cache variable */
    
    unsigned int C_S_updated;/**< Cache variable */
    gslpp::complex C_S_cache;/**< Cache variable */
    
    unsigned int C_P_updated;/**< Cache variable */
    gslpp::complex C_P_cache;/**< Cache variable */
    
    unsigned int C_Sp_updated;/**< Cache variable */
    gslpp::complex C_Sp_cache;/**< Cache variable */
    
    unsigned int C_Pp_updated;/**< Cache variable */
    gslpp::complex C_Pp_cache;/**< Cache variable */
    
    unsigned int Yupdated;/**< Cache variable */
    gslpp::vector<double> Ycache;/**< Cache variable */
    
    unsigned int H_V0updated;/**< Cache variable */
    gslpp::vector<double> H_V0cache;/**< Cache variable */
    gslpp::complex H_V0Ccache[2];/**< Cache variable */
    
    unsigned int H_A0updated;/**< Cache variable */
    
    unsigned int H_Supdated;/**< Cache variable */
    gslpp::vector<double> H_Scache;/**< Cache variable */
    
    unsigned int H_P_updated;/**< Cache variable */
    gslpp::vector<double> H_P_cache;/**< Cache variable */
    
    unsigned int I0_updated;/**< Cache variable */
    unsigned int I2_updated;/**< Cache variable */
    unsigned int I8_updated;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta2Cached;/**< Cache variable */
    
    double avaSigma0;/**< Gsl integral variable */
    double avaSigma2;/**< Gsl integral variable */
    
    double errSigma0;/**< Gsl integral variable */
    double errSigma2;/**< Gsl integral variable */
    
    double avaDelta0;/**< Gsl integral variable */
    double avaDelta2;/**< Gsl integral variable */
    
    double errDelta0;/**< Gsl integral variable */
    double errDelta2;/**< Gsl integral variable */
    
    gsl_function FS0;/**< Gsl integral variable */
    gsl_function FS2;/**< Gsl integral variable */
    
    gsl_function FD0;/**< Gsl integral variable */
    gsl_function FD2;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_sigma0;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma2;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_delta0;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta2;/**< Gsl integral variable */
    
    std::map<std::pair<double, double>, double > cacheSigma0;/**< Gsl integral variable */
    std::map<std::pair<double, double>, double > cacheSigma2;/**< Gsl integral variable */
    
    std::map<std::pair<double, double>, double > cacheDelta0;/**< Gsl integral variable */
    std::map<std::pair<double, double>, double > cacheDelta2;/**< Gsl integral variable */
    
};

#endif	/* MPLL_H */

    