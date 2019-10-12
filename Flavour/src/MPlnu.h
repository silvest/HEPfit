/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLNU_H
#define MPLNU_H

class StandardModel;
#include <gsl/gsl_integration.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResultPtr.h>
#include <gsl/gsl_spline.h>
#include <memory>

#define NBGL 3 /* ONLY 3 or 2*/

class MPlnu {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i);
    
     /**
     * @brief Destructor.
     */
    virtual ~MPlnu();
    std::vector<std::string> initializeMPlnuParameters();
    
    /**
    * @brief The integral of \f$ d\Gamma/dw \f$ from \f$w_{min}\f$ to \f$w_{max}\f$
    * @param[in] w_min maximum q^2 of the integral
    * @param[in] w_max minimum q^2 of the integral
    * @return \f$ <d\Gamma/dw> \f$ 
    */
    double getDeltaGammaDeltaw(double w_min, double w_max);

    /**
    * @brief Weak Unitarity constraint for BGL parameters related to 1- resonances
    * @return \f$ \Sum_i (af+_i^2) \f$ 
    */
    double get_unitarity_1min_BGL();
 
    /**
    * @brief Weak Unitarity constraint for BGL parameters related to 0+ resonances
    * @return \f$ \Sum_i (af0_i^2) \f$ 
    */
    double get_unitarity_0plus_BGL();
    
    /**
    * @brief Strong Unitarity constraint for BGL parameters using HQET
    * @return \f$ \Sum_i (af0_i^2) \f$ 
    */
    double get_strong_unitarity_BGL();

    /**
     * @brief return fplus form factor at \f$ q^2 \f$
     * @return \f$ f_{+}(q^2) \f$ 
     */
    double get_fplus(double q2);
    
    /**
    * @brief return f0 form factor at \f$ q^2 \f$
    * @return \f$ f_{0}(q^2) \f$ 
    */
    double get_f0(double q2);
    
    /**
     * @brief return fT form factor at \f$ q^2 \f$
     * @return \f$ q^2 \f$ 
     */
    double get_fT(double q2);

    /**
    * @brief The width of the meson M
    * @return \f$ \Gamma_M \f$ 
    */
    double getMwidth()
    {
        updateParameters();
        return width;
    }
    
    /**
    * @brief The BGL parameter \f$ a_0^{f_0}\f$
    * @return \f$ a_0^{f_0}\f$
    */
    double getaf0_0()
    {
        updateParameters();
        return af0_0;
    }
    
private:
    const StandardModel& mySM;/**< Model type */
    QCD::lepton lep;/**< Final leptons type */
    QCD::meson meson;/**< Initial meson type */
    QCD::meson pseudoscalarM;/**< Final vector meson type */
    std::vector<std::string> mplnuParameters;/**< The string of mandatory MPlnu parameters */
    bool CLNflag; /**< A flag for switching to CLN parameterization */
    bool btocNPpmflag; /**< A flag for switching to the +/- basis for NP Wilson coefficients */
    bool NPanalysis; /**< A flag to switch to BSM analysis */
    
    double GF;            /**<Fermi constant */
    double Mlep;          /**<Lepton mass */
    double Mnu;           /**<Neutrino mass */
    double MM;            /**<Initial meson mass */
    double MP;            /**<Final pseudoscalar meson mass */
    double w0;            /**<Kinematic variable w at q2=0 */
    double z0;            /**<Kinematic variable z at q2=0 */
    double RV;            /**<Dimensionless meson - vector mass ratio */
    double mu_b;          /**<b mass scale */
    double Mb;            /**<b quark mass */
    double Mc;            /**<charm quark mass */
    double width;         /**<Initial meson width */
    double ale_mub;   /**<@f\aplha_{em}(\mu_b)$@f$ */
    gslpp::complex Vcb;   /**<CKM factor of the decay*/
    double amplsq_factor;   /**< Overall helicity |A|^2 factor*/
    double q2min, q2max; /**< min and max lepton-neutrino invariant mass squared*/
    
    double eta_EW; /**<EW correction @f$\eta_{EW}@f$*/
    double CV_SM; /**<Wilson coeffients @f$C_{V}@f$*/
    
    double CS; /**<Wilson coeffients @f$C_{S}@f$*/
    double CP; /**<Wilson coeffients @f$C_{P}@f$*/
    double CSp; /**<Wilson coeffients @f$C_{S}'@f$*/
    double CPp; /**<Wilson coeffients @f$C_{P}'@f$*/
    double CV; /**<Wilson coeffients @f$C_{V}@f$*/
    double CA; /**<Wilson coeffients @f$C_{A}@f$*/
    double CVp; /**<Wilson coeffients @f$C_{V}'@f$*/
    double CAp; /**<Wilson coeffients @f$C_{A}'@f$*/
    double C7; /**<Wilson coeffients @f$C_{7}@f$*/
    double C7p; /**<Wilson coeffients @f$C_{7}'@f$*/
    double CT; /**<Wilson coeffients @f$C_{T}@f$*/
    double CTp; /**<Wilson coeffients @f$C_{Tp}@f$*/
    
    double fplusz0,rho1to2; /**<CLN form factor parameters*/
    double N_0, alpha_0, alpha_p, beta_0, beta_p, gamma_0, gamma_p; /**<CLN form factor parameter modification*/
    double af0_0,af0_1,af0_2,afplus_0,afplus_1,afplus_2; /**<BGL form factor parameters*/
#if NBGL == 3
    double af0_3, afplus_3;
#endif    
    double mBc1m_1,mBc1m_2,mBc1m_3,mBc1m_4,mBc0p_1,mBc0p_2,chitildeT,chiL,nI; /**<BGL form factor parameters*/
    double z1m_1,z1m_2,z1m_3,z0p_1,z0p_2;
    gslpp::complex z1m_4;
    double cached_intJ1_tau, cached_intJ2_tau, cached_intJ3_tau,
                cached_intJ1_mu, cached_intJ2_mu, cached_intJ3_mu,
                cached_intJ1_el, cached_intJ2_el, cached_intJ3_el; /**< caching Js integral btw q2min and q2mx*/
    double fplusz0_cache,rho1to2_cache;
    double N_0_cache, alpha_0_cache, alpha_p_cache, beta_0_cache, beta_p_cache, gamma_0_cache, gamma_p_cache;
    double af0_1_cache,af0_2_cache,afplus_0_cache,afplus_1_cache,afplus_2_cache;
#if NBGL == 3
    double af0_3_cache,afplus_3_cache;
#endif    
    double CS_cache,CSp_cache,CP_cache,CPp_cache,CV_cache,CVp_cache,CA_cache,CAp_cache,CT_cache,CTp_cache;
    bool checkcache_int_tau, checkcache_int_mu, checkcache_int_el;
    
    /**
     * @brief The update parameter method for MPlnu.
     */
    void updateParameters();

    /**
    * @brief kinematic function \f$ \lambda_{1/2} \f$.
    * @param[in] Mass \f$MM^2\f$ of the decay
    * @param[in] Mass \f$MP^2\f$ of the decay
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \sqrt{\lambda} \f$
    */   
    double lambda_half(double a, double b, double c);
    
    /**
    * @brief BGL outer function \f$ \phi_g \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_g \f$
    */
    double phi_fplus(double q2);

    /**
    * @brief Form factor function \f$ f_{+} \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ f_{+} \f$
    */
    double fplus(double q2);
    
    /**
    * @brief BGL outer function \f$ \phi_f_{0} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_f \f$
    */
    double phi_f0(double z);

    /**
    * @brief Form factor function \f$ f_{0} \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ f_{0} \f$
    */
    double f0(double q2);
    
    /**
    * @brief Form factor function \f$ f_{T} \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ f_{T} \f$
    */
    double fT(double q2);
 
    /**
    * @brief The helicity amplitude \f$ H_{V} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{V} \f$
    */
    gslpp::complex HV(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{A} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{A} \f$
    */
    gslpp::complex HA(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_P \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_P \f$
    */
    gslpp::complex HP(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_S \f$
    */
    gslpp::complex HS(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{T} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{T} \f$
    */
    gslpp::complex HT(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Tt} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Tt} \f$
    */
    gslpp::complex HTt(double q2);  
 
    /**
    * @brief The helicity amplitude \f$ G_{0} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{0} \f$
    */
    gslpp::complex G0(double q2);

    /**
    * @brief The helicity amplitude \f$ G_{1} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{1} \f$
    */    
    gslpp::complex G1(double q2); 
    
    /**
    * @brief The helicity amplitude \f$ G_{2} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{2} \f$
    */    
    gslpp::complex G2(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{1} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{1} \f$
    */
    double  J1(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{2} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{2} \f$
    */
    double  J2(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{3} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{3} \f$
    */
    double  J3(double q2);
    
    /**
    * @brief Differential Squared Amplitude \f$ d\Gamma/dq2 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ d\Gamma/dq2 \f$
    */
    double  dGammadq2(double q2);

    /**
    * @brief \f$ <J_{i}> \f$ 
    * @param[in] i, angular coefficient index (i = 1,2,3)
    * @param[in] q2_min, lower extreme \f$q^2\f$ of intgrated decay
    * @param[in] q2_max, upper extreme \f$q^2\f$ of intgrated decay
    * @return \f$ <J_{i}> \f$ 
    */
    double integrateJ(int i, double q2_min, double q2_max) ;
    
    /**
     * @brief \f$ d\Gamma/dw \f$ 
     * @param[in] w related to \f$q^2\f$ of the decay
     * @return \f$ d\Gamma/dw \f$ 
     */    
    double dGammadw(double w);
    
    gsl_error_handler_t * old_handler; /**< GSL error handler store */
    gsl_function FJ;/**< GSL integral variable */
    double J_res;/**< GSL integral variable */
    double J_err;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_J;/**< GSL integral variable */
 
};

#endif /* MPLNU_H */

