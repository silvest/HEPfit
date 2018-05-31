/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVLNU_H
#define MVLNU_H

class StandardModel;
#include <gsl/gsl_integration.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResultPtr.h>
#include <gsl/gsl_spline.h>
#include <memory>

class MVlnu {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i);
    
     /**
     * @brief Destructor.
     */
    virtual ~MVlnu();
    std::vector<std::string> initializeMVlnuParameters();
    
    /**
    * @brief The integral of \f$ d\Gamma/dw \f$ from \f$w_{min}\f$ to \f$w_{max}\f$
    * @param[in] w_min maximum q^2 of the integral
    * @param[in] w_max minimum q^2 of the integral
    * @return \f$ <d\Gamma/dw> \f$ 
    */
    double getDeltaGammaDeltaw(double w_min, double w_max);
    
    /**
    * @brief The integral of \f$ d\Gamma/dcos\theta_{l} \f$ from \f$cos\theta_{l,min}\f$ to \f$cos\theta_{l,max}\f$
    * @param[in] cl_min minmum cos angle value of the integral
    * @param[in] cl_max maximum cos angle value of the integral
    * @return \f$ <d\Gamma/dcos\theta_{l}> \f$ 
    */
    double getDeltaGammaDeltacl(double cl_min, double cl_max);
    
    /**
    * @brief The integral of \f$ d\Gamma/dcos\theta_{V} \f$ from \f$cos\theta_{V,min}\f$ to \f$cos\theta_{V,max}\f$
    * @param[in] cV_min minmum cos angle value of the integral
    * @param[in] cV_max maximum cos angle value of the integral
    * @return \f$ <d\Gamma/dcos\theta_{V}> \f$ 
    */
    double getDeltaGammaDeltacV(double cV_min, double cV_max);
    
    /**
    * @brief The integral of \f$ d\Gamma/dcos\chi \f$ from \f$cos\chi\f$ to \f$cos\chi\f$
    * @param[in] chi_min minmum cos angle value of the integral
    * @param[in] chi_max maximum cos angle value of the integral
    * @return \f$ <d\Gamma/dcos\chi> \f$ 
    */
    double getDeltaGammaDeltachi(double chi_min, double chi_max);

    /**
    * @brief The width of the meson M
    * @return \f$ \Gamma_M \f$ 
    */
    double getMwidth(){
        updateParameters();
        return width;
    }
    
private:
    const StandardModel& mySM;/**< Model type */
    QCD::lepton lep;/**< Final leptons type */
    QCD::meson meson;/**< Initial meson type */
    QCD::meson vectorM;/**< Final vector meson type */
    std::vector<std::string> mvlnuParameters;/**< The string of mandatory MVlnu parameters */
    bool CLNflag;
    
    double GF;            /**<Fermi constant */
    double Mlep;          /**<Lepton mass */
    double Mnu;           /**<Neutrino mass */
    double MM;            /**<Initial meson mass */
    double MV;            /**<Final vector meson mass */
    double w0;            /**<Kinematic variable w at q2=0 */
    double RV;            /**<Dimensionless meson - vector mass ratio */
    double mu_b;          /**<b mass scale */
    double Mb;            /**<b quark mass */
    double mb_pole;       /**<b quark pole mass */
    double Mc;            /**<charm quark mass */
    double mc_pole;       /**<charm quark pole mass */
    double width;         /**<Initial meson width */
    double ale_mub;   /**<@f\aplha_{em}(\mu_b)$@f$ */
    gslpp::complex Vcb;   /**<CKM factor of the decay*/
    double amplsq_factor;   /**< Overall helicity |A|^2 factor*/
    
    double CV_SM; /**<Wilson coeffients @f$C_{V}@f$*/
    double CA_SM; /**<Wilson coeffients @f$C_{A}@f$*/
    
    gslpp::complex CS; /**<Wilson coeffients @f$C_{S}@f$*/
    gslpp::complex CP; /**<Wilson coeffients @f$C_{P}@f$*/
    gslpp::complex CSp; /**<Wilson coeffients @f$C_{S}'@f$*/
    gslpp::complex CPp; /**<Wilson coeffients @f$C_{P}'@f$*/
    gslpp::complex CV; /**<Wilson coeffients @f$C_{V}@f$*/
    gslpp::complex CA; /**<Wilson coeffients @f$C_{A}@f$*/
    gslpp::complex CVp; /**<Wilson coeffients @f$C_{V}'@f$*/
    gslpp::complex CAp; /**<Wilson coeffients @f$C_{A}'@f$*/
    gslpp::complex C7; /**<Wilson coeffients @f$C_{7}@f$*/
    gslpp::complex C7p; /**<Wilson coeffients @f$C_{7}'@f$*/
    gslpp::complex CT; /**<Wilson coeffients @f$C_{T}@f$*/
    gslpp::complex CTp; /**<Wilson coeffients @f$C_{Tp}@f$*/
    
    double hA1w1,rho2,R1w1,R2w1; /**<CLN form factor parameters*/
    double af0,af1,af2,ag0,ag1,ag2,aF11,aF12; /**<BGL form factor parameters*/
    double mBcstV1,mBcstV2,mBcstV3,mBcstV4,mBcstA1,mBcstA2,mBcstA3,mBcstA4,chiTV,chiTA,nI; /**<BGL form factor parameters*/
    double zV1,zV2,zV3,zV4,zA1,zA2,zA3,zA4;
    
    /**
     * @brief The update parameter method for MVll.
     */
    void updateParameters();

    /**
    * @brief kinematic function \f$ \lambda_{1/2} \f$.
    * @param[in] Mass \f$MM^2\f$ of the decay
    * @param[in] Mass \f$MV^2\f$ of the decay
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \sqrt{\lambda} \f$
    */   
    double lambda_half(double a, double b, double c);
    
    /**
    * @brief BGL outer function \f$ \phi_f \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_f \f$
    */
    double phi_f(double z);

    /**
    * @brief BGL form factor function \f$ f \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ f \f$
    */
    double f_BGL(double q2);
    
    /**
    * @brief BGL outer function \f$ \phi_g \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_g \f$
    */
    double phi_g(double q2);

    /**
    * @brief BGL form factor function \f$ g \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ g \f$
    */
    double g_BGL(double q2);

    /**
    * @brief BGL outer function \f$ \phi_F1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \phi_F1 \f$
    */
    double phi_F1(double q2);
    
    /**
    * @brief BGL form factor function \f$ F_{1} \f$.
    * @param[in] z \f$z\f$ of the decay
    * @return \f$ F_{1} \f$
    */
    double F1_BGL(double q2);
    
    /**
    * @brief HQET form factor \f$ h_{A1} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_{A1} \f$
    */
    double hA1(double q2);
    
    /**
    * @brief Ratio of HQET form factors \f$ R_{1} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{1} \f$
    */
    double R1(double q2);
    
    /**
    * @brief Ratio of HQET form factors \f$ R_{2} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{2} \f$
    */
    double R2(double q2);
    
    /**
    * @brief Ratio of HQET form factors \f$ R_{0} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ R_{0} \f$
    */
    double R0(double q2);
    
    /**
    * @brief The transverse form factor \f$ V \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V \f$
    */
    double V(double q2);

    /**
    * @brief The transverse form factor \f$ A_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_0 \f$
    */
    double A0(double q2);

    
    /**
    * @brief The transverse form factor \f$ A_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_1 \f$
    */
    double A1(double q2);
    
    /**
    * @brief The transverse form factor \f$ A_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_2 \f$
    */
    double A2(double q2);
    
    /**
    * @brief The transverse form factor \f$ T_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_1 \f$
    */
    double T1(double q2);

    
    /**
    * @brief The transverse form factor \f$ T_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_2 \f$
    */
    double T2(double q2);
    
    /**
    * @brief The helicity form factor \f$ A_{12} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_{12} \f$
    */
    double A12(double q2); 
    
    /**
    * @brief The helicity form factor \f$ T_{23} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_{23} \f$
    */
    double T23(double q2);
 
    /**
    * @brief The helicity amplitude \f$ H_{V0} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{V0} \f$
    */
    gslpp::complex HV0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_{Vp} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Vp} \f$
    */
    gslpp::complex HVp(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Vm} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Vm} \f$
    */
    gslpp::complex HVm(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{A0} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{A0} \f$
    */
    gslpp::complex HA0(double q2);

    /**
    * @brief The helicity amplitude \f$ H_{Ap} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Ap} \f$
    */
    gslpp::complex HAp(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Am} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Am} \f$
    */
    gslpp::complex HAm(double q2);
    
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
    * @brief The helicity amplitude \f$ H_{T0} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{T0} \f$
    */
    gslpp::complex HT0(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{T0t} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{T0t} \f$
    */
    gslpp::complex HT0t(double q2);

    /**
    * @brief The helicity amplitude \f$ H_{Tp} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Tp} \f$
    */
    gslpp::complex HTp(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Tpt} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Tpt} \f$
    */
    gslpp::complex HTpt(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Tm} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Tm} \f$
    */
    gslpp::complex HTm(double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_{Tmt} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_{Tmt} \f$
    */
    gslpp::complex HTmt(double q2);  
 
    /**
    * @brief The helicity amplitude \f$ G_{000} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{000} \f$
    */
    gslpp::complex G000(double q2);

    /**
    * @brief The helicity amplitude \f$ G_{010} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{010} \f$
    */    
    gslpp::complex G010(double q2); 
    
    /**
    * @brief The helicity amplitude \f$ G_{020} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{020} \f$
    */    
    gslpp::complex G020(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{200} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{200} \f$
    */    
    gslpp::complex G200(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{210} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{210} \f$
    */    
    gslpp::complex G210(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{220} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{220} \f$
    */    
    gslpp::complex G220(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{211} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{221} \f$
    */    
    gslpp::complex G211(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{221} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{221} \f$
    */    
    gslpp::complex G221(double q2);
    
    /**
    * @brief The helicity amplitude \f$ G_{222} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ G_{222} \f$
    */       
    gslpp::complex G222(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{1s} \f$
    */
    double  J1s(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{1c} \f$
    */
    double  J1c(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{2s} \f$
    */
    double  J2s(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{2c} \f$
    */
    double  J2c(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{3} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{3} \f$
    */
    double  J3(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{4} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{4} \f$
    */
    double  J4(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{5} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{5} \f$
    */
    double  J5(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{6s} \f$
    */
    double  J6s(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{6c} \f$
    */
    double  J6c(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{7} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{7} \f$
    */
    double  J7(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{8} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{8} \f$
    */
    double  J8(double q2);
    
    /**
    * @brief The angular coefficient \f$ J_{9} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ J_{9} \f$
    */
    double  J9(double q2);

    /**
    * @brief \f$ <J_{i}> \f$ 
    * @param[in] i, angular coefficient index (i = 1,...,12)
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

    /**
     * @brief \f$ d\Gamma/dcos \theta_{l}/dq2 \f$ 
     * @param[in] q2 \f$q^2\f$ of the decay
     * @param[in] cl \f$ cos \theta_{l}\f$ of the decay
     * @return \f$ d\Gamma/(dcos \theta_{l} dq2) \f$ 
     */    
    double dGammadcldq2(double q2, double cl);

    /**
     * @brief \f$ d\Gamma/dcos \theta_{l} \f$ 
     * @param[in] cl \f$ cos \theta_{l}\f$ of the decay
     * @return \f$ d\Gamma/dcos \theta_{l} \f$ 
     */
    double dGammadcl(double cl);

    /**
     * @brief \f$ d\Gamma/dcos \theta_{V}/dq2 \f$ 
     * @param[in] q2 \f$q^2\f$ of the decay
     * @param[in] cV \f$ cos \theta_{V}\f$ of the decay
     * @return \f$ d\Gamma/(dcos \theta_{V} dq2) \f$ 
     */    
    double dGammadcVdq2(double q2, double cV);
    
    /**
     * @brief \f$ d\Gamma/dcos \theta_{V} \f$ 
     * @param[in] cV \f$ cos \theta_{V}\f$ of the decay
     * @return \f$ d\Gamma/dcos \theta_{V} \f$ 
     */
    double dGammadcV(double cV);
    
    /**
     * @brief \f$ d\Gamma/dcos \chi/dq2 \f$ 
     * @param[in] q2 \f$q^2\f$ of the decay
     * @param[in] chi \f$ chi\f$ of the decay
     * @return \f$ d\Gamma/( dq2 dcos \chi) \f$ 
     */    
    double dGammadchidq2(double q2, double chi);
    
    /**
     * @brief \f$ d\Gamma/dcos \chi \f$ 
     * @param[in] q2 \f$q^2\f$ of the decay
     * @param[in] cV \f$ cos \chi\f$ of the decay
     * @return \f$ d\Gamma/dcos \chi \f$ 
     */
    double dGammadchi(double chi);
    
    gsl_error_handler_t * old_handler; /**< GSL error handler store */
    gsl_function FJ;/**< GSL integral variable */
    double J_res;/**< GSL integral variable */
    double J_err;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_J;/**< GSL integral variable */
 
};

#endif /* MVLNU_H */

