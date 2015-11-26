/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMA_H
#define	BSGAMMA_H

#include <ThObservable.h>
#include <gsl/gsl_integration.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>

#define FOUR_BODY false

/**
 * @class Bsgamma
 * @ingroup Flavour
 * @brief A class for the @f$b \to s \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$b \to s \gamma@f$ decay, following
 * the prescriptions of @cite Misiak:2006zs and @cite Misiak:2015xwa. After the
 * Wilson coefficients are computed in computeCoeff() and the cache is checked in 
 * checkCache(), the parameters are updated in updateParameters(), in which part
 * of the parameters are taken from the fit in @cite Alberti:2014yda and used to compute
 * the ratio \f$C = | \frac{V_{ub}}{V_{cb}} |^2 \frac{\Gamma[\bar{B} \to X_c e \bar{\nu}]}
 * {\Gamma[\bar{B} \to X_u e \bar{\nu}]} \f$ in C_sem().
 * 
 * The perturbative part of the Branching Ratio is computed order by order:
 * 
 * @li at Leading Order it is computed in P0(), in which are
 * taken into account both the leading term due to @f$C_7@f$ @cite Misiak:2006ab and the subleading term
 * due to the 4-body contribution, computed in P0_4body() @cite Kaminski:2012eb ;
 * 
 * @li at Next to Leading Order it is computed in P11() and P21(),
 * where the latter is build from the Kij_1() function @cite Misiak:2006ab, which make use of the functions
 * Phi11_1(), Phi12_1(), Phi13_1(), Phi14_1(), Phi15_1(), Phi16_1(), Phi17_1(), Phi18_1(), 
 * Phi22_1(), Phi23_1(), Phi24_1(), Phi25_1(), Phi26_1(), Phi27_1(), Phi28_1(), Phi33_1(), 
 * Phi34_1(), Phi35_1(), Phi36_1(), Phi37_1(), Phi38_1(), Phi44_1(), Phi45_1(), Phi46_1(), 
 * Phi47_1(), Phi48_1(), Phi55_1(), Phi56_1(), Phi57_1(), Phi58_1(), Phi66_1(), Phi67_1(), 
 * Phi68_1(), Phi77_1(), Phi78_1() and Phi88_1() @cite Buras:2002tp, @cite Gambino:2001ew, 
 * @cite Pott:1995if . The subleading terms due to the
 * 4-body contributions @cite Huber:2014nna are currently hard-coded in Phi23_1_4body(), Phi24_1_4body(), 
 * Phi25_1_4body() and Phi26_1_4body(), and switched off due to setting the marco
 * FOUR_BODY to false.
 * 
 * The @f$V_{ub}@f$ corrections at LO are automatically taken into account in P0_4body() @cite Kaminski:2012eb,
 * while at NLO are computed in the function Vub_NLO() @cite Gambino:2001ew, @cite Huber:2014nna , which considers contributions
 * from 2-body, 3-body and 4-body decays, with the former switched off due to setting the marco
 * FOUR_BODY to false.
 * 
 * All the perturbative corrections are eventually added in the function P(). The 
 * non-perturbative corrections are computed in the function N(), which follows 
 * the prescription of @cite Benzke:2010js. The observables are finally computed in
 * the computeThValue() function.
 */
class Bsgamma : public ThObservable {
public: 
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    * @param[in] obsFlag flag to choose which observable to compute
    */
    Bsgamma(const StandardModel& SM_i, int obsFlag);
    
    
    /**
    * @brief The cutoff energy function \f$ \delta = 1 - \frac{2 E_0}{M_b^{\rm kin}} \f$.
    * @param[in] E0 cutoff energy
    * @return \f$ \delta(E0) \f$ 
    */
    double delta(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ \rho \f$ as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @return \f$ \rho(E0) \f$ 
    */
    double rho(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ \omega \f$ as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @return \f$ \omega(E0) \f$ 
    */
    double omega(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ T_1 \f$ as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_1(E0) \f$ 
    */
    double T1(double E0, double t);
    
    
    /**
    * @brief The cutoff energy function \f$ T_2 \f$ as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_2(E0) \f$ 
    */
    double T2(double E0, double t);
    
    
    /**
    * @brief The cutoff energy function \f$ T_3 \f$ as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_3(E0) \f$ 
    */
    double T3(double E0, double t);
    
    
    /**
    * @brief The 4-body LO contribution as defined in @cite Kaminski:2012eb .
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ P_{tree}^{(0)} \f$ 
    */
    double P0_4body(double E0, double t);
    
    
    /**
    * @brief The squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$, \f$ z \f$.
    * @return \f$ z \f$
    */
    double zeta();
    
    
    /**
    * @brief The funcion \f$ a(z) \f$ as defined in @cite Buras:2002tp .
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ a(z) \f$ 
    */
    gslpp::complex a(double z);
    
    
    /**
    * @brief The funcion \f$ b(z) \f$ as defined in @cite Buras:2002tp .
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ b(z) \f$ 
    */
    gslpp::complex b(double z);
    
    
    /**
    * @brief The funcion \f$ r_i^{(1)}(z) \f$ as defined in @cite Buras:2002tp .
    * @param[in] i function index
    * @param[in] z squared ratio between m_c and m_b^{\rm kin}
    * @return \f$ r_i(z)^{(1)} \f$
    */
    gslpp::complex r1(int i, double z);
    
    
    /**
    * @brief The function \f$ \Gamma \f$ as defined in @cite Gambino:2001ew .
    * @param[in] t dummy variable to be integrated out
    * @return \f$ \Gamma \f$ 
    */
    gslpp::complex Gamma_t(double t);
    
    
    /**
    * @brief The function \f$ k \f$ as defined in @cite Pott:1995if .
    * @param[in] Mq quark mass
    * @param[in] t dummy variable to be integrated out
    * @return \f$ k \f$ 
    */
    gslpp::complex kappa(double Mq, double t);
    
    
    /**
    * @brief The function \f$|k_c(t)|^2 t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_c(t)|^2t\f$
    */
    double getKc_abs2_t(double t)
    {
        return kappa(Mc,t).abs2() * t;
    };
    
    
    /**
    * @brief The function \f$|k_c(t)|^2(1 - t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_c(t)|^2(1 - t)\f$
    */
    double getKc_abs2_1mt(double t)
    {
        return kappa(Mc,t).abs2() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_c(t)|^2t(1 - t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_c(t)|^2t(1 - t)\f$
    */
    double getKc_abs2_t_1mt(double t)
    {
        return kappa(Mc,t).abs2() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$t|k_c(t)|^2(1 - t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_c(t)|^2(1 - t)^2\f$
    */
    double getKc_abs2_1mt2(double t)
    {
        return kappa(Mc,t).abs2() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_c(t))t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_c(t))t\f$ 
    */
    double getKc_re_t(double t)
    {
        return kappa(Mc,t).real() * t ;
    };
    
    
    /**
    * @brief The function \f$Re(k_c(t))t(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_c(t))t(1-t)\f$
    */
    double getKc_re_t_1mt(double t)
    {
        return kappa(Mc,t).real() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_c(t))t(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_c(t))t(1-t)^2\f$
    */
    double getKc_re_t_1mt2(double t)
    {
        return kappa(Mc,t).real() * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_c(t))(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_c(t))(1-t)\f$ 
    */
    double getKc_re_1mt(double t)
    {
        return kappa(Mc,t).real() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Im(k_c(t))(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Im(k_c(t))(1-t)\f$ 
    */
    double getKc_im_1mt(double t)
    {
        return kappa(Mc,t).imag() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_c(t))(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_c(t))(1-t)^2\f$
    */
    double getKc_re_1mt2(double t)
    {
        return kappa(Mc,t).real() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Im(k_c(t))(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Im(k_c(t))(1-t)^2\f$
    */
    double getKc_im_1mt2(double t)
    {
        return kappa(Mc,t).imag() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2(1 - t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2(1 - t)\f$
    */
    double getKb_abs2_1mt(double t)
    {
        return kappa(Mb_kin,t).abs2() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2(1 - t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2(1 - t)^2\f$
    */
    double getKb_abs2_1mt2(double t)
    {
        return kappa(Mb_kin,t).abs2() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2t(1 - t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2t(1 - t)\f$
    */
    double getKb_abs2_t_1mt(double t)
    {
        return kappa(Mb_kin,t).abs2() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2t(1 - t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2t(1 - t)^2\f$
    */
    double getKb_abs2_t_1mt2(double t)
    {
        return kappa(Mb_kin,t).abs2() * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2t^2(1 - t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2t^2(1 - t)\f$
    */
    double getKb_abs2_t2_1mt(double t)
    {
        return kappa(Mb_kin,t).abs2() * t * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$|k_b(t)|^2t^2(1 - t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_b(t)|^2t^2(1 - t)^2\f$
    */
    double getKb_abs2_t2_1mt2(double t)
    {
        return kappa(Mb_kin,t).abs2() * t * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))t\f$ 
    */
    double getKb_re_t(double t)
    {
        return kappa(Mb_kin,t).real() * t ;
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))t(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))t(1-t)\f$
    */
    double getKb_re_t_1mt(double t)
    {
        return kappa(Mb_kin,t).real() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))t^2(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))t^2(1-t)\f$
    */
    double getKb_re_t2_1mt(double t)
    {
        return kappa(Mb_kin,t).real() * t * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))t^2(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))t^2(1-t)^2\f$
    */
    double getKb_re_t2_1mt2(double t)
    {
        return kappa(Mb_kin,t).real() * t * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))t(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))t(1-t)^2\f$
    */
    double getKb_re_t_1mt2(double t)
    {
        return kappa(Mb_kin,t).real() * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))(1-t)\f$ 
    */
    double getKb_re_1mt(double t)
    {
        return kappa(Mb_kin,t).real() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))(1-t)^2\f$
    */
    double getKb_re_1mt2(double t)
    {
        return kappa(Mb_kin,t).real() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))Re(k_c(t))(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))Re(k_c(t))(1-t)\f$
    */
    double getKc_re_Kb_1mt(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))Re(k_c(t))(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))Re(k_c(t))(1-t)^2\f$
    */
    double getKc_re_Kb_1mt2(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))Re(k_c(t)t(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))Re(k_c(t)t(1-t)\f$
    */
    double getKc_re_Kb_t_1mt(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))Re(k_c(t)t(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))Re(k_c(t)t(1-t)^2\f$
    */
    double getKc_re_Kb_t_1mt2(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief Integral of the functions getKb_re_1mt() and getKb_re_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))(1-t) + \int_{1-\delta(E_0)}^1 Re(k_b(t))(1-t)^2\f$
    */
    double Int_b1(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_re_t_1mt() and getKb_re_t_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))t(1-t) + \int_{1-\delta(E_0)}^1 Re(k_b(t))t(1-t)^2\f$
    */
    double Int_b2(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_re_t() and getKb_re_t_1mt().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))t + \int_{1-\delta(E_0)}^1 Re(k_b(t))t(1-t)\f$
    */
    double Int_b3(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_re_t2_1mt() and getKb_re_t2_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))t^2(1-t) + \int_{1-\delta(E_0)}^1 Re(k_b(t))t^2(1-t)^2\f$
    */
    double Int_b4(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_abs2_1mt() and getKb_abs2_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |(k_b(t)|^2(1-t) + \int_{1-\delta(E_0)}^1 |(k_b(t)|^2(1-t)^2\f$
    */
    double Int_bb1(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_abs2_t_1mt() and getKb_abs2_t_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |(k_b(t)|^2t(1-t) + \int_{1-\delta(E_0)}^1 |(k_b(t)|^2t(1-t)^2\f$
    */
    double Int_bb2(double E0);
    
    
    /**
    * @brief Integral of the functions getKb_abs2_t2_1mt() and getKb_abs2_t2_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |(k_b(t)|^2t^2(1-t) + \int_{1-\delta(E_0)}^1 |(k_b(t)|^2t^2(1-t)^2\f$
    */
    double Int_bb4(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_re_Kb_1mt() and getKc_re_Kb_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))Re(k_c(t))(1-t) + \int_{1-\delta(E_0)}^1 Re(k_b(t))Re(k_c(t))(1-t)^2\f$
    */
    double Int_bc1(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_re_Kb_t_1mt() and getKc_re_Kb_t_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_b(t))Re(k_c(t))t(1-t) + \int_{1-\delta(E_0)}^1 Re(k_b(t))Re(k_c(t))t(1-t)^2\f$
    */
    double Int_bc2(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_re_1mt() and getKc_re_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_c(t))(1-t) + \int_{1-\delta(E_0)}^1 Re(k_c(t))(1-t)^2\f$
    */
    double Int_c1(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_im_1mt() and getKc_im_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Im(k_c(t))(1-t) + \int_{1-\delta(E_0)}^1 Im(k_c(t))(1-t)^2\f$
    */
    double Int_c1_im(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_re_t_1mt() and getKc_re_t_1mt2().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_c(t))t(1-t) + \int_{1-\delta(E_0)}^1 Re(k_c(t))t(1-t)^2\f$
    */
    double Int_c2(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_re_t() and getKc_re_t_1mt().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} Re(k_c(t))t + \int_{1-\delta(E_0)}^1 Re(k_c(t))t(1-t)\f$
    */
    double Int_c3(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_abs2_t() and getKc_abs2_t_1mt().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |k_c(t)|^2t + 2\int_{1-\delta(E_0)}^1 |k_c(t)|^2t(1-t)\f$
    */
    double Int_cc(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_abs2_1mt() and getKc_abs2_1mt^().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |k_c(t)|^2(1-t) + \int_{1-\delta(E_0)}^1 |k_c(t)|^2(1-t)^2\f$
    */
    double Int_cc1(double E0);
    
    
    /**
    * @brief Integral of the functions getKc_abs2_1mt().
    * @param[in] E0 energy cutoff
    * @return \f$\delta(E_0)\int_0^{1-\delta(E_0)} |k_c(t)|^2(1-t)\f$
    */
    double Int_cc1_part1(double E0);
    
    
    /**
    * @brief The 4-body NLO correction due to \f$Q_7\f$ and d, \f$ff^7_{d,MP}\f$, from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ff^7_{d,MP}\f$
    */
    double ff7_dMP(double E0);
    
    
    /**
    * @brief The 4-body NLO correction due to \f$Q_7\f$ and s, \f$ff^7_{s,MP}\f$, from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ff^7_{s,MP}\f$
    */
    double ff7_sMP(double E0);
    
    
    /**
    * @brief The 4-body NLO correction due to \f$Q_8\f$ and d, \f$ff^8_{d,MP}\f$, from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ff^8_{d,MP}\f$
    */
    double ff8_dMP(double E0);
    
    
    /**
    * @brief The 4-body NLO correction due to \f$Q_8\f$ and s, \f$ff^8_{s,MP}\f$, from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ff^8_{s,MP}\f$
    */
    double ff8_sMP(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{11}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{11}^{(1)} \f$
    */
    double Phi11_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{12}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{12}^{(1)} \f$
    */
    double Phi12_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{13}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{13}^{(1)} \f$
    */
    double Phi13_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{14}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{14}^{(1)} \f$
    */
    double Phi14_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{15}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{15}^{(1)} \f$
    */
    double Phi15_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{16}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{16}^{(1)} \f$
    */
    double Phi16_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{17}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ \Phi_{17}^{(1)} \f$
    */
    double Phi17_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{18}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ \Phi_{18}^{(1)} \f$
    */
    double Phi18_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{22}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{22}^{(1)} \f$
    */
    double Phi22_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{23}^{(1),{\rm 4-body}} \f$ function obtained from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{23}^{(1),{\rm 4-body}} \f$
    */
    double Phi23_1_4body(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{23}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{23}^{(1)} \f$
    */
    double Phi23_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{24}^{(1),{\rm 4-body}}  \f$ function obtained from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{24}^{(1),{\rm 4-body}} \f$
    */
    double Phi24_1_4body(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{24}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{24}^{(1)} \f$
    */
    double Phi24_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{25}^{(1),{\rm 4-body}}  \f$ function obtained from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{25}^{(1),{\rm 4-body}} \f$
    */
    double Phi25_1_4body(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{25}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{25}^{(1)} \f$
    */
    double Phi25_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{26}^{(1),{\rm 4-body}}  \f$ function obtained from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{26}^{(1),{\rm 4-body}} \f$
    */
    double Phi26_1_4body(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{26}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{26}^{(1)} \f$
    */
    double Phi26_1(double E0);
    
    
    /**
    * @brief The \f$ \Re \Phi_{27}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ \Re \Phi_{27}^{(1)} \f$
    */
    double Phi27_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Im\Phi_{27}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ \Im\Phi_{27}^{(1)} \f$
    */
    double Phi27_1_im(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{28}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{\rm kin}\f$
    * @return \f$ \Phi_{28}^{(1)} \f$
    */
    double Phi28_1(double E0, double z);
    
    /**
    * @brief The \f$ \Phi_{33}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{33}^{(1)} \f$
    */
    double Phi33_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{34}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{34}^{(1)} \f$
    */
    double Phi34_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{35}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{35}^{(1)} \f$
    */
    double Phi35_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{36}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{36}^{(1)} \f$
    */
    double Phi36_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{37}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{37}^{(1)} \f$
    */
    double Phi37_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{38}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{38}^{(1)} \f$
    */
    double Phi38_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{44}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{44}^{(1)} \f$
    */
    double Phi44_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{45}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{45}^{(1)} \f$
    */
    double Phi45_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{46}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{46}^{(1)} \f$
    */
    double Phi46_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{47}^{(1)} \f$ function from @cite Gambino:2001ew  and 
    * adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{47}^{(1)} \f$
    */
    double Phi47_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{48}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{48}^{(1)} \f$
    */
    double Phi48_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{55}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{55}^{(1)} \f$
    */
    double Phi55_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{56}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{56}^{(1)} \f$
    */
    double Phi56_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{57}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{57}^{(1)} \f$
    */
    double Phi57_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{58}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{58}^{(1)} \f$
    */
    double Phi58_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{66}^{(1)} \f$ function obtained using the prescription of @cite Chetyrkin:1996vx .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{66}^{(1)} \f$
    */
    double Phi66_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{67}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{67}^{(1)} \f$
    */
    double Phi67_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{68}^{(1)} \f$ function obtained using the prescription 
    * of @cite Chetyrkin:1996vx  and adding the 4-body contribution from @cite Huber:2014nna .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{68}^{(1)} \f$
    */
    double Phi68_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{77}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{77}^{(1)} \f$
    */
    double Phi77_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{78}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{78}^{(1)} \f$
    */
    double Phi78_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{88}^{(1)} \f$ function from @cite Gambino:2001ew .
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{88}^{(1)} \f$
    */
    double Phi88_1(double E0);
    
    
    /**
    * @brief The \f$ K_{ij}^{(1)} \f$ function from @cite Misiak:2006ab .
    * @param[in] i first index
    * @param[in] j second index
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ K_{ij}^{(1)} \f$
    */
    double Kij_1(int i, int j, double E0, double mu);
    
    
    /**
    * @brief Compute the Wilson Coefficient.
    * @param[in] mu low scale of the decay
    */
    void computeCoeff(double mu);
    
    
    /**
    * @brief The perturbative part \f$ P^{(0)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @param[in] E0 energy cutoff
    * @return \f$ P^{(0)} \f$
    */
    double P0(double E0);
    
    
    /**
    * @brief The perturbative part \f$ P_1^{(1)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @return \f$ P_1^{(1)} \f$
    */
    double P11();
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(1)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_2^{(1)} \f$
    */
    double P21(double E0, double mu);
    
    
    /**
    * @brief The perturbative part \f$ P_1^{(2)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @return \f$ P_1^{(2)} \f$
    */
    double P12();
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(2)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @param[in] E0 energy cutoff
    * @param[in] mu_b b quark scale
    * @param[in] mu_c c quark scale
    * @return \f$ P_2^{(2)} \f$
    */
    double P22(double E0, double mu_b, double mu_c);
    
    
    /**
    * @brief The perturbative part \f$ P_3^{(2)} \f$ of the BR as defined in @cite Misiak:2006ab .
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_3^{(2)} \f$
    */
    double P32(double E0, double mu);
    
    
    /**
    * @brief The 2 body NLO Vub part of the \f$BR\f$ as defined in @cite Gambino:2001ew , \f$Vub^{NLO}_{2b}\f$.
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{2b}\f$
    */
    double Vub_NLO_2body(bool CPodd);
    
    
    /**
    * @brief The 3 body NLO Vub part of the \f$BR\f$, \f$Vub^{NLO}_{3b}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{3b}\f$
    */
    double Vub_NLO_3body(double E0, bool CPodd);
    
    
    /**
    * @brief The 4 body NLO Vub part of the \f$BR\f$ obtained from @cite Huber:2014nna , \f$Vub^{NLO}_{4b}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{4b}\f$
    */
    double Vub_NLO_4body(double E0, bool CPodd);
    
    
    /**
    * @brief The total NLO Vub part of the \f$BR\f$, \f$Vub^{NLO}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}\f$
    */
    double Vub_NLO(double E0, bool CPodd);
    
    
    /**
    * @brief The perturbative part of the \f$BR\f$ as defined in @cite Misiak:2006ab , \f$P\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu_b b quark scale
    * @param[in] mu_c c quark scale
    * @param[in] order perturbation theory order
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$P\f$
    */
    double P(double E0, double mu_b, double mu_c, orders order, bool CPodd);
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ due to \f$Q_2-Q_7\f$ 
    * interference as defined in @cite Gambino:2001ew , \f$N_{27}\f$.
    * @return \f$N_{27}\f$
    */
    double N_27();
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ due to \f$Q_7-Q_7\f$ 
    * interference as defined in @cite Ewerth:2009yr , \f$N_{77}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu b quark scale
    * @return \f$N_{77}\f$
    */
    double N_77(double E0, double mu);
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ as defined in @cite Benzke:2010js , \f$N\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu b quark scale
    * @return \f$N\f$
    */
    double N(double E0, double mu);
    
    
    /**
    * @brief The ratio \f$C = | \frac{V_{ub}}{V_{cb}} |^2 \frac{\Gamma[\bar{B} \to X_c e \bar{\nu}]}{\Gamma[\bar{B} \to X_u e \bar{\nu}]} \f$ 
    * as defined in @cite Gambino:2013rza , but with coefficients 
    * slightly modified due to different imput parameters (obtained by private
    * conversation with Paolo Gambino).
    * @return \f$C\f$
    */
    double C_sem();
    
    
    /**
    * @brief The update parameter method for bsgamma.
    */
    void updateParameters();
    
    
    /**
    * @brief Computes the Branching Ratio for the \f$b \to s \gamma\f$ decay.
    * @return \f$BR\f$
    */
    double computeThValue();
    
    
private:
    double ale; /**<alpha electromagnetic */
    double alsUps; /**<alpha strong Upsilon */
    double Alstilde; /**<alpha strong divided by 4 pi */
    double E0; /**<energy cutoff */
    double mu_b; /**<b quark mass scale */
    double mu_c; /**<c quark mass scale */
    double mu_kin; /**<kinetic mass scale */
    double Mb_kin; /**<b quark mass in the kinetic scheme */
    double Mc; /**<c quark mass scale */
    double Ms;/**<s quark mass scale */
    double BRsl; /**<BR of the semileptonic decay \f$B \to X_c e \nu\f$ */
    double C; /**<The semileptonic phase space ratio */
    gslpp::complex lambda_t; /**<Vckm factor */
    gslpp::complex CKMu; /**<Vckm factor */
    gslpp::complex V_cb; /**<Vckm factor */
    double overall; /**<overall BR factor */
    double mu_pi2; /**<B meson expectation value of one of the relevant dim. 5 and 6 local operators*/
    double mu_G2; /**<B meson expectation value of one of the relevant dim. 5 and 6 local operators*/
    double rho_D3; /**<B meson expectation value of one of the relevant dim. 5 and 6 local operators*/
    double rho_LS3; /**<B meson expectation value of one of the relevant dim. 5 and 6 local operators*/
    double BLNPcorr; /**<non perturbative correction from @cite Benzke:2010js  */
    
    int obs; /**<observable type*/
    
    double BR; /**<BR of the decay */
    double BR_CPodd; /**<BR of the decay */
    
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients */
    
    gslpp::complex C1_0;/**<LO term of the Wilson coeffients @f$C_1@f$*/
    gslpp::complex C2_0;/**<LO term of the Wilson coeffients @f$C_2@f$*/
    gslpp::complex C3_0;/**<LO term of the Wilson coeffients @f$C_3@f$*/
    gslpp::complex C4_0;/**<LO term of the Wilson coeffients @f$C_4@f$*/
    gslpp::complex C5_0;/**<LO term of the Wilson coeffients @f$C_5@f$*/
    gslpp::complex C6_0;/**<LO term of the Wilson coeffients @f$C_6@f$*/
    gslpp::complex C7_0;/**<LO term of the Wilson coeffients @f$C_7@f$*/
    gslpp::complex C8_0;/**<LO term of the Wilson coeffients @f$C_8@f$*/
    
    gslpp::complex C1_1;/**<NLO term of the Wilson coeffients @f$C_1@f$*/
    gslpp::complex C2_1;/**<NLO term of the Wilson coeffients @f$C_2@f$*/
    gslpp::complex C3_1;/**<NLO term of the Wilson coeffients @f$C_3@f$*/
    gslpp::complex C4_1;/**<NLO term of the Wilson coeffients @f$C_4@f$*/
    gslpp::complex C5_1;/**<NLO term of the Wilson coeffients @f$C_5@f$*/
    gslpp::complex C6_1;/**<NLO term of the Wilson coeffients @f$C_6@f$*/
    gslpp::complex C7_1;/**<NLO term of the Wilson coeffients @f$C_7@f$*/
    gslpp::complex C8_1;/**<NLO term of the Wilson coeffients @f$C_8@f$*/
    
    gslpp::complex C7_2;/**<NNLO term of the Wilson coeffients @f$C_7@f$*/
    
    gslpp::complex C7p_0;/**<LO term of the Wilson coeffients @f$C'_7@f$*/
    gslpp::complex C7p_1;/**<NLO term of the Wilson coeffients @f$C_7@f$*/
    
    gsl_function INT;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_INT;/**< Gsl integral variable */
    double avaINT;/**< Gsl integral variable */    
    double errINT;/**< Gsl integral variable */
    
    unsigned int Intb1Cached;/**< Cache variable */
    unsigned int Intb2Cached;/**< Cache variable */
    unsigned int Intb3Cached;/**< Cache variable */
    unsigned int Intb4Cached;/**< Cache variable */
    unsigned int Intbb1Cached;/**< Cache variable */
    unsigned int Intbb2Cached;/**< Cache variable */
    unsigned int Intbb4Cached;/**< Cache variable */
    unsigned int Intbc1Cached;/**< Cache variable */
    unsigned int Intbc2Cached;/**< Cache variable */
    unsigned int Intc1Cached;/**< Cache variable */
    unsigned int Intc1imCached;/**< Cache variable */
    unsigned int Intc2Cached;/**< Cache variable */
    unsigned int Intc3Cached;/**< Cache variable */
    unsigned int IntccCached;/**< Cache variable */
    unsigned int Intcc1Cached;/**< Cache variable */
    unsigned int Intcc1p1Cached;/**< Cache variable */
    unsigned int IntPhi772rCached;/**< Cache variable */
    
    double CacheIntb1;/**< Cache variable */
    double CacheIntb2;/**< Cache variable */
    double CacheIntb3;/**< Cache variable */
    double CacheIntb4;/**< Cache variable */
    double CacheIntbb1;/**< Cache variable */
    double CacheIntbb2;/**< Cache variable */
    double CacheIntbb4;/**< Cache variable */
    double CacheIntbc1;/**< Cache variable */
    double CacheIntbc2;/**< Cache variable */
    double CacheIntc1;/**< Cache variable */
    double CacheIntc1im;/**< Cache variable */
    double CacheIntc2;/**< Cache variable */
    double CacheIntc3;/**< Cache variable */
    double CacheIntcc;/**< Cache variable */
    double CacheIntcc1;/**< Cache variable */
    double CacheIntcc1p1;/**< Cache variable */
    double CacheIntPhi772r;/**< Cache variable */
    
    unsigned int Intb_updated;/**< Cache variable */
    unsigned int Intbc_updated;/**< Cache variable */
    
    double Intb_cache;/**< Cache variable */
    gslpp::vector<double> Intbc_cache;/**< Cache variable */
    
    /**
     * @brief The caching method for bsgamma.
     */
    void checkCache();
    
};

#endif	/* BSGAMMA_H */