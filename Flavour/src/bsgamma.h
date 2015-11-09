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

/**
 * @class Bsgamma
 * @ingroup Flavour
 * @brief A class for the @f$b \to s \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$b \to s \gamma@f$ decay.
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
    * @brief The cutoff energy function \f$ \delta \f$.
    * @param[in] E0 cutoff energy
    * @return \f$ \delta(E0) \f$ 
    */
    double delta(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ \rho \f$ as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @return \f$ \rho(E0) \f$ 
    */
    double rho(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ \omega \f$ as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @return \f$ \omega(E0) \f$ 
    */
    double omega(double E0);
    
    
    /**
    * @brief The cutoff energy function \f$ T_1 \f$ as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_1(E0) \f$ 
    */
    double T1(double E0, double t);
    
    
    /**
    * @brief The cutoff energy function \f$ T_2 \f$ as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_2(E0) \f$ 
    */
    double T2(double E0, double t);
    
    
    /**
    * @brief The cutoff energy function \f$ T_3 \f$ as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ T_3(E0) \f$ 
    */
    double T3(double E0, double t);
    
    
    /**
    * @brief The tree level LO contribution as defined in arXiv:1209.0965v3.
    * @param[in] E0 cutoff energy
    * @param[in] t squared ratio between b quark and s quark masses
    * @return \f$ P_{tree}^{(0)} \f$ 
    */
    double P0tree(double E0, double t);
    
    
    /**
    * @brief The squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$, \f$ z \f$.
    * @return \f$ z \f$
    */
    double zeta();
    
    
    /**
    * @brief The funcion \f$ a(z) \f$ as defined in hep-ph/0203135.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ a(z) \f$ 
    */
    gslpp::complex a(double z);
    
    
    /**
    * @brief The funcion \f$ b(z) \f$ as defined in hep-ph/0203135.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ b(z) \f$ 
    */
    gslpp::complex b(double z);
    
    
    /**
    * @brief The funcion \f$ r_i^{(1)}(z) \f$ as defined in hep-ph/0203135.
    * @param[in] i function index
    * @param[in] z squared ratio between m_c and m_b^{1s}
    * @return \f$ r_i(z)^{(1)} \f$
    */
    gslpp::complex r1(int i, double z);
    
    
    /**
    * @brief The function \f$ \Gamma \f$ as defined in hep-ph/0104034v2.
    * @param[in] t dummy variable to be integrated out
    * @return \f$ \Gamma \f$ 
    */
    gslpp::complex Gamma_t(double t);
    
    
    /**
    * @brief The function \f$ k \f$ as defined in hep-ph/9512252v2.
    * @param[in] Mq quark mass
    * @param[in] t dummy variable to be integrated out
    * @return \f$ k \f$ 
    */
    gslpp::complex kappa(double Mq, double t);
    
    
    /**
    * @brief The function \f$|k_c(t)|^2 t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|k_c(t)|^2(1 - t)\f$
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
    * @return \f$|k_c(t)|^2(1 - t)\f$
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
    * @brief The function \f$Re(k_b(t))k_c(t)(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))k_c(t)(1-t)\f$
    */
    double getKc_re_Kb_1mt(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))k_c(t)(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))k_c(t)(1-t)^2\f$
    */
    double getKc_re_Kb_1mt2(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))k_c(t)t(1-t)\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))k_c(t)t(1-t)\f$
    */
    double getKc_re_Kb_t_1mt(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * t * (1. - t);
    };
    
    
    /**
    * @brief The function \f$Re(k_b(t))k_c(t)t(1-t)^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$Re(k_b(t))k_c(t)t(1-t)^2\f$
    */
    double getKc_re_Kb_t_1mt2(double t)
    {
        return kappa(Mc,t).real() * kappa(Mb_kin,t).real() * t * (1. - t) * (1. - t);
    };
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_b1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_b2(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_b3(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_b4(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_bb1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_bb2(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_bb4(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_bc1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_bc2(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_c1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_c1_im(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_c2(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_c3(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_cc(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_cc1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double Int_cc1_part1(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double ff7_dMP(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double ff7_sMP(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double ff8_dMP(double E0);
    
    
    /**
    * @brief xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxxxxxxxxx
    */
    double ff8_sMP(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{11}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{11}^{(1)} \f$
    */
    double Phi11_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{12}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{12}^{(1)} \f$
    */
    double Phi12_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{13}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{13}^{(1)} \f$
    */
    double Phi13_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{14}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{14}^{(1)} \f$
    */
    double Phi14_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{15}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{15}^{(1)} \f$
    */
    double Phi15_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{16}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{16}^{(1)} \f$
    */
    double Phi16_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{17}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Phi_{17}^{(1)} \f$
    */
    double Phi17_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{18}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Phi_{18}^{(1)} \f$
    */
    double Phi18_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{22}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{22}^{(1)} \f$
    */
    double Phi22_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{23}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{23}^{(1)} \f$
    */
    double Phi23_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{24}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{24}^{(1)} \f$
    */
    double Phi24_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{25}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{25}^{(1)} \f$
    */
    double Phi25_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{26}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{26}^{(1)} \f$
    */
    double Phi26_1(double E0);
    
    
    /**
    * @brief The \f$ \Re \Phi_{27}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Re \Phi_{27}^{(1)} \f$
    */
    double Phi27_1(double E0, double z);
    
    
    /**
    * @brief The \f$ \Im\Phi_{27}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Im\Phi_{27}^{(1)} \f$
    */
    double Phi27_1_im(double E0, double z);
    
    
    /**
    * @brief The \f$ \Phi_{28}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Phi_{28}^{(1)} \f$
    */
    double Phi28_1(double E0, double z);
    
    /**
    * @brief The \f$ \Phi_{33}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{33}^{(1)} \f$
    */
    double Phi33_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{34}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{34}^{(1)} \f$
    */
    double Phi34_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{35}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{35}^{(1)} \f$
    */
    double Phi35_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{36}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{36}^{(1)} \f$
    */
    double Phi36_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{37}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{37}^{(1)} \f$
    */
    double Phi37_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{38}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{38}^{(1)} \f$
    */
    double Phi38_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{44}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{44}^{(1)} \f$
    */
    double Phi44_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{45}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{45}^{(1)} \f$
    */
    double Phi45_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{46}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{46}^{(1)} \f$
    */
    double Phi46_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{47}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{47}^{(1)} \f$
    */
    double Phi47_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{48}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{48}^{(1)} \f$
    */
    double Phi48_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{55}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{55}^{(1)} \f$
    */
    double Phi55_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{56}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{56}^{(1)} \f$
    */
    double Phi56_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{57}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{57}^{(1)} \f$
    */
    double Phi57_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{58}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{58}^{(1)} \f$
    */
    double Phi58_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{66}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{66}^{(1)} \f$
    */
    double Phi66_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{67}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{67}^{(1)} \f$
    */
    double Phi67_1(double E0);
    
    /**
    * @brief The \f$ \Phi_{68}^{(1)} \f$ function from xxxxxxxxxxxxxxxx.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{68}^{(1)} \f$
    */
    double Phi68_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{77}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{77}^{(1)} \f$
    */
    double Phi77_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{78}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{78}^{(1)} \f$
    */
    double Phi78_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{88}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{88}^{(1)} \f$
    */
    double Phi88_1(double E0);
    
    
    /**
    * @brief The \f$ K_{ij}^{(1)} \f$ function from arXiv:1005.1173.
    * @param[in] i first index
    * @param[in] j second index
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ K_{ij}^{(1)} \f$
    */
    double Kij_1(int i, int j, double E0, double mu);
    
    
    /**
    * @brief The \f$ \Re r_2^{(2)} \f$ function from yyyyyyyyyyyyy.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ \Re r_2^{(2)} \f$
    */
    double Rer22(double z);
    
    
    /**
    * @brief The \f$ \Phi_{22}^{(2)\beta_0} \f$ function from yyyyyyyyyyyyy.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ \Phi_{22}^{(2)\beta_0} \f$
    */
    double Phi22_2beta0(double E0, double mu);
    
    
    /**
    * @brief The \f$ \Phi_{28}^{(2)\beta_0} \f$ function from yyyyyyyyyyyyy.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ \Phi_{28}^{(2)\beta_0} \f$
    */
    double Phi28_2beta0(double E0, double mu);
    
    
    /**
    * @brief The \f$ \Phi_{77}^{(2)\beta_0} \f$ function from yyyyyyyyyyyyy.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ \Phi_{77}^{(2)\beta_0} \f$
    */
    double Phi77_2beta0(double E0, double mu);
    
    
    /**
    * @brief The \f$ \Phi_{88}^{(2)\beta_0} \f$ function from yyyyyyyyyyyyy.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ \Phi_{88}^{(2)\beta_0} \f$
    */
    double Phi88_2beta0(double E0, double mu);
    
    
    /**
    * @brief The \f$ \delta Y^{(1)}(z_0,\mu) \f$ function from arXiv:0805.3911v2.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ \delta Y^{(1)}(z_0,\mu) \f$
    */
    double dY1(double E0);
    
    
    /**
    * @brief The \f$ Y^{(1)}(z_0,\mu) \f$ function from arXiv:0805.3911v2.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(1)}(z_0,\mu) \f$
    */
    double Y1(double E0, double mu);
    
    
    /**
    * @brief The \f$ Y^{(2,CF)}(z_0,\mu) \f$ function from arXiv:1005.5587v1.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2,CF)}(z_0,\mu) \f$
    */
    double Y2CF(double E0, double mu);
    
    
    /**
    * @brief The \f$ Y^{(2,CA)}(z_0,\mu) \f$ function from arXiv:1005.5587v1.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2,CA)}(z_0,\mu) \f$
    */
    double Y2CA(double E0, double mu);
    
    
    /**
    * @brief The \f$ Y^{(2,NL)}(z_0,\mu) \f$ function from arXiv:0805.3911v2.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2,NL)}(z_0,\mu) \f$
    */
    double Y2NL(double E0, double mu);
    
    
    /**
    * @brief The \f$ \Phi_1(\rho) \f$ function from arXiv:0805.3911v2.
    * @param[in] rho squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$\Phi_1(\rho) \f$
    */
    double Y2NV_PHI1(double rho);
    
    
    /**
    * @brief The \f$ \Phi_2(\rho) \f$ function from arXiv:0805.3911v2.
    * @param[in] rho squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$\Phi_2(\rho) \f$
    */
    double Y2NV_PHI2(double rho);
    
    
    /**
    * @brief The \f$ \Phi_3(\rho) \f$ function from arXiv:0805.3911v2.
    * @param[in] rho squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$\Phi_3(\rho) \f$
    */
    double Y2NV_PHI3(double rho);
    
    
    /**
    * @brief The \f$ \Phi_4(\rho) \f$ function from arXiv:0805.3911v2.
    * @param[in] rho squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$\Phi_4(\rho) \f$
    */
    double Y2NV_PHI4(double rho);
    
    
    /**
    * @brief The \f$ Y^{(2,NL)}(z_0,\mu) \f$ function from arXiv:0805.3911v2.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2,NL)}(z_0,\mu) \f$
    */
    double Y2NV(double E0, double mu);
    
    
    /**
    * @brief The \f$ Y^{(2,NL)}(z_0,\mu) \f$ function from arXiv:0805.3911v2.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2,NL)}(z_0,\mu) \f$
    */
    double Y2NH(double E0, double mu);
    
    
    /**
    * @brief The \f$ Y^{(2)}(z_0,\mu) \f$ function from arXiv:0805.3911v2 and arXiv:1005.5587v1.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ Y^{(2)}(z_0,\mu) \f$
    */
    double Y2(double E0, double mu);
    
    
    /**
    * @brief The \f$ f_{\rm NLO}(z,1) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ f_{\rm NLO}(z,1) \f$
    */
    double f_NLO_1(double z);
    
    
    /**
    * @brief The \f$ z \frac{d}{dz}f_{\rm NLO}(z,\delta) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return \f$ z \frac{d}{dz}f_{\rm NLO}(z,\delta)) \f$
    */
    double zdz_f_NLO(double z, double E0);
    
    
    /**
    * @brief The \f$ (1. - \delta)\frac{d}{d\delta}f_{\rm NLO}(z,\delta) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return \f$ (1. - \delta)\frac{d}{d\delta}f_{\rm NLO}(z,\delta) \f$
    */
    double mddel_f_NLO(double z, double E0);
    
    
    /**
    * @brief The \f$ h_{27}^{(2)}(z,\delta) \f$ function from arXiv:1009.5685 and arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return \f$ h_{27}^{(2)}(z,\delta) \f$
    */
    double h27_2(double z, double E0);
    
    
    /**
    * @brief The \f$ f_{q}(z,1) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return \f$ f_{q}(z,1) \f$
    */
    double f_q(double z, double E0);
    
    
    /**
    * @brief The \f$ f_{b}(z) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ f_{b}(z) \f$
    */
    double f_b(double z);
    
    
    /**
    * @brief The \f$ f_{c}(z) \f$ function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ f_{c}(z) \f$
    */
    double f_c(double z);
    
    
    /**
    * @brief The \f$ F_{1}(z,1) \f$ interpolated function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ F_{1}(z,1) \f$
    */
    double F_1(double z);
    
    
    /**
    * @brief The \f$ F_{2}(z,1) \f$ interpolated function from arXiv:1503.01791.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ F_{2}(z,1) \f$
    */
    double F_2(double z);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxx
    */
    double delddel_Phi22_1(double E0);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxx
    */
    double zdz_Phi22_1(double E0);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxx
    */
    double delddel_Phi28_1(double z, double E0);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxx
    */
    double zdz_Phi28_1(double z, double E0);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] E0 energy cutoff
    * @return xxxxxxxxx
    */
    double delddel_Phi88_1(double E0);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] r ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return xxxxxxxxx
    */
    double f_AEGG(double r);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] r ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return xxxxxxxxx
    */
    double delta_GBGS(double r);
    
    
    /**
    * @brief The xxxxxxxxx function from yyyyyyyy.
    * @param[in] r ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return xxxxxxxxx
    */
    double f_u(double r);
    
    
    /**
    * @brief The xxxxxx function from yyyyyyyyy
    * @param[in] z integration variable
    * @return xxxxxxxx
    */
    double omega77(double z);
    
    
    /**
    * @brief The xxxxxx function from yyyyyyyyy
    * @param[in] E0 energy cutoff
    * @return xxxxxxxx
    */
    double Int_Phi77_2rem(double E0);
    
    
    /**
    * @brief The xxxxxx function from yyyyyyyyy
    * @param[in] E0 energy cutoff
    * @return xxxxxxxx
    */
    double Phi77_2rem(double E0);
    
    
    /**
    * @brief The xxxxxx function from yyyyyyyyy
    * @param[in] E0 energy cutoff
    * @param[in] mu b quark scale
    * @return xxxxxxxx
    */
    double K77_2_z1(double E0, double mu);
    
    
    /**
    * @brief The \f$ K_{ij}^{(2)} \f$ function from arXiv:1503.01791.
    * @param[in] i first index
    * @param[in] j second index
    * @param[in] E0 energy cutoff
    * @param[in] mu_b b quark scale
    * @param[in] mu_c c quark scale
    * @return \f$ K_{ij}^{(2)} \f$
    */
    double Kij_2(int i, int j, double E0, double mu_b, double mu_c);
    
    
    /**
    * @brief Compute the Wilson Coefficient.
    * @param[in] mu low scale of the decay
    */
    void computeCoeff(double mu);
    
    
    /**
    * @brief The perturbative part \f$ P^{(0)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @return \f$ P^{(0)} \f$
    */
    double P0(double E0);
    
    
    /**
    * @brief The perturbative part \f$ P_1^{(1)} \f$ of the BR as defined in arXiv:1005.1173.
    * @return \f$ P_1^{(1)} \f$
    */
    double P11();
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(1)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_2^{(1)} \f$
    */
    double P21(double E0, double mu);
    
    
    /**
    * @brief The perturbative part \f$ P_1^{(2)} \f$ of the BR as defined in arXiv:1005.1173.
    * @return \f$ P_1^{(2)} \f$
    */
    double P12();
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(2)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu_b b quark scale
    * @param[in] mu_c c quark scale
    * @return \f$ P_2^{(2)} \f$
    */
    double P22(double E0, double mu_b, double mu_c);
    
    
    /**
    * @brief The perturbative part \f$ P_3^{(2)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_3^{(2)} \f$
    */
    double P32(double E0, double mu);
    
    
    /**
    * @brief The 2 body NLO Vub part of the \f$BR\f$ as defined in xxxxxxxxx, \f$Vub^{NLO}_{2b}\f$.
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{2b}\f$
    */
    double Vub_NLO_2body(bool CPodd);
    
    
    /**
    * @brief The 3 body NLO Vub part of the \f$BR\f$ as defined in xxxxxxxxx, \f$Vub^{NLO}_{3b}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{3b}\f$
    */
    double Vub_NLO_3body(double E0, bool CPodd);
    
    
    /**
    * @brief The 4 body NLO Vub part of the \f$BR\f$ as defined in xxxxxxxxx, \f$Vub^{NLO}_{4b}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}_{4b}\f$
    */
    double Vub_NLO_4body(double E0, bool CPodd);
    
    
    /**
    * @brief The NLO Vub part of the \f$BR\f$ as defined in xxxxxxxxx, \f$Vub^{NLO}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$Vub^{NLO}\f$
    */
    double Vub_NLO(double E0, bool CPodd);
    
    
    /**
    * @brief The perturbative part of the \f$BR\f$ as defined in arXiv:1005.1173, \f$P\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu_b b quark scale
    * @param[in] mu_c c quark scale
    * @param[in] order perturbation theory order
    * @param[in] CPodd switch to allow for CPodd terms
    * @return \f$P\f$
    */
    double P(double E0, double mu_b, double mu_c, orders order, bool CPodd);
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ as defined in xxxxxxxxx, \f$N_{27}\f$.
    * @return \f$N_{27}\f$
    */
    double N_27();
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ as defined in xxxxxxxxx, \f$N_{77}\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu b quark scale
    * @return \f$N_{77}\f$
    */
    double N_77(double E0, double mu);
    
    
    /**
    * @brief The non perturbative part of the \f$BR\f$ as defined in xxxxxxxxx, \f$N\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu b quark scale
    * @return \f$N\f$
    */
    double N(double E0, double mu);
    
    
    /**
    * @brief The \f$BR\f$ as defined in arXiv:1005.1173.
    * @param[in] order perturbation theory order
    */
    void computeBR(orders order);
    
    
    /**
    * @brief The \f$BR\f$ as defined in arXiv:1005.1173.
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
    double Mb_kin; /**<b quark mass in the 1s scheme */
    double Mc; /**<c quark mass scale */
    double Ms;/**<s quark mass scale */
    double BRsl; /**<BR of the semileptonic decay \f$B \to X_c e \nu\f$ */
    double C; /**<The semileptonic phase space ratio */
    gslpp::complex lambda_t; /**<Vckm factor */
    gslpp::complex CKMu; /**<Vckm factor */
    gslpp::complex V_cb; /**<Vckm factor */
    double overall; /**<overall BR factor */
    double mu_pi2;
    double mu_G2;
    double rho_D3;
    double rho_LS3;
    
    int obs; /**<observable type*/
    
    double BR; /**<BR of the decay */
    double BR_CPodd; /**<BR of the decay */
    
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    
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
    
    gsl_function INT;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_INT;/**< Gsl integral variable */
    double avaINT;/**< Gsl integral variable */    
    double errINT;/**< Gsl integral variable */
    
    unsigned int Intb1Cached;
    unsigned int Intb2Cached;
    unsigned int Intb3Cached;
    unsigned int Intb4Cached;
    unsigned int Intbb1Cached;
    unsigned int Intbb2Cached;
    unsigned int Intbb4Cached;
    unsigned int Intbc1Cached;
    unsigned int Intbc2Cached;
    unsigned int Intc1Cached;
    unsigned int Intc1imCached;
    unsigned int Intc2Cached;
    unsigned int Intc3Cached;
    unsigned int IntccCached;
    unsigned int Intcc1Cached;
    unsigned int Intcc1p1Cached;
    unsigned int IntPhi772rCached;
    
    double CacheIntb1;
    double CacheIntb2;
    double CacheIntb3;
    double CacheIntb4;
    double CacheIntbb1;
    double CacheIntbb2;
    double CacheIntbb4;
    double CacheIntbc1;
    double CacheIntbc2;
    double CacheIntc1;
    double CacheIntc1im;
    double CacheIntc2;
    double CacheIntc3;
    double CacheIntcc;
    double CacheIntcc1;
    double CacheIntcc1p1;
    double CacheIntPhi772r;
    
    /**
     * @brief The caching method for bsgamma.
     */
    void checkCache();
    
    unsigned int Intb_updated;/**< Cache variable */
    unsigned int Intbc_updated;/**< Cache variable */
    
    double Intb_cache;/**< Cache variable */
    gslpp::vector<double> Intbc_cache;/**< Cache variable */
    
    EvolDB1bsg myevol;
    
};

#endif	/* BSGAMMA_H */