/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVGAMMA_H
#define	MVGAMMA_H

class StandardModel;
#include "AmpDB2.h"
#include "BXqll.h"
#include "ThObservable.h"
#include "gslpp_function_adapter.h"
#include <gsl/gsl_integration.h>

#define NFPOLARBASIS_MVGAMMA true

/**
 * @class MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$M \to V \gamma@f$ decay.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to
 * compute the observables relative to the @f$M \to V \gamma@f$ decays, where
 * @f$M@f$ is a generic meson and @f$V@f$ is a vector meson. 
 * 
 * @anchor MVllParameters
 * <h3>%MVll parameters</h3>
 *
 * The mandatory parameters of %MVll are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T1</td>
 *   <td class="mod_symb">@f$a_0^{T_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_1@f$ of the @f$B\to K^*@f$ at @f$q^2=0@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a_0T1phi</td>
 *   <td class="mod_symb">@f$a_0^{T_1}@f$</td>
 *   <td class="mod_desc">The fit parameters for the form factor @f$T_1@f$ of the @f$B\to\phi@f$ at @f$q^2=0@f$.</td>
 * </tr>
  *   <td class="mod_name">%absh_p</td>
 *   <td class="mod_symb">@f$\mathrm{Abs}h_+^{(0)}@f$</td>
 *   <td class="mod_desc">The constant term of the absolute value of the hadronic parameter @f$h_+@f$ of the @f$B\to K^*@f$ at @f$q^2=0@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%argh_p</td>
 *   <td class="mod_symb">@f$\mathrm{Arg}h_+^{(0)}@f$</td>
 *   <td class="mod_desc">The constant term of the argument of the hadronic parameter @f$h_+@f$ of the @f$B\to K^*@f$ at @f$q^2=0@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%absh_m</td>
 *   <td class="mod_symb">@f$\mathrm{Abs}h_-^{(0)}@f$</td>
 *   <td class="mod_desc">The constant term of the absolute value of the hadronic parameter @f$h_-@f$ of the @f$B\to K^*@f$ at @f$q^2=0@f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%argh_m</td>
 *   <td class="mod_symb">@f$\mathrm{Arg}h_-^{(0)}@f$</td>
 *   <td class="mod_desc">The constant term of the argument of the hadronic parameter @f$h_-@f$ of the @f$B\to K^*@f$ at @f$q^2=0@f$.</td>
 * </tr>
 * </table>
 *  
 * This kind of decays can be described
 * by means of the @f$\Delta B = 1 @f$ weak effective Hamiltonian
 * @f[
 *   \mathcal{H}_\mathrm{eff}^{\Delta B = 1} = \mathcal{H}_\mathrm{eff}^\mathrm{had} +
 *   \mathcal{H}_\mathrm{eff}^\mathrm{\gamma},
 * @f]  
 * where the first term is the hadronic contribution 
 * @f[
 * \mathcal{H}_\mathrm{eff}^\mathrm{had} = \frac{4G_F}{\sqrt{2}}\Bigg[\sum_{p=u,c}\lambda_p\bigg(C_1 Q^{p}_1 
 * + C_2 Q^{p}_2\bigg) -\lambda_t \bigg(\sum_{i=3}^{6} C_i P_i + C_{8}Q_{8g} \bigg)\Bigg] \,,
 * @f]
 * involving current-current, chromodynamic penguin and chromomagnetic dipole operators, while the second one, given by
 * @f[
 * \mathcal{H}_\mathrm{eff}^\mathrm{\gamma} = - \frac{4G_F}{\sqrt{2}}\lambda_t C_7Q_{7\gamma} \,, 
 * @f]
 * includes the electromagnetic penguin operator.
 * 
 * Considering the matrix element of @f$\mathcal{H}_\mathrm{eff}^{\Delta B = 1}@f$
 * between the initial state @f$M@f$ and the final state @f$V \gamma@f$, only the contribution of 
 * @f$\mathcal{H}_\mathrm{eff}^\mathrm{\gamma}@f$ clearly factorizes into the 
 * product of hadronic form factors and leptonic tensors at all orders in strong interactions. 
 * Following @cite Jager:2012uw, we implemented the amplitude in the helicity basis; 
 * hence we made use of the helicity form factor @f$ T_-(0)@f$, which is related to the
 * ones in the transverse basis through the following relation:
 * @f[
 * T_{-}\left( q^{2}\right) = \frac{m_M^2 - m_V^2}{m_M^2}T_1\left( q^{2}\right)\,.
 * @f]
 * 
 * The effect of the operators of @f$\mathcal{H}_\mathrm{eff}^\mathrm{had}@f$ due to
 * exchange of soft gluon can be reabsorbed in the following parameterization,
 * @f[
 * h_\lambda(q^2) = \frac{\epsilon^*_\mu(\lambda)}{m_M^2} 
 * \int d^4x e^{iqx} \langle \bar V \vert T\{j^{\mu}_\mathrm{em} (x) 
 * \mathcal{H}_\mathrm{eff}^\mathrm{had} (0)\} \vert \bar M \rangle = 
 * h_\lambda^{(0)}\,,
 * @f]
 * while the effect due to exchange of hard gluons can be parametrized following 
 * the prescription of @cite Bosch:2001gv as a shift to the Wilson coefficient @f$C_7@f$ :
 * @f[
 * \Delta C_{7} = \frac{\alpha_s(\mu) C_F}{4\pi} \left( C_1(\mu) G_1(s_p)+ C_8(\mu) G_8\right)  
  + \frac{\alpha_s(\mu_h) C_F}{4\pi} \left( C_1(\mu_h) H_1(s_p)+ C_8(\mu_h) H_8\right)\,,
 * @f]
 * where the terms proportional to @f$G_i@f$ are the ones describing the 
 * corrections where the spectator quark is connected to the hard process only 
 * through soft interactions, while the ones proportional to @f$H_i@f$ 
 * (involving leading twist light-cone distributions) are the ones describing 
 * the corrections where the spectactor quark is involved in the hard process,
 * and @f$s=\frac{m_c^2}{m_b^2}@f$.
 * 
 * The amplitude can be therefore parametrized in terms of the following helicity amplitudes:
 * @f[
 * H_V^+ = \lambda_t \Big[- C_{7}' {T}_{-}
            - \frac{m_M}{m_b} 8 \pi^2 h_\lambda \Big]  \,,  \\
 * H_V^- = \lambda_t \Big[ C_{7} T_{-}
            - \frac{m_M}{m_b} 8 \pi^2 h_\lambda \Big]  \,.
 * @f]
 * 
 * Squaring the amplitude and summing over the spins it is possible to obtain 
 * the Branching Ratio, which is
 * @f[
 * BR = \frac {\alpha_e G_F^2 M_b^2 M_M \lambda}{(4\pi)^2 4 w_M} ( |H_V^+|^2 + |H_A^+|^2 +|\overline{H}_V^-|^2 + |\overline{H}_A^-|^2) \,.
 * @f]
 * 
 * The class is build as follows: after the
 * parameters are updated in updateParameters() and the form factor @f$ T_1 @f$
 * is computed in T_1() following @cite Straub:2015ica, the QCDF corrections to the Wilson coefficient @f$ C_7 @f$
 * is computed in the functions G1(), G8(), H1() and H8(). The helicity amplitudes
 * @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$ are build in H_V_p(), H_V_m(), H_V_p_bar()
 * and H_V_m_bar(), in order to be further used to compute the observables.
 */
class MVgamma {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~MVgamma();
    
    /**
     * @brief The update parameter method for MVgamma.
     */
    void updateParameters();
    
    /**
    * @brief A get method for the parameters necessary for MVgamma.
    * @return the vector of MVgamma specific parameters
    */
    std::vector<std::string> getMVgammaParameters()
    {
        return mVgammaParameters;
    }
    
    std::vector<std::string> parametersForMVgamma;
    
    double GF;            /**<Fermi constant */
    double ale;           /**<alpha electromagnetic */
    double MM;            /**<initial meson mass */
    double MM2;           /**<square of the initial meson mass */
    double MV;            /**<final vector meson mass */
    double Mb;            /**<b quark mass */
    double Mc;            /**<c quark mass */
    double mu_b;          /**<b mass scale */
    double mu_h;          /**<sqrt(mu_b*lambda_QCD) */
    double width;         /**<initial meson width */
    double fperp;         /**<vector meson perpendicular decay constant*/
    double fpara;         /**<vector meson decay constant*/
    double fB;            /**<B meson decay constant*/
    double Ms;            /**<s quark mass */
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor lambds_t*/
    gslpp::complex lambda_u;     /**<Vckm factor lambda_u*/
    gslpp::complex h[2];         /**<parameter that contains the contribution from the hadronic hamiltonian */
    double lambda;        /**<kinematic parameter */
    double spectator_charge; /**<charge of the spectator quark. */
    double alpha_s_mub; /**<@f\aplha_s(\mu_b)$@f$ */
    gslpp::complex DC7_QCDF;
    gslpp::complex DC7_QCDF_bar;
    
    double a_0T1;/**<LCSR fit parameter */
    double a_1T1;/**<LCSR fit parameter */
    double a_2T1;/**<LCSR fit parameter */
    double MRT1_2;/**<LCSR fit parameter */
    
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients at mub*/
    gslpp::vector<gslpp::complex> ** allcoeffh;/**<vector that contains the Wilson coeffients at muh*/
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients at mub*/
    
    gslpp::complex C_1;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_3;/**<Wilson coeffients @f$C_3@f$*/
    gslpp::complex C_4;/**<Wilson coeffients @f$C_4@f$*/
    gslpp::complex C_5;/**<Wilson coeffients @f$C_5@f$*/
    gslpp::complex C_6;/**<Wilson coeffients @f$C_6@f$*/
    gslpp::complex C_7;/**<Wilson coeffients @f$C_7@f$*/
    gslpp::complex C_7p;/**<Wilson coeffients @f$C_7'@f$*/
    gslpp::complex C_1_bar;/**<Wilson coeffients @f$C_1_bar(mu_b)@f$*/
    gslpp::complex C_2_bar;/**<Wilson coeffients @f$C_2_bar(mu_b)@f$*/
    gslpp::complex C_8;/**<Wilson coeffients @f$C_8(mu_b)@f$*/
    gslpp::complex C_2h_bar;/**<Wilson coeffients @f$C_2(mu_h)@f$*/
    gslpp::complex C_8h;/**<Wilson coeffients @f$C_8(mu_h)@f$*/
    
    
    /**
    * @brief The transverse form factor @f$ T_1 @f$.
    * @return @f$ T_1 @f$ 
    */
    double T_1();

    /**
    * @brief The helicity amplitude @f$ H_V^+ @f$.
    * @return @f$ H_V^+ @f$ 
    */
    gslpp::complex H_V_p();
    
    /**
    * @brief The helicity amplitude @f$ H_V^- @f$.
    * @return @f$ H_V^- @f$ 
    */
    gslpp::complex H_V_m();
    
    /**
    * @brief The helicity amplitude @f$ \bar{H}_V^+ @f$.
    * @return @f$ \bar{H}_V^+ @f$ 
    */
    gslpp::complex H_V_p_bar();
    
    /**
    * @brief The helicity amplitude @f$ \bar{H}_V^- @f$.
    * @return @f$ \bar{H}_V^- @f$ 
    */
    gslpp::complex H_V_m_bar();

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @param s @f$ m_c^2/m_b^2  @f$
     * @return @f$ G_1(s) @f$
     */
//    gslpp::complex G1(double s);
    
    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @return @f$ G_8 @f$
     */
//    gslpp::complex G8();

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @param s @f$ m_c^2/m_b^2  @f$
     * @return @f$ H_1(s) @f$
     */
//    gslpp::complex H1(double s); */

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @return  @f$ H_8 @f$
     */
//    gslpp::complex H8();
    
    /**
     * @brief QCDF Correction from various BFS papers (hep-ph/0403185, hep-ph/0412400) and Greub et. al (arXiv:0810.4077)..
     * @param conjugate a boolean to control conjugation
     * @return @f$ \Delta C_{7}^{QCDF} @f$
     */
    gslpp::complex deltaC7_QCDF(bool conjugate);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0412400). Part of Weak Annihilation.
     * @param conjugate a boolean to control conjugation
     * @return @f$ C_{34}^{q} @f$
     */
    gslpp::complex Cq34(bool conjugate);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0412400). Weak Annihilation.
     * @return @f$ T_{perp}^{ann,1} @f$
     */
    gslpp::complex T_perp_WA_1();
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0412400). Weak Annihilation.
     * @param conjugate a boolean to control conjugation
     * @return @f$ T_{perp}^{ann,2} @f$
     */
    gslpp::complex T_perp_WA_2(bool conjugate);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0106067). Part of non-factorizable spectator contribution
     * @param x complex argument
     * @return @f$ L_{1} @f$
     */
    gslpp::complex L1(gslpp::complex x);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0106067).Vector meson distribution amplitude
     * @param u integration variable in the range [0, 1]
     * @return @f$ \Delta L_{1} @f$
     */
    double phi_V(double u);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0106067). Part of 4 quark operator contribution.
     * @param u integration variable in the range [0, 1]
     * @param m mass of the quark 
     * @return @f$ t_{perp} @f$
     */
    gslpp::complex t_perp(double u, double m);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0106067). 4 quark operator contribution.
     * @param u integration variable in the range [0, 1]
     * @param conjugate a boolean to control conjugation
     * @return @f$ T_{perp,+}^{QSS} @f$
     */
    gslpp::complex T_perp_plus_QSS(double u, bool conjugate);
    
    /**
     * @brief QCDF Correction from various BFS paper (hep-ph/0106067). Chromomagnetic dipole contribution contribution.
     * @param u integration variable in the range [0, 1]
     * @return @f$ T_{perp,+}^{O8} @f$
     */
    gslpp::complex T_perp_plus_O8(double u);
    
    /**
     * @brief QCDF Correction from various BFS papers (hep-ph/0106067, hep-ph/0412400). Total.
     * @param u integration variable in the range [0, 1]
     * @param conjugate a boolean to control conjugation
     * @return @f$ T_{perp}^{total} @f$
     */
    gslpp::complex T_perp(double u, bool conjugate);

    /**
     * @brief Real part of the integrand for QCDF Correction
     * @param u integration variable in the range [0, 1]
     * @return the real part of the integrand
     */
    double getT_perp_integrand_real(double u) {
        return T_perp(u, false).real();
    };

    /**
     * @brief Imaginary part of the integrand for QCDF Correction
     * @param u integration variable in the range [0, 1]
     * @return the imaginary part of the integrand
     */
    double getT_perp_integrand_imag(double u) {
        return T_perp(u, false).imag();
    };

    /**
     * @brief Real part of the conjugate integrand for QCDF Correction
     * @param u integration variable in the range [0, 1]
     * @return the real part of the conjugate integrand
     */
    double getT_perp_bar_integrand_real(double u) {
        return T_perp(u, true).real();
    };

    /**
     * @brief Imaginary part of the conjugate integrand for QCDF Correction
     * @param u integration variable in the range [0, 1]
     * @return the imaginary part of the conjugate integrand
     */
    double getT_perp_bar_integrand_imag(double u) {
        return T_perp(u, true).imag();
    };
    
    /**
     * @brief QCDF Correction from various BFS papers (hep-ph/0106067, hep-ph/0412400). Total, integrated, in the helicity basis.
     * @param conjugate a boolean to control conjugation
     * @return @f$ T_{-}^{QCDF} @f$
     */
    gslpp::complex T_QCDF_minus(bool conjugate);
    
private:
    QCD::meson meson;
    QCD::meson vectorM;
    const StandardModel& SM;
    BXqll myBXqll;
    double T_perp_real;
    double T_perp_imag;
    double T_perp_bar_real;
    double T_perp_bar_imag;
    
    double average;/**< GSL integral variable */  
    double error;/**< GSL integral variable */    
    gsl_function f_GSL;/**< GSL integral variable */
    gsl_integration_cquad_workspace * w_GSL;/**< GSL integral variable */
    
    std::vector<std::string> mVgammaParameters;/**< The string of mandatory MVgamma parameters */
};



/**
 * @class BR_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$BR@f$ in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$BR@f$ in @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * BR = \frac {\alpha_e G_F^2 M_b^2 M_M \lambda}{(4\pi)^2 4 w_M} ( |H_V^+|^2 + |H_A^+|^2 +|\overline{H}_V^-|^2 + |\overline{H}_A^-|^2) \,.
 * @f]
 */
class BR_MVgamma : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    BR_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);
    
    /**
    * @brief The @f$BR@f$ in @f$M \to V \gamma@f$.
    * @return @f$BR@f$
    */
    double computeThValue ();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
    AmpDB2 myAmpDB2;
    double arg;
    double ADG; /**< @f$A_{\Delta\Gamma}@f$ */
    double ys; /** @f$\frac{\Delta\Gamma}{\Gamma}@f$ */
    double t_int; /** The factor that comes into CP averaged measurements due to finite lifetime differences. */
};



/**
 * @class ACP_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$C@f$ parameter of CPV in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$C@f$ parameter of CPV in @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * C = \frac {( |H_V^+|^2 + |H_V^-|^2 - |\overline{H}_V^+|^2 - |\overline{H}_V^-|^2)}{( |H_V^+|^2 + |H_V^+|^2 +|\overline{H}_V^-|^2 + |\overline{H}_V^-|^2)} \,.
 * @f]
 */
class C_MVgamma : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    C_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);
    
    /**
    * @brief The @f$C@f$ parameter of CPV in @f$M \to V \gamma@f$.
    * @return @f$C@f$
    */
    double computeThValue ();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class S_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$S@f$ parameter for CPV in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$S@f$ parameter of CPV in @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * S = \frac {2Im(q/p * ( H_V^{+*} \overline{H}_V^+ + H_V^{-*} \overline{H}_V^-))}{( |H_V^+|^2 + |H_V^-|^2 +|\overline{H}_V^+|^2 + |\overline{H}_V^-|^2)} \,.
 * @f]
 */
class S_MVgamma : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    S_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);
    
    /**
    * @brief The @f$S@f$ parameter for CPV in @f$M \to V \gamma@f$.
    * @return @f$S@f$
    */
    double computeThValue ();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
    AmpDB2 myAmpDB2;
    double arg;
};

/**
 * @class ADG_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$A_{\Delta\Gamma}@f$ parameter for @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$A_{\Delta\Gamma}@f$ parameter for @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * A_{\Delta\Gamma} = \frac {2Re(q/p * ( H_V^{+*} \overline{H}_V^+ + H_V^{-*} \overline{H}_V^-))}{( |H_V^+|^2 + |H_V^-|^2 +|\overline{H}_V^+|^2 + |\overline{H}_V^-|^2)} \,.
 * @f]
 */
class ADG_MVgamma : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ADG_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);
    
    /**
    * @brief The @f$S@f$ parameter for CPV in @f$M \to V \gamma@f$.
    * @return @f$S@f$
    */
    double computeThValue ();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
    AmpDB2 myAmpDB2;
    double arg;
};


/**
 * @class DC7_1
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^1@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^1 = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_- - h_+| \,.
 * @f]
 */
class DC7_1 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    DC7_1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^1@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class DC7_2
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^2@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^2 = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_- + h_+| \,.
 * @f]
 */
class DC7_2 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    DC7_2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class hp0_hm0
 * @ingroup Flavour
 * @brief A class for the absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in @f$B \to K^*@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in 
 * @f$B \to K^*@f$
 */
class hp0_hm0 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    hp0_hm0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The absolute value of the ratio @f$h_+^{(0)}/h_-^{(0)}@f$ in @f$B \to K^*@f$.
    * @return @f$h_+^{(0)}/h_-^{(0)}@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class AbsDC7_L
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^L@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^L = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_-| \,.
 * @f]
 */
class AbsDC7_L : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    AbsDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^L@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^L@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class AbsDC7_R
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^R@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^R = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_+| \,.
 * @f]
 */
class AbsDC7_R : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    AbsDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class ReDC7_L
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^L@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^L = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_-| \,.
 * @f]
 */
class ReDC7_L : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ReDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^L@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^L@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class ReDC7_R
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^R@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^R = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_+| \,.
 * @f]
 */
class ReDC7_R : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ReDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class ImDC7_L
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^L@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^L = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_-| \,.
 * @f]
 */
class ImDC7_L : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ImDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^L@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^L@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class ImDC7_R
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^R@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^R = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_+| \,.
 * @f]
 */
class ImDC7_R : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ImDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class AbsDC7_QCDF
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^R@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^R = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_+| \,.
 * @f]
 */
class AbsDC7_QCDF : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    AbsDC7_QCDF(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};

/**
 * @class AbsDC7_QCDF_bar
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^R@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^R = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_+| \,.
 * @f]
 */
class AbsDC7_QCDF_bar : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    AbsDC7_QCDF_bar(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    QCD::meson meson; /**< Initial meson type. */
    QCD::meson vectorM; /**< Final vector meson type. */
};
#endif	/* MVLL_H */

