/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVLL_H
#define	MVLL_H

#include "StandardModel.h"
#include "ThObservable.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>
#include <gsl/gsl_monte_plain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResultPtr.h>

#define SWITCH 8.2

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
 * @ingroup Flavour
 * @brief A class for the @f$M \to V l^+ l^-@f$ decay.  
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * build the observables relative to the @f$M \to V l^+ l^-@f$ decays, where
 * @f$M@f$ is a generic meson and @f$V@f$ is a vector meson. This kind of decays can be described
 * by means of the @f$\Delta B = 1 @f$ weak effective Hamiltonian
 * @f[
 *   \mathcal{H}_\mathrm{eff}^{\Delta B = 1} = \mathcal{H}_\mathrm{eff}^\mathrm{had} +
 *   \mathcal{H}_\mathrm{eff}^\mathrm{sl+\gamma},
 * @f]  
 * where the first term is the hadronic contribution 
 * @f[
 * \mathcal{H}_\mathrm{eff}^\mathrm{had} = \frac{4G_F}{\sqrt{2}}\Bigg[\sum_{p=u,c}\lambda_p\bigg(C_1 Q^{p}_1 
 * + C_2 Q^{p}_2\bigg) -\lambda_t \bigg(\sum_{i=3}^{6} C_i P_i + C_{8}Q_{8g} \bigg)\Bigg] \,,
 * @f]
 * involving current-current, chromodynamic penguin and chromomagnetic dipole operators, while the second one, given by
 * @f[
 * \mathcal{H}_\mathrm{eff}^\mathrm{sl+\gamma} = - \frac{4G_F}{\sqrt{2}}\lambda_t
 * \bigg( C_7Q_{7\gamma} + C_9Q_{9V} + C_{10}Q_{10A} \bigg) \,, 
 * @f]
 * includes the electromagnetic penguin plus the semileptonic operators.
 * 
 * Considering the matrix element of @f$\mathcal{H}_\mathrm{eff}^{\Delta B = 1}@f$
 * between the initial state @f$M@f$ and the final state @f$V l^+ l^-@f$, only the contribution of 
 * @f$\mathcal{H}_\mathrm{eff}^\mathrm{sl+\gamma}@f$ clearly factorizes into the 
 * product of hadronic form factors and leptonic tensors at all orders in strong interactions. 
 * Following @cite Jager:2012uw, we implemented the amplitude in the helicity basis; hence we made use of the helicity
 * form factors @f$ \tilde{V}_\lambda(q^2), \tilde{T}_\lambda(q^2)@f$ and @f$\tilde{S}(q^2) @f$ 
 * (where @f$\lambda=+,-,0@f$ represents the helicity), which are related to the
 * ones in the transverse basis through the following relations :
 * @f[
 * \tilde{V}_0(q^2) = \frac{4m_V}{\sqrt{q^2}}A_{12}(q^2)\,,\\
 * \tilde{V}_{\pm}\left( q^{2}\right) = \frac{1}{2} \bigg[ \Big( 1 + \frac{m_V}{m_M} \Big) A_1\left( q^{2}\right) 
 * \mp \frac{\lambda^{1/2}(q^2)}{m_M(m_M + m_V)} V\left( q^{2}\right) \bigg]\,, \\
 * \tilde{T}_0(q^2)=\frac{2\sqrt{q^2}m_V}{m_M(m_M + m_V)}T_{23}(q^2)\,,\\
 * \tilde{T}_{\pm}\left( q^{2}\right) = \frac{m_M^2 - m_V^2}{2m_M^2}T_2\left( q^{2}\right) 
 * \mp \frac{\lambda^{1/2}(q^2)}{2m_M^2}T_1\left( q^{2}\right)\,,\\
 * \tilde{S}\left( q^{2}\right) = -\frac{\lambda^{1/2}(q^2)}{2m_M(m_b+m_s)}A_0\left( q^{2}\right)\,,
 * @f]
 * where @f$\lambda(q^2) = 4m_M^2|\vec{k}|^2@f$, with @f$\vec{k}@f$ as the 3-momentum
 * of the meson @f$V@f$ in the @f$M@f$ rest frame.
 * 
 * The effect of the operators of @f$\mathcal{H}_\mathrm{eff}^\mathrm{had}@f$ due to
 * exchange of soft gluon can be reabsorbed in the following parameterization,
 * @f[
 * h_\lambda(q^2) = \frac{\epsilon^*_\mu(\lambda)}{m_M^2} 
 * \int d^4x e^{iqx} \langle \bar V \vert T\{j^{\mu}_\mathrm{em} (x) 
 * \mathcal{H}_\mathrm{eff}^\mathrm{had} (0)\} \vert \bar M \rangle = 
 * h_\lambda^{(0)} + \frac{q^2}{1\,\mathrm{GeV}^2} h_\lambda^{(1)} + \frac{q^4}{1\, \mathrm{GeV}^4} h_\lambda^{(2)} \,,
 * @f]
 * while the effect due to exchange of hard gluons can be parametrized following 
 * the prescription of @cite Beneke:2001at as a shift to the Wilson coefficient @f$C_9@f$ :
 * one first have to define the corrections
 * @f[
 * \Delta \mathcal{T}_a = \frac{\alpha_sC_F}{4\pi} C_a + \frac{\alpha_sC_F}{4}\frac{\pi}{N_c}\frac{f_Mf_{V,a}}{m_V F_a(q^2)}\Xi_a
 * \sum_{\pm}\int \frac{d\omega}{\omega}\Phi_{V,\pm}(\omega)\int_0^1du\Phi_{M,a}(u)T_{a,\pm}(u,\omega)\,,
 * @f]
 * where @f$a=\perp,\parallel@f$, @f$F_\perp(q^2) = T_1(q^2) @f$, @f$F_\parallel(q^2) = T_1(q^2) - T_3(q^2)@f$, 
 * @f$\Xi_\perp(q^2) = 1 @f$, @f$\Xi_\parallel(q^2) = \frac{2m_Vm_M}{m_M^2-q^2}@f$,
 * and @f$\Phi_X@f$ are leading twist light-cone distributions; the term proportional
 * to @f$C_a@f$ is the one describing the corrections where the spectator quark is
 * connected to the hard process only through soft interactions, while the one 
 * proportional to @f$T_{a,\pm}@f$ is the one describing the corrections where
 * the spectactor quark is involved in the hard process. Therefore, it is possible
 * to define the correction to the Wilson coefficient in the following way:
 * @f[
 * \Delta C_{9,\pm} = \frac{1}{q^2}\frac{m_b}{m_M} \left((m_M^2-m_V^2) \frac{m_M^2 - q^2}{m_M^2} 
 * \mp \sqrt{\lambda(q^2)}\right) \Delta T_{\perp}(q^2)\,,\\
 * \Delta C_{9,0} = \frac{1}{ 2 m_V m_M \sqrt{q^2} } \left(\left[(m_M^2-m_V^2) ( m_M^2-m_V^2 - q^2) - \lambda(q^2)\right]
 * (m_M^2 - q^2) \frac{m_b}{m_M^2q^2} \Delta T_{\perp}(q^2) - \lambda(q^2) 
 * \frac{m_b}{m_M^2-m_V^2}\left(\Delta T_{\parallel}(q^2) + \Delta T_\perp(q^2)\right)\right)\,.
 * @f]
 * 
 * The amplitude can be therefore parametrized in terms of the following helicity amplitudes:
 * @f[
 * H_V(\lambda) = -i\, N \Big\{C_{9} \tilde{V}_{L\lambda} +C_{9}'  \tilde{V}_{R\lambda}
 * + \frac{m_M^2}{q^2} \Big[\frac{2\, m_b}{m_M} (C_{7} \tilde{T}_{L\lambda} +  C_{7}'  \tilde{T}_{R\lambda})
            - 16 \pi^2 h_\lambda \Big] \Big\} \,,  \\
 * H_A(\lambda) = -i\, N (C_{10}  \tilde{V}_{L\lambda} + C_{10}'\tilde{V}_{R\lambda}) \,, \\
 * H_S = i\, N \frac{ m_b}{m_W} (C_S \tilde{S}_L + C_S' \tilde{S}_R)\,, \\
 * H_P = i\, N \Big\{ \frac{ m_b}{m_W} (C_P \tilde{S}_L + C_P' \tilde{S}_R)
 * + \frac{2\,m_\ell m_b}{q^2} \left[C_{10} \Big(\tilde{S}_L - \frac{m_s}{m_b} \tilde{S}_R \Big) 
 * + C_{10}' \Big(\tilde{S}_R - \frac{m_s}{m_b} \tilde{S}_L\Big) \right] \Big\} \,,
 * @f]
 * where @f$ N = - \frac{4 G_F m_M}{\sqrt{2}}\frac{e^2}{16\pi^2}\lambda_t@f$ and we have defined
 * @f[
 * \tilde{V}_{L\pm}(q^2) = -\tilde{V}_{R\mp}(q^2)=\tilde{V}_\pm(q^2)\,,\\
 * \tilde{T}_{L\pm}(q^2) = -\tilde{T}_{R\mp}(q^2)=\tilde{T}_\pm(q^2)\,,\\
 * \tilde{S}_L(q^2) = -\tilde{S}_R(q^2)=\tilde{S}(q^2)\,.
 * @f]
 * Squaring the amplitude and summing over the spins it is possible to obtain 
 * the fully differential decay rate, which is
 * @f[
 * \frac{d^{(4)} \Gamma}{dq^2\,d(\cos\theta_l)d(\cos\theta_k)d\phi} = \frac{9}{32\,\pi} 
 * \Big( I^s_1\sin^2\theta_k+I^c_1\cos^2\theta_k +(I^s_2\sin^2\theta_k+I^c_2\cos^2\theta_k)\cos2\theta_l \\
 * + I_3\sin^2\theta_k\sin^2\theta_l\cos2\phi +I_4\sin2\theta_k\sin2\theta_l\cos\phi 
 * + I_5\sin2\theta_k\sin\theta_l\cos\phi \\
 * +(I_6^s\sin^2\theta_k + I_6^c \cos^2\theta_K) \cos\theta_l
 * + I_7\sin2\theta_k\sin\theta_l\sin\phi+I_8\sin2\theta_k\sin2\theta_l\sin\phi +I_9\sin^2\theta_k\sin^2\theta_l\sin2\phi \Big) 
 * @f]
 * The angular coefficients involved in the differential decay rate are related to the
 * helicity amplitudes according to the following relations:
 * @f[
 * I_1^c = F \left\{ \frac{1}{2}\left(|H_V^0|^2+|H_A^0|^2\right)+ 
 * |H_P|^2+\frac{2m_\ell^2}{q^2}\left(|H_V^0|^2-|H_A^0|^2\right) + \beta^2 |H_S|^2 \right\}\,,\\
 * I_1^s = F \left\{\frac{\beta^2\!+\!2}{8}\left(|H_V^+|^2+|H_V^-|^2+(V\rightarrow A)\right)
 * +\frac{m_\ell^2}{q^2}\left(|H_V^+|^2+|H_V^-|^2-(V\rightarrow A)\right)\right\}\,,\\
 * I_2^c = -F\, \frac{\beta^2}{2}\left(|H_V^0|^2+|H_A^0|^2\right)\,,\\
 * I_2^s = F\, \frac{\beta^2}{8}\left(|H_V^+|^2+|H_V^-|^2\right)+(V\rightarrow A)\,,\\
 * I_3 = -\frac{F}{2}{\rm Re} \left[H_V^+(H_V^-)^*\right]+(V\rightarrow A)\,,\\
 * I_4 = F\, \frac{\beta^2}{4}{\rm Re}\left[(H_V^-+H_V^+)\left(H_V^0\right)^*\right]+(V\rightarrow A)\,,\\
 * I_5 = F\left\{ \frac{\beta}{2}{\rm Re}\left[(H_V^--H_V^+)\left(H_A^0\right)^*\right]
 * +(V\leftrightarrow A) - \frac{\beta\,m_\ell}{\sqrt{q^2}} {\rm Re} 
 * \left[H_S^* (H_V^+ + H_V^-)\right]\right\}\,,\\
 * I_6^s = F \beta\,{\rm Re}\left[H_V^-(H_A^-)^*-H_V^+(H_A^+)^*\right]\,,\\
 * I_6^c = 2 F \frac{\beta\, m_\ell}{\sqrt{q^2}} {\rm Re} \left[ H_S^* H_V^0 \right]\,,\\
 * I_7 = F \left\{ \frac{\beta}{2}\,{\rm Im}\left[\left(H_A^++H_A^-\right)
 * (H_V^0)^* \, +(V\leftrightarrow A) \right] - \frac{\beta\, m_\ell}{\sqrt{q^2} }\, {\rm Im} 
 * \left[ H_S^*(H_V^{-} -  H_V^{+}) \right] \right\}\,,\\
 * I_8 = F\, \frac{\beta^2}{4}{\rm Im}\left[(H_V^--H_V^+)(H_V^0)^*\right]+(V\rightarrow A)\,,\\
 * I_9 = F\, \frac{\beta^2}{2}{\rm Im}\left[H_V^+(H_V^-)^*\right]+(V\rightarrow A)\,,
 * @f]
 * where
 * @f[
 * F=\frac{ \lambda^{1/2}\beta\, q^2}{3 \times 2^{5} \,\pi^3\, m_M^3} BF(V \to {\rm final \, state})\,,
 * \qquad \beta = \sqrt{1 - \frac{4 m_\ell^2}{q^2} }\,.
 * @f]
 * The final observables are hence build employing CP-averages @f$\Sigma_i@f$ or CP-asymmetries @f$\Delta_i@f$ of
 * such angular coefficients; however, since on the experimental side the observables
 * are averaged over @f$ q^2 @f$ bins, an integration of the coeffiecients over such
 * bins has to be performed before they are combined in order to build the observables.
 * 
 * The class is organized as follows: after the parameters are updated in updateParameters() and the cache is checked in 
 * checkCache(), the form factor are build in the transverse basis in the functions
 * V(), A_0(), A_1(), A_2(), T_1() and  T_2() using the fit function FF_fit() from @cite Straub:2015ica .
 * The form factor are consequentely translated in the helicity basis through the
 * functions V_0t(), V_p(), V_m(), T_0t(), T_p(), T_m() and S_L() .
 * The basic elements required to compute the hard gluon corrections to the Wilson coefficient @f$C_9@f$
 * are build in the functions Tperpplus(), Tparplus(), Tparminus(), Cperp() and Cpar();
 * these corrections have to be integrated to be computed, so the final correction is
 * either obtaind through direct integration in the functions DeltaC9_p(), DeltaC9_m()
 * and DeltaC9_0(), or obtained through fitting in the functions fDeltaC9_p(), 
 * fDeltaC9_m() and fDeltaC9_0(). Form factors, Wilson coefficients and parameters 
 * are combined together in the functions H_V_0(), H_V_p(), H_V_m(), H_A_0(), 
 * H_A_p(), H_A_m(), H_S() and H_P() in order to build the helicity aplitudes, 
 * which are consequentely combined to create the angular coefficients in the 
 * function I_1c(), I_1s(), I_2c(), I_2s(), I_3(), I_4(), I_5(), I_6c(), I_6s(), 
 * I_7(), I_8(), I_9(). Those coefficients are used to create the CP averaged 
 * coefficients in the functions getSigma1c(), getSigma1s(), getSigma2c(), getSigma2s(), 
 * getSigma3(), getSigma4(), getSigma5(), getSigma6c(), getSigma6s(), getSigma7(), 
 * getSigma8(), getSigma9(), and the CP asymmetric coefficients in the function 
 * Delta(). The CP averaged and asymmetric coefficients are integrated over the 
 * @f$q^2@f$ bin in the functions integrateSigma() and integrateDelta(), in order
 * to be further used to build the observables.
 */
class MVll {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~MVll();  
    
    /**
    * @brief The integral of \f$ \Sigma_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Sigma_{i}> \f$ 
    */
    double integrateSigma(int i, double q_min, double q_max);
    
    /**
    * @brief The value of \f$ \Sigma_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] \f$ q^2 \f$ value of the function
    * @return \f$ <\Sigma_{i}> \f$ 
    */
    double getSigma(int i, double q_2);
    
    /**
    * @brief The integral of \f$ \Delta_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Delta_{i}> \f$ 
    */
    double integrateDelta(int i, double q_min, double q_max);
    
    /**
    * @brief The width of the meson M
    * @return \f$ \Gamma_M \f$ 
    */
    double getwidth(){
        updateParameters();
        return width;
    }
    
    /**
    * @brief The mass of the lepton l
    * @return \f$ m_l \f$ 
    */
    double getMlep(){
        updateParameters();
        return Mlep;
    }
    
    /**
    * @brief The factor \f$ \beta \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta \f$
    */
    double beta (double q2);
    
    /**
    * @brief The form factor \f$ V_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_0 \f$
    */
    double getV0(double q2)
    {
        updateParameters();
        return (2. * MM * sqrt(q2))/sqrt(lambda(q2)) * V_0t(q2);
    };
    
    /**
    * @brief The form factor \f$ V_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_+ \f$
    */
    double getVp(double q2)
    {
        updateParameters();
        return V_p(q2);
    };
    
    /**
    * @brief The form factor \f$ V_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_- \f$
    */
    double getVm(double q2)
    {
        return V_m(q2);
    };
    
    /**
    * @brief The form factor \f$ T_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_0 \f$
    */
    double getT0(double q2)
    {
        updateParameters();
        return twoMM3/sqrt(q2 * lambda(q2)) * T_0t(q2);
    };
    
    /**
    * @brief The form factor \f$ T_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_+ \f$
    */
    double getTp(double q2)
    {
        updateParameters();
        return T_p(q2);
    };
    
    /**
    * @brief The form factor \f$ T_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_- \f$
    */
    double getTm(double q2)
    {
        updateParameters();
        return T_m(q2);
    };
    
    /**
    * @brief The form factor \f$ S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S \f$
    */
    double getS(double q2)
    {
        updateParameters();
        return S_L_pre/sqrt(lambda(q2)) * S_L(q2);
    };
    
    /**
     * @brief The non-pertubative ccbar contributions to the helicity amplitudes
     * @param hel helicity
     * @param q2 \f$q^2\f$
     * @return \f$h_{hel}(q^2)\f$
     */
    gslpp::complex h_lambda(int hel, double q2);
    
    /**
    * @brief The helicity amplitude \f$ H_V^0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^0 \f$
    */
    gslpp::complex H_V_0(double q2, bool bar);
    
    /**
    * @brief The helicity amplitude \f$ H_V^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^+ \f$
    */
    gslpp::complex H_V_p(double q2, bool bar);
    
    /**
    * @brief The helicity amplitude \f$ H_V^- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^- \f$
    */
    gslpp::complex H_V_m(double q2, bool bar);

    /**
    * @brief The helicity amplitude \f$ H_A^0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^0 \f$
    */
    gslpp::complex H_A_0(double q2, bool bar);

    /**
    * @brief The helicity amplitude \f$ H_A^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^+ \f$
    */
    gslpp::complex H_A_p(double q2, bool bar);

    /**
    * @brief The helicity amplitude \f$ H_A^- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^- \f$
    */
    gslpp::complex H_A_m(double q2, bool bar);

    /**
    * @brief The helicity amplitude \f$ H_S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_S \f$
    */
    gslpp::complex H_S(double q2, bool bar);

    /**
    * @brief The helicity amplitude \f$ H_P \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_P \f$
    */
    gslpp::complex H_P(double q2, bool bar);
    
    /**
    * @brief The real part of \f$ \tilde{g}^1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^1) \f$ 
    */
    double getgtilde_1_re(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_1_pre/(sqrt(lambda(q2)) * V(q2)) * (h_lambda(2,q2)-h_lambda(1,q2))).real()/q2;
    }
    
    /**
    * @brief The immaginary part of \f$ \tilde{g}^1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^1) \f$ 
    */
    double getgtilde_1_im(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_1_pre/(sqrt(lambda(q2)) * V(q2)) * (h_lambda(2,q2)-h_lambda(1,q2))).imag()/q2;
    }
    
    double getQCDf_1(double q2)
    {
        /* NOTE THE ZEROS */
        return (gtilde_1_pre/(sqrt(lambda(q2)) * V(q2)) * 1./(16. * M_PI * M_PI * MM*MM) * ((0.*C_9 + fDeltaC9_m(q2) + 0.*Y(q2)) * V_m(q2) - (C_9 + fDeltaC9_p(q2) + Y(q2)) * V_p(q2))).abs();
    }
    
    double getQCDf_2(double q2)
    {
        /* NOTE THE ZEROS */
        return (gtilde_2_pre/A_1(q2) * 1./(16. * M_PI * M_PI * MM*MM) * ((0.*C_9 + fDeltaC9_m(q2) + 0.*Y(q2)) * V_m(q2) + (0.*C_9 + fDeltaC9_p(q2) + 0.*Y(q2)) * V_p(q2))).abs();
    }
    
    double getQCDf_3(double q2)
    {
        /* NOTE THE ZEROS */
        return (gtilde_3_pre/(lambda(q2) * A_2(q2)) * (sqrt(q2) * 1./(16. * M_PI * M_PI * MM*MM) * (0.*C_9 + fDeltaC9_0(q2) + 0.*Y(q2)) * V_0t(q2) - (MM2mMV2 - q2)/(4.*MV) * 1./(16. * M_PI * M_PI * MM*MM) * ((0.*C_9 + fDeltaC9_m(q2) + 0.*Y(q2)) * V_m(q2) + (0.*C_9 + fDeltaC9_p(q2) + 0.*Y(q2)) * V_p(q2)))).abs();
    }
    
    /**
    * @brief The real part of \f$ \tilde{g}^2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^2) \f$ 
    */
    double getgtilde_2_re(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_2_pre/A_1(q2) * (h_lambda(1,q2)+h_lambda(2,q2))).real()/q2;
    }
    
    /**
    * @brief The immaginary part of \f$ \tilde{g}^2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^2) \f$ 
    */
    double getgtilde_2_im(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_2_pre/A_1(q2) * (h_lambda(1,q2)+h_lambda(2,q2))).imag()/q2;
    }
    
    /**
    * @brief The real part of \f$ \tilde{g}^3 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^3) \f$ 
    */
    double getgtilde_3_re(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_3_pre/(lambda(q2) * A_2(q2)) * (sqrt(q2)*h_lambda(0,q2)/q2-
                (MM2mMV2 - q2)/(4.*MV) * (h_lambda(1,q2)+h_lambda(2,q2))/q2)).real();
    }

    /**
    * @brief The imaginary part of \f$ \tilde{g}^3 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^3) \f$ 
    */
    double getgtilde_3_im(double q2)
    {
        updateParameters();
        return C2_inv * (gtilde_3_pre/(lambda(q2) * A_2(q2)) * (sqrt(q2)*h_lambda(0,q2)/q2-
                (MM2mMV2 - q2)/(4.*MV) * (h_lambda(1,q2)+h_lambda(2,q2))/q2)).imag();
    }
    
    /**
    * @brief The real part of \f$ h_0 \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_0) \f$
    */
    double geth_0_re(double q2)
    {
        return (sixteenM_PI2MM2 * h_lambda(0,q2)/q2).real();
    }
    
    /**
    * @brief The imaginary part of \f$ h_0 \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_0) \f$
    */
    double geth_0_im(double q2)
    {
        return (sixteenM_PI2MM2 * h_lambda(0,q2)/q2).imag();
    }
    
    /**
    * @brief The real part of \f$ h_+ \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_+) \f$
    */
    double geth_p_re(double q2)
    {
        return (sixteenM_PI2MM2 * h_lambda(1,q2)/q2).real();
    }
    
    /**
    * @brief The imaginary part of \f$ h_+ \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_+) \f$
    */
    double geth_p_im(double q2)
    {
        return (sixteenM_PI2MM2 * h_lambda(1,q2)/q2).imag();
    }
    
    /**
    * @brief The real part of \f$ h_- \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_-) \f$
    */
    double geth_m_re(double q2)
    {
        return (sixteenM_PI2MM2 * h_lambda(2,q2)/q2).real();
    }

    /**
    * @brief The imaginary part of \f$ h_- \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_-) \f$
    */
    double geth_m_im(double q2)
    {
        return (sixteenM_PI2MM2 *  h_lambda(2,q2)/q2).imag();
    }

    /**
    * @brief The absolute value of the ratio \f$ h_+^{(0)}/h_0^{(0)} \f$.  
    * @return \f$ h_+^{(0)}/h_0^{(0)} \f$
    */
    double gethp0_hm0_abs()
    {
        updateParameters();
        return (h_0[1]/h_0[2]).abs();
    }
    
    /**
    * @brief The absolute value of the ratio \f$ h_-^{(0)}/h_0^{(0)} \f$.
    * @return \f$ h_-^{(0)}/h_0^{(0)} \f$
    */
    double gethm0_h00_abs()
    {
        updateParameters();
        return (h_0[2]/h_0[0]).abs();
    }
    
private:
    const StandardModel& mySM;/**< Model type */
    StandardModel::lepton lep;/**< Final leptons type */
    StandardModel::meson meson;/**< Initial meson type */
    StandardModel::meson vectorM;/**< Final vector meson type */
    
    double GF;            /**<Fermi constant */
    double ale;           /**<Alpha electromagnetic */
    double Mlep;          /**<Muon mass */
    double MM;            /**<Initial meson mass */
    double MV;            /**<Final vector meson mass */
    double Mb;            /**<b quark mass */
    double mu_b;          /**<b mass scale */
    double mu_h;          /**<\f$\sqrt{\mu_b*\lambda_{QCD}}\f$ */
    double Mc;            /**<c quark mass */
    double Ms;            /**<s quark mass */
    double spectator_charge;  /**<spectator quark charge */
    double width;         /**<Initial meson width */
    double fperp;         /**<vector meson perpendicular decay constant*/
    double ys;            /**<CP-violation factor \f$\frac{\Delta \Gamma}{2\Gamma}\f$*/
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor */
    double b;             /**<BF of the decay V -> final states */
    gslpp::complex h_0[3];         /**<Parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_1[3];         /**<Parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_2[3];         /**<Parameter that contains the contribution from the hadronic hamiltonian */
    
    double t_p;/**< Cache variable */
    double t_m;/**< Cache variable */
    double t_0;/**< Cache variable */
    double z_0;/**< Cache variable */
    double MMpMV;/**< Cache variable */
    double MMpMV2;/**< Cache variable */
    double MMmMV;/**< Cache variable */
    double MMmMV2;/**< Cache variable */
    double MM2;/**< Cache variable */
    double MM4;/**< Cache variable */
    double MV2;/**< Cache variable */
    double MV4;/**< Cache variable */
    double MMMV;/**< Cache variable */
    double MM2mMV2;/**< Cache variable */
    double fourMV;/**< Cache variable */
    double onepMMoMV;/**< Cache variable */
    double MM_MMpMV;/**< Cache variable */
    double twoMM2;/**< Cache variable */
    double twoMV2;/**< Cache variable */
    double twoMM_mbpms;/**< Cache variable */
    double fourMM2;/**< Cache variable */
    double Mlep2;/**< Cache variable */
    double twoMlepMb;/**< Cache variable */
    double MboMW;/**< Cache variable */
    double MsoMb;/**< Cache variable */
    double M_PI2osix;/**< Cache variable */
    double twoMM;/**< Cache variable */
    double m4downcharge;/**< Cache variable */
    double threeGegen0;/**< Cache variable */
    double threeGegen1otwo;/**< Cache variable */
    double twoMc2;/**< Cache variable */
    double ninetysixM_PI3MM3;/**< Cache variable */
    double sixteenM_PI2;/**< Cache variable */
    double sixteenM_PI2MM2;/**< Cache variable */
    double twoMboMM;/**< Cache variable */
    gslpp::complex H_0_pre;/**< Cache variable */
    double mu_b2;/**< Cache variable */
    double Mc2;/**< Cache variable */
    double Mb2;/**< Cache variable */
    double fourMc2;/**< Cache variable */
    double fourMb2;/**< Cache variable */
    double logMc;/**< Cache variable */
    double logMb;/**< Cache variable */
    gslpp::complex H_0_WC;/**< Cache variable */
    gslpp::complex H_c_WC;/**< Cache variable */
    gslpp::complex H_b_WC;/**< Cache variable */
    double fournineth;/**< Cache variable */
    double half;/**< Cache variable */
    double twothird;/**< Cache variable */
    gslpp::complex ihalfMPI;/**< Cache variable */
    double twoMM3;/**< Cache variable */
    double gtilde_1_pre;/**< Cache variable */
    double gtilde_2_pre;/**< Cache variable */
    double gtilde_3_pre;/**< Cache variable */
    double C2_inv;/**< Cache variable */
    double S_L_pre;/**< Cache variable */
    gslpp::complex NN;/**< coupling including the CKM element */
    gslpp::complex NN_conjugate;/**< conjugate of the coupling including the CKM element */
    double sixMMoMb;/**< Cache variable */
    double CF;/**< Cache variable */
    double deltaT_0;/**< Cache variable */
    double deltaT_1par;/**< Cache variable */
    double deltaT_1perp;/**< Cache variable */
    double Ee;/**< Cache variable */
    
    bool h_pole;
    
    gslpp::complex ubar;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex arg1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex B01;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex B00;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex xp;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex xm;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex yp;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex ym;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex L1xp;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex L1xm;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex L1yp;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex L1ym;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F87_0;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F87_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F87_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F87_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_0;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_L1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_L1_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_L1_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F29_L1_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_0;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_L1_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_L1_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    gslpp::complex F27_L1_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    double F89_0;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    double F89_1;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    double F89_2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    double F89_3;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    
    double a_0V;/**<LCSR fit parameter */
    double a_1V;/**<LCSR fit parameter */
    double a_2V;/**<LCSR fit parameter */
    double MRV_2;/**<LCSR fit parameter */
    double a_0A0;/**<LCSR fit parameter */
    double a_1A0;/**<LCSR fit parameter */
    double a_2A0;/**<LCSR fit parameter */
    double MRA0_2;/**<LCSR fit parameter */
    double a_0A1;/**<LCSR fit parameter */
    double a_1A1;/**<LCSR fit parameter */
    double a_2A1;/**<LCSR fit parameter */
    double MRA1_2;/**<LCSR fit parameter */
    double a_0A12;/**<LCSR fit parameter */
    double a_1A12;/**<LCSR fit parameter */
    double a_2A12;/**<LCSR fit parameter */
    double MRA12_2;/**<LCSR fit parameter */
    double a_0T1;/**<LCSR fit parameter */
    double a_1T1;/**<LCSR fit parameter */
    double a_2T1;/**<LCSR fit parameter */
    double MRT1_2;/**<LCSR fit parameter */
    double a_0T2;/**<LCSR fit parameter */
    double a_1T2;/**<LCSR fit parameter */
    double a_2T2;/**<LCSR fit parameter */
    double MRT2_2;/**<LCSR fit parameter */
    double a_0T23;/**<LCSR fit parameter */
    double a_1T23;/**<LCSR fit parameter */
    double a_2T23;/**<LCSR fit parameter */
    double MRT23_2;/**<LCSR fit parameter */

    gslpp::vector<gslpp::complex> ** allcoeff;/**<Vector that contains the Wilson coeffients */
    gslpp::vector<gslpp::complex> ** allcoeffh;/**<Vector that contains the Wilson coeffients at scale @f$\mu_h@f$ */
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<Vector that contains the primed Wilson coeffients */
    
    gslpp::complex C_1;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_1L_bar;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_1Lh_bar;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_2L_bar;/**<Leading order Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_2Lh_bar;/**<Leading order Wilson coeffients @f$C_2@f$ at scale @f$\mu_h@f$*/
    gslpp::complex C_3;/**<Wilson coeffients @f$C_3@f$*/
    gslpp::complex C_4;/**<Wilson coeffients @f$C_4@f$*/
    gslpp::complex C_5;/**<Wilson coeffients @f$C_5@f$*/
    gslpp::complex C_6;/**<Wilson coeffients @f$C_6@f$*/
    gslpp::complex C_7;/**<Wilson coeffients @f$C_7@f$*/
    gslpp::complex C_8L;/**<Leading order Wilson coeffients @f$C_8@f$*/
    gslpp::complex C_8Lh;/**<Leading order Wilson coeffients @f$C_8@f$ at scale @f$\mu_h@f$*/
    gslpp::complex C_9;/**<Wilson coeffients @f$C_9@f$*/
    gslpp::complex C_10;/**<Wilson coeffients @f$C_{10}@f$*/
    gslpp::complex C_S;/**<Wilson coeffients @f$C_S@f$*/
    gslpp::complex C_P;/**<Wilson coeffients @f$C_P@f$*/
    
    gslpp::complex C_7p;/**<Wilson coeffients @f$C_7'@f$*/
    gslpp::complex C_9p;/**<Wilson coeffients @f$C_9'@f$*/
    gslpp::complex C_10p;/**<Wilson coeffients @f$C_{10}'@f$*/
    gslpp::complex C_Sp;/**<Wilson coeffients @f$C_S'@f$*/
    gslpp::complex C_Pp;/**<Wilson coeffients @f$C_P'@f$*/
    
    std::vector<double> ReDeltaC9_p_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_p_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ReDeltaC9_m_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_m_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ReDeltaC9_0_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_0_mumu;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ReDeltaC9_p_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_p_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ReDeltaC9_m_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_m_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ReDeltaC9_0_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9_0_ee;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> myq2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    
    TFitResultPtr refres_p_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_p_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr refres_m_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_m_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr refres_0_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_0_mumu;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    
    TFitResultPtr refres_p_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_p_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr refres_m_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_m_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr refres_0_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres_0_ee;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
        
    TGraph gr1;/**<Tgraph to be used for fitting the QCDF @f$\Delta C_9@f$ */
    TGraph gr2;/**<Tgraph to be used for fitting the QCDF @f$\Delta C_9@f$ */
    
    TF1 reffit;/**<TF1 to be used for fitting the QCDF @f$\Delta C_9@f$ */
    TF1 imffit;/**<TF1 to be used for fitting the QCDF @f$\Delta C_9@f$ */
    
    double tmpq2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    
    double avaSigma;/**< Gsl integral variable */
    double avaDelta;/**< Gsl integral variable */
    double avaDTPPR;/**< Gsl integral variable */    
    
    double errSigma;/**< Gsl integral variable */
    double errDelta;/**< Gsl integral variable */
    double errDTPPR;/**< Gsl integral variable */
    
    gsl_function FS;/**< Gsl integral variable */
    gsl_function FD;/**< Gsl integral variable */
    gsl_function DTPPR;/**< Gsl integral variable */
    
    gsl_integration_cquad_workspace * w_DTPPR;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_sigma;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_delta;/**< Gsl integral variable */
    
    gsl_error_handler_t * old_handler; /**< GSL error handler store */
    
    std::map<std::pair<double, double>, gslpp::complex > cacheI1;/**< Cache variable */
    
    std::map<std::pair<double, double>, double > cacheSigma0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma1;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma2;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma3;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma4;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma5;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma6;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma7;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma9;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma10;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma11;/**< Cache variable */
    
    std::map<std::pair<double, double>, double > cacheDelta0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta1;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta2;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta3;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta7;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta11;/**< Cache variable */
    
    unsigned int N_updated;/**< Cache variable */
    gslpp::vector<double> N_cache;/**< Cache variable */
    gslpp::complex Nc_cache;/**< Cache variable */
    
    unsigned int V_updated;/**< Cache variable */
    gslpp::vector<double> V_cache;/**< Cache variable */
    
    unsigned int A0_updated;/**< Cache variable */
    gslpp::vector<double> A0_cache;/**< Cache variable */
    
    unsigned int A1_updated;/**< Cache variable */
    gslpp::vector<double> A1_cache;/**< Cache variable */
    
    unsigned int T1_updated;/**< Cache variable */
    gslpp::vector<double> T1_cache;/**< Cache variable */
    
    unsigned int T2_updated;/**< Cache variable */
    gslpp::vector<double> T2_cache;/**< Cache variable */
    
    unsigned int k2_updated;/**< Cache variable */
    gslpp::vector<double> k2_cache;/**< Cache variable */
    
    unsigned int z_updated;/**< Cache variable */
    
    unsigned int lambda_updated;/**< Cache variable */
    
    unsigned int beta_updated;/**< Cache variable */
    double beta_cache;/**< Cache variable */
    
    unsigned int F_updated;/**< Cache variable */
    
    unsigned int VL1_updated;/**< Cache variable */
    unsigned int VL2_updated;/**< Cache variable */
    
    unsigned int TL1_updated;/**< Cache variable */
    unsigned int TL2_updated;/**< Cache variable */

    unsigned int VR1_updated;/**< Cache variable */
    unsigned int VR2_updated;/**< Cache variable */
    
    unsigned int TR1_updated;/**< Cache variable */
    unsigned int TR2_updated;/**< Cache variable */
    
    unsigned int VL0_updated;/**< Cache variable */
    gslpp::vector<double> VL0_cache;/**< Cache variable */
    
    unsigned int TL0_updated;/**< Cache variable */
    gslpp::vector<double> TL0_cache;/**< Cache variable */
    
    unsigned int VR0_updated;/**< Cache variable */
    
    unsigned int TR0_updated;/**< Cache variable */
    
    unsigned int Mb_Ms_updated;/**< Cache variable */
    
    unsigned int SL_updated;/**< Cache variable */
    gslpp::vector<double> SL_cache;/**< Cache variable */
    
    unsigned int SR_updated;/**< Cache variable */
    
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
    
    unsigned int C_2Lh_updated;/**< Cache variable */
    gslpp::complex C_2Lh_cache;/**< Cache variable */
    
    unsigned int C_8Lh_updated;/**< Cache variable */
    gslpp::complex C_8Lh_cache;/**< Cache variable */
    
    unsigned int Yupdated;/**< Cache variable */
    gslpp::vector<double> Ycache;/**< Cache variable */
    
    gslpp::complex h0Ccache[3];/**< Cache variable */
    gslpp::complex h1Ccache[3];/**< Cache variable */
    gslpp::complex h2Ccache[3];/**< Cache variable */
    
    unsigned int h0_updated;/**< Cache variable */
    unsigned int h1_updated;/**< Cache variable */
    unsigned int h2_updated;/**< Cache variable */
    
    unsigned int H_V0updated;/**< Cache variable */
    gslpp::vector<double> H_V0cache;/**< Cache variable */
    
    unsigned int H_V1updated;/**< Cache variable */
    gslpp::vector<double> H_V1cache;/**< Cache variable */
    
    unsigned int H_V2updated;/**< Cache variable */
    gslpp::vector<double> H_V2cache;/**< Cache variable */
    
    unsigned int H_A0updated;/**< Cache variable */
    unsigned int H_A1updated;/**< Cache variable */
    unsigned int H_A2updated;/**< Cache variable */
    
    unsigned int H_Supdated;/**< Cache variable */
    gslpp::vector<double> H_Scache;/**< Cache variable */
    
    unsigned int H_Pupdated;/**< Cache variable */
    gslpp::vector<double> H_Pcache;/**< Cache variable */
    
    unsigned int I0_updated;/**< Cache variable */
    unsigned int I1_updated;/**< Cache variable */
    unsigned int I2_updated;/**< Cache variable */
    unsigned int I3_updated;/**< Cache variable */
    unsigned int I4_updated;/**< Cache variable */
    unsigned int I5_updated;/**< Cache variable */
    unsigned int I6_updated;/**< Cache variable */
    unsigned int I7_updated;/**< Cache variable */
    unsigned int I8_updated;/**< Cache variable */
    unsigned int I9_updated;/**< Cache variable */
    unsigned int I10_updated;/**< Cache variable */
    unsigned int I11_updated;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > I1Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma1Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma3Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma4Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma5Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma6Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma7Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma9Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma10Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma11Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta1Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta2Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta3Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta7Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta11Cached;/**< Cache variable */
    
    std::map<double, unsigned int> deltaTparpCached;/**< Cache variable */
    std::map<double, unsigned int> deltaTparmCached;/**< Cache variable */
    std::map<double, unsigned int> deltaTperpCached;/**< Cache variable */
    
    std::map<double, gslpp::complex> cacheDeltaTparp;/**< Cache variable */
    std::map<double, gslpp::complex> cacheDeltaTparm;/**< Cache variable */
    std::map<double, gslpp::complex> cacheDeltaTperp;/**< Cache variable */
    
    unsigned int deltaTparpupdated;/**< Cache variable */
    unsigned int deltaTparmupdated;/**< Cache variable */
    unsigned int deltaTperpupdated;/**< Cache variable */
    
    unsigned int T_updated;/**< Cache variable */
    gslpp::vector<double> T_cache;/**< Cache variable */
    
    /**
     * @brief The update parameter method for MVll.
     */
    void updateParameters();
    
    /**
     * @brief The caching method for MVll.
     */
    void checkCache();
    
    /**
    * @brief The lattice parameter \f$ z \f$ from arXiv:1310.3722v3.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double z(double q2);
    
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
    double A_0(double q2);

    
    /**
    * @brief The transverse form factor \f$ A_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_1 \f$
    */
    double A_1(double q2);
    
    /**
    * @brief The transverse form factor \f$ A_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_2 \f$
    */
    double A_2(double q2);
    
    /**
    * @brief The transverse form factor \f$ T_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_1 \f$
    */
    double T_1(double q2);

    
    /**
    * @brief The transverse form factor \f$ T_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_2 \f$
    */
    double T_2(double q2);
    
    /**
    * @brief The helicity form factor \f$ \tilde{V}_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \tilde{V}_0 \f$
    */
    double V_0t(double q2);
    
    /**
    * @brief The helicity form factor \f$ V_+ \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_+ \f$
    */
    double V_p(double q2);
    
    /**
    * @brief The helicity form factor \f$ V_- \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_- \f$
    */
    double V_m(double q2);

    /**
    * @brief The helicity form factor \f$ \tilde{T}_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \tilde{T}_0 \f$
    */
    double T_0t(double q2);
    
    /**
    * @brief The helicity form factor \f$ T_+ \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_+ \f$
    */
    double T_p(double q2);
    
    /**
    * @brief The helicity form factor \f$ T_- \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_- \f$
    */
    double T_m(double q2);

    /**
    * @brief The helicity form factor \f$ S_L \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S_L \f$
    */
    double S_L(double q2);
    
    /**
    * @brief The \f$ h(q^2,0) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h(q^2,0) \f$
    */
    gslpp::complex H_0(double q2);
    
    /**
    * @brief The \f$ h(q^2,m_c) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] mu mass scale
    * @return \f$ h(q^2,m_c) \f$
    */
    gslpp::complex H_c(double q2, double mu);
    
    /**
    * @brief The \f$ h(q^2,m_b) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h(q^2,m_b) \f$
    */
    gslpp::complex H_b(double q2);
    
    /**
    * @brief The \f$ Y(q^2) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ Y(q^2) \f$
    */
    gslpp::complex Y(double q2);   
        
    /**
    * @brief The square of the 3-momentum of the recoiling meson in the M rest frame, \f$ k^2 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ k^2 \f$ 
    */
    double k2 (double q2);
        
    /**
    * @brief The factor \f$ \beta^2 \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta^2 \f$
    */
    double beta2 (double q2);
    
    /**
    * @brief The factor \f$ \lambda \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \lambda \f$
    */
    double lambda(double q2);

    /**
    * @brief The factor \f$ F \f$ used in the angular coefficients I_i. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] b_i BF of the decay \f$ V \to M_1 M_2\f$ 
    * @return \f$ F \f$
    */
    double F(double q2, double b_i);
    
    /**
    * @brief The angular coefficient \f$ I_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1c} \f$
    */
    double  I_1c(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1s} \f$
    */
    double  I_1s(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{2c} \f$
    */
    double  I_2c(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{2s} \f$
    */
    double  I_2s(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_3 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_3 \f$
    */
    double  I_3(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_4 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_4 \f$
    */
    double  I_4(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_5 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_5 \f$
    */
    double  I_5(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{6c} \f$
    */
    double  I_6c(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{6s} \f$
    */
    double  I_6s(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_7 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_7 \f$
    */
    double  I_7(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_8 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_8 \f$
    */
    double  I_8(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ I_9 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_9 \f$
    */
    double  I_9(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_1s \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_1s \f$
    */
    double  h_1s(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_1c \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_1c \f$
    */
    double  h_1c(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_2s \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_2s \f$
    */
    double  h_2s(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_2c \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_2c \f$
    */
    double  h_2c(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_3 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_3 \f$
    */
    double  h_3(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_4 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_4 \f$
    */
    double  h_4(double q2, bool bar);
    
    /**
    * @brief The angular coefficient \f$ h_7 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ h_7 \f$
    */
    double  h_7(double q2, bool bar);
    
    /**
    * @brief The CP average \f$ \Sigma_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{1s} \f$
    */
    double getSigma1c(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_1c(q2, 0) + I_1c(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_1c(q2, 0) + I_1c(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_1c(q2, 0) - ys * h_1c(q2, 0) + I_1c(q2, 1) - ys * h_1c(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma1c : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{1c} \f$
    */
    double getSigma1s(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_1s(q2, 0) + I_1s(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_1s(q2, 0) + I_1s(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_1s(q2, 0) - ys * h_1s(q2, 0) + I_1s(q2, 1) - ys * h_1s(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma1s : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2s} \f$
    */
    double getSigma2c(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_2c(q2, 0) + I_2c(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_2c(q2, 0) + I_2c(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_2c(q2, 0) - ys * h_2c(q2, 0) + I_2c(q2, 1) - ys * h_2c(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma2c : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2c} \f$
    */
    double getSigma2s(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_2s(q2, 0) + I_2s(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_2s(q2, 0) + I_2s(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_2s(q2, 0) - ys * h_2s(q2, 0) + I_2s(q2, 1) - ys * h_2s(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma2s : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{3} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{3} \f$
    */ 
    double getSigma3(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_3(q2, 0) + I_3(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_3(q2, 0) + I_3(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_3(q2, 0) - ys * h_3(q2, 0) + I_3(q2, 1) - ys * h_3(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma3 : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{4} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{4} \f$
    */ 
    double getSigma4(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_4(q2, 0) + I_4(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_4(q2, 0) + I_4(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_4(q2, 0) - ys * h_4(q2, 0) + I_4(q2, 1) - ys * h_4(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma4 : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{5} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{5} \f$
    */ 
    double getSigma5(double q2)
    {
        return (I_5(q2, 0) + I_5(q2, 1))/2.;
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{6s} \f$
    */ 
    double getSigma6s(double q2)
    {
        return (I_6s(q2, 0) + I_6s(q2, 1))/2.;
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{6c} \f$
    */ 
    double getSigma6c(double q2)
    {
        return (I_6c(q2, 0) + I_6c(q2, 1))/2.;
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{7} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{7} \f$
    */ 
    double getSigma7(double q2)
    {
        switch(vectorM){
            case StandardModel::K_star:
                return (I_7(q2, 0) + I_7(q2, 1))/2.;
                break;
            case StandardModel::K_star_P:
                return (I_7(q2, 0) + I_7(q2, 1))/2.;
                break;
            case StandardModel::PHI:
                return (I_7(q2, 0) - ys * h_7(q2, 0) + I_7(q2, 1) - ys * h_7(q2, 1))/2.;
                break;
            default:
                std::stringstream out;
                out << vectorM;
                throw std::runtime_error("MVll::getSigma7 : vector " + out.str() + " not implemented");
        }
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{8} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{8} \f$
    */ 
    double getSigma8(double q2)
    {
        return (I_8(q2, 0) + I_8(q2, 1))/2.;
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{9} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{9} \f$
    */ 
    double getSigma9(double q2)
    {
        return (I_9(q2, 0) + I_9(q2, 1))/2.;
    };
    
        /**
    * @brief The CP difference \f$ \Delta_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1s} \f$
    */
    double getDelta1c(double q2)
    {
        return (I_1c(q2, 0) - I_1c(q2, 1)) / 2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1c} \f$
    */
    double getDelta1s(double q2)
    {
        return (I_1s(q2, 0) - I_1s(q2, 1)) / 2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2s} \f$
    */
    double getDelta2c(double q2)
    {
        return (I_2c(q2, 0) - I_2c(q2, 1)) / 2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2c} \f$
    */
    double getDelta2s(double q2)
    {
        return (I_2s(q2, 0) - I_2s(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{3} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{3} \f$
    */ 
    double getDelta3(double q2)
    {
        return (I_3(q2, 0) - I_3(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{4} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{4} \f$
    */ 
    double getDelta4(double q2)
    {
        return (I_4(q2, 0) - I_4(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{5} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{5} \f$
    */ 
    double getDelta5(double q2)
    {
        return (I_5(q2, 0) - I_5(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{6s} \f$
    */ 
    double getDelta6s(double q2)
    {
        return (I_6s(q2, 0) - I_6s(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{6c} \f$
    */ 
    double getDelta6c(double q2)
    {
        return (I_6c(q2, 0) - I_6c(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{7} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{7} \f$
    */ 
    double getDelta7(double q2)
    {
        return (I_7(q2, 0) - I_7(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{8} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{8} \f$
    */ 
    double getDelta8(double q2)
    {
        return (I_8(q2, 0) - I_8(q2, 1))/2.;
    };
    
    /**
    * @brief The CP difference \f$ \Delta_{9} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{9} \f$
    */ 
    double getDelta9(double q2)
    {
        return (I_9(q2, 0) - I_9(q2, 1))/2.;
    };
    
    /**
    * @brief The fit function from @cite Straub:2015ica, \f$ FF^{\rm fit} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] a_2 fit parameter
    * @param[in] MR2 square of the nearest resonance mass
    * @return \f$ FF^{\rm fit} \f$
    */
    double FF_fit(double q2, double a_0, double a_1, double a_2, double MR2);
    
    /**
    * @brief The \f$ I_1 \f$ function from @cite Beneke:2001at .
    * @param[in] u dummy variable to be integrated out
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_1 \f$
    */
    gslpp::complex I1(double u, double q2);
    
    /**
    * @brief The \f$ T^{\perp}_+ \f$ function from @cite Beneke:2001at .
    * @param[in] u dummy variable to be integrated out
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T^{\perp}_+ \f$
    */
    gslpp::complex Tperpplus(double u, double q2);
    
    /**
    * @brief The \f$ T^{\parallel}_+ \f$ function from @cite Beneke:2001at .
    * @param[in] u dummy variable to be integrated out
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T^{\parallel}_+ \f$
    */
    gslpp::complex Tparplus(double u, double q2);
    
    /**
    * @brief The \f$ T^{\parallel}_- \f$ function from @cite Beneke:2001at .
    * @param[in] u dummy variable to be integrated out
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T^{\parallel}_- \f$
    */
    gslpp::complex Tparminus(double u, double q2);
    
    /**
    * @brief The real part of the integral involving \f$ T^{\perp}_+ \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Re T^{\perp}_+ \Phi^{\perp}\f$
    */
    double Integrand_ReTperpplus(double up);
    
    /**
    * @brief The imaginary part of the integral involving \f$ T^{\perp}_+ \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Im T^{\perp}_+ \Phi^{\perp}\f$
    */
    double Integrand_ImTperpplus(double up);
    
    /**
    * @brief The real part of the integral involving \f$ T^{\parallel}_+ \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Re T^{\parallel}_+ \Phi^{\parallel}\f$
    */
    double Integrand_ReTparplus(double up);
    
    /**
    * @brief The imaginary part of the integral involving \f$ T^{\parallel}_+ \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Im T^{\parallel}_+ \Phi^{\parallel}\f$
    */
    double Integrand_ImTparplus(double up);
    
    /**
    * @brief The real part of the integral involving \f$ T^{\parallel}_- \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Re T^{\parallel}_- \Phi^{\parallel}\f$
    */
    double Integrand_ReTparminus(double up);
    
    /**
    * @brief The imaginary part of the integral involving \f$ T^{\parallel}_- \f$ at fixed \f$ q^2 \f$, according to @cite Beneke:2001at .
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Im T^{\parallel}_- \Phi^{\parallel}\f$
    */
    double Integrand_ImTparminus(double up);
    
    /**
    * @brief The sum of Integrand_ReTparplus() and Integrand_ReTparminus().
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Re T^{\parallel}_+ \Phi^{\parallel} + Re T^{\parallel}_- \Phi^{\parallel}\f$
    */
    double Integrand_ReTpar_pm(double up);
    
    /**
    * @brief The sum of Integrand_ImTparplus() and Integrand_ImTparminus().
    * @param[in] up dummy variable to be integrated out
    * @return \f$ Im T^{\parallel}_+ \Phi^{\parallel} + Im T^{\parallel}_- \Phi^{\parallel}\f$
    */
    double Integrand_ImTpar_pm(double up);

    /**
    * @brief The correction \f$ F_{19} \f$ from @cite Asatrian:2001de.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F_{19} \f$
    */
    gslpp::complex F19(double q2);

    /**
    * @brief The correction \f$ F_{27} \f$ from @cite Asatrian:2001de.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F_{27} \f$
    */
    gslpp::complex F27(double q2);

    /**
    * @brief The correction \f$ F_{29} \f$ from @cite Asatrian:2001de.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F_{29} \f$
    */
    gslpp::complex F29(double q2);

    /**
    * @brief The correction \f$ F_{87} \f$ from @cite Asatrian:2001de.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F_{87} \f$
    */
    gslpp::complex F87(double q2);

    /**
    * @brief The correction \f$ F_{89} \f$ from @cite Asatrian:2001de.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F_{89} \f$
    */
    double F89(double q2);
    
    /**
    * @brief The correction \f$ C_{\perp} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ C_{\perp} \f$
    */
    gslpp::complex Cperp(double q2);
    
    /**
    * @brief The correction \f$ C_{\parallel} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ C_{\parallel} \f$
    */
    gslpp::complex Cpar(double q2);
    
    /**
    * @brief The total correction \f$ \Delta \mathcal{T}^{\perp} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta \mathcal{T}^{\perp} \f$
    */
    gslpp::complex deltaTperp(double q2);
    
    /**
    * @brief The total correction \f$ \Delta \mathcal{T}^{\parallel} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta \mathcal{T}^{\parallel} \f$
    */
    gslpp::complex deltaTpar(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^0 \f$ computed integrating over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^0 \f$
    */
    gslpp::complex DeltaC9_0(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^+ \f$ computed integrating over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^+ \f$
    */
    gslpp::complex DeltaC9_p(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^- \f$ computed integrating over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^- \f$
    */
    gslpp::complex DeltaC9_m(double q2);
    
    /**
    * @brief The fit function for the real part of the QCDF correction \f$ \Delta C_9^{\lambda} \f$.
    * @param[in] x fit variable
    * @param[in] p fit parameters vector
    * @return \f$ f_{Re \Delta C_9^{\lambda}} \f$
    */
    double reDC9fit(double* x, double* p);
    
    /**
    * @brief The fit function for the imaginary part of the QCDF correction \f$ \Delta C_9^{\lambda} \f$.
    * @param[in] x fit variable
    * @param[in] p fit parameters vector
    * @return \f$ f_{Im \Delta C_9^{\lambda}} \f$
    */
    double imDC9fit(double* x, double* p);
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^+ \f$ in the muon channel.
    */
    void fit_DeltaC9_p_mumu();
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^- \f$ in the muon channel.
    */
    void fit_DeltaC9_m_mumu();
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^0 \f$ in the muon channel.
    */
    void fit_DeltaC9_0_mumu();
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^+ \f$ in the electron channel.
    */
    void fit_DeltaC9_p_ee();
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^- \f$ in the electron channel.
    */
    void fit_DeltaC9_m_ee();
    
    /**
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^0 \f$ in the electron channel.
    */
    void fit_DeltaC9_0_ee();
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^+ \f$ computed fitting over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^+ \f$
    */
    gslpp::complex fDeltaC9_p(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^- \f$ computed fitting over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^- \f$
    */
    gslpp::complex fDeltaC9_m(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9^0 \f$ computed fitting over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9^0 \f$
    */
    gslpp::complex fDeltaC9_0(double q2);
    
};

#endif	/* MVLL_H */

