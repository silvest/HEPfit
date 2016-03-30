/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLL_H
#define	MPLL_H

#include <StandardModel.h>
#include <ThObservable.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>
#include <gsl/gsl_monte_plain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFitResultPtr.h>

#define SWITCH 8.2

/**
 * @class MPll
 * @ingroup Flavour
 * @brief A class for the @f$M \to P l^+ l^-@f$ decay.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to
 * build the observables relative to the @f$M \to P l^+ l^-@f$ decays, where
 * @f$M@f$ is a generic meson and @f$P@f$ is a pseudoscalar meson. This kind of decays can be described
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
 * between the initial state @f$M@f$ and the final state @f$P l^+ l^-@f$, only the contribution of 
 * @f$\mathcal{H}_\mathrm{eff}^\mathrm{sl+\gamma}@f$ clearly factorizes into the 
 * product of hadronic form factors and leptonic tensors at all orders in strong interactions. 
 * Following @cite Jager:2012uw, we implemented the amplitude in the helicity basis; 
 * hence we made use of the helicity form factors @f$ \tilde{V}_0(q^2), 
 * \tilde{T}_0(q^2)@f$ and @f$\tilde{S}(q^2) @f$, which are related to the
 * ones in the transverse basis through the following relations :
 * @f[
 * \tilde{V}_0(q^2) = i \frac{\sqrt{\lambda(q^2)}}{2m_M\sqrt{q^2}}f_+(q^2)\,,\\
 * \tilde{T}_0(q^2) = i \frac{\sqrt{\lambda(q^2)q^2}}{2m_M^2(m_M+m_P)}f_T(q^2)\,,\\
 * \tilde{S}(q^2) = -\frac{m_M^2-m_P^2}{2m_M(m_b+m_s)}\frac{1+m_s/m_b}{1-m_s/m_b}f_0(q^2)\,,
 * @f]
 * where @f$\lambda(q^2) = 4m_M^2|\vec{k}|^2@f$, with @f$\vec{k}@f$ as the 3-momentum
 * of the meson @f$P@f$ in the @f$M@f$ rest frame.
 * 
 * The effect of the operators of @f$\mathcal{H}_\mathrm{eff}^\mathrm{had}@f$ due to
 * exchange of soft gluon can be reabsorbed in the following parameterization,
 * @f[
 * h_0(q^2) = \frac{\epsilon^*_\mu(\lambda)}{m_M^2} 
 * \int d^4x e^{iqx} \langle \bar P \vert T\{j^{\mu}_\mathrm{em} (x) 
 * \mathcal{H}_\mathrm{eff}^\mathrm{had} (0)\} \vert \bar M \rangle = 
 * h_0^{(0)} + \frac{q^2}{1\,\mathrm{GeV}^2} h_0^{(1)}\,.
 * @f]
 * 
 * The amplitude can be therefore parametrized in terms of the following helicity amplitudes:
 * @f[
 * H_V = -i\, N \Big\{C_{9} \tilde{V}_{L,0} +C_{9}'  \tilde{V}_{R,0}
 * + \frac{m_M^2}{q^2} \Big[\frac{2\, m_b}{m_M} (C_{7} \tilde{T}_{L,0} +  C_{7}'  \tilde{T}_{R,0})
            - 16 \pi^2 h_0 \Big] \Big\} \,,  \\
 * H_A = -i\, N (C_{10}  \tilde{V}_{L,0} + C_{10}'\tilde{V}_{R,0}) \,, \\
 * H_S = i\, N \frac{ m_b}{m_W} (C_S \tilde{S}_L + C_S' \tilde{S}_R)\,, \\
 * H_P = i\, N \Big\{ \frac{ m_b}{m_W} (C_P \tilde{S}_L + C_P' \tilde{S}_R)
 * + \frac{2\,m_\ell m_b}{q^2} \left[C_{10} \Big(\tilde{S}_L - \frac{m_s}{m_b} \tilde{S}_R \Big) 
 * + C_{10}' \Big(\tilde{S}_R - \frac{m_s}{m_b} \tilde{S}_L\Big) \right] \Big\} \,,
 * @f]
 * where @f$ N = - \frac{4 G_F m_M}{\sqrt{2}}\frac{e^2}{16\pi^2}\lambda_t@f$ and we have defined
 * @f[
 * \tilde{V}_{L,0}(q^2) = -\tilde{V}_{R,0}(q^2)=\tilde{V}_0(q^2)\,,\\
 * \tilde{T}_{L,0}(q^2) = -\tilde{T}_{R,0}(q^2)=\tilde{V}_0(q^2)\,,\\
 * \tilde{S}_L(q^2) = -\tilde{S}_R(q^2)=\tilde{S}(q^2)\,.
 * @f]
 * Squaring the amplitude and summing over the spins it is possible to obtain 
 * the fully differential decay rate, which is
 * @f[
 * \frac{d^{(4)} \Gamma}{dq^2\,d(\cos\theta_l)} = \frac{9}{32\,\pi} 
 * \Big( I^c_1 +I^c_2\cos2\theta_l + I_6^c \cos\theta_l \Big) 
 * @f]
 * The angular coefficients involved in the differential decay rate are related to the
 * helicity amplitudes according to the following relations:
 * @f[
 * I_1^c = F \left\{ \frac{1}{2}\left(|H_V^0|^2+|H_A^0|^2\right)+ 
 * |H_P|^2+\frac{2m_\ell^2}{q^2}\left(|H_V^0|^2-|H_A^0|^2\right) + \beta^2 |H_S|^2 \right\}\,,\\
 * I_2^c = -F\, \frac{\beta^2}{2}\left(|H_V^0|^2+|H_A^0|^2\right)\,,\\
 * I_6^c = 2 F \frac{\beta\, m_\ell}{\sqrt{q^2}} {\rm Re} \left[ H_S^* H_V^0 \right]\,,\\
 * @f]
 * where
 * @f[
 * F=\frac{ \lambda^{1/2}\beta\, q^2}{3 \times 2^{5} \,\pi^3\, m_M^3}\,,
 * \qquad \beta = \sqrt{1 - \frac{4 m_\ell^2}{q^2} }\,.
 * @f]
 * The final observables are hence build employing CP-averages @f$\Sigma_i@f$ or CP-asymmetries @f$\Delta_i@f$ of
 * such angular coefficients; however, since on the experimental side the observables
 * are averaged over @f$ q^2 @f$ bins, an integration of the coeffiecients over such
 * bins has to be performed before they are combined in order to build the observables.
 * 
 * The class is organized as follows: after the parameters are updated in 
 * updateParameters() and the cache is checked in
 * checkCache(), the form factor are build in the transverse basis in the functions
 * f_plus(), f_0() and f_T() @cite Ball:2004ye. They are consequentely translated in the helicity basis through the
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
    * @brief The helicity amplitude \f$ H_V^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_V^{\lambda} \f$
    */
    gslpp::complex H_V(double q2);


    /**
    * @brief The helicity amplitude \f$ H_A^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_A^{\lambda} \f$
    */
    gslpp::complex H_A(double q2);


    /**
    * @brief The helicity amplitude \f$ H_S^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_S^{\lambda} \f$
    */
    gslpp::complex H_S(double q2);


    /**
    * @brief The helicity amplitude \f$ H_P^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ H_P^{\lambda} \f$
    */
    gslpp::complex H_P(double q2);
    
    
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
    

private:
    const StandardModel& mySM;/**< Model type */
    StandardModel::lepton lep;/**< Final leptons type */
    StandardModel::meson meson;/**< Initial meson type */
    StandardModel::meson pseudoscalar;/**< Final pseudoscalar meson type */
    
    double GF;            /**<Fermi constant */
    double ale;           /**<alpha electromagnetic */
    double Mlep;          /**<muon mass */
    double MM;            /**<initial meson mass */
    double MP;            /**<final pseudoscalar meson mass */
    double Mb;            /**<b quark mass */
    double mu_b;          /**<b mass scale */
    double mu_h;          /**<\f$\sqrt{\mu_b*\lambda_{QCD}}\f$ */
    double Mc;            /**<c quark mass */
    double Ms;            /**<s quark mass */
    double width;         /**<initial meson width */
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor */
    gslpp::complex h_0;          /**<parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_0_1;        /**<parameter that contains the contribution from the hadronic hamiltonian */
//    double q2;            /**<\f$q^2\f$ of the decay */
    
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
    gslpp::vector<gslpp::complex> ** allcoeffh;/**<Vector that contains the Wilson coeffients at scale @f$\mu_h@f$ */
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients */
    
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
    
    
    std::vector<double> ReDeltaC9;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> ImDeltaC9;/**<Vector that samples the QCDF @f$\Delta C_9@f$ */
    std::vector<double> myq2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    TFitResultPtr refres;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
    TFitResultPtr imfres;/**<Vector that contains the fitted QCDF @f$\Delta C_9@f$ */
        
    TGraph gr1;/**<Tgraph to be used for fitting the QCDF @f$\Delta C_9@f$ */
    TGraph gr2;/**<Tgraph to be used for fitting the QCDF @f$\Delta C_9@f$ */
    
    TF1 reffit;/**<TF1 to be used for fitting the QCDF @f$\Delta C_9@f$ */
    TF1 imffit;/**<TF1 to be used for fitting the QCDF @f$\Delta C_9@f$ */
    
    double tmpq2;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    
    gslpp::complex H_0_pre;/**< Cache variable */
    gslpp::complex H_0_WC;/**< Cache variable */
    gslpp::complex H_c_WC;/**< Cache variable */
    gslpp::complex H_b_WC;/**< Cache variable */
    
    gslpp::complex ihalfMPI;/**< Cache variable */
    double fournineth;/**< Cache variable */
    double half;/**< Cache variable */
    double twothird;/**< Cache variable */
    double Mc2;/**< Cache variable */
    double Mb2;/**< Cache variable */
    double logMc;/**< Cache variable */
    double logMb;/**< Cache variable */
    double mu_b2;/**< Cache variable */
    double fourMc2;/**< Cache variable */
    double fourMb2;/**< Cache variable */
    double Mlep2;/**< Cache variable */
    double NN;/**< Cache variable */
    double MM2;/**< Cache variable */
    double MM4;/**< Cache variable */
    double MP2;/**< Cache variable */
    double MP4;/**< Cache variable */
    double MM2mMP2;/**< Cache variable */
    double twoMP2;/**< Cache variable */
    double twoMM;/**< Cache variable */
    double twoMM2;/**< Cache variable */
    double twoMM2_MMpMP;/**< Cache variable */
    double twoMM_MbpMs;/**< Cache variable */
    double S_L_pre;/**< Cache variable */
    double fourMM2;/**< Cache variable */
    double twoMboMM;/**< Cache variable */
    double sixteenM_PI2;/**< Cache variable */
    double ninetysixM_PI3MM3;/**< Cache variable */
    double MboMW;/**< Cache variable */
    double MboMM;/**< Cache variable */
    double MsoMb;/**< Cache variable */
    double twoMlepMb;/**< Cache variable */
    double DC9pre;/**< Cache variable */
    double threeGegen0;/**< Cache variable */
    double threeGegen1otwo;/**< Cache variable */
    double M_PI2osix;/**< Cache variable */
    double twoMc2;/**< Cache variable */
    double sixMMoMb;/**< Cache variable */
    double CF;/**< Cache variable */
    double deltaT_0;/**< Cache variable */
    double deltaT_1par;/**< Cache variable */
    
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
    double Ee;/**<Variable used to compute the QCDF @f$\Delta C_9@f$ */
    
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
    
    unsigned int F_updated;/**< Cache variable */
    
    unsigned int VL_updated;/**< Cache variable */
    
    unsigned int TL_updated;/**< Cache variable */
    
    unsigned int SL_updated;/**< Cache variable */
    gslpp::vector<double> SL_cache;/**< Cache variable */
    
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
    
    unsigned int C_2Lh_updated;/**< Cache variable */
    gslpp::complex C_2Lh_cache;/**< Cache variable */
    
    unsigned int C_8Lh_updated;/**< Cache variable */
    gslpp::complex C_8Lh_cache;/**< Cache variable */
    
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
    
    std::map<std::pair<double, double>, unsigned int > I1Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta2Cached;/**< Cache variable */
    
    double avaSigma;/**< Gsl integral variable */
    double avaDelta;/**< Gsl integral variable */
    double avaDTPPR;/**< Gsl integral variable */ 
    
    double errSigma;/**< Gsl integral variable */
    double errDelta;/**< Gsl integral variable */
    double errDTPPR;/**< Gsl integral variable */
    
    gsl_function FS;/**< Gsl integral variable */
    gsl_function FD;/**< Gsl integral variable */
    gsl_function DTPPR;/**< Gsl integral variable */
    
    gsl_integration_cquad_workspace * w_sigma;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_delta;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_DTPPR;/**< Gsl integral variable */
    
    gsl_error_handler_t * old_handler; /**< GSL error handler store */
    
    std::map<std::pair<double, double>, gslpp::complex > cacheI1;/**< Cache variable */
    
    std::map<std::pair<double, double>, double > cacheSigma0;/**< Gsl integral variable */
    std::map<std::pair<double, double>, double > cacheSigma2;/**< Gsl integral variable */
    
    std::map<std::pair<double, double>, double > cacheDelta0;/**< Gsl integral variable */
    std::map<std::pair<double, double>, double > cacheDelta2;/**< Gsl integral variable */
    
    std::map<double, unsigned int> deltaTparpCached;/**< Cache variable */
    std::map<double, unsigned int> deltaTparmCached;/**< Cache variable */
    
    std::map<double, gslpp::complex> cacheDeltaTparp;/**< Cache variable */
    std::map<double, gslpp::complex> cacheDeltaTparm;/**< Cache variable */
    
    unsigned int deltaTparpupdated;/**< Cache variable */
    unsigned int deltaTparmupdated;/**< Cache variable */
    
    unsigned int T_updated;/**< Cache variable */
    gslpp::vector<double> T_cache;/**< Cache variable */
    
    
    
    /**
     * @brief The update parameter method for MPll.
     */
    void updateParameters();
    
    /**
     * @brief The caching method for MPll.
     */
    void checkCache();
    
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
    * @brief The helicity form factor \f$ T_L^0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_L^{\lambda} \f$
    */
    gslpp::complex T_L(double q2);


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
    * @param[in] mu2 mass scale
    * @return \f$ h(q^2,m_c) \f$
    */
    gslpp::complex H_c(double q2, double mu2);
    
    /**
    * @brief The \f$ h(q^2,m_b) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] mu2 mass scale
    * @return \f$ h(q^2,m_b) \f$
    */
    gslpp::complex H_b(double q2, double mu2);
    
    
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
    * @brief The factor \f$ \beta \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta \f$
    */
    double beta (double q2);
        
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
    * @brief The factor \f$ F \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ F \f$
    */
    double F(double q2);
    
    /**
    * @brief The angular coefficient \f$ I_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{1c} \f$
    */
    double  I_1c(double q2);
    
    /**
    * @brief The angular coefficient \f$ I_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{2c} \f$
    */
    double  I_2c(double q2);
    
    /**
    * @brief The angular coefficient \f$ I_{6c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_{6c} \f$
    */
    double  I_6c(double q2);
    
    
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
    * @return \f$ \Sigma_{1c} \f$
    */
    double getSigma1c(double q2)
    {
        return I_1c(q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2c} \f$
    */
    double getSigma2c(double q2)
    {
        return I_2c(q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{6c} \f$
    */
    double getSigma6c(double q2)
    {
        return I_6c(q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1c} \f$
    */
    double getDelta1c(double q2)
    {
        return Delta(0, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2s} \f$
    */
    double getDelta2c(double q2)
    {
        return Delta(2, q2);
    };
    
    /**
    * @brief The \f$ I_1 \f$ function from @cite Beneke:2001at .
    * @param[in] u dummy variable to be integrated out
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ I_1 \f$
    */
    gslpp::complex I1(double u, double q2);
    
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
    * @brief The correction \f$ C_{\parallel} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ C_{\parallel} \f$
    */
    gslpp::complex Cpar(double q2);
    
    /**
    * @brief The total correction \f$ \Delta \mathcal{T}^{\parallel} \f$ from @cite Beneke:2001at .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta \mathcal{T}^{\parallel} \f$
    */
    gslpp::complex deltaTpar(double q2);
    
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
    * @brief The fitting routine for the QCDF correction \f$ \Delta C_9^0 \f$ in the muon channel.
    */
    void fit_DeltaC9_mumu();
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9 \f$ computed fitting over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9 \f$
    */
    gslpp::complex fDeltaC9(double q2);
    
    /**
    * @brief The total QCDF correction \f$ \Delta C_9 \f$ computed integrating over \f$ u \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta C_9 \f$
    */
    gslpp::complex DeltaC9(double q2);
    
};

#endif	/* MPLL_H */

    