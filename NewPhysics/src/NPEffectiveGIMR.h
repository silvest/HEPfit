/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVEGIMR_H
#define	NPEFFECTIVEGIMR_H

#include "gslpp.h"
#include "NPbase.h"
#include <string.h>
#include <stdexcept>


/**
 * @class NPEffectiveGIMR
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian, with lepton-flavour universality and
 * without quark-flavour universality.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the dimension-six effective Lagrangian introduced in NPEffective
 * class, where lepton-flavour universality is assumed, but quark-flavour universality
 * is not.
 *
*
 * @anchor NPEffectiveGIMRInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPEffectiveGIMR(), it is required to call the initialization method
 * InitializeModel(), inherited from NPEffective class.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEffectiveGIMRParameters
 * <h3>%Model parameters</h3>
 *
 * *** This table has to be modified!! ***
 *
 *
 * The model parameters of %NPEffectiveGIMR are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cWB_NP </td>
 *   <td class="mod_symb">\f$C_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{WB}=\big(H^\dagger\tau^a H\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cH_NP </td>
 *   <td class="mod_symb">\f$C_H\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{H}=\big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cLL_NP </td>
 *   <td class="mod_symb">\f$C_{LL}  \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{LL}=\frac 12\big(\overline{L}\,\gamma^\mu\tau^a L\big)
 *   \big(\overline{L}\,\gamma_\mu\tau^a L\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHLp_NP </td>
 *   <td class="mod_symb">\f$ C_{HL}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HL}^\prime =i\big(H^\dagger D_\mu\tau^a H\big)
 *   \big(\overline{L}\,\gamma^\mu\tau^a L\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1p_NP, %cHQ2p_NP, %cHQ3p_NP </td>
 *   <td class="mod_symb">\f$ C_{HQ_i}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HQ_i}^\prime =i\big(H^\dagger D_\mu\tau^a H\big)
 *   \big(\overline{Q_i}\,\gamma^\mu\tau^a Q_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL_NP </td>
 *   <td class="mod_symb">\f$ C_{HL} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HL} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{L}\,\gamma^\mu L\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1_NP, %cHQ2_NP, %cHQ3_NP </td>
 *   <td class="mod_symb">\f$ C_{HQ_i} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HQ_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{Q_i}\,\gamma^\mu Q_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE_NP </td>
 *   <td class="mod_symb">\f$ C_{HE} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HE} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{E}\,\gamma^\mu E\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU1_NP, %cHU2_NP, %cHU3_NP </td>
 *   <td class="mod_symb">\f$ C_{HU_i}  \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HU_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{U_i}\,\gamma^\mu U_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD1_NP, %cHD2_NP, %cHD3_NP </td>
 *   <td class="mod_symb">\f$ C_{HD_i} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HD_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{D_i}\,\gamma^\mu D_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPEffectiveGIMRFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPEffectiveGIMRFunctions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 *
 */
class NPEffectiveGIMR : public NPbase {
public:

    /**
     *　@brief The number of the model parameters in %NPEffectiveGIMR. 
     */
    static const int NNPEffectiveGIMRVars = 210;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveGIMR.
     */
    static const std::string NPEffectiveGIMRVars[NNPEffectiveGIMRVars];

    /**
     *　@brief The number of the model parameters in %NPEffectiveGIMR
     * with lepton and quark flavour universalities.
     */
    static const int NNPEffectiveGIMRVars_LFU_QFU = 114;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveGIMR with lepton and quark flavour universalities.
     */
    static const std::string NPEffectiveGIMRVars_LFU_QFU[NNPEffectiveGIMRVars_LFU_QFU];

    /**
     * @brief Constructor.
     * @param[in] FlagLeptonUniversal_in
     * @param[in] FlagQuarkUniversal_in
     */
    NPEffectiveGIMR(const bool FlagLeptonUniversal_in = false, const bool FlagQuarkUniversal_in = false);

    /**
     * @brief The post-update method for %NPEffectiveGIMR.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief A method to check if all the mandatory parameters for %NPEffectiveGIMR
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPEffectiveGIMR.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the Fermi constant.
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimental value measured through muon decays, 
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const;

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the value of @f$U@f$
     */
    virtual double obliqueU() const;

    /**
     * @brief @copybrief StandardModel::Mw()
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /**
     * @brief @copybrief StandardModel::GammaW()
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;

    /**
     * @brief New physics contribution to @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual double deltaGV_f(const Particle p) const;

    /**
     * @brief New physics contribution to @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual double deltaGA_f(const Particle p) const;

    /**
     * @brief New physics contribution to @f$g_L^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_L^f@f$
     */
        // no generation mixing
    double deltaGL_f(const Particle p) const;

    /**
     * @brief New physics contribution to @f$g_R^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_R^f@f$
     */
        // no generation mixing
    double deltaGR_f(const Particle p) const;


        ////////////////////////////////////////////////////////////////////////
    
    virtual double sigma_eeTOffbar(const Particle p, const double sqrt_s) const;
    
    virtual gslpp::complex sigma_eeTOffbarF(const Particle p, const double sqrt_s) const;
    virtual gslpp::complex sigma_eeTOffbarB(const Particle p, const double sqrt_s) const;
    
    double sigma_eeTOmumu(const double sqrt_s) const;
    double sigma_eeTOqq(const double sqrt_s) const;
    
    double AFB_mu(const double sqrt_s) const;

    ////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{L}@f$
     */
    // no generation mixing
    gslpp::complex deltaGL_Wff(const Particle pbar, const Particle p) const;
    /**
     * @brief New physics contribution to the charged current coupling @f$W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Wff}^{R}@f$
     */
    // no generation mixing
    gslpp::complex deltaGR_Wff(const Particle pbar, const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H G_{\mu\nu}^AG^{A \mu\nu}@f$.
     * @return @f$\delta g_{HGG}@f$
     */
    double deltaG_hgg() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu\nu}^\dagger W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(1)}@f$
     */
    double deltaG1_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\nu}^\dagger \partial^\mu W^{\mu\nu}@f$.
     * @return @f$\delta g_{HWW}^{(2)}@f$
     */
    double deltaG2_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_{\mu}^\dagger W^{\mu}@f$.
     * @return @f$\delta g_{HWW}^{(3)}@f$
     */
    double deltaG3_hWW() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(1)}@f$
     */
    double deltaG1_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu Z^{\mu\nu}@f$.
     * @return @f$\delta g_{HZZ}^{(2)}@f$
     */
    double deltaG2_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu} Z^{\mu}@f$.
     * @return @f$\delta g_{HZZ}^{(3)}@f$
     */
    double deltaG3_hZZ() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(1)}@f$
     */
    double deltaG1_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_{\nu} \partial^\mu F^{\mu\nu}@f$.
     * @return @f$\delta g_{HZA}^{(2)}@f$
     */
    double deltaG2_hZA() const;
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H F_{\mu\nu} F^{\mu\nu}@f$.
     * @return @f$\delta g_{HAA}@f$
     */
    double deltaG_hAA() const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H f\bar{f}@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{Hff}@f$
     */
    // no generation mixing
    gslpp::complex deltaG_hff(const Particle p) const;

    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{WffH}^{L}@f$
     */
    // no generation mixing
    gslpp::complex deltaGL_Wffh(const Particle pbar, const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H W_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] pbar a lepton or quark
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{WffH}^{R}@f$
     */
    // no generation mixing
    gslpp::complex deltaGR_Wffh(const Particle pbar, const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_\mu \bar{f_L}\gamma^mu f_L@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{ZffH}^{L}@f$
     */
    // no generation mixing
    double deltaGL_Zffh(const Particle p) const;
    
    /**
     * @brief The new physics contribution to the coupling of the effective interaction @f$H Z_\mu \bar{f_R}\gamma^mu f_R@f$.
     * @param[in] p a lepton or quark
     * @return @f$\delta g_{ZffH}^{R}@f$
     */
    // no generation mixing 
    double deltaGR_Zffh(const Particle p) const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Loop function entering in the calculation of the effecive @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the particle in the loop.
     * @return @f$f(\tau)@f$
     */
    gslpp::complex f_triangle(const double tau) const;
    /**
     * @brief Fermionic loop function entering in the calculation of the effecive @f$Hgg@f$ and @f$H\gamma\gamma@f$ couplings.
     * @details
     * @f$A^H_f(\tau)=2\tau [1+(1-\tau)f(\tau)]@f$
     * @param[in] @f$\tau=4 M^2/m_h^2@f$, with @f$M@f$ the mass of the fermion in the loop.
     * @return @f$A^H_f(\tau)@f$
     */
    gslpp::complex AH_f(const double tau) const;

    /**
     * @brief The ratio @f$\mu_{ggH}@f$ between the gluon-gluon fusion Higgs
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH}@f$
     */
    virtual double muggH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF}@f$ between the vector-boson fusion Higgs
     * production cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF}@f$
     */
    virtual double muVBF(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{WH}@f$ between the W-Higgs associated production
     * cross-section in the current model and in the Standard Model. 
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{WH}@f$
     */
    virtual double muWH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ZH}@f$ between the Z-Higgs associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ZH}@f$
     */
    virtual double muZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{eeZH}@f$ between the 
     * @f$e^{+}e^{-}\to ZH}@f$ associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{eeZH}@f$
     */
    virtual double mueeZH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VH}@f$ between the WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VH}@f$
     */
    virtual double muVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{VBF+VH}@f$ between the sum of VBF and WH+ZH associated production
     * cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{VBF+VH}@f$
     */
    virtual double muVBFpVH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ttH}@f$ between the t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ttH}@f$
     */
    virtual double muttH(const double sqrt_s) const;
    /**
     * @brief The ratio @f$\mu_{ggH+ttH}@f$ between the sum of gluon-gluon fusion
     * and t-tbar-Higgs associated 
     * production cross-section in the current model and in the Standard Model.
     * @param[in] sqrt_s the center-of-mass energy in TeV
     * @return @f$\mu_{ggH+ttH}@f$
     */
    virtual double muggHpttH(const double sqrt_s) const;
    /**
     * @brief The ratio of the Br@f$(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to gg)@f$/Br@f$(H\to gg)_{\mathrm{SM}}@f$
     */
    virtual double BrHggRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to WW)@f$/Br@f$(H\to WW)_{\mathrm{SM}}@f$
     */
    virtual double BrHWWRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to ZZ)@f$/Br@f$(H\to ZZ)_{\mathrm{SM}}@f$
     */
    virtual double BrHZZRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to Z\gamma)@f$/Br@f$(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHZgaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \gamma\gamma)@f$/Br@f$(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    virtual double BrHgagaRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to \tau^+\tau^-)@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to \tau^+\tau^-)@f$/Br@f$(H\to \tau^+\tau^-)_{\mathrm{SM}}@f$
     */
    virtual double BrHtautauRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to c\bar{c})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to c\bar{c})@f$/Br@f$(H\to c\bar{c})_{\mathrm{SM}}@f$
     */
    virtual double BrHccRatio() const;
    /**
     * @brief The ratio of the Br@f$(H\to b\bar{b})@f$ in the current model
     * and in the Standard Model.
     * @return Br@f$(H\to b\bar{b})@f$/Br@f$(H\to b\bar{b})_{\mathrm{SM}}@f$
     */
    virtual double BrHbbRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double computeGammaTotalRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double deltaGammaTotalRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H)@f$ 
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H)@f$/@f$\Gamma(H)_{\mathrm{SM}}@f$
     */
    virtual double deltaGammaTotalRatio2() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to gg)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    double GammaHggRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to WW)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double GammaHWWRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to ZZ)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double GammaHZZRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to Z\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHZgaRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to \gamma\gamma)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double GammaHgagaRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to \tau\tau)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double GammaHtautauRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to cc)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double GammaHccRatio() const;
    /**
     * @brief The ratio of the @f$\Gamma(H\to bb)@f$ in the current model
     * and in the Standard Model.
     * @return @f$\Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double GammaHbbRatio() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */
    double deltaGammaHggRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. Only terms that are
     * linear in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio1() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to gg)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to gg)@f$/@f$\Gamma(H\to gg)_{\mathrm{SM}}@f$
     */   
    double deltaGammaHggRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to WW)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to WW)@f$/@f$\Gamma(H\to WW)_{\mathrm{SM}}@f$
     */
    double deltaGammaHWWRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to ZZ)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to ZZ)@f$/@f$\Gamma(H\to ZZ)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZZRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to Z\gamma)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to Z\gamma)@f$/@f$\Gamma(H\to Z\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHZgaRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \gamma\gamma)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to \gamma\gamma)@f$/@f$\Gamma(H\to \gamma\gamma)_{\mathrm{SM}}@f$
     */
    double deltaGammaHgagaRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to \tau\tau)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to \tau\tau)@f$/@f$\Gamma(H\to \tau\tau)_{\mathrm{SM}}@f$
     */
    double deltaGammaHtautauRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to cc)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to cc)@f$/@f$\Gamma(H\to cc)_{\mathrm{SM}}@f$
     */
    double deltaGammaHccRatio2() const;
    /**
     * @brief The new physics contribution to the ratio of the @f$\Gamma(H\to bb)@f$
     * in the current model and in the Standard Model. Only terms that are
     * quadratic in the effective Lagrangian coefficients.
     * @return @f$\delta \Gamma(H\to bb)@f$/@f$\Gamma(H\to bb)_{\mathrm{SM}}@f$
     */
    double deltaGammaHbbRatio2() const;
    
    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    double CW; ///< The dimension-6 operator coefficient \f$C_{W}\f$.
    double CHG;
    double CHW;
    double CHB;
    double CHWB;
    double CHD;
    double CHbox;
    double CH;
    double CHL1_11, CHL1_12r, CHL1_13r, CHL1_22, CHL1_23r, CHL1_33;
    double CHL1_12i, CHL1_13i, CHL1_23i;
    double CHL3_11, CHL3_12r, CHL3_13r, CHL3_22, CHL3_23r, CHL3_33;
    double CHL3_12i, CHL3_13i, CHL3_23i;
    double CHe_11, CHe_12r, CHe_13r, CHe_22, CHe_23r, CHe_33;
    double CHe_12i, CHe_13i, CHe_23i;
    double CHQ1_11, CHQ1_12r, CHQ1_13r, CHQ1_22, CHQ1_23r, CHQ1_33;
    double CHQ1_12i, CHQ1_13i, CHQ1_23i;
    double CHQ3_11, CHQ3_12r, CHQ3_13r, CHQ3_22, CHQ3_23r, CHQ3_33;
    double CHQ3_12i, CHQ3_13i, CHQ3_23i;
    double CHu_11, CHu_12r, CHu_13r, CHu_22, CHu_23r, CHu_33;
    double CHu_12i, CHu_13i, CHu_23i;
    double CHd_11, CHd_12r, CHd_13r, CHd_22, CHd_23r, CHd_33;
    double CHd_12i, CHd_13i, CHd_23i;
    double CHud_11r, CHud_12r, CHud_13r, CHud_22r, CHud_23r, CHud_33r;
    double CHud_11i, CHud_12i, CHud_13i, CHud_22i, CHud_23i, CHud_33i;
    double CeH_11r, CeH_12r, CeH_13r, CeH_22r, CeH_23r, CeH_33r;
    double CeH_11i, CeH_12i, CeH_13i, CeH_22i, CeH_23i, CeH_33i;
    double CuH_11r, CuH_12r, CuH_13r, CuH_22r, CuH_23r, CuH_33r;
    double CuH_11i, CuH_12i, CuH_13i, CuH_22i, CuH_23i, CuH_33i;
    double CdH_11r, CdH_12r, CdH_13r, CdH_22r, CdH_23r, CdH_33r;
    double CdH_11i, CdH_12i, CdH_13i, CdH_22i, CdH_23i, CdH_33i;
    double CLL_1221, CLL_2112;
    double CLQ1;
    double CLQ3;
    double Cee;
    double Ceu;
    double Ced;
    double CLe;
    double CLu;
    double CLd;
    double CQe;
    double Lambda_NP;
// The error in the parameters multiplying the dimension-6 operator coefficients in the production cross sections.
    double eVBF2_HZZ1, eVBF2_HZZ2, eVBF2_HZZ3, eVBF2_HZA1, eVBF2_HZA2, eVBF2_HAA;///< VBF (Tevatron 1.96 TeV).
    double eVBF2_HWW1, eVBF2_HWW2, eVBF2_HWW3, eVBF2_Hgg, eVBF2_HZuL, eVBF2_HZuR;
    double eVBF2_HZdL, eVBF2_HZdR, eVBF2_HWud, eVBF2_ZuL, eVBF2_ZuR, eVBF2_ZdL;
    double eVBF2_ZdR, eVBF2_Wud;
    double eVBF78_HZZ1, eVBF78_HZZ2, eVBF78_HZZ3, eVBF78_HZA1, eVBF78_HZA2, eVBF78_HAA;///< VBF (LHC 7 & 8 TeV).
    double eVBF78_HWW1, eVBF78_HWW2, eVBF78_HWW3, eVBF78_Hgg, eVBF78_HZuL, eVBF78_HZuR;
    double eVBF78_HZdL, eVBF78_HZdR, eVBF78_HWud, eVBF78_ZuL, eVBF78_ZuR, eVBF78_ZdL;
    double eVBF78_ZdR, eVBF78_Wud;
    
    double eWH2_HWW1, eWH2_HWW2, eWH2_HWW3, eWH2_HWud, eWH2_Wud;///< WH (Tevatron 1.96 TeV).
    double eWH78_HWW1, eWH78_HWW2, eWH78_HWW3, eWH78_HWud, eWH78_Wud;///< WH (LHC 7 & 8 TeV).
    
    double eZH2_HZZ1, eZH2_HZZ2, eZH2_HZZ3, eZH2_HZA1, eZH2_HZA2, eZH2_HZuL, eZH2_HZuR;///< ZH (Tevatron 1.96 TeV).
    double eZH2_HZdL, eZH2_HZdR, eZH2_ZuL, eZH2_ZuR, eZH2_ZdL, eZH2_ZdR;
    double eZH78_HZZ1, eZH78_HZZ2, eZH78_HZZ3, eZH78_HZA1, eZH78_HZA2, eZH78_HZuL, eZH78_HZuR;///< ZH (LHC 7 & 8 TeV).
    double eZH78_HZdL, eZH78_HZdR, eZH78_ZuL, eZH78_ZuR, eZH78_ZdL, eZH78_ZdR;
    
    double ettH2_Htt, ettH2_Hgg;///< ttH (Tevatron 1.96 TeV).
    double ettH78_Htt, ettH78_Hgg;///< ttH (LHC 7 & 8 TeV).
    
    double MwInput;

    double LambdaNP2;
    double v2_over_LambdaNP2;
    double cW_tree, sW_tree;
    double cW2_tree, sW2_tree;
    double delta_ZZ, delta_AZ, delta_AA;
    double delta_h;

    double CHF1_diag(const Particle F) const;
    double CHF3_diag(const Particle F) const;
    double CHf_diag(const Particle f) const;
    gslpp::complex CHud_diag(const Particle u) const;
    gslpp::complex CfH_diag(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
private:

    bool FlagMwInput; ///< A boolean flag that is true if the W mass is taken as an input parameter.
    bool FlagQuadraticTerms; ///< A boolean flag that is true if the quadratic terms in cross sections and widths are switched on.

    /**
     * @brief An internal boolean flag that is true if assuming lepton flavour
     * universality.
     */
    const bool FlagLeptonUniversal;

    /**
     * @brief An internal boolean flag that is true if assuming quark flavour
     * universality.
     */
    const bool FlagQuarkUniversal;

};

#endif	/* NPEFFECTIVEGIMR_H */

