/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVEBS_H
#define	NPEFFECTIVEBS_H

#include "NPbase.h"
#include <string.h>
#include <stdexcept>


/**
 * @addtogroup NewPhysics
 * @brief A module for model-independent studies of new physics.
 * @{
 */

/**
 * @class NPEffectiveBS
 * @brief A model class for new physics in the form the
 * dimension-six effective Lagrangian in the basis of Barbieri and Strumia \cite Barbieri:1999tm.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the dimension-six effective Lagrangian:
 * @f[
 * \mathcal{L}_\mathrm{eff}
 * = \mathcal{L}_\mathrm{SM}
 *   + \sum_i \frac{C_i}{\Lambda^2} \mathcal{O}_i
 * @f]
 * where only the SM gauge-invariant and flavour-diagonal dimension-six operators
 * that are relevant to the electroweak precision observables are considered:
 * @f{align}{
 * \mathcal{O}_{WB} &= (H^\dagger\tau^a H) W^a_{\mu\nu}B^{\mu\nu},
 * &
 * \mathcal{O}_{H} &= | H^\dagger D_\mu H|^2\,,
 * \nonumber\\
 * \mathcal{O}_{LL} &= \frac{1}{2}(\overline{L}\gamma_\mu\tau^a L)^2\,,
 * &
 * \mathcal{O}^\prime_{HL} &= i(H^\dagger D_\mu\tau^a H)(\overline{L}\gamma^\mu\tau^a L)\,,
 * \nonumber\\
 * \mathcal{O}^\prime_{HQ} &= i(H^\dagger D_\mu\tau^a H)(\overline{Q}\gamma^\mu\tau^a Q)\,,
 * &
 * \mathcal{O}_{HL} &= i(H^\dagger D_\mu H)(\overline{L}\gamma^\mu L)\,,
 * \nonumber\\
 * \mathcal{O}_{HQ} &= i(H^\dagger D_\mu H)(\overline{Q}\gamma^\mu Q)\,,
 * &
 * \mathcal{O}_{HE} &= i(H^\dagger D_\mu H)(\overline{E}\gamma^\mu E)\,,
 * \nonumber\\
 * \mathcal{O}_{HU} &= i(H^\dagger D_\mu H)(\overline{U}\gamma^\mu U)\,,
 * &
 * \mathcal{O}_{HD} &= i(H^\dagger D_\mu H)(\overline{D}\gamma^\mu D)\,,
 * @f}
 * as well as the Hermitian conjugate for operators @f$\mathcal{O}^\prime_{HL}@f$
 * to @f$\mathcal{O}_{HD}@f$. 
 *
 * 
 * @anchor NPEffectiveBSInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPEffectiveBS(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEffectiveBSParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEffectiveBS are summarized below:
 *
 *
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
 *   \f$\mathcal{O}_{WB}= \big(H^\dagger\tau^aH\big)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cH_NP </td>
 *   <td class="mod_symb">\f$C_H\f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{H}= \big|H^\dagger D_\mu H\big|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cLL_NP </td>
 *   <td class="mod_symb">\f$C_{LL}  \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{LL}=\frac 12 \big(\overline{L}\,\gamma^\mu\tau^a L\big)
 *    \big(\overline{L}\,\gamma_\mu\tau^a L\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHLp_NP </td>
 *   <td class="mod_symb">\f$C_{HL}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HL}^\prime=i \big(H^\dagger D_\mu\tau^a H\big)
 *    \big(\overline{L}\,\gamma^\mu\tau^a L\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQp_NP </td>
 *   <td class="mod_symb">\f$C_{HQ}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HQ}^\prime=i \big(H^\dagger D_\mu\tau^a H\big)
 *    \big(\overline{Q}\,\gamma^\mu\tau^a Q\big)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL_NP </td>
 *   <td class="mod_symb">\f$C_{HL} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HL}=i \big(H^\dagger D_\mu H\big)
 *    \big(\overline{L}\,\gamma^\mu L\big) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ_NP </td>
 *   <td class="mod_symb">\f$C_{HQ} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HQ}=i \big(H^\dagger D_\mu H\big)
 *    \big(\overline{Q}\,\gamma^\mu Q\big) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE_NP </td>
 *   <td class="mod_symb">\f$C_{HE} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HE}=i \big(H^\dagger D_\mu H\big)
 *    \big(\overline{E}\,\gamma^\mu E\big) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU_NP </td>
 *   <td class="mod_symb">\f$C_{HU} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HU}=i \big(H^\dagger D_\mu H\big)
 *    \big(\overline{U}\,\gamma^\mu U\big) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD_NP </td>
 *   <td class="mod_symb">\f$C_{HD} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f$\mathcal{O}_{HD}=i \big(H^\dagger D_\mu H\big)
 *    \big(\overline{D}\,\gamma^\mu D\big) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 *
 * The parameters to be used as inputs when lepton or quark universality
 * are not assumed are the following:
 *
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
 *   <td class="mod_name">%cHL1p_NP, cHL2p_NP, cHL3p_NP </td>
 *   <td class="mod_symb">\f$ C_{HL_i}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HL_i}^\prime =i\big(H^\dagger D_\mu\tau^a H\big)
 *   \big(\overline{L_i}\,\gamma^\mu\tau^a L_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1p_NP, %cHQ2p_NP, %cHQ3p_NP </td>
 *   <td class="mod_symb">\f$ C_{HQ_i}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HQ_i}^\prime =i\big(H^\dagger D_\mu\tau^a H\big)
 *   \big(\overline{Q_i}\,\gamma^\mu\tau^a Q_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL1_NP, cHL2_NP, cHL3_NP </td>
 *   <td class="mod_symb">\f$ C_{HL_i} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HL_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{L_i}\,\gamma^\mu L_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1_NP, %cHQ2_NP, %cHQ3_NP </td>
 *   <td class="mod_symb">\f$ C_{HQ_i} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HQ_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{Q_i}\,\gamma^\mu Q_i\big)\f$ for @f$i=1,2,3@f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE1_NP, cHE2_NP, cHE3_NP </td>
 *   <td class="mod_symb">\f$ C_{HE_i} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator
 *   \f${\cal O}_{HE_i} =i\big(H^\dagger D_\mu H\big)
 *   \big(\overline{E_i}\,\gamma^\mu E_i\big)\f$ for @f$i=1,2,3@f$. </td>
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
 * @anchor NPEffectiveBSFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPEffectiveBSFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$v@f$&nbsp;&nbsp;(with v()),
 * @li @f$M_{W}^{\mathrm{tree}}@f$&nbsp;&nbsp;(with Mw_tree()),
 * @li @f$\Delta G@f$&nbsp;&nbsp;(with DeltaGF()),
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()),
 *
 * and also the get methods getCoeff() to retrieve the value of each coefficient.
 * 
 */
class NPEffectiveBS : public NPbase {
public:

    /**
     *　@brief The number of the model parameters in %NPEffectiveBS.
     */
    static const int NNPEffectiveBSVars = 30;

    /**
     *　@brief The number of the model parameters in %NPEffectiveBS with lepton
     * flavour universality.
     */
    static const int NNPEffectiveBSVars_LFU = 19;

    /**
     *　@brief The number of the model parameters in %NPEffectiveBS with quark
     * flavour universality.
     */
    static const int NNPEffectiveBSVars_QFU = 22;

    /**
     *　@brief The number of the model parameters in %NPEffectiveBS with lepton
     * and quark flavour universalities.
     */
    static const int NNPEffectiveBSVars_LFU_QFU = 11;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveBS.
     */
    static const std::string NPEffectiveBSVars[NNPEffectiveBSVars];

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveBS with lepton flavour universality.
     */
    static const std::string NPEffectiveBSVars_LFU[NNPEffectiveBSVars_LFU];

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveBS with quark flavour universality.
     */
    static const std::string NPEffectiveBSVars_QFU[NNPEffectiveBSVars_QFU];

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveBS with lepton and quark flavour universalities.
     */
    static const std::string NPEffectiveBSVars_LFU_QFU[NNPEffectiveBSVars_LFU_QFU];

    /**
     * @brief Constructor.
     * @param[in] FlagLeptonUniversal_in an internal boolean flag that is true if assuming lepton flavour universality
     * @param[in] FlagQuarkUniversal_in an internal boolean flag that is true if assuming quark flavour universality
     */
    NPEffectiveBS(const bool FlagLeptonUniversal_in = false, const bool FlagQuarkUniversal_in = false);

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The SM Higgs vacuum expectation value \f$v\f$.
     * @details
     * @f[
     *    v = \left(\frac{1}{\sqrt{2} G_F}\right)^{1/2}.
     * @f]
     * @return the SM value of the Higgs VEV in GeV
     * 
     */
    virtual double v() const;

    /**
     * @brief @copybrief StandardModel::Mw_tree()
     * @details The tree-level mass is obtained by subtracting the effect
     * of @f$\Delta G@f$. 
     * @return @f$M_W^{\mathrm{tree}}@f$ in GeV.
     */
    virtual double Mw_tree() const;

    /**
     * @brief @copybrief NPbase::DeltaGF()
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimentl value measured through muon decays,
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * The dimension-six operators yield the following corrections:
     * @f[
     *  \Delta G_F
     *  = - C_{LL}\left(\frac{v}{\Lambda}\right)^2
     *    + C^\prime_{HL_1}\left(\frac{v}{\Lambda}\right)^2
     *    + C^\prime_{HL_2}\left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @brief @copybrief NPbase::obliqueS()
     * @details The operator @f$\mathcal{O}_{WB}@f$ contributes to the @f$S@f$
     * parameter:
     * @f[
     * S = \frac{4s_Wc_W\, C_{WB} }{\alpha(M_Z^2)} \left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief @copybrief NPbase::obliqueT()
     * @details The operator @f$\mathcal{O}_{H}@f$ contributes to the @f$T@f$
     * parameter:
     * @f[
     * T = - \frac{C_H}{2\alpha(M_Z^2)}\left(\frac{v}{\Lambda}\right)^2.
     * @f]
     * 
     * See @cite Ciuchini:2013pca and references therein.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief @copybrief NPbase::obliqueU()
     * @return @f$U=0@f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief NPbase::GammaW()
     * @details The width of the @f$W@f$ boson receives the corrections from the
     * operators @f$\mathcal{O}^\prime_{HL}@f$ and @f$\mathcal{O}^\prime_{HQ}@f$,
     * in addition to the corrections via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$:
     * @f[
     * \Gamma_W
     * =
     * \Gamma_{W,\mathrm{SM}}
     * \left[
     * 1 - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *     \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     *   - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     *   + \left( C_{HL_1}^\prime + C_{HL_2}^\prime  + C_{HL_3}^\prime
     *          + C_{HQ_1}^\prime + C_{HQ_2}^\prime\right)
     * \left(\frac{v}{\Lambda}\right)^2
     * \right].
     * @f]
     * 
     * See @cite Ciuchini:2013pca and references therein.
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief @copybrief NPbase::deltaGV_f()
     * @details New physics contribution to the neutral-current vector
     * coupling @f$g_V^f@f$ is given by
     * @f[
     * \delta g_V^f = \delta g_L^f + \delta g_R^f.
     * @f]
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$
     */
    virtual double deltaGV_f(const Particle f) const;

    /**
     * @brief @copybrief NPbase::deltaGA_f()
     * @details New physics contribution to the neutral-current axial-vector
     * coupling @f$g_A^f@f$ is given by
     * @f[
     * \delta g_A^f = \delta g_L^f - \delta g_R^f.
     * @f]
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$
     */
    virtual double deltaGA_f(const Particle f) const;

    ////////////////////////////////////////////////////////////////////////    
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    double cWB; ///< The dimension-6 operator coefficient \f$C_{WB}\f$.
    double cH; ///< The dimension-6 operator coefficient \f$C_{H}\f$.
    double cL1L1; ///< The dimension-6 operator coefficient \f$C_{L_1L_1}\f$.
    double cL1L2; ///< The dimension-6 operator coefficient \f$C_{L_1L_2}\f$.
    double cL1L3; ///< The dimension-6 operator coefficient \f$C_{L_1L_3}\f$.
    double cL2L2; ///< The dimension-6 operator coefficient \f$C_{L_2L_2}\f$.
    double cL2L3; ///< The dimension-6 operator coefficient \f$C_{L_2L_3}\f$.
    double cL3L3; ///< The dimension-6 operator coefficient \f$C_{L_3L_3}\f$.
    double cHL1p; ///< The dimension-6 operator coefficient \f$C_{HL_1}^\prime\f$.
    double cHL2p; ///< The dimension-6 operator coefficient \f$C_{HL_2}^\prime\f$.
    double cHL3p; ///< The dimension-6 operator coefficient \f$C_{HL_3}^\prime\f$.
    double cHQ1p; ///< The dimension-6 operator coefficient \f$C_{HQ_1}^\prime\f$.
    double cHQ2p; ///< The dimension-6 operator coefficient \f$C_{HQ_2}^\prime\f$.
    double cHQ3p; ///< The dimension-6 operator coefficient \f$C_{HQ_3}^\prime\f$.
    double cHL1; ///< The dimension-6 operator coefficient \f$C_{HL_1}\f$.
    double cHL2; ///< The dimension-6 operator coefficient \f$C_{HL_2}\f$.
    double cHL3; ///< The dimension-6 operator coefficient \f$C_{HL_3}\f$.
    double cHQ1; ///< The dimension-6 operator coefficient \f$C_{HQ_1}\f$.
    double cHQ2; ///< The dimension-6 operator coefficient \f$C_{HQ_2}\f$.
    double cHQ3; ///< The dimension-6 operator coefficient \f$C_{HQ_3}\f$.
    double cHE1; ///< The dimension-6 operator coefficient \f$C_{HE_1}\f$.
    double cHE2; ///< The dimension-6 operator coefficient \f$C_{HE_2}\f$.
    double cHE3; ///< The dimension-6 operator coefficient \f$C_{HE_3}\f$.
    double cHU1; ///< The dimension-6 operator coefficient \f$C_{HU_1}\f$.
    double cHU2; ///< The dimension-6 operator coefficient \f$C_{HU_2}\f$.
    double cHU3; ///< The dimension-6 operator coefficient \f$C_{HU_3}\f$.
    double cHD1; ///< The dimension-6 operator coefficient \f$C_{HD_1}\f$.
    double cHD2; ///< The dimension-6 operator coefficient \f$C_{HD_2}\f$.
    double cHD3; ///< The dimension-6 operator coefficient \f$C_{HD_3}\f$.
    double Lambda_NP; ///< The new physics scale \f$\Lambda\f$.

    ////////////////////////////////////////////////////////////////////////
private:

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

    /**
     * @brief New physics contribution to @f$g_L^f@f$.
     * @details New physics contributions to the neutral-current left-handed
     * coupling @f$g_L^f@f$: e.g., 
     * @f[
     * \delta g_L^{\nu_i}
     * = \frac{C_{HL_i}^\prime-C_{HL_i}}{2} \left(\frac{v}{\Lambda}\right)^2,
     * \qquad
     * \delta g_L^{e_i}
     * = -\frac{C_{HL_i}^\prime+C_{HL_i}}{2} \left(\frac{v}{\Lambda}\right)^2,
     * @f]
     *
     * where @f$L\to Q@f$, @f$\nu_i\to u_i@f$ and @f$e_i\to d_i@f$ for quarks.
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_L^f@f$
     *
     * @attention The new physics contribution via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$ are not included in this function.
     */
    double deltaGL_f_tmp(const Particle f) const;

    /**
     * @brief New physics contribution to @f$g_R^f@f$.
     * @details New physics contributions to the neutral-current right-handed
     * coupling @f$g_R^f@f$: e.g., 
     * @f[
     * \delta g_R^{e_i} = -\frac{C_{HE_i}}{2} \left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_R^f@f$
     *
     * @attention The new physics contribution via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$ are not included in this function.
     */
    double deltaGR_f_tmp(const Particle f) const;

};

/**
 * @}
 */

#endif	/* NPEFFECTIVEBS_H */

