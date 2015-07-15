/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVEGIMR_H
#define	NPEFFECTIVEGIMR_H

#include <string.h>
#include <stdexcept>
#include <gslpp.h>
#include "NPbase.h"

using namespace gslpp;


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
    static const int NNPEffectiveGIMRVars = 121;

    /**
     * @brief A string array containing the labels of the model parameters in
     * %NPEffectiveGIMR.
     */
    static const std::string NPEffectiveGIMRVars[NNPEffectiveGIMRVars];

    /**
     *　@brief The number of the model parameters in %NPEffectiveGIMR
     * with lepton and quark flavour universalities.
     */
    static const int NNPEffectiveGIMRVars_LFU_QFU = 25;

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
     * @brief 
     * @return
     */
    virtual bool PostUpdate();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief @copybrief Model::setFlag()
     * @copydetails Model::setFlag()
     */
    virtual bool setFlag(const std::string name, const bool value);


    ////////////////////////////////////////////////////////////////////////

    virtual double DeltaGF() const;

    virtual double obliqueS() const;

    virtual double obliqueT() const;

    virtual double obliqueU() const;

    virtual double Mw() const;

    virtual double GammaW() const;

    virtual double deltaGV_f(const Particle p) const;

    virtual double deltaGA_f(const Particle p) const;

    // no generation mixing
    double deltaGL_f(const Particle p) const;

    // no generation mixing
    double deltaGR_f(const Particle p) const;


    ////////////////////////////////////////////////////////////////////////

    // no generation mixing
    complex deltaGL_Wff(const Particle pbar, const Particle p) const;
    // no generation mixing
    complex deltaGR_Wff(const Particle pbar, const Particle p) const;

    double deltaG_hgg() const;
    double deltaG1_hWW() const;
    double deltaG2_hWW() const;
    double deltaG3_hWW() const;
    double deltaG1_hZZ() const;
    double deltaG2_hZZ() const;
    double deltaG3_hZZ() const;
    double deltaG1_hZA() const;
    double deltaG2_hZA() const;
    double deltaG_hAA() const;

    // no generation mixing
    complex deltaG_hff(const Particle p) const;

    // no generation mixing
    complex deltaGL_Wffh(const Particle pbar, const Particle p) const;
    // no generation mixing
    complex deltaGR_Wffh(const Particle pbar, const Particle p) const;
    // no generation mixing
    double deltaGL_Zffh(const Particle p) const;
    // no generation mixing 
    double deltaGR_Zffh(const Particle p) const;


    ////////////////////////////////////////////////////////////////////////

    complex f_triangle(const double tau) const;
    complex AH_f(const double tau) const;

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double mueeZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
    virtual double muggHpttH(const double sqrt_s) const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double computeGammaTotalRatio() const;

    double GammaHggRatio() const;
    double GammaHWWRatio() const;
    double GammaHZZRatio() const;
    double GammaHZgaRatio() const;
    double GammaHgagaRatio() const;
    double GammaHtautauRatio() const;
    double GammaHccRatio() const;
    double GammaHbbRatio() const;


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
    double Lambda_NP;
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
    complex CHud_diag(const Particle u) const;
    complex CfH_diag(const Particle f) const;


    ////////////////////////////////////////////////////////////////////////
private:

    bool FlagMwInput; ///< A boolean flag that is true if the W mass is taken as an input parameter.

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

