/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE2_H
#define	NPEFFECTIVE2_H

#include "NPEffective.h"

/**
 * @class NPEffective2
 * @brief A model class for new physics in the form of the
 * dimension-six effective Lagrangian, with lepton-flavour universality and
 * without quark-flavour universality.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the dimension-six effective Lagrangian introduced in NPEffective
 * class, where lepton-flavour universality is assumed, but quark-flavour universality
 * is not.
 *
 *
 * @anchor NPEffective2Initialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPEffective2(), it is required to call the initialization method
 * InitializeModel(), inherited from NPEffective class.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEffective2Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEffective2 are summarized below:
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
 * @anchor NPEffective2Flags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPEffective2Functions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 *
 */
class NPEffective2 : public NPEffective {
public:

    /**
     * @brief The number of the model parameters in %NPEffective2.
     */
    static const int NNPEffectiveVars = 19;

    /**
     * @brief A string array containing the labels of the model parameters in %NPEffective2.
     */
    static const std::string NPEffectiveVars[NNPEffectiveVars];

    /**
     * @brief The default constructor.
     */
    NPEffective2();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);


};

#endif	/* NPEFFECTIVE2_H */

