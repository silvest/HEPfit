/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE1_H
#define	NPEFFECTIVE1_H

#include "NPEffective.h"

/**
 * @class NPEffective1
 * @brief A model class for new physics in the form the
 * dimension-six effective Lagrangian, with lepton-flavour and quark-flavour
 * universalities. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions
 * associated with the dimension-six effective Lagrangian introduced in NPEffective
 * class, where both lepton-flavour and quark-flavour universalities are assumed. 
 * 
 *
 * @anchor NPEffective1Initialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPEffective1(), it is required to call the initialization method
 * InitializeModel(), inherited from NPEffective class.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPEffective1Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEffective1 are summarized below:
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
 * 
 * @anchor NPEffective1Flags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPEffective1Functions
 * <h3>Important member functions</h3>
 *
 * See the base classes of the current class.
 * 
 */
class NPEffective1 : public NPEffective {
public:

    /**
     * @brief The number of the model parameters in %NPEffective1.
     */
    static const int NNPEffectiveVars = 11;

    /**
     * @brief A string array containing the labels of the model parameters in %NPEffective1.
     */
    static const std::string NPEffectiveVars[NNPEffectiveVars];

    /**
     * @brief The default constructor.
     */
    NPEffective1();

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

#endif	/* NPEFFECTIVE1_H */

