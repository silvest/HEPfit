/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE1_H
#define	NPEFFECTIVE1_H

#include "NPEffective.h"

/**
 * @class NPEffective1
 * @brief A class for new physics in the form of contributions to the 
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details Base class for an effective Lagrangian extension in the basis of 
 * \cite Barbieri:1999tm, assuming diagonal and family-universal \f${\cal O}_{HQ,HL,HE,HU,HD}\f$
 * and \f${\cal O}_{HQ,HL}^\prime\f$ interactions.
 * 
 * @anchor NPEffective1Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of NPEffective1 are summarized below: 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cWB_NP </td>
 *   <td class="mod_symb">\f$c_{WB} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{WB}=\left(H^\dagger\sigma_aH\right)W_{\mu\nu}^a B^{\mu\nu}\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cH_NP </td>
 *   <td class="mod_symb">\f$c_H\f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{H}=\left|H^\dagger D_\mu H\right|^2\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cLL_NP </td>
 *   <td class="mod_symb">\f$c_{LL}  \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{LL}=\frac 12\left(\overline{L}\gamma^\mu\sigma_a L\right)\left(\overline{L}\gamma_\mu\sigma_a L\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHLp_NP </td>
 *   <td class="mod_symb">\f$c_{HL}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HL}^\prime=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L}\sigma_a\gamma^\mu L\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQp_NP </td>
 *   <td class="mod_symb">\f$c_{HQ}^\prime \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HQ}^\prime=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q}\sigma_a\gamma^\mu Q\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL_NP </td>
 *   <td class="mod_symb">\f$c_{HL} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HL}=i\left(H^\dagger D_\mu H\right)\left(\overline{L}\gamma^\mu L\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ_NP </td>
 *   <td class="mod_symb">\f$c_{HQ} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HQ}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q}\gamma^\mu Q\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE_NP </td>
 *   <td class="mod_symb">\f$c_{HE} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HE}=i\left(H^\dagger D_\mu H\right)\left(\overline{E}\gamma^\mu E\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU_NP </td>
 *   <td class="mod_symb">\f$c_{HU} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HU}=i\left(H^\dagger D_\mu H\right)\left(\overline{U}\gamma^\mu U\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD_NP </td>
 *   <td class="mod_symb">\f$c_{HD} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f${\cal O}_{HD}=i\left(H^\dagger D_\mu H\right)\left(\overline{D}\gamma^\mu D\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda_{NP} \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 * 
 */
class NPEffective1 : public NPEffective {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NNPEffectiveVars = 11;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string NPEffectiveVars[NNPEffectiveVars];

   /**
     * @brief Constructor.
     */
    NPEffective1()
    : NPEffective()
    {
    }

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPEffective1";
    }

    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


protected:
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
      
};

#endif	/* NPEFFECTIVE1_H */

