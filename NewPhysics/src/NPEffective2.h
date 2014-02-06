/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE2_H
#define	NPEFFECTIVE2_H

#include "NPEffective.h"

/**
 * @class NPEffective2
 * @brief A class for new physics in the form of contributions to the 
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details Base class for an effective Lagrangian extension in the basis of 
 * \cite Barbieri:1999tm, assuming diagonal and non-universal \f${\cal O}_{HQ,HL,HE,HU,HD}\f$
 * and \f${\cal O}_{HQ,HL}^\prime\f$ interactions.
 * 
 * @anchor NPEffective2Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of NPEffective2 are summarized below: 
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
 *   <td class="mod_name">%cHL1p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}^\prime\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^1}\sigma_a\gamma^\mu L^1\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL2p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}^\prime\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^2}\sigma_a\gamma^\mu L^2\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL3p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}^\prime\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^3}\sigma_a\gamma^\mu L^3\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}^\prime\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^1}\sigma_a\gamma^\mu Q^1\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ2p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}^\prime\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^2}\sigma_a\gamma^\mu Q^2\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ3p_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}^\prime\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^3}\sigma_a\gamma^\mu Q^3\right)\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL1_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^1}\gamma^\mu L^1\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL2_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^2}\gamma^\mu L^2\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHL3_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HL}\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HL}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^3}\gamma^\mu L^3\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ1_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^1}\gamma^\mu Q^1\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ2_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^2}\gamma^\mu Q^2\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHQ3_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HQ}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^3}\gamma^\mu Q^3\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE1_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HE}\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HE}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^1}\gamma^\mu E^1\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE2_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HE}\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HE}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^2}\gamma^\mu E^2\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHE3_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HE}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^3}\gamma^\mu E^3\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU1_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HU}\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HU}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^1}\gamma^\mu U^1\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU2_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HU}\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HU}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^2}\gamma^\mu U^2\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHU3_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HU}\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HU}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^3}\gamma^\mu U^3\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD1_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HD}\right)_{11} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HD}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^1}\gamma^\mu D^1\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD2_NP </td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{22} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HD}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^2}\gamma^\mu D^2\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cHD3_NP</td>
 *   <td class="mod_symb">\f$\left(c_{HQ}\right)_{33} \f$</td>
 *   <td class="mod_desc">The coefficient of the operator \f$\left({\cal O}_{HD}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^3}\gamma^\mu D^3\right) \f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%Lambda_NP </td>
 *   <td class="mod_symb">\f$\Lambda_{NP} \f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 * 
 */
class NPEffective2 : public NPEffective {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NNPEffectiveVars = 19;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string NPEffectiveVars[NNPEffectiveVars];

    /**
     * @brief Constructor.
     */
    NPEffective2()
    : NPEffective()
    {
    }

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPEffective2";
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

#endif	/* NPEFFECTIVE2_H */

