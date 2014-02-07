/* 
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPZBBBAR_H
#define	NPZBBBAR_H

#include "NPbase.h"

/**
 * @class NPZbbbar
 * @brief A model class for new physics with non-standard @f$Zb\bar{b}@f$
 * couplings.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is a Model class containing the necessary functions to 
 * work with new physics contributions to electroweak precision observables,
 * in the form of contributions to the neutral current couplings of the bottom
 * quark:
 * @f{eqnarray}{
 * \mathcal{L}
 * &=&
 * \frac{e}{2 s_W c_W}\,
 * Z_\mu\, \bar{b}
 * \left[
 *   \big(g_{V,\mathrm{SM}}^b + \delta g_V^b\big) \gamma_\mu
 *    - \big(g_{A,\mathrm{SM}}^b + \delta g_A^b\big) \gamma_\mu\gamma_5
 * \right] b\,,
 * \\
 * &=&
 * \frac{e}{2 s_W c_W}\sqrt{\rho_{Z,\mathrm{SM}}^b + \delta \rho_Z^b}\,
 * Z_\mu\, \bar{b}
 * \Big\{ \big[ I_3^b - 2Q_b\big(\kappa_{Z,\mathrm{SM}}^b + \delta \kappa_Z^b\big) s_W^2\big]\gamma^\mu
 *   - I_3^b\gamma^\mu\gamma_5 \Big\} b\,.
 * @f}
 * Equivalently, the shifts to the left-handed and right-handed couplings are
 * given by
 * @f[
 * \mathcal{L} =
 * \frac{e}{2s_W c_W}\,
 * Z_\mu\, \bar{b}
 * \left[ \big(g_{R,\mathrm{SM}}^b+\delta g_R^b\big) \gamma_\mu (1 + \gamma_5)
 * + \big(g_{L,\mathrm{SM}}^b+\delta g_L^b\big) \gamma_\mu (1 - \gamma_5)
 * \right] b
 * @f]
 * with
 * @f[
 *  \delta g_R^b = \frac{\delta g_V^b - \delta g_A^b}{2}\,,\qquad
 *  \delta g_L^b = \frac{\delta g_V^b + \delta g_A^b}{2}\,.
 * @f]
 *
 * 
 * @anchor NPZbbbarInitialization
 * <h3>Initialization</h3>
 *
 * The constructor NPZbbbar() initializes the model flags explained below to their
 * default values. After creating an instance of the current class,
 * it is required to call the initialization method InitializeModel(), which
 * allocates memory to the pointer #myEWSM inherited from StndardModel, by using
 * type EWNPZbbbar. 
 * This pointer is then used in computing the fermionic neutral-current
 * couplings in the presence of new physics contribution to @f$Zb\bar{b}@f$.
 * In a Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPZbbbarParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPZbbbar are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGVb</td>
 *   <td class="mod_symb">\f$\delta g_{V}^b\f$</td>
 *   <td class="mod_desc">New physics New physics contribution to \f$g_{V}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGAb</td>
 *   <td class="mod_symb">\f$\delta g_{A}^b\f$</td>
 *   <td class="mod_desc">New physics New physics contribution to \f$g_{A}^b\f$.</td>
 * </tr>
 * </table>
 * 
 * @anchor NPZbbbarFlags
 * <h3>%Model Flags</h3>
 *
 * The flags of %NPZbbbar are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%NPZbbbarLR</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if \f$\delta g_{L}^b\f$ and
 *   \f$\delta g_{R}^b\f$ are used for the model parameters, instead of
 *   \f$\delta g_{V}^b\f$ and \f$\delta g_{A}^b\f$, respectively. Note that
 *   the label %deltaGVb (%deltaGAb) represents \f$\delta g_{L}^b\f$
 *   (\f$\delta g_{R}^b\f$) in a model configuration file when the flag is TRUE.
 *   The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%NotLinearizedNP</td>
 *   <td class="mod_valu">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if new physics contributions
 *   to electroweak observables are not linearized in \f$\delta g^b\f$.
 *   The default value is FALSE; new physics contributions are linearized.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPZbbbarReimplemented
 * <h3>Reimplemented quantities</h3>
 * 
 * Compared to the base class NPbase, the methods for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()).
 *
 */
class NPZbbbar : public NPbase  {
public:

    /**
     * @brief The number of the model parameters in %NPZbbbar.
     */
    static const int NZbbbarVars = 2;

    /**
     * @brief A string array containing the labels of the model parameters in %NPZbbbar.
     */
    static const std::string ZbbbarVars[NZbbbarVars];

    /**
     * @brief Constructor.
     */
    NPZbbbar();

    /**
     * @brief A method to fetch the name of %NPZbbbar.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPZbbbar";
    }

    /**
     * @brief A method to initialize %NPZbbbar.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers #myEWSM and #myStandardModelMatching, which are used
     * for %EW precision and flavour observables, respectively.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);  
    
    /**
     * @brief The update method for %NPZbbbar.
     * @details This method updates all the model parameters with giving DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for %NPZbbbar
     * have been provided in the model configuration file.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPZbbbar.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to check the sanity of the set of model flags.
     * @return a boolean that is true if the set of model flags is sane
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to retrieve the model flag NotLinearizedNP.
     * @details See @ref NPZbbbarFlags for detail.
     * @return a boolean that is true if new physics contributions
     * to electroweak observables are not linearized in \f$\delta g^b\f$.
     */
    bool IsFlagNotLinearizedNP() const
    {
        return FlagNotLinearizedNP;
    }

    
    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief New physics contribution to @f$g_V^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_V^l@f$ (zero in the current model)
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @brief New physics contribution to @f$g_V^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\delta g_V^q@f$ (non-zero only for \f$q=b\f$)
     */
    virtual double deltaGVq(StandardModel::quark q) const;
    
    /**
     * @brief New physics contribution to @f$g_A^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_A^l@f$ (zero in the current model)
     */
    virtual double deltaGAl(StandardModel::lepton l) const;
    
    /**
     * @brief New physics contribution to @f$g_A^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\delta g_A^q@f$ (non-zero only for \f$q=b\f$)
     */
    virtual double deltaGAq(StandardModel::quark q) const;
        
    
    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief A method to set the value of a parameter of %NPZbbbar.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
    

    ////////////////////////////////////////////////////////////////////////
private:

    /* These variables may be used as the deviations in the left-handed 
     * and right-handed couplings if the flag "NPZbbbarLR" is set to true.
     * Therefore, they should not be used directly. Instead, the functions
     * deltaGVq() and deltaGAq() have to be called. */
    double myDeltaGVb;///< New physics contribution to \f$g_{V}^b\f$.
    double myDeltaGAb;///< New physics contribution to \f$g_{A}^b\f$.

    /**
     * @brief A boolean flag that is true if \f$\delta g_{L,R}^b\f$ are used
     * instead of \f$\delta g_{V,A}^b\f$.
     * @details If this flag is set to true,
     *    @li the label myDeltaGVb denotes @f$\delta g_L^b@f$;
     *    @li the label myDeltaGAb denotes @f$\delta g_R^b@f$.
     */
    bool FlagNPZbbbarLR;

    /**
     * @brief A boolean flag that is true if new physics contributions
     * to electroweak observables are not linearized in \f$\delta g^b\f$.
     */
    bool FlagNotLinearizedNP;

};

#endif	/* NPZBBBAR_H */

