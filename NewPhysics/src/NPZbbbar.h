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
 * @brief A class for new physics with non-standard @f$Zb\bar{b}@f$ couplings. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to work with new physics 
 * corrections to electroweak precision observables, in the form of contributions
 * to the neutral current couplings of the bottom quark, \f$\delta g_{V,A}^b\f$ or 
 * \f$\delta g_{L,R}^b\f$.
 *
 * @anchor NPZbbbarParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of NPZbbbar are summarized below: 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myDeltaGVb</td>
 *   <td class="mod_symb">\f$\delta g_{V}^b\f$</td>
 *   <td class="mod_desc">New physics correction to \f$g_{V}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myDeltaGAb</td>
 *   <td class="mod_symb">\f$\delta g_{A}^b\f$</td>
 *   <td class="mod_desc">New physics correction to \f$g_{A}^b\f$.</td>
 * </tr>
 * </table>
 * 
 * @anchor NPZbbbarFlags
 * <h3>%Model Flags</h3>
 *
 * The flags of NPZbbbar are summarized below: 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FlagNPZbbbarLR</td>
 *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">New physics parameterized in terms of \f$\delta g_{L,R}^b\f$. </td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">New physics parameterized in terms of \f$\delta g_{V,A}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FlagNotLinearizedNP</td>
  *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">New physics contributions to electroweak
 * observables not linearized in \f$\delta g^b\f$.</td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">New physics contributions to electroweak
 * observables linear in \f$\delta g^b\f$.</td>
 * </tr>
 * </table>
 *
 */
class NPZbbbar : public NPbase  {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NZbbbarVars = 2;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string ZbbbarVars[NZbbbarVars];

    /**
     * @brief Constructor.
     */
    NPZbbbar();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPZbbbar";
    }

    /**
     * @brief A method to initialize the model.
     * @return true is model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A method to initialize the model.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return true is model initialization is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);  
    
    /**
     * @brief The update method for the model class.
     * @details This method updates all the parameters of the model every time a
     * new set of parameters is generated.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful.
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A set method to fix the flags for the model.
     * @param[in] name the name of the flag
     * @param[in] value the value of the flag that can be true or false
     * @return a boolean to designate the success or failure of this procedure
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to check the sanity of the set of flags.
     * @return true if the set of flags is sane.
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to check the value of the flag \b FlagNotLinearizedNP
     * @return the value of \b FlagNotLinearizedNP
     */
    bool IsFlagNotLinearizedNP() const
    {
        return FlagNotLinearizedNP;
    }

    
    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief The new physics correction to @f$g_V^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^l@f$ 
     * (zero in this model) 
     */
    virtual double deltaGVl(StandardModel::lepton l) const;
    
    /**
     * @brief The new physics correction to @f$g_V^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^q@f$
     * (non-zero only for \f$q=b\f$)
     */
    virtual double deltaGVq(StandardModel::quark q) const;
    
    /**
     * @brief The new physics correction to @f$g_A^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^l@f$
     * (zero in this model)
     */
    virtual double deltaGAl(StandardModel::lepton l) const;
    
    /**
     * @brief The new physics correction to @f$g_A^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^q@f$
     * (non-zero only for \f$q=b\f$)
     */
    virtual double deltaGAq(StandardModel::quark q) const;
        
    
    ////////////////////////////////////////////////////////////////////////
protected:
    /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
    

    ////////////////////////////////////////////////////////////////////////
private:
    /* These variables may be used as the deviations in the left-handed 
     * and right-handed couplings if the flag "NPZbbbarLR" is set to true.
     * Therefore, they should not be used directly. Instead, the functions
     * deltaGVq() and deltaGAq() have to be called. */
    double myDeltaGVb;///< New physics correction to \f$g_{V}^b\f$.
    double myDeltaGAb;///< New physics correction to \f$g_{A}^b\f$.

    /*
     * If true,
     *    myDeltaGVb --> delta g_L^b
     *    myDeltaGAb --> delta g_R^b
     */
    bool FlagNPZbbbarLR;///< Flag: if true uses \f$\delta g_{L,R}^b\f$ instead of \f$\delta g_{V,A}^b\f$. 

    bool FlagNotLinearizedNP;///< Flag: if false uses only correction linear in \f$\delta g^b\f$.

};

#endif	/* NPZBBBAR_H */

