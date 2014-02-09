/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_H
#define	NPEPSILONS_H

#include <EWSM.h>
#include "NPbase.h"

/**
 * @class NPEpsilons
 * @brief A class for new physics in the form of contributions to the
 * \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to work with new physics 
 * corrections to electroweak precision observables, in the form of contributions
 * to the \f$\varepsilon_{1,2,3,b}\f$ parameters \cite Altarelli:1990zd, \cite Altarelli:1991fk,\cite Altarelli:1993sz. 
 * Both SM and new physics contributions to \f$\varepsilon_i\f$ are parameterized. 
 *
 * @anchor NPEpsilonsParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEpsilons are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_1 </td>
 *   <td class="mod_symb">\f$\varepsilon_1\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_1\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_2 </td>
 *   <td class="mod_symb">\f$\varepsilon_2\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_3 </td>
 *   <td class="mod_symb">\f$\varepsilon_3\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_3\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon_b </td>
 *   <td class="mod_symb">\f$\varepsilon_b\f$</td>
 *   <td class="mod_desc">The parameter \f$\varepsilon_b\f$.</td>
 * </tr>
 * </table>
 * 
 * @anchor NPEpsilonsFlags
 * <h3>%Model Flags</h3>
 *
 * The flags of %NPEpsilons are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>Value</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon1SM</td>
 *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">Only the SM value of \f$\varepsilon_1\f$ is considered.</td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">Includes both, the SM and new physics, contributions to \f$\varepsilon_1\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon2SM</td>
 *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">Only the SM value of \f$\varepsilon_2\f$ is considered.</td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">Includes both, the SM and new physics, contributions to \f$\varepsilon_2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilon3SM</td>
 *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">Only the SM value of \f$\varepsilon_3\f$ is considered.</td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">Includes both, the SM and new physics, contributions to \f$\varepsilon_3\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%epsilonbSM</td>
 *   <td class="mod_desc">True</td>
 *   <td class="mod_desc">Only the SM value of \f$\varepsilon_b\f$ is considered.</td>
 * <tr>
 *   <td class="mod_name"> </td>
 *   <td class="mod_desc">False</td>
 *   <td class="mod_desc">Includes both, the SM and new physics, contributions to \f$\varepsilon_b\f$.</td>
 * </tr>
 * </table>
 *
 * 
 */
class NPEpsilons : public NPbase  {
public:

    /**
     * @brief The number of the model parameters in %NPEpsilons.
     */
    static const int NEPSILONvars = 4;

    /**
     * @brief A string array containing the labels of the model parameters in %NPEpsilons.
     */
    static const std::string EPSILONvars[NEPSILONvars];
    
    /**
     * @brief The default constructor.
     */
    NPEpsilons();

    /**
     * @brief @copybrief Model::ModelName()
     * @copydetails Model::ModelName()
     */
    virtual std::string ModelName() const 
    {
        return "NPEpsilons";
    }

    /**
     * @brief @copybrief StandardModel::InitializeModel()
     * @details This method allocates memory to the pointer #myEWSM, inherited
     * from StndardModel, with type EWNPEpsilons.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief @copybrief Model::Init()
     * @copydetails Model::Init()
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief @copybrief Model::Update()
     * @copydetails Model::Update()
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
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
    
    /**
     * @brief @copybrief Model::CheckFlags()
     * @copydetails Model::CheckFlags()
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the SM value (FlagEpsilon1SM=True)
     * or the SM plus new physics value (FlagEpsilon1SM=False) of \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the SM value (FlagEpsilon2SM=True)
     * or the SM plus new physics value (FlagEpsilon2SM=False) of \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const;

    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the SM value (FlagEpsilon3SM=True)
     * or the SM plus new physics value (FlagEpsilon3SM=False) of \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const;
 
    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value (FlagEpsilonbSM=True)
     * or the SM plus new physics value (FlagEpsilonbSM=False) of \f$\varepsilon_b\f$
     */
    virtual double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The \f$W\f$ boson mass.
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /** 
     * @brief The \f$W\f$ decay width \f$\Gamma_W\f$.
     * @return the total width of the \f$W\f$ boson in GeV [NOT IMPLEMENTED YET]
     */
    virtual double GammaW() const;

    
    ////////////////////////////////////////////////////////////////////////   
protected:    

    double myEpsilon_1;///< The parameter \f$\varepsilon_1\f$.
    double myEpsilon_2;///< The parameter \f$\varepsilon_2\f$.
    double myEpsilon_3;///< The parameter \f$\varepsilon_3\f$.
    double myEpsilon_b;///< The parameter \f$\varepsilon_b\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);
    
    
    ////////////////////////////////////////////////////////////////////////         
private:

    bool FlagEpsilon1SM;///< Flag: if true only the SM value of \f$\varepsilon_1\f$ is considered.
    bool FlagEpsilon2SM;///< Flag: if true only the SM value of \f$\varepsilon_2\f$ is considered.
    bool FlagEpsilon3SM;///< Flag: if true only the SM value of \f$\varepsilon_3\f$ is considered. 
    bool FlagEpsilonbSM;///< Flag: if true only the SM value of \f$\varepsilon_b\f$ is considered. 

    
};

#endif	/* NPEPSILONS_H */

