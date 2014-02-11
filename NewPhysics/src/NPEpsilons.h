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
 * @brief A model class for new physics in the form of contributions to the
 * \f$\varepsilon_{1,2,3,b}\f$ parameters.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing parameters and functions to work
 * with the epsilon paramterization 
 * \cite Altarelli:1990zd, \cite Altarelli:1991fk,\cite Altarelli:1993sz
 * for new physics contributions to the electroweak precision observables.
 * It is noted that the \f$\varepsilon_{1,2,3,b}\f$ parameters include both the
 * SM and new physics contributions.
 *
 *
 * @anchor NPEpsilonsInitialization
 * <h3>Initialization</h3>
 *
 * The constructor NPEpsilons() initializes the model flags explained below to
 * their default values. After creating an instance of the current class,
 * it is required to call the initialization method InitializeModel(), which
 * allocates memory to the pointer #myEWSM, inherited from StndardModel, with
 * type EWNPEpsilons.
 * This pointer is then used in computing the \f$W\f$-boson mass and the
 * fermionic neutral-current couplings with the epsilon parameters.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
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
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_1\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilon2SM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_2\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilon3SM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_3\f$. The default value is FALSE.</td>
 * <tr>
 *   <td class="mod_name">%epsilonbSM</td>
 *   <td class="mod_desc">TRUE&nbsp;/&nbsp;<b>FALSE</b></td>
 *   <td class="mod_desc">This flag is set to TRUE if only the SM contribution
 *   is considered for \f$\varepsilon_b\f$. The default value is FALSE.</td>
 * </table>
 *
 *
 * @anchor NPEpsilonsFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$M_W@f$&nbsp;&nbsp; (with Mw()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 *
 * where the latter is not available, since it cannot be simply parameterized
 * by the epsilon parameters.
 * In addition, the functions for the epsilon parameters are also provided:
 * 
 * @li @f$\varepsilon_1@f$, @f$\varepsilon_2@f$, @f$\varepsilon_3@f$ and
 * @f$\varepsilon_b@f$&nbsp;&nbsp;
 * (with epsilon1(), epsilon2(), epsilon3() and epsilonb()).
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
     * @return the SM value (FlagEpsilon1SM=true) or the SM plus new physics
     * value (FlagEpsilon1SM=false) of \f$\varepsilon_1\f$
     */
    double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the SM value (FlagEpsilon2SM=true) or the SM plus new physics
     * value (FlagEpsilon2SM=false) of \f$\varepsilon_2\f$
     */
    double epsilon2() const;

    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the SM value (FlagEpsilon3SM=true) or the SM plus new physics
     * value (FlagEpsilon3SM=false) of \f$\varepsilon_3\f$
     */
    double epsilon3() const;
 
    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value (FlagEpsilonbSM=true) or the SM plus new physics
     * value (FlagEpsilonbSM=false) of \f$\varepsilon_b\f$
     */
    double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief @copybrief StandardModel::Mw()
     * @details This function calls EWNPEpsilons::Mw() via
     * EWNPEpsilons::Mw_NPEpsilons().
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /** 
     * @brief @copybrief StandardModel::GammaW()
     * 
     * @warning This function is not available.
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

    bool FlagEpsilon1SM;///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_1\f$.
    bool FlagEpsilon2SM;///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_2\f$.
    bool FlagEpsilon3SM;///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_3\f$.
    bool FlagEpsilonbSM;///< A boolean flag that is true if only the SM contribution is considered for \f$\varepsilon_b\f$.

    
};

#endif	/* NPEPSILONS_H */

