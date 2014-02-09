/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEPSILONS_PURENP_H
#define	NPEPSILONS_PURENP_H

#include "NPbase.h"

/**
 * @class NPEpsilons_pureNP
 * @brief A class for new physics in the form of contributions to the
 * \f$\delta\varepsilon_{1,2,3,b}\f$ parameters.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to work with new physics 
 * corrections to electroweak precision observables, in the form of contributions
 * to the \f$\varepsilon_{1,2,3,b}\f$ parameters \cite Altarelli:1990zd, \cite Altarelli:1991fk,\cite Altarelli:1993sz. 
 * Only new physics contributions to \f$\varepsilon_i\f$ are parameterized.
 *
 * @anchor NPEpsilons_pureNPParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPEpsilons_pureNP are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delEps_1 </td>
 *   <td class="mod_symb">\f$\delta_{NP}\varepsilon_1\f$</td>
 *   <td class="mod_desc">The new physics contribution to \f$\varepsilon_1\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delEps_2 </td>
 *   <td class="mod_symb">\f$\delta_{NP}\varepsilon_2\f$</td>
 *   <td class="mod_desc">The new physics contribution to \f$\varepsilon_2\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delEps_3 </td>
 *   <td class="mod_symb">\f$\delta_{NP}\varepsilon_3\f$</td>
 *   <td class="mod_desc">The new physics contribution to \f$\varepsilon_3\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%delEps_b </td>
 *   <td class="mod_symb">\f$\delta_{NP}\varepsilon_b\f$</td>
 *   <td class="mod_desc">The new physics contribution to \f$\varepsilon_b\f$.</td>
 * </tr>
 * </table>
 * 
 */
class NPEpsilons_pureNP : public NPbase {
public:
    
    /**
     * @brief The number of the model parameters in %NPEpsilons_pureNP.
     */
    static const int NEPSILONpureNPvars = 4;
    
    /**
     * @brief A string array containing the labels of the model parameters in %NPEpsilons_pureNP.
     */
    static const std::string EPSILONpureNPvars[NEPSILONpureNPvars];

    /**
     * @brief The default constructor.
     */
    NPEpsilons_pureNP();

    /**
     * @brief @copybrief Model::ModelName()
     * @copydetails Model::ModelName()
     */
    virtual std::string ModelName() const
    {
        return "NPEpsilons_pureNP";
    }

    /**
     * @brief @copybrief StandardModel::InitializeModel()
     * @copydetails NPbase::InitializeModel()
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
     * @return the pure new physics contribution to \f$\varepsilon_1\f$
     */
    virtual double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the pure new physics contribution to \f$\varepsilon_2\f$
     */
    virtual double epsilon2() const;
    
    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the pure new physics contribution to \f$\varepsilon_3\f$
     */
    virtual double epsilon3() const;
    
    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the pure new physics contribution to \f$\varepsilon_b\f$
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

    /**
     * @brief The new physics correction to @f$g_V^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current vector coupling @f$g_A^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_V^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^q@f$
     */  
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @brief The new physics correction to @f$g_A^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^l@f$
     */  
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_A^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^q@f$
     */
    virtual double deltaGAq(StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////
protected:

    double deltaEps_1;///< The new physics contribution to \f$\varepsilon_1\f$.
    double deltaEps_2;///< The new physics contribution to \f$\varepsilon_2\f$.
    double deltaEps_3;///< The new physics contribution to \f$\varepsilon_3\f$.
    double deltaEps_b;///< The new physics contribution to \f$\varepsilon_b\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    
};

#endif	/* NPEPSILONS_PURENP_H */

