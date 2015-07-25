/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTU_H
#define	NPSTU_H

#include "NPbase.h"

/**
 * @class NPSTU
 * @brief A model class for new physics in the form of contributions to the oblique
 * parameters \f$S,~T\f$ and \f$U\f$. 
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing the necessary functions to compute
 * new physics contributions to the electroweak precision observables with the
 * Peskin-Takeuchi oblique parameters \cite Peskin:1990zt, \cite Peskin:1991sw.
 *
 *
 * @anchor NPSTUInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPSTU(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 * 
 * @anchor NPSTUParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPSTU are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueS </td>
 *   <td class="mod_symb">\f$S \f$</td>
 *   <td class="mod_desc">The oblique parameter \f$S\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueT </td>
 *   <td class="mod_symb">\f$T \f$</td>
 *   <td class="mod_desc">The oblique parameter \f$T\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueU </td>
 *   <td class="mod_symb">\f$U \f$</td>
 *   <td class="mod_desc">The oblique parameter \f$U\f$.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPSTUFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPSTUFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the
 * following quantities are reimplemented in the current class:
 *
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 *
 * In addition, the functions for the epsilon parameters are also provided:
 *
 * @li @f$\varepsilon_1@f$, @f$\varepsilon_2@f$, @f$\varepsilon_3@f$ and
 * @f$\varepsilon_b@f$&nbsp;&nbsp;
 * (with epsilon1(), epsilon2(), epsilon3() and epsilonb()).
 * 
 */
class NPSTU : public NPbase {
public:

    /**
     * @brief The number of the model parameters in %NPSTU.
     */
    static const int NSTUvars = 3;

    /**
     * @brief A string array containing the labels of the model parameters in %NPSTU.
     */
    static const std::string STUvars[NSTUvars];

    /**
     * @brief The default constructor.
     */
    NPSTU();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief NPbase::obliqueS()
     * @copydetails NPbase::obliqueS()
     */
    virtual double obliqueS() const
    {
        return myObliqueS;
    }

    /**
     * @brief @copybrief NPbase::obliqueT()
     * @copydetails NPbase::obliqueT()
     */
    virtual double obliqueT() const
    {
        return myObliqueT;
    }

    /**
     * @brief @copybrief NPbase::obliqueU()
     * @copydetails NPbase::obliqueU()
     */
    virtual double obliqueU() const
    {
        return myObliqueU;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the value of the @f$\varepsilon_1@f$ parameter
     * (SM plus new physics contributions)
     */
    double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the value of the @f$\varepsilon_2@f$ parameter
     * (SM plus new physics contributions)
     */
    double epsilon2() const;

    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the value of the @f$\varepsilon_3@f$ parameter
     *  (SM plus new physics contributions)
     */
    double epsilon3() const;

    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value of the @f$\varepsilon_b@f$ parameter
     */
    double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////
protected:

    double myObliqueS; ///< The oblique parameter \f$S\f$.
    double myObliqueT; ///< The oblique parameter \f$T\f$.
    double myObliqueU; ///< The oblique parameter \f$U\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);


};

#endif	/* NPSTU_H */

