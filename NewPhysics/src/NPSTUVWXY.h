/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTUVWXY_H
#define	NPSTUVWXY_H

#include "NPbase.h"
#include <cmath>

/**
 * @class NPSTUVWXY
 * @brief A model class for new physics in the form of contributions to the extended
 * oblique parameters \f$\hat{S},~\hat{T},~\hat{U},~V,~W,~X\f$ and \f$Y\f$.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing the necessary functions to compute
 * new physics contributions to the electroweak precision observables and so forth,
 * with the extended set of oblique parameters \f$\hat{S},~\hat{T},~\hat{U},~V,~W,~X\f$
 * and \f$Y\f$ \cite Barbieri:2004qk.
 *
 * The hatted parameters are related to the original Peskin-Takeuchi parameters
 * \cite Peskin:1990zt, \cite Peskin:1991sw through the relations:
 * @f[
 * \widehat{S} = \frac{\alpha}{4s_W^2}\,S\,,\qquad
 * \widehat{T} = \alpha\,T\,,\qquad
 * \widehat{U} = -\frac{\alpha}{4s_W^2}\,U\,.
 * @f]
 * Moreover, the following combinations of the parameters are relevant to the 
 * electroweak precision observables: 
 * @f{eqnarray}{
 * \frac{\alpha}{4s_W^2}\, S'
 * &=&
 * \widehat{S} - W + \frac{X}{s_Wc_W} - Y\,,
 * \\
 * \alpha\,T'
 * &=&
 * \widehat{T} - W + \frac{2s_W}{c_W}\,X - \frac{s_W^2}{c_W^2}\,Y\,,
 * \\
 * \frac{\alpha}{4s_W^2}\, U'
 * &=& - \widehat{U} + V + W - \frac{2s_W}{c_W}\,X\,.
 * @f}
 *
 *
 * @anchor NPSTUVWXYInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPSTUVWXY(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPSTUVWXYParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPSTUVWXY are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueShat </td>
 *   <td class="mod_symb">\f$\hat{S}\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$\hat{S}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueThat </td>
 *   <td class="mod_symb">\f$\hat{T} \f$</td>
 *   <td class="mod_desc">The oblique parameter \f$\hat{T}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueUhat </td>
 *   <td class="mod_symb">\f$\hat{U}\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$\hat{U}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueV </td>
 *   <td class="mod_symb">\f$V\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$V\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueW </td>
 *   <td class="mod_symb">\f$W\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$W\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueX </td>
 *   <td class="mod_symb">\f$X\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$X\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%obliqueY </td>
 *   <td class="mod_symb">\f$Y\f$</td>
 *   <td class="mod_desc">The oblique parameter \f$Y\f$.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPSTUVWXYFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPSTUVWXYFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the
 * following quantities are reimplemented in the current class:
 *
 * @li @f$S'@f$, @f$T'@f$ and @f$U'@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()).
 *
 * It is noted that these functions do not represent the original Peskin-Takeuchi
 * parameters @f$S@f$, @f$T@f$ and @f$U@f$.
 * In addition, the functions for the extended oblique parameters as well as
 * the epsilon parameters are also provided:
 *
 * @li @f$\hat{S}@f$, @f$\hat{T}@f$, @f$\hat{U}@f$,
 * @f$V@f$, @f$W@f$, @f$X@f$ and @f$Y@f$&nbsp;&nbsp;
 * (with obliqueShat(), obliqueThat(), obliqueUhat(),
 * obliqueV(), obliqueW(), obliqueX() and obliqueY()),
 * @li @f$\varepsilon_1@f$, @f$\varepsilon_2@f$, @f$\varepsilon_3@f$ and
 * @f$\varepsilon_b@f$&nbsp;&nbsp;
 * (with epsilon1(), epsilon2(), epsilon3() and epsilonb()).
 *
 */
class NPSTUVWXY : public NPbase {
public:

    /**
     * @brief The number of the model parameters in %NPSTUVWXY.
     */
    static const int NSTUVWXYvars = 7;

    /**
     * @brief A string array containing the labels of the model parameters in %NPSTUVWXY.
     */
    static const std::string STUVWXYvars[NSTUVWXYvars];

    /**
     * @brief The default constructor.
     */
    NPSTUVWXY();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The oblique parameter \f$\hat{S}\f$.
     * @return the value of \f$\displaystyle\hat{S}\f$
     */
    virtual double obliqueShat() const
    {
        return myObliqueShat;
    }

    /**
     * @brief The oblique parameter \f$\hat{T}\f$.
     * @return the value of \f$\displaystyle\hat{T}\f$
     */
    virtual double obliqueThat() const
    {
        return myObliqueThat;
    }

    /**
     * @brief The oblique parameter \f$\hat{U}\f$.
     * @return the value of \f$\displaystyle\hat{U}\f$
     */
    virtual double obliqueUhat() const
    {
        return myObliqueUhat;
    }

    /**
     * @brief The oblique parameter \f$V\f$.
     * @return the value of \f$V\f$
     */
    virtual double obliqueV() const
    {
        return myObliqueV;
    }

    /**
     * @brief The oblique parameter \f$W\f$.
     * @return the value of \f$W\f$
     */
    virtual double obliqueW() const
    {
        return myObliqueW;
    }

    /**
     * @brief The oblique parameter \f$X\f$.
     * @return the value of \f$X\f$
     */
    virtual double obliqueX() const
    {
        return myObliqueX;
    }

    /**
     * @brief The oblique parameter \f$Y\f$.
     * @return the value of \f$Y\f$
     */
    virtual double obliqueY() const
    {
        return myObliqueY;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief NPSTU::epsilon1()
     * @copydetails NPSTU::epsilon1()
     */
    double epsilon1() const;

    /**
     * @brief @copybrief NPSTU::epsilon2()
     * @copydetails NPSTU::epsilon2()

     */
    double epsilon2() const;

    /**
     * @brief @copybrief NPSTU::epsilon3()
     * @copydetails NPSTU::epsilon3()

     */
    double epsilon3() const;

    /**
     * @brief @copybrief NPSTU::epsilonb()
     * @copydetails NPSTU::epsilonb()
     */
    double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////     
    // Combinations of the extended oblique parameters

    /**
     * @brief The combination of extended oblique parameters, \f$S'\f$.
     * @return the value of \f$S'\f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The combination of extended oblique parameters, \f$T'\f$.
     * @return the value of \f$T'\f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The combination of extended oblique parameters, \f$U'\f$.
     * @return the value of \f$U'\f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief NPbase::GammaW()
     * @details
     * @f[
     * \Gamma_W = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S' - 2c_W^2\,T' - \frac{c_W^2-s_W^2}{2s_W^2}\,U'
     *         - 2 (c_W^2 - s_W^2)\, \overline{W} \right)
     * - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     * \right],
     * @f]
     * where @f$\alpha(M_Z^2)\,\overline{W} = V - W@f$.
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////
protected:

    double myObliqueShat; ///< The oblique parameter \f$\hat{S}\f$.
    double myObliqueThat; ///< The oblique parameter \f$\hat{T}\f$.
    double myObliqueUhat; ///< The oblique parameter \f$\hat{U}\f$.
    double myObliqueV; ///< The oblique parameter \f$V\f$.
    double myObliqueW; ///< The oblique parameter \f$W\f$.
    double myObliqueX; ///< The oblique parameter \f$X\f$.
    double myObliqueY; ///< The oblique parameter \f$Y\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);


};

#endif	/* NPSTUVWXY_H */

