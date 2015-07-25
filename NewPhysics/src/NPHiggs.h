/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGS_H
#define	NPHIGGS_H

#include "NPbase.h"

/**
 * @class NPHiggs
 * @brief A model class for new physics in the form of an electroweak chiral
 * Lagrangian with a light Higgs-like scalar.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class describing new physics models based on the
 * electroweak chiral Lagrangian with a light neutral scalar \f$h\f$. 
 * The effective Lagrangian expansion is truncated at the level of two-derivatives
 * \cite Contino:2010mh,
 * 
 * \f{eqnarray}{
 *   {\cal L} &=&
 *   \frac {1}{2} (\partial_\mu h)^2 - V(h)
 *   + \frac{v^2}{4}\mathrm{Tr}(D_\mu \Sigma^\dagger D^\mu \Sigma)
 *     \left(1+2a\frac{h}{v}+b\frac{h^2}{v^2}+\cdots\right)
 *  \\
 *  &&
 *   - \left[
 *     m_{u,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}u_{R,i}\\0\end{array}\right)\left(1+c_u \frac{h}{v}+\cdots\right)
 *     + m_{d,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}0\\d_{R,i}\end{array}\right)\left(1+c_d \frac{h}{v}+\cdots\right)
 *     + m_{e,i}\overline{L}_{L,i}\Sigma \left(\begin{array}{c}0\\e_{R,i}\end{array}\right)\left(1+c_e \frac{h}{v}+\cdots\right)
 *     + \mathrm{h.c.}
 * \right],
 * \f}
 * where @f$V(h)@f$ is given by 
 * \f[ V(h)
 *   = \frac{m_h}{2}h^2
 *   + \frac{d_3}{6}\left(\frac{3m_h^2}{v}\right)\,h^3
 *   + \frac{d_4}{24}\left(\frac{3m_h^2}{v^2}\right)\,h^4 + \cdots\,.
 * \f]
 *
 * The SM corresponds to the choice @f$a=b=c_u=c_d=c_e=d_3=d_4=1@f$.
 * The dominant deviations from the SM in the electroweak precision observables
 * are induced by the non-standard @f$HVV@f$ coupling @f$a\neq 1@f$. 
 * This generates extra contributions to the oblique parameters
 * @cite Barbieri:2007bh ,
 * @f{eqnarray}{
 * S &=& \frac{1}{12\pi} (1 - a^2)
 *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,,\\
 * T &=& - \frac{3}{16\pi c_W^2} (1 - a^2)
 *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,,
 * @f}
 * where @f$\Lambda@f$ is the cutoff scale of the model.
 *
 *
 * @anchor NPHiggsInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor 
 * NPHiggs(), it is required to call the initialization method
 * InitializeModel(). 
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 * 
 * @anchor NPHiggsParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPHiggs are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a </td>
 *   <td class="mod_symb">\f$a \f$</td>
 *   <td class="mod_desc">The \f$hVV\f$ coupling \f$a\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%b </td>
 *   <td class="mod_symb">\f$b \f$</td>
 *   <td class="mod_desc">The \f$hhVV\f$ coupling \f$b\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_u </td>
 *   <td class="mod_symb">\f$c_u \f$</td>
 *   <td class="mod_desc">The \f$h\overline{u}u\f$ coupling \f$c_u\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_d</td>
 *   <td class="mod_symb">\f$c_d \f$</td>
 *   <td class="mod_desc">The \f$h\overline{d}d\f$ coupling \f$c_d\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_e </td>
 *   <td class="mod_symb">\f$c_e \f$</td>
 *   <td class="mod_desc">The \f$h\overline{e}e\f$ coupling \f$c_e\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d_3 </td>
 *   <td class="mod_symb">\f$d_3\f$</td>
 *   <td class="mod_desc">The \f$hhh\f$ coupling \f$d_3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d_4 </td>
 *   <td class="mod_symb">\f$d_4\f$</td>
 *   <td class="mod_desc">The \f$hhhh\f$ coupling \f$d_4\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%LambdaNP</td>
 *   <td class="mod_symb">\f$\Lambda\f$</td>
 *   <td class="mod_desc">The cutoff scale of the effective Lagrangian. When
 *   the value of this parameter is set to 0 in model initialization,
 *   the cutoff scale is taken to be @f$\Lambda = 4\pi v/\sqrt{|1-a^2|}@f$,
 *   where the @f$W_LW_L@f$ scattering becomes non-perturbative. </td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPHiggsFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPHiggsFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the
 * following quantities are reimplemented in the current class:
 *
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 *
 */
class NPHiggs : public NPbase {
public:

    /**
     * @brief The number of the model parameters in %NPHiggs.
     */
    static const int NNPHIGGSvars = 8;

    /**
     * @brief A string array containing the labels of the model parameters in %NPHiggs.
     */
    static const std::string NPHIGGSvars[NNPHIGGSvars];

    /**
     * @brief The default constructor.
     */
    NPHiggs();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////     

    /**
     * @brief The oblique parameter \f$S\f$.
     * @details
     * @f[ 
     *  S = \frac{1}{12\pi} (1 - a^2) \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,.
     * @f]
     *
     * See @cite Barbieri:2007bh.
     * @return \f$S\f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * @details
     * @f[
     *  T = - \frac{3}{16\pi c_W^2} (1 - a^2)
     *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,.
     * @f]
     *
     * See @cite Barbieri:2007bh.
     * @return \f$T\f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return \f$U=0\f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////    
protected:

    double a; ///< The \f$hVV\f$ coupling \f$a\f$.
    double b; ///< The \f$hhVV\f$ coupling \f$b\f$.
    double c_u; ///< The \f$h\overline{u}u\f$ coupling \f$c_u\f$.
    double c_d; ///< The \f$h\overline{d}d\f$ coupling \f$c_d\f$.
    double c_e; ///< The \f$h\overline{e}e\f$ coupling \f$c_e\f$.
    double d_3; ///< The \f$hhh\f$ coupling \f$d_3\f$.
    double d_4; ///< The \f$hhhh\f$ coupling \f$d_4\f$.
    double LambdaNP_in; ///< The new physics scale \f$\Lambda\f$

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);


};

#endif	/* NPHIGGS_H */

