/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPZBBBAR_H
#define	NPZBBBAR_H

#include "NPbase.h"

/**
 * @class NPZbbbar
 * @brief A model class for new physics in the form of contributions to the
 * @f$Zb\bar{b}@f$ couplings.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is a Model class containing the necessary functions to
 * work with new physics contributions to electroweak precision observables,
 * in the form of contributions to the neutral current couplings of the bottom
 * quark:
 * @f[
 * \mathcal{L}
 * =
 * \frac{e}{2 s_W c_W}\,
 * Z_\mu\, \bar{b}
 * \left[
 *   \big(g_{V,\mathrm{SM}}^b + \delta g_V^b\big) \gamma_\mu
 *    - \big(g_{A,\mathrm{SM}}^b + \delta g_A^b\big) \gamma_\mu\gamma_5
 * \right] b\,.
 * @f]
 * Equivalently, the shifts in the left-handed and right-handed couplings are
 * defined as
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
 * it is required to call the initialization method InitializeModel(). 
 * This pointer is then used in computing the fermionic neutral-current
 * couplings in the presence of new physics contribution to @f$Zb\bar{b}@f$.
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPZbbbarParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPZbbbar are summarized below:
 *
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGVb</td>
 *   <td class="mod_symb">\f$\delta g_{V}^b\f$</td>
 *   <td class="mod_desc">New physics contribution to \f$g_{V}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGAb</td>
 *   <td class="mod_symb">\f$\delta g_{A}^b\f$</td>
 *   <td class="mod_desc">New physics contribution to \f$g_{A}^b\f$.</td>
 * </tr>
 * </table>
 *
 * which are applicable when the instance of the current class is created with
 * the constructor NPZbbbar() or NPZbbbar(false). 
 * Alternatively, when using NPZbbbar(true), which corresponds to the use of the 
 * model name "NPZbbbarLR" instead of "NPZbbbar" in InputParser::ModelFactory(),
 * the following parameters have to be used as inputs:
 *
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGLb</td>
 *   <td class="mod_symb">\f$\delta g_{L}^b\f$</td>
 *   <td class="mod_desc">New physics contribution to \f$g_{L}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGRb</td>
 *   <td class="mod_symb">\f$\delta g_{R}^b\f$</td>
 *   <td class="mod_desc">New physics contribution to \f$g_{R}^b\f$.</td>
 * </tr>
 * </table>
 *
 * @anchor NPZbbbarFlags
 * <h3>%Model Flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPZbbbarFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()).
 *
 */
class NPZbbbar : public NPbase {
public:

    /**
     * @brief The number of the model parameters in %NPZbbbar.
     */
    static const int NZbbbarVars = 2;

    /**
     * @brief A string array containing the labels of the model parameters in %NPZbbbar,
     * used for NPZbbbarLR=FALSE.
     */
    static const std::string ZbbbarVAVars[NZbbbarVars];

    /**
     * @brief A string array containing the labels of the model parameters in %NPZbbbar,
     * used for NPZbbbarLR=TRUE.
     */
    static const std::string ZbbbarLRVars[NZbbbarVars];

    /**
     * @brief Constructor.
     * @param[in] FlagNPZbbbarLR_in a flag that is true if \f$\delta g_{L}^b\f$
     * and  \f$\delta g_{R}^b\f$ are used for the model parameters, instead of
     * \f$\delta g_{V}^b\f$ and \f$\delta g_{A}^b\f$. For the model parameters, 
     * the labels %deltaGVb and %deltaGAb are applicable for FlagNPZbbbarLR_in=fALSE,
     * while %deltaGLb and %deltaGRb are applicable for FlagNPZbbbarLR_in=TRUE.
     */
    NPZbbbar(const bool FlagNPZbbbarLR_in = false);

    /**
     * @brief The post-update method for %NPZbbbar.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_V^f@f$ (zero in the current model, except for @f$f=b@f$)
     */
    virtual double deltaGV_f(const Particle f) const;

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f a lepton or quark
     * @return @f$\delta g_A^f@f$ (zero in the current model, except for @f$f=b@f$)
     */
    virtual double deltaGA_f(const Particle f) const;

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief StandardModel::Mw()
     * @return @f$M_W@f$ in GeV
     */
    virtual double Mw() const;

    /**
     * @brief The total width of the @f$W@f$ boson, @f$\Gamma_W@f$.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;

    /** 
     * @brief @copybrief StandardModel::Gamma_Z()
     * @return @f$\Gamma_Z@f$ in GeV
     */
    virtual double Gamma_Z() const;

    /** 
     * @brief @copybrief StandardModel::sigma0_had()
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    virtual double sigma0_had() const;

    /** 
     * @brief @copybrief StandardModel::sin2thetaEff()
     * @param[in] f a lepton or quark
     * @return @f$\sin^2\theta_{\rm eff}^{\,\ell}@f$
     */
    virtual double sin2thetaEff(const Particle f) const;

    /** 
     * @brief @copybrief StandardModel::A_f()
     * @param[in] f a lepton or quark
     * @return @f$\mathcal{A}_\ell@f$
     */
    virtual double A_f(const Particle f) const;

    /** 
     * @brief The forward-backward asymmetry in @f$e^+e^-\to Z\to f \bar{f}@f$ at the
     * @f$Z@f$-pole, @f$A^f_{FB}@f$.
     * @param[in] f a lepton or quark
     * @return @f$A^f_{FB}@f$
     */
    virtual double AFB(const Particle f) const;

    /** 
     * @brief The ratio @f$R_\ell^0=\Gamma_{\mathrm{had}}/\Gamma_\ell@f$
     * or @f$R_q^0=\Gamma_q/\Gamma_{\mathrm{had}}@f$, for leptons or quarks, respectively.
     * @param[in] f a lepton or quark
     * @return @f$R_f^0@f$
     */
    virtual double R0_f(const Particle f) const;

    ////////////////////////////////////////////////////////////////////////
protected:

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

    ////////////////////////////////////////////////////////////////////////
private:

    double myDeltaGVb; ///< New physics contribution to \f$g_{V}^b\f$.
    double myDeltaGAb; ///< New physics contribution to \f$g_{A}^b\f$.
    double myDeltaGLb; ///< New physics contribution to \f$g_{L}^b\f$.
    double myDeltaGRb; ///< New physics contribution to \f$g_{R}^b\f$.

    /**
     * @brief An internal boolean flag. This flag is initialized in the
     * constructor. 
     */
    const bool FlagNPZbbbarLR;
};

#endif	/* NPZBBBAR_H */

