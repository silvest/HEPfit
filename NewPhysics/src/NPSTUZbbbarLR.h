/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTUZBBBARLR_H
#define	NPSTUZBBBARLR_H

#include "NPSTU.h"

/**
 * @class NPSTUZbbbarLR
 * @brief A model class for new physics in the form of oblique and 
 * @f$Zb\bar{b}@f$-vertex corrections.
 * @ingroup NewPhysics
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class containing the necessary functions to compute
 * new physics contributions to the electroweak precision observables with the
 * Peskin-Takeuchi oblique parameters \cite Peskin:1990zt, \cite Peskin:1991sw
 * and with additional shifts to the left-handed and right-handed
 * @f$Zb\bar{b}@f$ vertices. The effective @f$Zb\bar{b}@f$ couplings involve
 * new physics contributions both from the oblique corrections and from the 
 * vertex corrections:
 * @f[
 * \mathcal{L}_{\mathrm{eff}} =
 * \frac{e}{2s_W c_W}\,
 * Z_\mu\, \bar{b}
 * \left[ g_{R,\mathrm{eff}}^b\, \gamma_\mu (1 + \gamma_5)
 * + g_{L,\mathrm{eff}}^b\, \gamma_\mu (1 - \gamma_5)
 * \right] b
 * @f]
 * where the effective couplings are given by 
 * @f[
 * g_{R,\mathrm{eff}}^b = g_{R,\mathrm{SM}}^b
 * + g_{R,\mathrm{SM}}^b \frac{\alpha(M_Z^2)\, T}{2}
 * + g_{R,\mathrm{SM}}^b
 *   \frac{\alpha(M_Z^2)\left( S - 4\,c_W^2s_W^2\, T \right)}
 *   {4s_W^2\,(c_W^2-s_W^2)}
 * + \delta g_{R,\mathrm{vert}}^b\,,
 * \\
 * g_{L,\mathrm{eff}}^b = g_{L,\mathrm{SM}}^b
 * + g_{L,\mathrm{SM}}^b \frac{\alpha(M_Z^2)\, T}{2}
 * + g_{R,\mathrm{SM}}^b
 *   \frac{\alpha(M_Z^2)\left( S - 4\,c_W^2s_W^2\, T \right)}
 *   {4s_W^2\,(c_W^2-s_W^2)}
 * + \delta g_{L,,\mathrm{vert}}^b\,.
 * @f]
 *
 *
 * @anchor NPSTUZbbbarLRInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor
 * NPSTUZbbbarLR(), it is required to call the initialization method
 * InitializeModel().
 * In the Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 *
 * @anchor NPSTUZbbbarLRParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPSTUZbbbarLR are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGLb</td>
 *   <td class="mod_symb">\f$\delta g_{L,\mathrm{vert}}^b\f$</td>
 *   <td class="mod_desc">New physics vertex corrections to \f$g_{L}^b\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%deltaGRb</td>
 *   <td class="mod_symb">\f$\delta g_{R,\mathrm{vert}}^b\f$</td>
 *   <td class="mod_desc">New physics vertex corrections to \f$g_{R}^b\f$.</td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPSTUZbbbarLRFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPSTUZbbbarLRFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPSTU, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVf()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAf()).
 * 
 */
class NPSTUZbbbarLR : public NPSTU {
public:

    /**
     * @brief The number of the model parameters in %NPSTUZbbbarLR.
     */
    static const int NSTUZbbbarLRvars = 2;

    /**
     * @brief A string array containing the labels of the model parameters in %NPSTUZbbbarLR.
     */
    static const std::string STUZbbbarLRvars[NSTUZbbbarLRvars];

    /**
     * @brief Constructor.
     */
    NPSTUZbbbarLR();

    /**
     * @brief @copybrief Model::CheckParameters()
     * @copydetails Model::CheckParameters()
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief New physics contribution to the neutral-current vector coupling @f$g_V^f@f$.
     * @param[in] f name of a fermion
     * @return @f$\delta g_V^f@f$ (non-zero only for @f$f=b@f$)
     */
    virtual double deltaGV_f(const Particle f) const;

    /**
     * @brief New physics contribution to the neutral-current axial-vector coupling @f$g_A^f@f$.
     * @param[in] f name of a fermion
     * @return @f$\delta g_A^f@f$ (non-zero only for @f$f=b@f$)
     */
    virtual double deltaGA_f(const Particle f) const;

    ////////////////////////////////////////////////////////////////////////
protected:

    double myDeltaGLb; ///< New physics contribution to \f$g_{L}^b\f$.
    double myDeltaGRb; ///< New physics contribution to \f$g_{R}^b\f$.

    /**
     * @brief @copybrief Model::setParameter()
     * @copydetails Model::setParameter()
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPSTUZBBBARLR_H */

