/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NEWPHYSICSPARAMS_H
#define	NEWPHYSICSPARAMS_H

#include <string>
#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class NewPhysicsParams
 * @brief A class for retrieving parameters associated with the models in
 * NewPhysics module.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is responsible for retrieving the values of parameters
 * in NPEffective1, NPEffective2, NPEpsilons, NPEpsilons_pureNP, NPSTU,
 * NPSTUVWXY and NPZbbbar classes.
 * The available parameters are listed in the description of the constructor
 * NewPhysicsParams().
 */
class NewPhysicsParams : public ThObservable {
public:

    /**
     * @brief Constructor. 
     * @param[in] ObsType a reference to an object of type ThObsType
     * @param[in] name_i the name of the parameter to be retrieved: 
     * @li "epsilon1", "epsilon2", "epsilon3", "epsilonb":&nbsp;
     * the epsilon parameters @f$\varepsilon_1@f$, @f$\varepsilon_2@f$,
     * @f$\varepsilon_3@f$ and @f$\varepsilon_b@f$,
     * @li "deltaGVb", "deltaGAb":&nbsp;
     * the new physics corrections to the vector and axial-vector @f$Zb\bar{b}@f$
     * couplings, denoted by @f$\delta g_V^b@f$ and @f$\delta g_A^b@f$,
     * @li "deltaGRb", "deltaGLb":&nbsp;
     * the new physics corrections to the right-handed and left-handed @f$Zb\bar{b}@f$
     * couplings, denoted by @f$\delta g_R^b@f$ and @f$g\delta _L^b@f$,
     * @li "deltaRhoZb", "deltaKappaZb":&nbsp;
     * the new physics corrections to the @f$Zb\bar{b}@f$ couplings, 
     * denoted by @f$\delta\rho_Z^b@f$ and @f$\delta\kappa_Z^b@f$,
     * @li "cHQ1pPLUScHQ2p_NP", "cHQ2pMINUScHQ2_NP", "cHQ3pPLUScHQ3_NP":&nbsp;
     * @f$C'_{HQ_1} + C'_{HQ_2}@f$, @f$C'_{HQ_2} - C_{HQ_2}@f$ and 
     * @f$C'_{HQ_3} + C_{HQ_3}@f$
     * @li "c_Ae_NP":&nbsp;
     * @f$C[\mathcal{A}_\ell]@f$ defined in @cite Ciuchini:2013pca,
     * @li "c_GammaZ_uds_NP":&nbsp;
     * @f$C[\Gamma_{uds}]@f$ defined in @cite Ciuchini:2013pca.
     */
    NewPhysicsParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    /**
     * @brief A method to retrieve the value of the desired new physics parameter
     * specified with the constructor.
     * @return the value of the parameter
     */
    double computeThValue();

private:
    const std::string name;///< The name of the parmeter to be retrieved.

};

#endif	/* NEWPHYSICSPARAMS_H */

