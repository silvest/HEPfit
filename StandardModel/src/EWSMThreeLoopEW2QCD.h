/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTHREELOOPEW2QCD_H
#define	EWSMTHREELOOPEW2QCD_H

#include "EWSMcache.h"
using namespace gslpp;

/**
 * @class EWSMThreeLoopEW2QCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha^2\alpha_s)@f$ three-loop corrections to the
 * %EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class handles three-loop mixed %EW-%QCD contributions of
 * @f$O(\alpha^2\alpha_s)@f$ to the following quantities: 
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$,
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$,
 * @li @f$\Delta\rho@f$,
 * @li @f$\Delta r_{\mathrm{rem}}@f$,
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$,
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$,
 *
 * where only @f$\Delta\rho@f$ is non-zero in the current class.
 * See also the description of EWSM class for details on the above quantities.
 */
class EWSMThreeLoopEW2QCD {

public:
 
    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMThreeLoopEW2QCD(const EWSMcache& cache_i);

    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief leptonic contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{lept}^{alpha^2 alpha_s}
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief top-quark contribution to alpha
     * @param[in] s invariant mass squared 
     * @return Delta alpha_{top}^{alpha^2 alpha_s}
     */
    double DeltaAlpha_t(const double s) const;
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{alpha^2 alpha_s}
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha^2 alpha_s}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha^2 alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha^2 alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha^2 alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha^2 alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;    
    
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;///< A reference to an object of type EWSMcache.
    

    ////////////////////////////////////////////////////////////////////////    
     
    
};

#endif	/* EWSMTHREELOOPEW2QCD_H */

