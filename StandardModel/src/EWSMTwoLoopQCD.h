/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOLOOPQCD_H
#define	EWSMTWOLOOPQCD_H

#include "EWSMcache.h"
using namespace gslpp;

/**
 * @class EWSMTwoLoopQCD
 * @ingroup StandardModel
 * @brief A class for @f$O(\alpha\alpha_s)@f$ two-loop corrections to the %EW
 * precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class handles two-loop %QCD contributions of
 * @f$O(\alpha\alpha_s)@f$ to the following quantities:
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_l()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha_t()),
 * @li @f$\Delta\rho@f$&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (with DeltaRho()),
 * @li @f$\Delta r_{\mathrm{rem}}@f$&nbsp;&nbsp; (with DeltaR_rem()),
 * @li @f$\delta\rho_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaRho_rem_l() and deltaRho_rem_q()),
 * @li @f$\delta\kappa_{\mathrm{rem}}^{f}@f$&nbsp;&nbsp; (with deltaKappa_rem_l() and deltaKappa_rem_q()).
 *
 * See also the description of EWSM class for details on the above quantities.
 * The @f$O(\alpha\alpha_s)@f$ two-loop %QCD contributions to the vacuum
 * polarization amplitudes of the weak gauge bosons were calculated in
 * @cite Djouadi:1987gn, @cite Djouadi:1987di, @cite Kniehl:1989yc, 
 * @cite Halzen:1990je, @cite Kniehl:1991gu, @cite Kniehl:1992dx and @cite Djouadi:1993ss. 
 */
class EWSMTwoLoopQCD {

public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMTwoLoopQCD(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Leptonic contribution of @f$O(\alpha\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$.
     * @details This contribution vanishes at @f$O(\alpha\alpha_s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{lept}}^{\alpha\alpha_s}=0@f$
     */
    double DeltaAlpha_l(const double s) const;

    /**
     * @brief Top-quark contribution of @f$O(\alpha\alpha_s)@f$
     * to the electromagnetic coupling @f$\alpha@f$.
     * @details A simple numerical formula presented in @cite Kuhn:1998ze is
     * employed. 
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}^{\alpha\alpha_s}@f$
     */
    double DeltaAlpha_t(const double s) const;
    
    /**
     * @brief Leading two-loop %QCD contribution of @f$O(\alpha\alpha_s)@f$
     * to @f$\Delta\rho@f$.
     * @details The formula used here is given by
     * @f[
     * \Delta\rho^{\alpha\alpha_s}
     * = 3\,X_t^\alpha\frac{\alpha_s(m_t^2)}{\pi} \delta^{\mathrm{QCD}}_2,
     * @f]
     * where @f$X_t^\alpha = \alpha\, m_t^2/(16\pi s_W^2 M_W^2)@f$, and
     * @f$\delta^{\mathrm{QCD}}_2@f$ is computed via deltaQCD_2().
     * See Chapter 8 of @cite Bardin:1999ak.
     * This quantity contributes to @f$\Delta r@f$ and the @f$Zf\bar{f}@f$
     * effective couplings @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$.
     * See also the description of EWSM class.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\Delta\rho^{\alpha\alpha_s}@f$
     */
    double DeltaRho(const double Mw_i) const;

    /**
     * @brief remainder contribution to Delta r
     * @param[in] Mw_i the W-boson mass
     * @return Delta r_{rem}^{alpha alpha_s}
     */
    double DeltaR_rem(const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{l, alpha alpha_s}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta rho_{rem}^{q, alpha alpha_s}
     */
    complex deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const;

    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{l, alpha alpha_s}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const;
                                                  
    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of quark
     * @param[in] Mw_i the W-boson mass
     * @return delta kappa_{rem}^{q, alpha alpha_s}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const;  
    

    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @return delta^QCD_2 for the leading contribution to Delta rho
     */
    double deltaQCD_2() const;

    /**
     * @param[in] x x=s/M_t^2
     * @param[in] Mw_i the W-boson mass
     * @return F_1(x)
     * @attention valid for 0<=x<1
     */
    double F1(const double x, const double Mw_i) const;
    
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return V_1(r)
     * @attention valid for 0<=r<1
     */
    double V1(const double r) const;
        
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return FA_1(r)
     * @attention valid for 0<=r<1
     */
    double A1(const double r) const;

    /**
     * @param[in] r r=s/(4M_t^2)
     * @return the derivative of V_1(r)
     * @attention valid for r << 1
     */
    double V1prime(const double r) const;
    
    /**
     * @param[in] r r=s/(4M_t^2)
     * @return the derivative of A_1(r)
     * @attention valid for r << 1
     */
    double A1prime(const double r) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta r^{ud}, contribution to Delta r from the light-quark doublets
     */
    double DeltaR_ud(const double Mw_i) const;

    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta r^{tb}, contribution to Delta r from the top-bottom doublet
     */
    double DeltaR_tb(const double Mw_i) const;
        
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{ud}, contribution to rho_Z^f from the light-quark doublets
     */
    double DeltaRho_ud(const double Mw_i) const;

    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta rho^{tb}, contribution to rho_Z^f from the top-bottom doublet
     */
    double DeltaRho_tb(const double Mw_i) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta kappa^{ud}, contribution to kappa_Z^f from the light-quark doublets
     */
    complex DeltaKappa_ud(const double Mw_i) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return Delta kappa^{tb}, contribution to kappa_Z^f from the top-bottom doublet
     */
    complex DeltaKappa_tb(const double Mw_i) const;
     
    
    ////////////////////////////////////////////////////////////////////////        
    
private:
    const EWSMcache& cache;///< A reference to an object of type EWSMcache.
    
    
};

#endif	/* EWSMTWOLOOPQCD_H */

