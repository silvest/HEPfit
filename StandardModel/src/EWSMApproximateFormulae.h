/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMAPPROXIMATEFORMULAE_H
#define	EWSMAPPROXIMATEFORMULAE_H

#include "StandardModel.h"

/**
 * @class EWSMApproximateFormulae
 * @ingroup StandardModel
 * @brief A class for approximate formulae of @f$M_W@f$, @f$\sin\theta_{\rm eff}^f@f$ and @f$R_b^0@f$.  
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWSMApproximateFormulae {
public:
      
    /**
     * @brief An EWSMApproximateFormulae constructor.
     * @param[in] SM_i A reference to a StandardModel object.
     */
    EWSMApproximateFormulae(const StandardModel& SM_i);    

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The W-boson mass with the full two-loop EW corrections. 
     * @param[in] DeltaAlpha_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return The W-boson mass obtained from an approximate two-loop formula.
     */
    double Mw(const double DeltaAlphaL5q_i) const;
    
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^\ell@f$ with the full two-loop EW corrections. 
     * @param[in] l Name of lepton.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return The effective weak mixing angle for @f$Z\to\ell\bar{\ell}@f$ obtained from an approximate two-loop formula.
     */
    double sin2thetaEff_l(const StandardModel::lepton l, const double DeltaAlphaL5q_i) const;

    /**
     * @brief @f$\sin^2\theta_{\rm eff}^q@f$ with the full two-loop EW corrections.
     * @param[in] q Name of quark.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return The effective weak mixing angle for @f$Z\to q\bar{q}@f$ obtained from an approximate two-loop formula.
     * @attention EW two-loop bosonic contribution is missing for @f$q=b@f$. 
     */
    double sin2thetaEff_q(const StandardModel::quark q, const double DeltaAlphaL5q_i) const;    
    
    /**
     * @brief @f$\Delta r_{\rm rem}^{(\alpha^2)}@f$
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$.
     * @param[in] Mw_i The W-boson mass.
     * @return Irreducible two-loop EW contribution to @f$\Delta r@f$ obtained from an approximate formula.
     */
    double DeltaR_TwoLoopEW_rem(const double DeltaAlphaL5q_i, const double Mw_i) const;

    /**
     * @brief EW two-loop contribution to @f$\Delta\kappa_Z^\ell = \kappa_Z^\ell - 1@f$.
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$.
     * @param[in] Mw_i The W-boson mass.
     * @return EW two-loop contribution to @f$\Delta\kappa^{\ell, \alpha^2}@f$ obtained from an approximate formula.
     */
    double DeltaKappa_l_TwoLoopEW_rem(const double DeltaAlphaL5q_i, const double Mw_i) const;

    /**
     * @brief EW two-loop fermionic contribution to @f$\Delta\kappa_Z^b = \kappa_Z^b - 1@f$.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$.
     * @param[in] Mw_i The W-boson mass.
     * @return EW two-loop fermionic contribution to @f$\Delta\kappa^{b, \alpha^2}@f$ obtained from an approximate formula.
     */
    double DeltaKappa_b_TwoLoopEW_rem(const double DeltaAlphaL5q_i, const double Mw_i) const;
    
    /**
     * @brief @f$R_b^0@f$ with the complete fermionic EW two-loop corrections.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return @f$R_b^0@f$ obtained from an approximate two-loop formula.
     */
    double R0_bottom(const double DeltaAlphaL5q_i) const;

    /**
     * @brief @f$\Gamma_u/\Gamma_b@f$.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return @f$\Gamma_u/\Gamma_b@f$.
     */
    double Gu_over_Gb(const double DeltaAlphaL5q_i) const;
    
    /**
     * @brief @f$\Gamma_u/\Gamma_b@f$.
     * @param[in] DeltaAlphaL5q_i The sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return @f$\Gamma_d/\Gamma_b@f$.
     */
    double Gd_over_Gb(const double DeltaAlphaL5q_i) const;

    double GammaZ(const double DeltaAlphaL5q_i) const;

    double sigmaHadron(const double DeltaAlphaL5q_i) const;

    
    ////////////////////////////////////////////////////////////////////////
    
private:
    const StandardModel& SM;
    
};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

