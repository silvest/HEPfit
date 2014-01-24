/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
 * @brief A class for approximate formulae of %EW precision observables
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 * @f$M_W@f$, @f$\sin\theta_{\rm eff}^f@f$, and partial decay widths of @f$Z@f$
 * boson.
 * @cite Freitas:2013dpa
 * @cite Freitas:2012sy
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
     * @brief The @f$W@f$-boson mass in the on-shell scheme with the full
     * two-loop %EW corrections.
     * @details This function is based on the approximate formula for @f$M_W@f$
     * presented in @cite Awramik:2003rn, which includes the complete two-loop
     * %EW corrections as well as leading three-loop corrections. The formula
     * used here approximates the full result to be better than 0.5 (0.2) MeV
     * over the range of 10 GeV @f$\leq m_h\leq@f$ 1 TeV (100 GeV @f$\leq m_h
     * \leq@f$ 1 TeV), if all other inputs vary within their @f$2\sigma@f$
     * ranges listed in @cite Awramik:2003rn.
     * @param[in] DeltaAlpha_i the sum of the leptonic and hadronic corrections
     * to @f$\alpha@f$ at @f$q^2=M_Z^2@f$
     * @return the @f$W@f$-boson mass in units of GeV
     */
    double Mw(const double DeltaAlphaL5q_i) const;
    
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^\ell@f$ in the on-shell scheme with the
     * full two-loop EW corrections. 
     * @details This function is based on the approximate formula for the
     * leptonic weak mixing angle presented in @cite Awramik:2006uz (see also 
     * @cite Awramik:2004ge). 
     * @param[in] l name of lepton
     * @param[in] DeltaAlphaL5q_i he sum of the leptonic and hadronic corrections to @f$\alpha@f$ at @f$M_Z@f$..
     * @return the effective weak mixing angle for @f$Z\to\ell\bar{\ell}@f$ obtained from an approximate two-loop formula.
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
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    
};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

