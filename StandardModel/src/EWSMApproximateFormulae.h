/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMAPPROXIMATEFORMULAE_H
#define	EWSMAPPROXIMATEFORMULAE_H

#include <string>
#include <math.h>
#include "EWSMcache.h"

/**
 * @class EWSMApproximateFormulae
 * @ingroup StandardModel
 * @brief A class for approximate formulae of the %EW precision observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The member functions in the current class compute the %EW precision 
 * observables @f$M_W@f$, @f$\sin\theta_{\rm eff}^f@f$, @f$\Gamma_f@f$,
 * @f$\Gamma_Z@f$, @f$\sigma^0_h@f$, @f$R^0_\ell@f$ @f$R^0_c@f$ and @f$R^0_b@f$,
 * based on the approximate formulae given in
 * @cite Awramik:2003rn, @cite Awramik:2004ge, @cite Awramik:2006uz,
 * @cite Awramik:2008gi, @cite Freitas:2012sy, @cite Freitas:2013dpa and
 * @cite Freitas:2014hra. 
 * (The actual implementation for @f$M_W@f$ corresponds to arXiv:hep-ph/0311148v2,
 * which updates the results presented in the journal version of @cite Awramik:2003rn.)
 * The maximal deviations to the full results and the valid ranges of input
 * parameters are summarized in the description of each function.
 */
class EWSMApproximateFormulae {
public:

    /**
     * @brief Constructor.
     * @param[in] cache_i a reference to an object of type EWSMcache
     */
    EWSMApproximateFormulae(const EWSMcache& cache_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The @f$W@f$-boson mass with the full two-loop %EW corrections.
     * @details This function is based on the approximate formula for @f$M_W@f$
     * presented in @cite Awramik:2003rn, which includes the complete two-loop
     * %EW corrections as well as leading three-loop corrections, and the 
     * four-loop corrections to the rho parameter.
     * (The four-loop effects are not included in the results presented in the 
     * journal version of @cite Awramik:2003rn. The parametrization used here 
     * corresponds to the results in arXiv:hep-ph/0311148v2, which updates the
     * the ones presented in the published version.)
     * The approximate formula reproduces the full result to be better than
     * 0.5 (0.25) MeV over the range of 10 GeV @f$\leq m_h\leq@f$ 1 TeV
     * (100 GeV @f$\leq m_h \leq@f$ 1 TeV), if other inputs vary within
     * their @f$2\sigma@f$ ranges of the 2003 data, where their @f$1\sigma@f$
     * ranges are given by
     * @f$\alpha_s = 0.1190\pm 0.0027@f$,
     * @f$\Delta\alpha^{\ell+5q} = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1875\pm 0.0021@f$ GeV, and
     * @f$m_t = 174.3\pm 5.1@f$ GeV.
     * @return the @f$W@f$-boson mass in units of GeV
     */
    double Mw() const;
    
    /**
     * @brief The value of @f$\Delta\alpha^{5}_{had}(M_Z^2)@f$ obtained from the @f$W@f$-boson mass, 
     * using the full two-loop %EW corrections.
     * @details This function is based on the approximate formula for @f$M_W@f$
     * presented in @cite Awramik:2003rn. See notes in @f$Mw()@f$ function.
     * @return the value of @f$\Delta\alpha^{5}_{had}(M_Z^2)@f$
     */
    double dAlpha5hMw() const;

    
    /**
     * @brief The value of the effective weak mixing anlge for a given fermion.
     * @param[in] p name of a particle
     * @return the value of @f$\sin^^{\theta_{Eff}^p}@f$
     */
    double sin2thetaEff(const Particle p) const
    {
        if (p.is("QUARK"))
            return sin2thetaEff_q((QCD::quark) (p.getIndex() - 6));
        else if (p.is("LEPTON"))
            return sin2thetaEff_l((QCD::lepton) p.getIndex());
        else
            throw std::runtime_error("EWSMApproximateFormulae::sin2thetaEff() called with wrong argument");
    }

    /**
     * @brief @f$\Delta r_{\rm rem}^{(\alpha^2)}@f$. 
     * @details This function is based on the approximate formula for the irreducible 
     * %EW two-loop contribution to @f$\Delta r@f$ presented in @cite Awramik:2006uz,
     * which includes the complete two-loop %EW corrections as well as leading
     * three-loop corrections.
     * The approximate formula reproduces the full result to be better than
     * @f$2.7\times 10^{-5}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.119\pm 0.002@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 172.5\pm 2.3@f$ GeV.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return irreducible two-loop %EW contribution to @f$\Delta r@f$
     */
    double DeltaR_TwoLoopEW_rem(const double Mw_i) const;

    /**
     * @brief @f$\Delta\kappa_Z^{\ell, (\alpha^2)}@f$. 
     * @details This function is based on the approximate formula for the irreducible
     * %EW two-loop contribution to @f$\Delta\kappa_Z^\ell = \kappa_Z^\ell - 1@f$ 
     * presented in @cite Awramik:2006uz, which includes the complete two-loop
     * %EW corrections as well as leading three-loop corrections.
     * The approximate formula reproduces the full result to be better than
     * @f$1.8\times 10^{-5}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.119\pm 0.002@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 172.5\pm 2.3@f$ GeV.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return irreducible two-loop %EW contribution to @f$\Delta\kappa_Z^\ell@f$
     */
    double DeltaKappa_l_TwoLoopEW_rem(const double Mw_i) const;

    /**
     * @brief @f$\Delta\kappa_Z^{b, (\alpha^2)}@f$. 
     * @details This function is based on the approximate formula for the irreducible
     * %EW two-loop contribution to @f$\Delta\kappa_Z^b = \kappa_Z^b - 1@f$
     * presented in @cite Awramik:2008gi, which includes the complete fermionic
     * two-loop %EW corrections as well as leading three-loop corrections.
     * The bosonic two-loop %EW corrections are not included.
     * The approximate formula reproduces the full result to be better than
     * @f$1.4\times 10^{-5}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.119\pm 0.002@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 172.5\pm 2.3@f$ GeV.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return irreducible fermionic two-loop %EW contribution to @f$\Delta\kappa_Z^b@f$
     */
    double DeltaKappa_b_TwoLoopEW_rem(const double Mw_i) const;

    /**
     * @brief @f$R_b^0@f$. 
     * @details This function is based on the approximate formula for 
     * @f$R_b^0=\Gamma_b/\Gamma_h@f$ presented in @cite Freitas:2012sy, which
     * includes the complete fermionic two-loop %EW corrections as well as
     * leading three-loop corrections.
     * The bosonic two-loop %EW corrections are not included.
     * The approximate formula reproduces the full result to be better than
     * @f$10^{-6}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0007@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05900\pm 0.00033@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 173.2\pm 0.9@f$ GeV.
     * @return @f$R_b^0=\Gamma_b/\Gamma_h@f$
     */
    double R0_bottom_OLD() const;

    /**
     * @brief @f$\Gamma_u/\Gamma_b@f$. 
     * @details This function is based on the approximate formula for
     * the ratio @f$\Gamma_u/\Gamma_b@f$ obtained from A. Freitas in private
     * communication on Sep. 21, 2013, which includes the complete fermionic
     * two-loop %EW corrections as well as leading three-loop corrections.
     * The bosonic two-loop %EW corrections are not included.
     * The approximate formula reproduces the full result to be better than
     * @f$3.3\times 10^{-6}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0007@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05900\pm 0.00033@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 173.2\pm 0.9@f$ GeV.
     * @return @f$\Gamma_u/\Gamma_b@f$
     */
    double Gu_over_Gb_OLD() const;

    /**
     * @brief @f$\Gamma_d/\Gamma_b@f$.
     * @details This function is based on the approximate formula for
     * the ratio @f$\Gamma_d/\Gamma_b@f$ obtained from A. Freitas in private
     * communication on Sep. 21, 2013, which includes the complete fermionic
     * two-loop %EW corrections as well as leading three-loop corrections.
     * The bosonic two-loop %EW corrections are not included.
     * The approximate formula reproduces the full result to be better than
     * @f$3.0\times 10^{-6}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0007@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05900\pm 0.00033@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 173.2\pm 0.9@f$ GeV.
     * @return @f$\Gamma_d/\Gamma_b@f$
     */
    double Gd_over_Gb_OLD() const;

    /**
     * @brief @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$R^0_\ell@f$, @f$R^0_c@f$, @f$R^0_b@f$, or @f$\sigma^0_h@f$.
     * @details This function is based on the approximate formulae for partial
     * and total widths of the @f$Z@f$ boson and hadronic @f$Z@f$-pole cross
     * section presented in @cite Freitas:2014hra, which include the complete
     * fermionic two-loop %EW corrections as well as leading three-loop
     * corrections. The bosonic two-loop %EW corrections are not included.
     * The approximate formulae reproduce the full results to be better than
     * 0.001 MeV, 0.01 MeV, 0.1 pb, @f$0.1\times 10^{-3}@f$ and
     * @f$0.01\times 10^{-3}@f$ for @f$\Gamma_f@f$, @f$\Gamma_Z@f$,
     * @f$\sigma^0_h@f$, @f$R^0_\ell@f$ and @f$R^0_{c,b}@f$, respectively,
     * if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0042@f$ GeV,
     * @f$m_t = 173.2\pm 2.0@f$ GeV and
     * @f$m_h = 125.7\pm 2.5@f$ GeV.
     * @param[in] observable name of the observable to be computed:
     * "Gamma_nu", "Gamma_e_mu", "Gamma_tau", "Gamma_u", "Gamma_c", "Gamma_d_s",
     * "Gamma_b", "GammaZ", "sigmaHadron", "R0_lepton", "R0_charm", "R0_bottom"
     * @return @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$\sigma^0_h@f$, @f$R^0_\ell@f$,  @f$R^0_c@f$, or @f$R^0_b@f$
     *
     * @attention The function EWSMApproximateFormulae::X_extended is applicable to
     * larger ranges of @f$M_Z@f$, @f$m_t@f$ and @f$m_h@f$.
     */
    double X(const std::string observable) const;

    /**
     * @brief @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$R^0_\ell@f$, @f$R^0_c@f$, @f$R^0_b@f$, or @f$\sigma^0_h@f$.
     * @details This function is based on the approximate formulae for partial
     * and total widths of the @f$Z@f$ boson and hadronic @f$Z@f$-pole cross
     * section presented in @cite Freitas:2014hra, which include the complete
     * fermionic two-loop %EW corrections as well as leading three-loop
     * corrections. The bosonic two-loop %EW corrections are not included.
     * The approximate formulae reproduce the full results to be better than
     * 0.001 MeV, 0.01 MeV, 0.1 pb, @f$0.1\times 10^{-3}@f$ and
     * @f$0.01\times 10^{-3}@f$ for @f$\Gamma_f@f$, @f$\Gamma_Z@f$,
     * @f$\sigma^0_h@f$, @f$R^0_\ell@f$ and @f$R^0_{c,b}@f$, respectively,
     * if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0084@f$ GeV,
     * @f$165 < m_t < 190@f$ GeV and
     * @f$70 < m_h < 1000@f$ GeV.
     * @param[in] observable name of the observable to be computed:
     * "Gamma_nu", "Gamma_e_mu", "Gamma_tau", "Gamma_u", "Gamma_c", "Gamma_d_s",
     * "Gamma_b", "GammaZ", "sigmaHadron", "R0_lepton", "R0_charm", "R0_bottom"
     * @return @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$\sigma^0_h@f$, @f$R^0_\ell@f$,  @f$R^0_c@f$, or @f$R^0_b@f$
     *
     * @attention The function EWSMApproximateFormulae::X is applicable to
     * smaller ranges of @f$M_Z@f$, @f$m_t@f$ and @f$m_h@f$.
     */
    double X_extended(const std::string observable) const;
    
    
    /**
     * @brief @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$R^0_\ell@f$, @f$R^0_c@f$, @f$R^0_b@f$, or @f$\sigma^0_h@f$.
     * @details This function is based on the approximate formulae for partial
     * and total widths of the @f$Z@f$ boson and hadronic @f$Z@f$-pole cross
     * section presented in arXiv: 1804.10236, which include the complete
     * two-loop %EW corrections as well as leading three-loop
     * corrections.
     * The approximate formulae reproduce the full results to be better than
     * 0.001 MeV, 0.002 MeV, 0.006 MeV, 0.012 MeV, @f$ 0.1\times 10^{-3}@f$,
     * @f$ 0.01\times 10^{-3}@f$ and 0.1 pb,
     * for @f$\Gamma_{e,\mu,\tau,\nu}@f$, @f$\Gamma_{q\not = b}@f$,
     * @f$\Gamma_{b}@f$, @f$\Gamma_Z@f$, @f$R^0_{l}@f$, @f$R^0_{c,b}@f$ and
     * @f$\sigma^0_h@f$, respectively,
     * if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0042@f$ GeV,
     * @f$169.2 < m_t < 177.2@f$ GeV and
     * @f$120.1 < m_h < 130.1@f$ GeV.
     * For @f$m_h@f$ beyond [85,165] GeV there are significant differences with some predicions
     * of X_extended, which go well beyond the expected size of the bosonic corrections (>~2x).
     * The function redirects to X_extended in that case.
     * @param[in] observable name of the observable to be computed:
     * "Gamma_nu", "Gamma_e_mu", "Gamma_tau", "Gamma_u", "Gamma_c", "Gamma_d_s",
     * "Gamma_b", "GammaZ", "sigmaHadron", "R0_lepton", "R0_charm", "R0_bottom"
     * @return @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$\sigma^0_h@f$, @f$R^0_\ell@f$,  @f$R^0_c@f$, or @f$R^0_b@f$
     *
     */
    double X_full_2_loop(const std::string observable) const;
    
    
    /**
     * @brief @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$, @f$\Gamma_\nu@f$, 
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$R^0_\ell@f$, @f$R^0_c@f$, @f$R^0_b@f$, or @f$\sigma^0_h@f$.
     * @details This function is based on the approximate formulae for partial
     * and total widths of the @f$Z@f$ boson and hadronic @f$Z@f$-pole cross
     * section presented in arXiv: 1906.08815, which include the complete
     * two-loop %EW corrections as well as leading three-loop
     * corrections.
     * The approximate formulae reproduce the full results to be better than
     * 0.0015 MeV, 0.0015 MeV, 0.002 MeV, 0.006 MeV, 0.006 MeV, 0.007 MeV, 0.007 MeV,
     * 0.04 MeV, @f$ 0.12\times 10^{-3}@f$, @f$ 0.1\times 10^{-3}@f$, @f$ 0.12\times 10^{-3}@f$,
     * and 0.15 pb,
     * for @f$\Gamma_{e,\mu,\tau,\nu}@f$, @f$\Gamma_{q\not = b}@f$,
     * @f$\Gamma_{b}@f$, @f$\Gamma_Z@f$, @f$R^0_{l}@f$, @f$R^0_{c,b}@f$ and
     * @f$\sigma^0_h@f$, respectively,
     * if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0084@f$ GeV,
     * @f$155 < m_t < 192@f$ GeV and
     * @f$25 < m_h < 225@f$ GeV.
     * @param[in] observable name of the observable to be computed:
     * "Gamma_nu", "Gamma_e_mu", "Gamma_tau", "Gamma_u", "Gamma_c", "Gamma_d_s",
     * "Gamma_b", "GammaZ", "sigmaHadron", "R0_lepton", "R0_charm", "R0_bottom"
     * @return @f$\Gamma_\nu@f$, @f$\Gamma_{e,\mu}@f$, @f$\Gamma_\tau@f$,
     * @f$\Gamma_u@f$, @f$\Gamma_c@f$, @f$\Gamma_{d,s}@f$, @f$\Gamma_b@f$,
     * @f$\Gamma_Z@f$, @f$\sigma^0_h@f$, @f$R^0_\ell@f$,  @f$R^0_c@f$, or @f$R^0_b@f$
     *
     */
    double X_full(const std::string observable) const;
    
    
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^b@f$ with the full two-loop %EW
     * corrections.
     * @details This function is based on the approximate formulae for 
     * presented in arXiv: 1906.08815, which include the complete
     * two-loop %EW corrections as well as leading three-loop
     * corrections.
     * The approximate formulae reproduce the full results to be better than
     * 0.0025 * 10^-4, if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0084@f$ GeV,
     * @f$155 < m_t < 192@f$ GeV and
     * @f$25 < m_h < 225@f$ GeV.
     * @return @f$\sin^2\theta_{\rm eff}^b@f$
     *
     */
    double sin2thetaEff_b_full() const;
    
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^l@f$ with the full two-loop %EW
     * corrections.
     * @details This function is based on the approximate formulae for 
     * presented in arXiv: 1906.08815, which include the complete
     * two-loop %EW corrections as well as leading three-loop
     * corrections.
     * The approximate formulae reproduce the full results to be better than
     * 0.0056 * 10^-4, if inputs vary within the ranges
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.0050@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.0590\pm 0.0005@f$,
     * @f$M_Z = 91.1876\pm 0.0084@f$ GeV,
     * @f$155 < m_t < 192@f$ GeV and
     * @f$25 < m_h < 225@f$ GeV.
     * @return @f$\sin^2\theta_{\rm eff}^l@f$
     *
     */
    double sin2thetaEff_l_full() const;
    
    
    //LEP2 Observables    
    
    /**
     * @brief The @f$e^+e^- \to \mu^+\mu^-@f$ cross section at LEP2
     * @details This function is based on the approximate formula for the
     * cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \mu^+\mu^-@f$ cross section in units of pb
     */
    double LEP2sigmaMuApprox(const double s) const;
    
    /**
     * @brief The @f$e^+e^- \to \mu^+\mu^-@f$ forward-backward asymmetry at LEP2
     * @details This function is based on the approximate formula for the
     * forward-backward asymmetry, obtained fitting the ZFitter predictions to 
     * a semi-analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \mu^+\mu^-@f$ forward-backward asymmetry in units of pb
     */
    double LEP2AFBmuApprox(const double s) const;
    
    /**
     * @brief The @f$e^+e^- \to \tau^+\tau^-@f$ cross section at LEP2
     * @details This function is based on the approximate formula for the
     * cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \tau^+\tau^-@f$ cross section in units of pb
     */
    double LEP2sigmaTauApprox(const double s) const;
    
    /**
     * @brief The @f$e^+e^- \to \tau^+\tau^-@f$ forward-backward asymmetry at LEP2
     * @details This function is based on the approximate formula for the
     * forward-backward asymmetry, obtained fitting the ZFitter predictions to 
     * a semi-analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \tau^+\tau^-@f$ forward-backward asymmetry in units of pb
     */
    double LEP2AFBtauApprox(const double s) const;
    
    /**
     * @brief The @f$e^+e^- \to hadrons@f$ cross section at LEP2
     * @details This function is based on the approximate formula for the
     * cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to hadrons@f$ cross section in units of pb
     */
    double LEP2sigmaHadronApprox(const double s) const;
    
       
    //LEP2 Differential Observables
    
    /**
     * @brief The @f$e^+e^- \to e^+ e^-@f$ differential cross section at LEP2
     * @details This function is based on the approximate formula for the
     * differential cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to e^+ e^-@f$ cross section in units of pb
     */
    double LEP2dsigmadcosEApprox(const double s, const double cos) const;
    
    /**
     * @brief The @f$e^+e^- \to \mu^+\mu^-@f$ differential cross section at LEP2
     * @details This function is based on the approximate formula for the
     * differential cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \mu^+\mu^-@f$ cross section in units of pb
     */
    double LEP2dsigmadcosMuApprox(const double s, const double cos) const;
    
    /**
     * @brief The @f$e^+e^- \to \tau^+\tau^-@f$ differential cross section at LEP2
     * @details This function is based on the approximate formula for the
     * differential cross section, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return the @f$e^+e^- \to \tau^+\tau^-@f$ cross section in units of pb
     */
    double LEP2dsigmadcosTauApprox(const double s, const double cos) const;

    // EW low-energy observables: neutrino-scattering
    
    /**
     * @brief The effective neutrino nucleon LH coupling: gLnuN2
     * @details This function is based on the approximate formula for the
     * coupling, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$g_L^2(\nu N)@f$
     */
    double LEgLnuN2Approx() const;
    
    /**
     * @brief The effective neutrino nucleon RH coupling: gRnuN2
     * @details This function is based on the approximate formula for the
     * coupling, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$g_R^2(\nu N)@f$
     */    
    double LEgRnuN2Approx() const;
    
    /**
     * @brief The effective neutrino nucleon LH parameter: ThetaLnuN
     * @details This function is based on the approximate formula for the
     * parameter, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$\theta_L(\nu N)@f$
     */     
    double LEThetaLnuNApprox() const;

    /**
     * @brief The effective neutrino nucleon RH parameter: ThetaRnuN
     * @details This function is based on the approximate formula for the
     * parameter, obtained fitting the ZFitter predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$\theta_R(\nu N)@f$
     */     
    double LEThetaRnuNApprox() const;
    
    /**
     * @brief The effective (muon) neutrino-electron vector coupling: gVnue
     * @details This function is based on the approximate formula for the
     * coupling, obtained fitting the SM predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$g_V^{\nu_\mu e}@f$
     */
    double LEgVnueApprox() const;
    
    /**
     * @brief The effective (muon) neutrino-electron axial-vector coupling: gAnue
     * @details This function is based on the approximate formula for the
     * coupling, obtained fitting the SM predictions to a semi-
     * analytical expression as function of the SM parameters.
     * @return @f$g_A^{\nu_\mu e}@f$
     */
    double LEgAnueApprox() const;
    
    ////////////////////////////////////////////////////////////////////////

private:
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^\ell@f$ with the full two-loop %EW corrections.
     * @details This function is based on the approximate formulae for the
     * leptonic weak mixing angles presented in @cite Awramik:2006uz
     * (see also @cite Awramik:2004ge), which include the complete two-loop
     * %EW corrections as well as leading three-loop corrections. 
     * The approximate formulae reproduce the full results to be better than
     * @f$4.5\times 10^{-6}@f$ for the Higgs mass 10 GeV @f$\leq m_h\leq@f$ 1 TeV,
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.119\pm 0.002@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 172.5\pm 2.3@f$ GeV.
     * @param[in] l name of a lepton (see QCD::lepton)
     * @return the effective weak mixing angle for @f$Z\to\ell\bar{\ell}@f$
     */
    double sin2thetaEff_l(const QCD::lepton l) const;

    /**
     * @brief @f$\sin^2\theta_{\rm eff}^q@f$ with the full two-loop %EW
     * corrections (bosonic two-loop %EW corrections are missing for @f$q=b@f$).
     * @details This function is based on the approximate formulae for the
     * weak mixing angles presented in @cite Awramik:2006uz and @cite Awramik:2008gi,
     * which include the complete two-loop %EW corrections as well as leading
     * three-loop corrections. It is noted that bosonic two-loop %EW corrections
     * are missing for @f$q=b@f$.
     * The approximate formulae reproduce the full results to be better than
     * @f$4.5\times 10^{-6}@f$ (@f$4.3\times 10^{-6}@f$) for the Higgs mass
     * 10 GeV @f$\leq m_h\leq@f$ 1 TeV in the case of @f$q=u,d,s,c@f$ (@f$q=b@f$),
     * if other inputs vary within their @f$2\sigma@f$ ranges of the
     * following outdated data:
     * @f$\alpha_s(M_Z^2) = 0.119\pm 0.002@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.05907\pm 0.00036@f$,
     * @f$M_Z = 91.1876\pm 0.0021@f$ GeV and
     * @f$m_t = 172.5\pm 2.3@f$ GeV.
     * @param[in] q name of a quark (see QCD::quark)
     * @return the effective weak mixing angle for @f$Z\to q\bar{q}@f$
     */
    double sin2thetaEff_q(const QCD::quark q) const;
    
    /**
     * @brief @f$\sin^2\theta_{\rm eff}^b@f$ with the full two-loop %EW
     * corrections.
     * @details This function is based on the approximate formulae for the
     * weak mixing angle presented in [arXiv:1607.08375 hep-ph] (add citation),
     * which include the complete two-loop %EW corrections as well as leading
     * three-loop corrections.
     * The approximate formulae reproduce the full results with average and 
     * maximal deviations of @f$2\times 10^{-7}@f$ and
     * @f$1.3\times 10^{-6}@f$, respectively, for the input parameters in the
     * following ranges:
     * @f$m_h = 125.1 \pm 5@f$ GeV, 
     * @f$\alpha_s(M_Z^2) = 0.1184\pm 0.005@f$,
     * @f$\Delta\alpha^{\ell+5q}(M_Z^2) = 0.059\pm 0.0005@f$, 
     * @f$M_Z = 91.1876\pm 0.0042@f$ GeV and
     * @f$m_t = 173.2\pm 4.0@f$ GeV.
     * @return the effective weak mixing angle for @f$Z\to b\bar{b}@f$
     */
    double sin2thetaEff_b() const;



    const EWSMcache& mycache; ///< A reference to an object of type StandardModel.

};

#endif	/* EWSMAPPROXIMATEFORMULAE_H */

