/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_H
#define	EW_H

#include <stdexcept>
#include <ThObsType.h>
#include <StandardModel.h>
#include "EW_NPZff.h"

using namespace gslpp;

/**
 * @class EW
 * @ingroup EW 
 * @brief The base class for the electroweak precision observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains basic functions for the computation of the
 * electroweak precision observables at the @f$Z@f$ pole, such as the left-right
 * and forward-backward asymmetries in @f$e^+e^-\to Z\to f\bar{f}@f$, the 
 * partial and total decay widths of the @f$Z@f$ boson, and the pole hadronic
 * cross section.
 *
 * The quantities calculated in the current class may contain new physics (NP)
 * contribution in addition to the standard model (SM) contribution, depending
 * on the model under consideration. In some models, the NP contribution to
 * an observable @f$\mathcal{O}@f$ is linearized and added to the SM contribution:
 * @f[
 * \mathcal{O}
 * = \mathcal{O}(g_{V,\mathrm{SM}}^f,\, g_{A,\mathrm{SM}}^f)
 *   + \mathcal{O}(\delta g_{V,\mathrm{NP}}^f,\, \delta g_{A,\mathrm{NP}}^f)\,,
 * @f]
 * where the observable depends on the effective @f$Zf\bar{f}@f$ couplings
 * @f$g_{V}^f@f$ and @f$g_{A}^f@f$, and the second term is a linear function
 * in terms of @f$\delta g_{V,\mathrm{NP}}^f@f$ and @f$\delta g_{A,\mathrm{NP}}^f@f$.
 * In this case, the NP contributions are not taken into account in the current
 * class, and the outputs for the observables correspond purely to their SM
 * predictions. The NP conributions are then added in each ThObservable class,
 * such as Alepton, GammaZ, Rlepton, etc., via the functions in EW_NPZff class.
 *
 * In the other models, the NP contribution is not linearized:
 * @f[
 * \mathcal{O}_i = \mathcal{O}(g_V^f,\, g_A^f)
 * = \mathcal{O}(g_{V,\mathrm{SM}}^f + \delta g_{V,\mathrm{NP}}^f,\,
 *               g_{A,\mathrm{SM}}^f + \delta g_{A,\mathrm{NP}}^f)\,,
 * @f]
 * and the outputs of the current class carry the NP contributions.
 *
 * The choice of the above two cases is determined by the function 
 * checkNPZff_linearized(). In the most models, implemented into the codes,
 * the NP contributions are linearized. Only in the case
 * of NPEpsilons model as well as NPZbbbar model with the model flag
 * @ref NPZbbbarFlags "NotLinearizedNP=TRUE", the NP contributions are not
 * linearized.
 *
 * See @cite Ciuchini:2013pca and references therein. 
 */
class EW : public ThObsType {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EW(const StandardModel& SM_i);

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief A method to check whether NP contributions to the @f$Z@f$-pole 
     * observables have to be linearized in the model under consideration.
     * @details If this method returns true, NP contributions to the @f$Z@f$-pole
     * observables have to be linearized, and the functions defined in EW_NPZff
     * class have to be employed.
     * @return a boolean that is true if the NP contributions have to be linearized
     *
     * @attention This function returns true in the case where the model is StandardModel.
     */
    bool checkNPZff_linearized() const;
    
    /**
     * @brief A get method to access the member reference of type StandardModel.
     * @return the member reference #SM
     */
    const StandardModel& getSM() const 
    {
        return SM;
    } 

    /**
     * @brief A get method to access the member object of type EW_NPZff.
     * @return the member object EW_NPZff
     */
    const EW_NPZff getMyEW_NPZff() const
    {
        return myEW_NPZff;
    }

    
    ////////////////////////////////////////////////////////////////////////
    // Final-state corrections to Z-decay widths

    /**
     * @brief The non-factorizable EW-%QCD corrections to the partial widths
     * for @f$Z\to q\bar{q}@f$, denoted as @f$\Delta_{\mathrm{EW/QCD}}@f$.
     * @details
     * See @cite Czarnecki:1996ei and @cite Harlander:1997zb.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Delta_{\mathrm{EW/QCD}}@f$ in GeV
     */
    double Delta_EWQCD(const QCD::quark q) const;

    /**
     * @brief The radiator factor associated with the final-state QED and %QCD
     * corrections to the the vector-current interactions, @f$R_V^q(M_Z^2)@f$.
     * @details
     * See @cite Chetyrkin:1994js, @cite Bardin:1999ak, @cite Bardin:1999yd
     * and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$R_V^q(M_Z^2)@f$
     */
    double RVq(const QCD::quark q) const;

    /**
     * @brief The radiator factor associated with the final-state QED and %QCD
     * corrections to the the axial-vector-current interactions, @f$R_A^q(M_Z^2)@f$.
     * @details
     * See @cite Chetyrkin:1994js, @cite Bardin:1999ak, @cite Bardin:1999yd,
     * @cite Baikov:2012er and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$R_A^q(M_Z^2)@f$
     */
    double RAq(const QCD::quark q) const;

    /**
     * @brief The singlet vector corrections to the hadronic @f$Z@f$-boson width,
     * denoted as @f$R_V^h@f$.
     * @details In addition to the final-state corrections represented by
     * the radiator factors @f$R_V^q(M_Z^2)@f$ and @f$R_A^q(M_Z^2)@f$,
     * there exist singlet vector corrections to the total hadronic width
     * @cite Chetyrkin:1994js, @cite Baikov:2012er, which is much smaller than
     * the other corrections.
     * 
     * The assignment of the singlet vector corrections to the partial widths
     * is ambiguous @cite Bardin:1997xq. See Gamma_had() for our prescription.
     * @return @f$R_V^h@f$
     */
    double RVh() const;


    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to l^+ l^-@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_l@f$.
     * @details The asymmetry  @f$\mathcal{A}_l@f$ is given by
     * @f[
     * \mathcal{A}_l =
     * \frac{2\, {\rm Re}\left(g_{V}^l/g_{A}^l\right)}
     * {1+\left[{\rm Re}\left(g_{V}^l/g_{A}^l\right)\right]^2}\,,
     * @f]
     * where the ratio of the effective couplings @f$g_{V}^l/g_{A}^l@f$ is
     * computed via the two-loop approximate formula of 
     * @f$\sin^2\theta_{\rm eff}^{\,l}@f$, EWSMApproximateFormulae::sin2thetaEff_l(),
     * when checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel
     * is set to APPROXIMATEFORMULA.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\mathcal{A}_l@f$
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @brief The left-right asymmetry in @f$e^+e^-\to Z\to q^+ q^-@f$ at the
     * @f$Z@f$-pole, @f$\mathcal{A}_q@f$.
     * @details The asymmetry @f$\mathcal{A}_q@f$ is given by
     * @f[
     * \mathcal{A}_q =
     * \frac{2\, {\rm Re}\left(g_{V}^q/g_{A}^q\right)}
     * {1+\left[{\rm Re}\left(g_{V}^q/g_{A}^q\right)\right]^2}\,,
     * @f]
     * where the ratio of the effective couplings @f$g_{V}^q/g_{A}^q@f$ is
     * computed via the two-loop approximate formula of
     * @f$\sin^2\theta_{\rm eff}^{\,q}@f$, EWSMApproximateFormulae::sin2thetaEff_q(),
     * when checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel
     * is set to APPROXIMATEFORMULA.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\mathcal{A}_q@f$
     */
    double A_q(const QCD::quark q) const;

    /**
     * @brief The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\,l}@f$
     * for @f$Zl\bar{l}@f$ at the the @f$Z@f$-mass scale.
     * @details
     * When checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel 
     * is set to APPROXIMATEFORMULA, this function uses the two-loop approximate
     * formula of @f$\sin^2\theta_{\rm eff}^{\,l}@f$ via
     * EWSMApproximateFormulae::sin2thetaEff_l().
     * Otherwise, the effective weak mixing angle is calculated from the coupling
     * @f$\kappa_Z^l@f$:
     * @f[
     * \sin^2\theta_{\rm eff}^{\,l} = {\rm Re}(\kappa_Z^l)\,s_W^2\,.
     * @f]
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\sin^2\theta_{\rm eff}^{\,l}@f$
     */
    double sin2thetaEff(const StandardModel::lepton l) const;
    
    /**
     * @brief The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\,q}@f$
     * for @f$Zq\bar{q}@f$ at the the @f$Z@f$-mass scale.
     * @details
     * When checkNPZff_linearized() returns true and
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel 
     * is set to APPROXIMATEFORMULA, this function uses the two-loop approximate
     * formula of @f$\sin^2\theta_{\rm eff}^{\,q}@f$ via
     * EWSMApproximateFormulae::sin2thetaEff_q().
     * Otherwise, the effective weak mixing angle is calculated from the coupling
     * @f$\kappa_Z^q@f$:
     * @f[
     * \sin^2\theta_{\rm eff}^{\,q} = {\rm Re}(\kappa_Z^q)\,s_W^2\,.
     * @f]
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\sin^2\theta_{\rm eff}^{\,q}@f$
     */
    double sin2thetaEff(const QCD::quark q) const;   
    
    /**
     * @brief The @f$Z\to l^+ l^-@f$ partial decay width, @f$\Gamma_l@f$.
     * @details
     * When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\Gamma_l@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the partial width is calculated with
     * @f$\rho_Z^l@f$ and @f$g_{V}^l/g_{A}^l@f$ @cite Bardin:1999ak :
     * @f[
     * \Gamma_l =
     * \Gamma_0 \big|\rho_Z^f\big|
     * \sqrt{1-\frac{4m_l^2}{M_Z^2}}
     * \left[ \left(1+\frac{2m_l^2}{M_Z^2}\right)
     *   \left(\left|\frac{g_{V}^l}{g_{A}^l}\right|^2 + 1 \right)
     *   - \frac{6m_l^2}{M_Z^2}
     * \right]
     * \left( 1 + \frac{3}{4}\frac{\alpha(M_Z^2)}{\pi}\, Q_l^2 \right)
     * @f]
     * with @f$\Gamma_0=G_\mu M_Z^3/(24\sqrt{2}\pi)@f$. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\Gamma_l@f$ in GeV
     */
    double Gamma_l(const StandardModel::lepton l) const;
    
    /**
     * @brief The @f$Z\to q\bar{q}@f$ partial decay width, @f$\Gamma_q@f$.
     * @details
     * When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\Gamma_q@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the partial width is calculated with
     * @f$\rho_Z^q@f$ and @f$g_{V}^q/g_{A}^q@f$ @cite Bardin:1999ak :
     * @f[
     * \Gamma_q = N_c\, \Gamma_0 \big|\rho_Z^q\big|\,
     * \left[ \left|\frac{g_{V}^q}{g_{A}^q}\right|^2 R_V^q(M_Z^2)
     *   + R_A^q(M_Z^2) \right]
     * + \Delta_{\rm EW/QCD}
     * @f]
     * with @f$\Gamma_0=G_\mu M_Z^3/(24\sqrt{2}\pi)@f$. 
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Gamma_q@f$ in GeV
     */
    double Gamma_q(const QCD::quark q) const;
    
    /**
     * @brief The invisible partial decay width of the @f$Z@f$ boson,
     * @f$\Gamma_{\mathrm{inv}}@f$.
     * @details
     * @f[
     * \Gamma_{\mathrm{inv}} = 3\,\Gamma_\nu\,,
     * @f]
     * where @f$\Gamma_{\nu}@f$ is the partial width for @f$Z\to\nu\bar{\nu}@f$.
     * @return @f$\Gamma_{\mathrm{inv}}@f$ in GeV
     */
    double Gamma_inv() const;

    /**
     * @brief The hadronic decay width of the @f$Z@f$ boson, @f$\Gamma_{h}@f$.
     * @details
     * The hadronic width is given by the sum, 
     * @f[
     *  \Gamma_h = \Gamma_u + \Gamma_d + \Gamma_c + \Gamma_s + \Gamma_b\,.
     * @f]
     * Furthermore, the singlet vector corrections are added, following the
     * prescription in @cite Bardin:1997xq :
     * @f[
     * \Gamma_h = \sum_q \Gamma_q + 4N_c\Gamma_0 R_V^h\,.
     * @f]
     * @return @f$\Gamma_{h}@f$ in GeV
     */
    double Gamma_had() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\Gamma_Z@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the total decay width is calculated with
     * @f[
     * \Gamma_Z = \Gamma_{e} + \Gamma_{\mu} + \Gamma_{\tau} 
     * + \Gamma_{\mathrm{inv}} + \Gamma_h\,.
     * @f]
     * @return @f$\Gamma_Z@f$ in GeV
     */
    double Gamma_Z() const;

    /**
     * @brief The hadronic cross section for @f$e^+e^- \to Z \to \mathrm{hadrons}@f$
     * at the @f$Z@f$-pole, @f$\sigma_h^0@f$.
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$\sigma_h^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, the hadronic cross section is calculated with
     * @f[
     * \sigma_h^0 = \frac{12\pi}{M_Z^2}\frac{\Gamma_e\Gamma_h}{\Gamma_Z^2}\,.
     * @f]
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     */
    double sigma0_had() const; 

    /**
     * @brief @copybrief Rlepton::computeThValue()
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$R_\ell^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, @f$R_\ell^0@f$ is calculated with
     * @f[
     * R_\ell^0 = \frac{\Gamma_h}{\Gamma_\ell}\,.
     * @f]
     * @return @f$R_\ell^0 @f$
     */
    double R0_l() const;

    /**
     * @brief @copybrief Rcharm::computeThValue()
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$R_c^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, @f$R_c^0@f$ is calculated with
     * @f[
     * R_c^0 = \frac{\Gamma_c}{\Gamma_h}\,.
     * @f]
     * @return @f$R_c^0 @f$
     */
    double R0_c() const;

    /**
     * @brief @copybrief Rbottom::computeThValue()
     * @details When checkNPZff_linearized() returns true and the model flag
     * @ref StandardModelFlags "NoApproximateGammaZ" of StandardModel is set
     * to false, this function uses the two-loop approximate formula of
     * @f$R_b^0@f$ via EWSMApproximateFormulae::X_extended().
     * Otherwise, @f$R_b^0@f$ is calculated with
     * @f[
     * R_b^0 = \frac{\Gamma_b}{\Gamma_h}\,.
     * @f]
     * @return @f$R_b^0 @f$
     */
    double R0_b() const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    const EW_NPZff myEW_NPZff;///< An object of type EW_NPZff.

};

#endif	/* EW_H */

