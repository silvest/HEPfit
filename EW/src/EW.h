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
 * @details This class includes basic functions for the computation of the SM
 * predictions for electroweak precision pseudo observables, such as partial
 * decay widths of the @f$Z@f$ boson, left-right asymmetries and cross sections
 * at the @f$Z@f$ pole.
 *
 *
 * new physics contributions are incorporated  into the effective couplings
 *
 * @f[
 * \mathcal{O}_i = \mathcal{O}(g_V^f, g_A^f)
 * @f]
 *
 *
 * @f[
 * \mathcal{O}_i 
 * = \mathcal{O}(g_{V,\mathrm{SM}}^f, g_{A,\mathrm{SM}}^f)
 *   + \mathcal{O}(\delta g_{V,\mathrm{NP}}^f, \delta g_{A,\mathrm{NP}}^f)
 * @f]
 *
 *
 * linealized
 * check by checkNPZff()
 *
 *
 *
 *
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

    /*
     * @brief The non-factorizable EW-QCD corrections to the partial widths
     * for @f$Z\to q\bar{q}@f$, denoted as @f$\Delta_{\mathrm{EW/QCD}}^q@f$.
     * @details
     * See @cite Czarnecki:1996ei and @cite Harlander:1997zb (and also
     * @ref Bardin:1999yd).
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Delta_{\mathrm{EW/QCD}}^q@f$ in GeV
     */
    double Delta_EWQCD(const QCD::quark q) const;

    /**
     * @param[in] q name of a quark (see QCD::quark)
     * @return Radiator functions to the vector current due to the
     * final-state QED and QCD corrections.
     */
    double RVq(const QCD::quark q) const;

    /**
     * @param[in] q name of a quark (see QCD::quark)
     * @return Radiator functions to the axial-vector current due to the
     * final-state QED and QCD corrections.
     */
    double RAq(const QCD::quark q) const;

    /**
     * @return Singlet vector corrections to the width of Z to hadrons.
     */
    double RVh() const;


    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The @f$Z@f$-pole leptonic left-right asymmetry, @f$A_l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return the asymmetry parameter for @f$Z\to l\bar{l}@f$, @f$A_l@f$
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @brief The @f$Z@f$-pole quark left-right asymmetry, @f$A_q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return the asymmetry parameter for @f$Z\to q\bar{q}@f$, @f$A_q@f$
     */
    double A_q(const QCD::quark q) const;

    /**
     * @brief The effective weak mixing angle @f$\sin^2\theta_{\rm eff}^{\,l}@f$
     * for @f$Zl\bar{l}@f$ at the the @f$Z@f$-mass scale.
     * @details
     * When checkNPZff() returns true and the model flag @ref StandardModelFlags
     * "KappaZ" is set to APPROXIMATEFORMULA, this functions uses the two-loop
     * approximate via EWSMApproximateFormulae::sin2thetaEff_l(). Otherwise,
     * the effective weak mixing angle is calculated from the coupling
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
     * When checkNPZff() returns true and the model flag @ref StandardModelFlags
     * "KappaZ" is set to APPROXIMATEFORMULA, this functions uses the two-loop
     * approximate via EWSMApproximateFormulae::sin2thetaEff_q(). Otherwise,
     * the effective weak mixing angle is calculated from the coupling
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
     * @details The partial width of @f$Z@f$ decaying into a charged-lepton pair,
     * including contribution from final-state QED interactions, is given in 
     * terms of the effective couplings @f$g_{V}^l@f$, @f$g_{A}^l@f$ and @f$\rho_Z^f@f$ 
     * @cite Bardin:1999ak :
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
     *
     *
     *
     *
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\Gamma_l@f$ in GeV
     */
    double Gamma_l(const StandardModel::lepton l) const;
    
    /**
     * @brief The @f$Z\to q\bar{q}@f$ partial decay width, @f$\Gamma_q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Gamma_q@f$ in GeV
     */
    double Gamma_q(const QCD::quark q) const;
    
    /**
     * @brief The @f$Z@f$-boson invisible partial decay width, @f$\Gamma_{inv}@f$.
     * @return the invisible decay width of the @f$Z@f$ boson in GeV
     */
    double Gamma_inv() const;

    /**
     * @brief The @f$Z\to\mbox{hadrons}@f$ partial decay width, @f$\Gamma_{had}@f$.
     * @return the hadronic decay width of the @f$Z@f$ boson in GeV
     */
    double Gamma_had() const;

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @return the total decay width of the @f$Z@f$ boson in GeV
     */
    double Gamma_Z() const;

    /**
     * @brief The @f$Z@f$-pole hadronic cross section, @f$\sigma_h^0@f$.
     * @return the cross section for @f$e^+e^- \to Z \to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole in GeV@f$^{-2}@f$
     */
    double sigma0_had() const; 

    double R0_l() const;
    double R0_c() const;
    double R0_b() const;



    
    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    const EW_NPZff myEW_NPZff;///< An object of type EW_NPZff.

};

#endif	/* EW_H */

