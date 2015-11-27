/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMTWOFERMIONSLEP2_HOLLIK_H
#define	EWSMTWOFERMIONSLEP2_HOLLIK_H

#include <gslpp.h>
#include <Polylogarithms.h>
#include <PVfunctions.h>
#include "EWSMOneLoopEW_HV.h"

/**
 * @class EWSMTwoFermionsLEP2_Hollik
 * @ingroup StandardModel
 * @brief A class for the cross sections and forward-backward asymmetries of @f$e^+e^-\to f\bar{f}@f$ at LEP-II. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The formulae used in the current class are referred to Hollik's pape, Fortschr. 
 * Phys 38 (1990), 3, 165, and calculated in the 't Hooft-Feynman gauge. 
 */
class EWSMTwoFermionsLEP2_Hollik {
public:

    /**
     * @brief Constructor. 
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSMTwoFermionsLEP2_Hollik(const StandardModel& SM_i);


    ////////////////////////////////////////////////////////////////////////  

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;

    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const QCD::quark q, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bQED with/without QED corrections
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const StandardModel::lepton l, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;

    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bQED with/without QED corrections
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const QCD::quark q, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;


    ////////////////////////////////////////////////////////////////////////  

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l_old(const StandardModel::lepton l, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bQED) const;

    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q_old(const QCD::quark q, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bQED) const;


    ////////////////////////////////////////////////////////////////////////  

    void TEST(const double s, const double Mw) const;


    ////////////////////////////////////////////////////////////////////////  
private:
    const StandardModel& SM;
    const EWSMOneLoopEW_HV myOneLoopEW_HV;
    const Polylogarithms Polylog;
    const PVfunctions PV;

    bool bUseHollik;


    ////////////////////////////////////////////////////////////////////////  

    /**
     * @param[in] l name of lepton
     * @param[in] Mw the W-boson mass
     * @return the tree-level vector coupling for Z->l lbar
     */
    double vl(const StandardModel::lepton l, const double Mw) const;

    /**
     * @param[in] q name of quark
     * @param[in] Mw the W-boson mass
     * @return the tree-level vector coupling for Z->q qbar
     */
    double vq(const QCD::quark q, const double Mw) const;

    /**
     * @param[in] l name of lepton
     * @param[in] Mw the W-boson mass
     * @return the tree-level axial-vector coupling for Z->l lbar
     */
    double al(const StandardModel::lepton l, const double Mw) const;

    /**
     * @param[in] q name of quark
     * @param[in] Mw the W-boson mass
     * @return the tree-level axial-vector coupling for Z->q qbar
     */
    double aq(const QCD::quark q, const double Mw) const;


    ////////////////////////////////////////////////////////////////////////  

    // Renormalized self-energies
    gslpp::complex Sigma_hat_ZZ(const double mu, const double s, const double Mw) const;
    gslpp::complex Sigma_hat_gZ(const double mu, const double s, const double Mw) const;
    gslpp::complex Sigma_hat_gg(const double mu, const double s, const double Mw) const;

    // Dressed gauge-boson propagators
    gslpp::complex chi_Z(const double mu, const double s, const double Mw,
            const bool bDP) const;
    gslpp::complex chi_gamma(const double mu, const double s, const double Mw,
            const bool bDP) const;
    gslpp::complex chi_gammaZ(const double mu, const double s, const double Mw,
            const bool bDP) const;

    // Renormalized vertex form factors for the Z-f-f vertex (non-QED part)
    gslpp::complex FVZ_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex FAZ_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex FVZ_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex FAZ_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex FL_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex FL_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex FL_u(const double s, const double Mw) const;
    gslpp::complex FL_d(const double s, const double Mw) const;
    gslpp::complex Lambda2(const double s, const double M) const;
    gslpp::complex Lambda3(const double s, const double M) const;
    // for f=b
    gslpp::complex deltaZL_fin(const double s, const double Mw) const;
    gslpp::complex Fb(const double s, const double Mw) const;
    gslpp::complex Fc(const double s, const double Mw) const;
    gslpp::complex Fd(const double s, const double Mw) const;
    gslpp::complex Fe(const double s, const double Mw) const;
    gslpp::complex Ff(const double s, const double Mw) const;
    gslpp::complex Fg(const double s, const double Mw) const;

    // Renormalized vertex form factors for the gamma-f-f vertex (non-QED part)
    gslpp::complex FVgamma_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex FAgamma_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex FVgamma_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex FAgamma_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex GL_l(const StandardModel::lepton l, const double s, const double Mw) const;
    gslpp::complex GL_q(const QCD::quark q, const double s, const double Mw) const;
    gslpp::complex GL_u(const double s, const double Mw) const;
    gslpp::complex GL_d(const double s, const double Mw) const;
    // for f=b
    gslpp::complex Gb(const double s, const double Mw) const;
    gslpp::complex Gc(const double s, const double Mw) const;
    gslpp::complex Gd(const double s, const double Mw) const;
    gslpp::complex Ge(const double s, const double Mw) const;
    gslpp::complex Gf(const double s, const double Mw) const;
    gslpp::complex Gg(const double s, const double Mw) const;

    // Born + dressed propagators + non-QED vertex corrections 
    gslpp::complex V_e(const int j, const double s, const double Mw, const bool bWEAK) const;
    gslpp::complex A_e(const int j, const double s, const double Mw, const bool bWEAK) const;
    gslpp::complex V_l(const int j, const StandardModel::lepton l, const double s,
            const double Mw, const bool bWEAK) const;
    gslpp::complex V_q(const int j, const QCD::quark q, const double s,
            const double Mw, const bool bWEAK) const;
    gslpp::complex A_l(const int j, const StandardModel::lepton l, const double s,
            const double Mw, const bool bWEAK) const;
    gslpp::complex A_q(const int j, const QCD::quark q, const double s,
            const double Mw, const bool bWEAK) const;
    gslpp::complex chi(const int j, const double s, const double Mw, const bool bDP) const;

    // QED corrections to e^+e^- -> f fbar cross sections
    double delta() const;
    double Bf(const double s, const double mf) const;
    double gamma_delta(const double s, const double mf, const double Qf) const;
    gslpp::complex gamma_delta_int(const double s, const double GammaZ, const double mf, const double Qf) const;
    double gamma_delta_res(const double s, const double GammaZ, const double mf, const double Qf) const;
    double gamma_tail(const double s, const double GammaZ) const;
    double gamma_fin(const double s, const double mf, const double Qf) const;
    double C11V(const double s, const double mf, const double Qf) const;
    double C11A(const double s, const double mf, const double Qf) const;
    gslpp::complex C12V(const double s, const double GammaZ, const double mf, const double Qf) const;
    gslpp::complex C12A(const double s, const double mf, const double Qf) const;
    double C22V(const double s, const double GammaZ, const double mf, const double Qf) const;
    double C22A(const double s, const double mf, const double Qf) const;

    double G1_l(const StandardModel::lepton l, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;
    double G1_q(const QCD::quark q, const double s, const
            double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;
    double G2_l(const StandardModel::lepton l, const double s, const double Mw,
            const bool bDP) const;
    double G2_q(const QCD::quark q, const double s, const double Mw,
            const bool bDP) const;
    double G3_l(const StandardModel::lepton l, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;
    double G3_q(const QCD::quark q, const double s,
            const double Mw, const double GammaZ,
            const bool bDP, const bool bWEAK, const bool bQED) const;

    // Loop functions
    gslpp::complex B0bar_Hollik(const double s, const double m1, const double m2) const;
    gslpp::complex B1bar_Hollik(const double s, const double m1, const double m2) const;
    gslpp::complex B1barPrime_Hollik(const double s, const double m1, const double m2) const;
    gslpp::complex C0_Hollik(const double s, const double M, const double Mprime) const;
    gslpp::complex C1plus_Hollik(const double s, const double M, const double Mprime) const;
    gslpp::complex C2zero_Hollik(const double s, const double M, const double Mprime) const;
    gslpp::complex C2plus_Hollik(const double s, const double M, const double Mprime) const;
    gslpp::complex C2minus_Hollik(const double s, const double M, const double Mprime) const;

    double sigma_f_old(const double s, const double Mw, const double GammaZ,
            const double mf, const double Qf, const double I3f, const double Ncf,
            const bool bDP = true, const bool bQED = true) const;

};

#endif	/* EWSMTWOFERMIONSLEP2_HOLLIK_H */

