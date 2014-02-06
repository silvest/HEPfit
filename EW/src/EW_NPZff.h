/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_NPZFF_H
#define	EW_NPZFF_H

#include <StandardModel.h>

/**
 * @class EW_NPZff
 * @ingroup EW
 * @brief A class for new physics contributions to \f$Z\f$-pole pseudo observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains functions to incorporate new physics
 * contributions to \f$Z\f$-pole pseudo observables, parameterized in terms
 * of corrections to the \f$Z f\bar{f}\f$ vertices.\n
 * New physics contributions are linearized in the corrections to the neutral-current
 * couplings, \f$g_V^f+\delta g_V^f\f$,\f$g_A^f+\delta g_A^f\f$.
 */
class EW_NPZff {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of StandardModel class
     */
    EW_NPZff(const StandardModel& SM_i);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The total decay width of the \f$Z\f$ boson, \f$\Gamma_Z\f$.
     * @param[in] GammaZ_SM the SM prediction for \f$\Gamma_Z\f$ [GeV] 
     * @return the prediction for \f$\Gamma_Z\f$ [GeV] including SM plus new 
     * physics effects
     */
    double GammaZ(const double GammaZ_SM) const;

    /**
     * @brief The \f$Z\f$-pole hadronic cross section, \f$\sigma_h^0\f$.
     * @param[in] sigmaHadron_SM the SM prediction for \f$\sigma_h^0\f$ [GeV\f$^{-2}\f$]
     * @return the prediction for \f$\sigma_h^0\f$ [GeV\f$^{-2}\f$] including SM plus new 
     * physics effects
     */
    double sigmaHadron(const double sigmaHadron_SM) const;
 
    /**
     * @brief The effective weak mixing angle, \f$\sin^2{\theta_{Eff}^\ell}\f$.
     * @param[in] sin2thetaEff_SM the SM prediction for \f$\sin^2{\theta_{\mathrm{Eff}}^\ell}\f$
     * @return the prediction for \f$\sin^2{\theta_{\mathrm{Eff}}^\ell}\f$ including SM plus new 
     * physics effects
     */
    double sin2thetaEff(const double sin2thetaEff_SM) const;

    /**
     * @brief The longitudinal polarization in @f$Z\to \tau^+ \tau^-@f$, \f$P_{\tau}^{\mathrm{pol}}\f$.
     * @param[in] PtauPol_SM the SM prediction for \f$P_{\tau}^{\mathrm{pol}}=A_\tau\f$
     * @return the prediction for \f$P_{\tau}^\mathrm{pol}=A_\tau\f$ including SM plus new 
     * physics effects
     */
    double PtauPol(const double PtauPol_SM) const;

    /**
     * @brief The \f$Z\f$-pole leptonic left-right asymmetry, \f$A_l\f$.
     * @param[in] Alepton_SM the SM prediction for \f$A_\ell\f$
     * @return the prediction for \f$A_\ell\f$ including SM plus new 
     * physics effects
     */
    double Alepton(const double Alepton_SM) const;

    /**
     * @brief The \f$Z\f$-pole charm left-right asymmetry, \f$A_c\f$.
     * @param[in] Acharm_SM the SM prediction for \f$A_c\f$
     * @return the prediction for \f$A_c\f$ including SM plus new 
     * physics effects
     */
    double Acharm(const double Acharm_SM) const;

    /**
     * @brief The \f$Z\f$-pole bottom left-right asymmetry, \f$A_b\f$.
     * @param[in] Abottom_SM the SM prediction for \f$A_b\f$
     * @return the prediction for \f$A_b\f$ including SM plus new 
     * physics effects
     */
    double Abottom(const double Abottom_SM) const;

    /**
     * @brief The \f$Z\f$-pole leptonic forward-backward asymmetry, \f$A_{FB}^{0,\ell}\f$.
     * @param[in] AFBlepton_SM the SM prediction for \f$A_{FB}^{0,\ell}\f$
     * @return the prediction for \f$A_{FB}^{0,\ell}\f$ including SM plus new 
     * physics effects
     */
    double AFBlepton(const double AFBlepton_SM) const;
    
    /**
     * @brief The \f$Z\f$-pole charm forward-backward asymmetry, \f$A_{FB}^{0,c}\f$.
     * @param[in] AFBcharm_SM the SM prediction for \f$A_{FB}^{0,c}\f$
     * @return the prediction for \f$A_{FB}^{0,c}\f$ including SM plus new 
     * physics effects
     */
    double AFBcharm(const double AFBcharm_SM) const;

    /**
     * @brief The \f$Z\f$-pole bottom forward-backward asymmetry, \f$A_{FB}^{0,b}\f$.
     * @param[in] AFBbottom_SM the SM prediction for \f$A_{FB}^{0,b}\f$
     * @return the prediction for \f$A_{FB}^{0,b}\f$ including SM plus new 
     * physics effects
     */
    double AFBbottom(const double AFBbottom_SM) const;

    /**
     * @brief The ratio \f$R_\ell^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \ell^+ \ell^-)\f$.
     * @param[in] Rlepton_SM the SM prediction for \f$R_\ell^0\f$
     * @return the prediction for \f$R_\ell^0\f$ including SM plus new 
     * physics effects
     */
    double Rlepton(const double Rlepton_SM) const;

    /**
     * @brief The ratio \f$R_c^0=\Gamma(Z\to c\bar{c})/\Gamma(Z\to {\rm hadrons})\f$.
     * @param[in] Rcharm_SM the SM prediction for \f$R_c^0\f$
     * @return the prediction for \f$R_c^0\f$ including SM plus new 
     * physics effects
     */
    double Rcharm(const double Rcharm_SM) const;

    /**
     * @brief The ratio \f$R_b^0=\Gamma(Z\to b\bar{b})/\Gamma(Z\to {\rm hadrons})\f$.
     * @param[in] Rbottom_SM the SM prediction for \f$R_b^0\f$
     * @return the prediction for \f$R_b^0\f$ including SM plus new 
     * physics effects
     */
    double Rbottom(const double Rbottom_SM) const;
    
private:
    const StandardModel& SM;///< A reference to an object of the StandardModel class.
    
};

#endif	/* EW_NPZFF_H */

