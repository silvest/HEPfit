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
 * of corrections to the \f$Z\rightarrow f\bar{f}\f$ vertices.\n
 * New physics contributions are linearized in the corrections to the neutral-current
 * couplings, \f$g_V^f+\delta g_V^\f$,\f$g_A^f+\delta g_A^\f$.
 */
class EW_NPZff {
public:

    EW_NPZff(const StandardModel& SM_i);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] GammaZ_SM the SM prediction for \f$\Gamma_Z\f$ [GeV] 
     * @return the prediction for \f$\Gamma_Z\f$ [GeV] including SM plus new 
     * physics effects
     */
    double GammaZ(const double GammaZ_SM) const;

    /**
     * @param[in] sigmaHadron_SM the SM prediction for \f$\sigma_h^0\f$ [nb]
     * @return the prediction for \f$\sigma_h^0\f$ [nb] including SM plus new 
     * physics effects
     */
    double sigmaHadron(const double sigmaHadron_SM) const;
 
    /**
     * @param[in] sin2thetaEff_SM the SM prediction for \f$\sin^2{\theta_{\mathrm{Eff}}^\ell}\f$
     * @return the prediction for \f$\sin^2{\theta_{\mathrm{Eff}}^\ell}\f$ including SM plus new 
     * physics effects
     */
    double sin2thetaEff(const double sin2thetaEff_SM) const;

    /**
     * @param[in] PtauPol_SM the SM prediction for \f$P_{\tau}^\mathrm{pol}=A_\tau\f$
     * @return the prediction for \f$P_{\tau}^\mathrm{pol}=A_\tau\f$ including SM plus new 
     * physics effects
     */
    double PtauPol(const double PtauPol_SM) const;

    /**
     * @param[in] Alepton_SM the SM prediction for \f$A_\ell\f$
     * @return the prediction for \f$A_\ell\f$ including SM plus new 
     * physics effects
     */
    double Alepton(const double Alepton_SM) const;

    /**
     * @param[in] Acharm_SM the SM prediction for \f$A_c\f$
     * @return the prediction for \f$A_c\f$ including SM plus new 
     * physics effects
     */
    double Acharm(const double Acharm_SM) const;

    /**
     * @param[in] Abottom_SM the SM prediction for \f$A_b\f$
     * @return the prediction for \f$A_b\f$ including SM plus new 
     * physics effects
     */
    double Abottom(const double Abottom_SM) const;

    /**
     * @param[in] AFBlepton_SM the SM prediction for \f$A_{FB}^{0,\ell}\f$
     * @return the prediction for \f$A_{FB}^{0,\ell}\f$ including SM plus new 
     * physics effects
     */
    double AFBlepton(const double AFBlepton_SM) const;

    
    /**
     * @param[in] AFBcharm_SM the SM prediction for \f$A_{FB}^{0,c}\f$
     * @return the prediction for \f$A_{FB}^{0,c}\f$ including SM plus new 
     * physics effects
     */
    double AFBcharm(const double AFBcharm_SM) const;

    /**
     * @param[in] AFBbottom_SM the SM prediction for \f$A_{FB}^{0,b}\f$
     * @return the prediction for \f$A_{FB}^{0,b}\f$ including SM plus new 
     * physics effects
     */
    double AFBbottom(const double AFBbottom_SM) const;

    /**
     * @param[in] Rlepton_SM the SM prediction for \f$R_\ell\f$
     * @return the prediction for \f$R_\ell\f$ including SM plus new 
     * physics effects
     */
    double Rlepton(const double Rlepton_SM) const;

    /**
     * @param[in] Rcharm_SM the SM prediction for \f$R_c\f$
     * @return the prediction for \f$R_c\f$ including SM plus new 
     * physics effects
     */
    double Rcharm(const double Rcharm_SM) const;

    /**
     * @param[in] Rbottom_SM the SM prediction for \f$R_b\f$
     * @return the prediction for \f$R_b\f$ including SM plus new 
     * physics effects
     */
    double Rbottom(const double Rbottom_SM) const;
    
private:
    const StandardModel& SM;
    
};

#endif	/* EW_NPZFF_H */

