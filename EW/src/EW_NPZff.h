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
 * @brief A class for new physics contributions to @f$Z@f$-pole pseudo observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains functions to incorporate new physics (NP)
 * contributions to the @f$Z@f$-pole pseudo observables, parameterized in terms
 * of the @f$Z f\bar{f}@f$ effective couplings, @f$g_V^f@f$ and @f$g_A^f@f$. 
 *
 * NP contributions are linearized in the corrections @f$\delta g_V^f@f$ and 
 * @f$\delta g_A^f@f$, as explained in the detailed description of EW class.
 * 
 */
class EW_NPZff {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EW_NPZff(const StandardModel& SM_i);

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$.
     * @param[in] GammaZ_SM the SM prediction for @f$\Gamma_Z@f$ in GeV
     * @return @f$\Gamma_Z@f$ in GeV, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double GammaZ(const double GammaZ_SM) const;

    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$.
     * @param[in] sigmaHadron_SM the SM prediction for @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$
     * @return @f$\sigma_h^0@f$ in GeV@f$^{-2}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double sigmaHadron(const double sigmaHadron_SM) const;
 
    /**
     * @brief @copybrief sin2thetaEff::computeThValue()
     * @param[in] sin2thetaEff_SM the SM prediction for @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$
     * @return @f$\sin^2\theta_{\rm eff}^{\rm lept}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double sin2thetaEff(const double sin2thetaEff_SM) const;

    /**
     * @brief @copybrief PtauPol::computeThValue()
     * @param[in] PtauPol_SM the SM prediction for @f$P_\tau^{\mathrm{pol}}@f$
     * @return @f$P_\tau^{\mathrm{pol}}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double PtauPol(const double PtauPol_SM) const;

    /**
     * @brief @copybrief Alepton::computeThValue()
     * @param[in] Alepton_SM the SM prediction for @f$\mathcal{A}_\ell@f$
     * @return @f$\mathcal{A}_\ell@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Alepton(const double Alepton_SM) const;

    /**
     * @brief @copybrief Acharm::computeThValue()
     * @param[in] Acharm_SM the SM prediction for @f$\mathcal{A}_c@f$
     * @return @f$\mathcal{A}_c@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Acharm(const double Acharm_SM) const;

    /**
     * @brief @copybrief Abottom::computeThValue()
     * @param[in] Abottom_SM the SM prediction for @f$\mathcal{A}_b@f$
     * @return @f$\mathcal{A}_b@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Abottom(const double Abottom_SM) const;

    /**
     * @brief @copybrief AFBlepton::computeThValue()
     * @param[in] AFBlepton_SM the SM prediction for @f$A^{0,\ell}_{\mathrm{FB}}@f$
     * @return @f$A^{0,\ell}_{\mathrm{FB}}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double AFBlepton(const double AFBlepton_SM) const;
    
    /**
     * @brief @copybrief AFBcharm::computeThValue()
     * @param[in] AFBcharm_SM the SM prediction for @f$A^{0,c}_{\mathrm{FB}}@f$
     * @return @f$A^{0,c}_{\mathrm{FB}}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double AFBcharm(const double AFBcharm_SM) const;

    /**
     * @brief @copybrief AFBbottom::computeThValue()
     * @param[in] AFBbottom_SM the SM prediction for @f$A^{0,b}_{\mathrm{FB}}@f$
     * @return @f$A^{0,b}_{\mathrm{FB}}@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double AFBbottom(const double AFBbottom_SM) const;

    /**
     * @brief @copybrief Rlepton::computeThValue()
     * @param[in] Rlepton_SM the SM prediction for @f$R_\ell^0@f$
     * @return @f$R_\ell^0@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Rlepton(const double Rlepton_SM) const;

    /**
     * @brief @copybrief Rcharm::computeThValue()
     * @param[in] Rcharm_SM the SM prediction for @f$R_c^0@f$
     * @return @f$R_c^0@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Rcharm(const double Rcharm_SM) const;

    /**
     * @brief @copybrief Rbottom::computeThValue()
     * @param[in] Rbottom_SM the SM prediction for @f$R_b^0@f$
     * @return @f$R_b^0@f$, including SM plus NP contributions
     *
     * @attention This function is applicable only to the NP model classes that
     * are inherited from NPbase.
     */
    double Rbottom(const double Rbottom_SM) const;

    
    ////////////////////////////////////////////////////////////////////////
private:
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    
};

#endif	/* EW_NPZFF_H */

