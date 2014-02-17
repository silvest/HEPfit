/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THFACTORY_H
#define	THFACTORY_H

#include <ThObservable.h>
#include <ModelObservable.h>
#include <StandardModel.h>
#include <EW.h>
#include <Flavour.h>
#include <LeptonFlavour.h>

/**
 * @class ThFactory
 * @ingroup InputParser 
 * @brief A class for an interface to theoretical predictions for observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The list of available observables are as follows:
 *
 * @anchor ListOfObservables
 *
 * @li Electroweak Precision Observables:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mw</td>
 *   <td class="obs_symb">@f$ M_W @f$</td>
 *   <td class="obs_desc">@copybrief Mw::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%GammaW</td>
 *   <td class="obs_symb">@f$ \Gamma_W @f$</td>
 *   <td class="obs_desc">@copybrief GammaW::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%GammaZ</td>
 *   <td class="obs_symb">@f$ \Gamma_Z @f$</td>
 *   <td class="obs_desc">@copybrief GammaZ::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmaHadron</td>
 *   <td class="obs_symb"> @f$\sigma_h^0@f$</td>
 *   <td class="obs_desc">@copybrief sigmaHadron::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sin2thetaEff</td>
 *   <td class="obs_symb">@f$\sin^2\theta_{\rm eff}^{\rm lept}@f$</td>
 *   <td class="obs_desc">@copybrief sin2thetaEff::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%PtauPol</td>
 *   <td class="obs_symb">@f$P_\tau^\mathrm{pol}@f$</td>
 *   <td class="obs_desc">@copybrief PtauPol::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Alepton</td>
 *   <td class="obs_symb">@f$\mathcal{A}_\ell@f$</td>
 *   <td class="obs_desc">@copybrief Alepton::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Acharm</td>
 *   <td class="obs_symb">@f$\mathcal{A}_c@f$</td>
 *   <td class="obs_desc">@copybrief Acharm::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Abottom</td>
 *   <td class="obs_symb">@f$\mathcal{A}_b@f$</td>
 *   <td class="obs_desc">@copybrief Abottom::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBlepton</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^{0,\ell}@f$</td>
 *   <td class="obs_desc">@copybrief AFBlepton::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBcharm</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^{0,c}@f$</td>
 *   <td class="obs_desc">@copybrief AFBcharm::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBbottom</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^{0,b}@f$</td>
 *   <td class="obs_desc">@copybrief AFBbottom::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rlepton</td>
 *   <td class="obs_symb">@f$R_\ell^0@f$</td>
 *   <td class="obs_desc">@copybrief Rlepton::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rcharm</td>
 *   <td class="obs_symb">@f$R_c^0@f$</td>
 *   <td class="obs_desc">@copybrief Rcharm::computeThValue()</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Rbottom</td>
 *   <td class="obs_symb">@f$R_b^0@f$</td>
 *   <td class="obs_desc">@copybrief Rbottom::computeThValue()</td>
 * </tr>
 * </table>
 *
 * @li Epsilon %Parameters:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon1</td>
 *   <td class="obs_symb">@f$\varepsilon_1@f$</td>
 *   <td class="obs_desc">The epsilon parameter @f$\varepsilon_1@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon2</td>
 *   <td class="obs_symb">@f$\varepsilon_2@f$</td>
 *   <td class="obs_desc">The epsilon parameter @f$\varepsilon_2@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilon3</td>
 *   <td class="obs_symb">@f$\varepsilon_3@f$</td>
 *   <td class="obs_desc">The epsilon parameter @f$\varepsilon_3@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%epsilonb</td>
 *   <td class="obs_symb">@f$\varepsilon_b@f$</td>
 *   <td class="obs_desc">The epsilon parameter @f$\varepsilon_b@f$.</td>
 * </tr>
 * </table>
 *
 * @li  LEP-II Two-Fermion Processes, where [sqrt_s] denotes the center of mass
 * energy:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmaqLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$\sigma(e^+e^-\to q\bar{q})@f$</td>
 *   <td class="obs_desc">The @f$e^+ e^-\rightarrow\mathrm{hadrons}@f$ cross section in pb measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmamuLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$\sigma(e^+e^-\to \mu^+\mu^-)@f$</td>
 *   <td class="obs_desc">The @f$e^+ e^-\rightarrow \mu^+ \mu^-@f$ cross section in pb measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%sigmatauLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$\sigma(e^+e^-\to \tau^+\tau^-)@f$</td>
 *   <td class="obs_desc">The @f$e^+ e^-\rightarrow \tau^+ \tau^-@f$ cross section in pb measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBmuLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^\mu@f$</td>
 *   <td class="obs_desc">The forward-backward asymmetry in @f$e^+ e^-\rightarrow \mu^+ \mu^-@f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBtauLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^\tau@f$</td>
 *   <td class="obs_desc">The forward-backward asymmetry in @f$e^+ e^-\rightarrow \tau^+ \tau^-@f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBbottomLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^b@f$</td>
 *   <td class="obs_desc">The forward-backward asymmetry in @f$e^+ e^-\rightarrow b \bar{b}@f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%AFBcharmLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$A_{\mathrm{FB}}^c@f$</td>
 *   <td class="obs_desc">The forward-backward asymmetry in @f$e^+ e^-\rightarrow c \bar{c}@f$ measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%RbottomLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$R_b@f$</td>
 *   <td class="obs_desc">The ratio of the @f$e^+ e^-\rightarrow b \bar{b}@f$ and @f$e^+ e^-\rightarrow\mathrm{hadrons}@f$ cross sections measured at LEP2.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%RcharmLEP2_[sqrt s]</td>
 *   <td class="obs_symb">@f$R_c@f$</td>
 *   <td class="obs_desc">The ratio of the @f$e^+ e^-\rightarrow c \bar{c}@f$ and @f$e^+ e^-\rightarrow\mathrm{hadrons}@f$ cross sections measured at LEP2.</td>
 * </tr>
 * </table>
 *
 * @li Flavour Observables:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Dmd1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Dms1</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%M12D</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%ArgD</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%EpsilonK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%EpsiloP_o_Epsilon</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%DmK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vud</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vus</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vub</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Vcb</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%alpha</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%alpha_2a</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%gamma</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%SJPsiK</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%SJPsiPhi</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BR_Bdmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BRbar_Bdmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Amumu_Bd</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Smumu_Bd</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BR_Bsmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%BRbar_Bsmumu</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Amumu_Bs</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Smumu_Bs</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * @li %Model %Parameters of StandardModel class:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mz</td>
 *   <td class="obs_symb">@f$M_Z@f$</td>
 *   <td class="obs_desc">The mass of the @f$Z@f$ boson in GeV.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%mHl</td>
 *   <td class="obs_symb">@f$m_h@f$</td>
 *   <td class="obs_desc">The Higgs mass in GeV.</td>
 * </tr>
 * </table>
 *
 * @li %Model %Parameters, etc. of NPEffective1 and NPEffective2 classes:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ1pPLUScHQ2p_NP</td>
 *   <td class="obs_symb">@f$C'_{HQ_1} + C'_{HQ_2}@f$</td>
 *   <td class="obs_desc">The sum of the coefficients of the operators
 *   @f${\cal O}'_{HQ_1}@f$ and @f${\cal O}'_{HQ_2}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ2pMINUScHQ2_NP</td>
 *   <td class="obs_symb">@f$C'_{HQ_2} - C_{HQ_2}@f$</td>
 *   <td class="obs_desc">The difference of the coefficients of the operators
 *   @f${\cal O}'_{HQ_2}@f$ and @f${\cal O}_{HQ_2}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%cHQ3pPLUScHQ3_NP</td>
 *   <td class="obs_symb">@f$C'_{HQ_3} + C_{HQ_3}@f$</td>
 *   <td class="obs_desc">The sum of the coefficients of the operators
 *   @f${\cal O}'_{HQ_3}@f$ and @f${\cal O}_{HQ_3}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%c_Ae_NP</td>
 *   <td class="obs_symb">@f$C[\mathcal{A}_\ell]@f$</td>
 *   <td class="obs_desc">A combination of the coefficients of the dim. 6
 *   operators, associated with @f$\mathcal{A}_\ell@f$.
 *   See @cite Ciuchini:2013pca.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%c_GammaZ_uds_NP</td>
 *   <td class="obs_symb">@f$C[\Gamma_{uds}]@f$</td>
 *   <td class="obs_desc">A combination of the coefficients of the dim. 6
 *   operators, associated with the light-quark contribution to @f$\Gamma_Z@f$. 
 *   See @cite Ciuchini:2013pca.</td>
 * </tr>
 * </table>
 *
 * @li %Model %Parameters, etc. of NPZbbbar class:
 * 
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGVb</td>
 *   <td class="obs_symb">@f$\delta g_V^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the vector coupling @f$g_V^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGAb</td>
 *   <td class="obs_symb">@f$\delta g_A^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the axial-vector coupling @f$g_A^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGLb</td>
 *   <td class="obs_symb">@f$\delta g_L^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the left-handed coupling @f$g_L^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaGRb</td>
 *   <td class="obs_symb">@f$\delta g_R^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the right-handed coupling @f$g_R^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaRhoZb</td>
 *   <td class="obs_symb">@f$\delta \rho_Z^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the coupling @f$\rho_Z^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%deltaKappaZb</td>
 *   <td class="obs_symb">@f$\delta \kappa_Z^b@f$</td>
 *   <td class="obs_desc">New physics contribution to the coupling @f$\kappa_Z^b@f$ for @f$Zb\bar{b}@f$.</td>
 * </tr>
 * </table>
 *
 * @li SUSY Spectra and Observables:
 *
 * <table class="observables">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHl</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHh</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHa</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%MHp</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msu1, %Msu2, %Msu3, %Msu4, %Msu5, %Msu6</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Msd1, %Msd2, %Msd3, %Msd4, %Msd5, %Msd6</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mch1, %Mch2</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mneu1, %Mneu2, %Mneu3, %Mneu4</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%Mw_dRho</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * <tr>
 *   <td class="obs_name">%OutputSLHAfromFH</td>
 *   <td class="obs_symb"></td>
 *   <td class="obs_desc"></td>
 * </tr>
 * </table>
 *
 * 
 */
class ThFactory {
public:
    
    /**
     * @brief Constructor.
     * @details All the theory observables are built by the constructor. These
     * are implemented as a map between the names of the observables and a
     * pointer to an object of type ThObservable.
     * @param[in] myModel a reference to an object of the type StandardModel
     */
    ThFactory(const StandardModel& myModel);
    
    /**
     * @brief Destructor for deleting all the non-NULL pointers
     * in the thobs map.
     */
    virtual ~ThFactory();
    
    /**
     * @brief This method checks for the existence of an observable of a specific name in the map
     * thobs and returns a pointer to the particular ThObservable.
     * @param[in] name the name of the observable
     * @return a pointer to an object of the type ThObservable
     */
    ThObservable* getThMethod(const std::string& name);

private:
    std::map<std::string, ThObservable *> thobs; ///< A map between a pointer of type ThObservable and the name of the observable.
    EW myEW; ///< An object of type EW.
    Flavour myFlavour; ///< An object of type Flavour.
    LeptonFlavour myLeptonFlavour;///< An object of type LeptonFlavour.
    ModelObservable myMO; ///< An object of type ModelObservable
};

#endif	/* THFACTORY_H */