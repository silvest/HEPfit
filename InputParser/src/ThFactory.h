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

/**
 * @class ThFactory
 * @ingroup InputParser 
 * @brief A class for an interface to theoretical predictions for observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 *
 *
 *
 * @anchor ListOfObservables
 * <h3>%Observables</h3>
 *
 * The list of available observables are as follows:
 * \li \b Mw: The mass of the \f$ W \f$ boson, \f$ M_W \f$.
 * \li \b GammaW: The total decay width of the \f$W\f$ boson, \f$ \Gamma_W \f$.
 * \li \b GammaZ: The total decay width of the \f$Z\f$ boson, \f$ \Gamma_Z \f$.
 * \li \b sigmaHadron: The \f$e^+ e^- \rightarrow Z \rightarrow \mathrm{hadrons}\f$
 * cross section at the \f$Z\f$ pole, @f$\sigma_h^0@f$.
 * \li \b sin2thetaEff: The leptonic effective weak mixing angle, \f$\sin^2{\theta_\mathrm{Eff}^\ell}\f$. 
 * \li \b PtauPol: The longitudinal polarization of the final state in 
 * \f$e^+ e^- \rightarrow Z \rightarrow \tau^+ \tau^-\f$, \f$P_\tau^\mathrm{pol}\f$.
 * \li \b Alepton: The leptonic left-right asymmetry parameter at the \f$Z\f$
 * pole, \f$A_\ell\f$ (asummes lepton universality).
 * \li \b Acharm: The charm quark left-right asymmetry parameter, \f$A_c\f$.
 * \li \b Abottom: The bottom quark left-right asymmetry parameter, \f$A_b\f$.
 * \li \b AFBlepton: The leptonic forward-backward asymmetry at the \f$Z\f$ pole, 
 * \f$A_{FB}^{0,\ell}\f$ (asummes lepton universality).
 * \li \b AFBcharm: The charm quark forward-backward asymmetry at the \f$Z\f$ pole, 
 * \f$A_{FB}^{0,c}\f$.
 * \li \b AFBbottom: The bottom quark forward-backward asymmetry at the \f$Z\f$ pole, 
 * \f$A_{FB}^{0,b}\f$.
 * \li \b Rlepton: The ratio between the \f$Z\rightarrow \mathrm{hadrons}\f$ and 
 * \f$Z\rightarrow \ell^+ \ell^-\f$ decay widths, \f$R_\ell^0\f$ (one lepton flavour, asummes
 * lepton universality).
 * \li \b Rcharm: The ratio between the \f$Z\rightarrow c\bar{c}\f$ and 
 * \f$Z\rightarrow \mathrm{hadrons}\f$ decay widths, \f$R_c^0\f$.
 * \li \b Rbottom: The ratio between the \f$Z\rightarrow b\bar{b}\f$ and 
 * \f$Z\rightarrow \mathrm{hadrons}\f$ decay widths, \f$R_b^0\f$.
 * \li \b epsilon1:
 * \li \b epsilon2:
 * \li \b epsilon3:
 * \li \b epsilonb:
 * \li \b sigmaqLEP2_: The \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross section measured
 * at LEP2
 * \li \b sigmamuLEP2_: The \f$e^+ e^-\rightarrow \mu^+ \mu^-\f$ cross section measured
 * at LEP2
 * \li \b sigmatauLEP2_: The \f$e^+ e^-\rightarrow \tau^+ \tau^-\f$ cross section measured
 * at LEP2
 * \li \b AFBmuLEP2_: The forward-backward asymmetry in \f$e^+ e^-\rightarrow \mu^+ \mu^-\f$ measured
 * at LEP2
 * \li \b AFBtauLEP2_:The forward-backward asymmetry in \f$e^+ e^-\rightarrow \tau^+ \tau^-\f$ measured
 * at LEP2
 * \li \b AFBbottomLEP2_: The forward-backward asymmetry in \f$e^+ e^-\rightarrow b \bar{b}\f$ measured
 * at LEP2
 * \li \b AFBcharmLEP2_: The forward-backward asymmetry in \f$e^+ e^-\rightarrow c \bar{c}\f$ measured
 * at LEP2
 * \li \b RbottomLEP2_: The ratio of the \f$e^+ e^-\rightarrow b \bar{b}\f$ and 
 * \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross sections measured at LEP2
 * \li \b RcharmLEP2_: The ratio of the \f$e^+ e^-\rightarrow c \bar{c}\f$ and 
 * \f$e^+ e^-\rightarrow\mathrm{hadrons}\f$ cross sections measured at LEP2
 * \li \b Dmd1:
 * \li \b Dms1:
 * \li \b M12D:
 * \li \b ArgD:
 * \li \b EpsilonK:
 * \li \b EpsiloP_o_Epsilon:
 * \li \b DmK:
 * \li \b Vud:
 * \li \b Vus:
 * \li \b Vub:
 * \li \b Vcb:
 * \li \b alpha:
 * \li \b alpha_2a:
 * \li \b gamma:
 * \li \b SJPsiK:
 * \li \b SJPsiPhi:
 * \li \b BR_Bdmumu:
 * \li \b BRbar_Bdmumu:
 * \li \b Amumu_Bd:
 * \li \b Smumu_Bd:
 * \li \b BR_Bsmumu:
 * \li \b BRbar_Bsmumu:
 * \li \b Amumu_Bs:
 * \li \b Smumu_Bs:
 * \li \b AlsMz: The strong coupling constat at the \f$ Z \f$ pole, \f$\alpha_s(M_Z) \f$
 * \li \b dAle5Mz:
 * \li \b Mz: The mass of the \f$ Z \f$ boson, \f$ M_Z \f$
 * \li \b mtop: The mass of the top quark, \f$ M_{top} \f$
 * \li \b mHl:
 * \li \b delMw:
 * \li \b delSin2th_l:
 * \li \b delGammaZ:
 * \li \b delRhoZ_nu:
 * \li \b delRhoZ_e:
 * \li \b delRhoZ_u:
 * \li \b delRhoZ_d:
 * \li \b delRhoZ_b:
 * \li \b cHLp_NP:
 * \li \b cHQp_NP:
 * \li \b cHQ_NP:
 * \li \b cHL_NP:
 * \li \b cHE_NP:
 * \li \b c_Ae_NP:
 * \li \b c_GammaZ_uds_NP:
 * \li \b cHU2_NP
 * \li \b cHD3_NP
 * \li \b cHQ1pPLUScHQ2p_NP:
 * \li \b cHQ2pMINUScHQ2_NP:
 * \li \b cHQ3pPLUScHQ3_NP:
 * \li \b deltaGVb:
 * \li \b deltaGAb:
 * \li \b deltaGLb:
 * \li \b deltaGRb:
 * \li \b deltaRhoZb:
 * \li \b deltaKappaZb:
 * \li \b MHl:
 * \li \b MHh:
 * \li \b MHa:
 * \li \b MHp:
 * \li \b Msu1:
 * \li \b Msu2:
 * \li \b Msu3:
 * \li \b Msu4:
 * \li \b Msu5:
 * \li \b Msu6:
 * \li \b Msd1:
 * \li \b Msd2:
 * \li \b Msd3:
 * \li \b Msd4:
 * \li \b Msd5:
 * \li \b Msd6:
 * \li \b Mch1:
 * \li \b Mch2:
 * \li \b Mneu1:
 * \li \b Mneu2:
 * \li \b Mneu3:
 * \li \b Mneu4:
 * \li \b Mw_dRho:
 *
 *  
 * Masses and widhts are computed in GeV. Cross sections are given in nb.
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
     * @brief The default destructor. It destroys all the non-NULL pointers
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
    EW myEW; ///< An object of type EW().
    Flavour myFlavour; ///< An object of type Flavour.
    ModelObservable myMO; ///< An object of type ModelObservable
};

#endif	/* THFACTORY_H */