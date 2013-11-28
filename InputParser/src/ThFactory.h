/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @details The list of available observables are as follows:
 * \li \b Mw: The mass of the \f$ W \f$ boson, \f$ M_W \f$
 * \li \b GammaW:
 * \li \b GammaZ:
 * \li \b sigmaHadron:
 * \li \b sin2thetaEff:
 * \li \b PtauPol:
 * \li \b Alepton:
 * \li \b Acharm:
 * \li \b Abottom:
 * \li \b AFBlepton:
 * \li \b AFBcharm:
 * \li \b AFBbottom:
 * \li \b Rlepton:
 * \li \b Rcharm:
 * \li \b Rbottom:
 * \li \b epsilon1:
 * \li \b epsilon2:
 * \li \b epsilon3:
 * \li \b epsilonb:
 * \li \b sigmaqLEP2_:
 * \li \b sigmamuLEP2_:
 * \li \b sigmatauLEP2_:
 * \li \b AFBmuLEP2_:
 * \li \b AFBtauLEP2_:
 * \li \b AFBbottomLEP2_:
 * \li \b AFBcharmLEP2_:
 * \li \b RbottomLEP2_:
 * \li \b RcharmLEP2_:
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
 */
class ThFactory {
public:
    
    /**
     * @brief The default constructor.\n
     */
    ThFactory(const StandardModel& myModel);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ThFactory();
    
    /**
     * @brief
     * @param[in]
     */
    ThObservable* getThMethod(const std::string& name);

private:
    std::map<std::string, ThObservable *> thobs; /**< */
    EW myEW; /**< */
    Flavour myFlavour; /**< */
    ModelObservable myMO; /**< */
};

#endif	/* THFACTORY_H */