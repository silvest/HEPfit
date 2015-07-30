/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWPO.h"
#include "StandardModel.h"

EWPO::EWPO(const StandardModel& SM_i, int obsFlag) 
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    if (obsFlag > 0 and obsFlag < 16) obs = obsFlag;
    else throw std::runtime_error("obsFlag in EWPO() called from "
            "ThFactory::ThFactory() can only be 1 (Al) or 2 (Ppoltau) "
            "or 3 (Ac) or 4 (Ab) or 5 (AFBl0) or 6 (AFBc0) or 7 (AFBb0) "
            "or 8 (GammaZ) or 9 (Rl0) or 10 (Rc0) or 11 (Rb0) "
            "or 12 (Sigmahad) or 13 (GammaW) or 14 (sinthetaeffl_2) or 15 (MW)");
};

double EWPO::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA=myTHDM->getMA();
    double mHh=myTHDM->getMHh();
    double mHp=myTHDM->getMHp();
    double vev=myTHDM->v();
    double sina=myTHDM->computeSina();
    double cosa=myTHDM->computeCosa();
    double tanb=myTHDM->getTanb();
    double sinb=myTHDM->getSinb();
    double cosb=myTHDM->getCosb();
    double sin_ba=myTHDM->getsin_ba();
    double m12_2=myTHDM->getM12_2();

    double AlSM=myTHDM->A_f(SM.getLeptons(SM.ELECTRON));
    double PpoltauSM=myTHDM->A_f(SM.getLeptons(SM.TAU));
    double AcSM=myTHDM->A_f(SM.getQuarks(SM.CHARM));
    double AbSM=myTHDM->A_f(SM.getQuarks(SM.BOTTOM));
    double AFBl0SM=myTHDM->AFB(SM.getLeptons(SM.ELECTRON));
    double AFBc0SM=myTHDM->AFB(SM.getQuarks(SM.CHARM));
    double AFBb0SM=myTHDM->AFB(SM.getQuarks(SM.BOTTOM));
    double GammaZSM=myTHDM->Gamma_Z();
    double Rl0SM=myTHDM->R0_f(SM.getLeptons(SM.ELECTRON));
    double Rc0SM=myTHDM->R0_f(SM.getQuarks(SM.CHARM));
    double Rb0SM=myTHDM->R0_f(SM.getQuarks(SM.BOTTOM));
    double SigmahadSM=myTHDM->sigma0_had();
    double GammaWSM=myTHDM->GammaW();
    double sinthetaeffl_2SM=myTHDM->sin2thetaEff(SM.getLeptons(SM.ELECTRON));
    double MWSM=myTHDM->Mw();

    double Al = AlSM;
    double Ppoltau = PpoltauSM;
    double Ac = AcSM;
    double Ab = AbSM;

    double AFBl0 = AFBl0SM;
    double AFBc0 = AFBc0SM;
    double AFBb0 = AFBb0SM;

    double GammaZ = GammaZSM;
    double Rl0 = Rl0SM;
    double Rc0 = Rc0SM;
    double Rb0 = Rb0SM;
    double Sigmahad = SigmahadSM;

    double GammaW = GammaWSM;

    double sinthetaeffl_2 = sinthetaeffl_2SM;
    double MW = MWSM;

    if (obs == 1) return(Al);
    if (obs == 2) return(Ppoltau);
    if (obs == 3) return(Ac);
    if (obs == 4) return(Ab);
    if (obs == 5) return(AFBl0);
    if (obs == 6) return(AFBc0);
    if (obs == 7) return(AFBb0);
    if (obs == 8) return(GammaZ);
    if (obs == 9) return(Rl0);
    if (obs == 10) return(Rc0);
    if (obs == 11) return(Rb0);
    if (obs == 12) return(Sigmahad);
    if (obs == 13) return(GammaW);
    if (obs == 14) return(sinthetaeffl_2);
    if (obs == 15) return(MW);

    throw std::runtime_error("EWPO::computeThValue(): Observable type not defined. Can be only any of (1) till (15)");
    return (EXIT_FAILURE);

}
