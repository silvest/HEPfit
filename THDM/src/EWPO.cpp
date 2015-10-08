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
            "ThFactory::ThFactory() can only be 1 (AlTHDM) or 2 (PpoltauTHDM) "
            "or 3 (AcTHDM) or 4 (AbTHDM) or 5 (AFBl0THDM) or 6 (AFBc0THDM) or 7 (AFBb0THDM) "
            "or 8 (GammaZTHDM) or 9 (Rl0THDM) or 10 (Rc0THDM) or 11 (Rb0THDM) "
            "or 12 (SigmahadTHDM) or 13 (GammaWTHDM) or 14 (sinthetaeffl_2THDM) or 15 (MWTHDM)");
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

    double AlTHDM = AlSM;
    double PpoltauTHDM = PpoltauSM;
    double AcTHDM = AcSM;
    double AbTHDM = AbSM;

    double AFBl0THDM = AFBl0SM;
    double AFBc0THDM = AFBc0SM;
    double AFBb0THDM = AFBb0SM;

    double GammaZTHDM = GammaZSM;
    double Rl0THDM = Rl0SM;
    double Rc0THDM = Rc0SM;
    double Rb0THDM = Rb0SM;
    double SigmahadTHDM = SigmahadSM;

    double GammaWTHDM = GammaWSM;

    double sinthetaeffl_2THDM = sinthetaeffl_2SM;
    double MWTHDM = MWSM;

    if (obs == 1) return(AlTHDM);
    if (obs == 2) return(PpoltauTHDM);
    if (obs == 3) return(AcTHDM);
    if (obs == 4) return(AbTHDM);
    if (obs == 5) return(AFBl0THDM);
    if (obs == 6) return(AFBc0THDM);
    if (obs == 7) return(AFBb0THDM);
    if (obs == 8) return(GammaZTHDM);
    if (obs == 9) return(Rl0THDM);
    if (obs == 10) return(Rc0THDM);
    if (obs == 11) return(Rb0THDM);
    if (obs == 12) return(SigmahadTHDM);
    if (obs == 13) return(GammaWTHDM);
    if (obs == 14) return(sinthetaeffl_2THDM);
    if (obs == 15) return(MWTHDM);

    throw std::runtime_error("EWPO::computeThValue(): Observable type not defined. Can be only any of (1) till (15)");
    return (EXIT_FAILURE);

}
