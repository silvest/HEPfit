/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWPO.h"
#include "StandardModel.h"

EWPO::EWPO(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
};

double EWPO::computeThValue()
{
    return 0.0;
}

void EWPO::computeTHDMcouplings()
{
//    double mHl=myTHDM->getMHl();
//    double mA=myTHDM->getMA();
//    double mHh=myTHDM->getMHh();
//    double mHp=myTHDM->getMHp();
//    double vev=myTHDM->v();
//    double sina=myTHDM->computeSina();
//    double cosa=myTHDM->computeCosa();
//    double tanb=myTHDM->getTanb();
//    double sinb=myTHDM->getSinb();
//    double cosb=myTHDM->getCosb();
//    double sin_ba=myTHDM->getsin_ba();
//    double m12_2=myTHDM->getM12_2();
    
//    get MW shift
//    calculate deltagV and deltagA = g(MW_corr)-g(MW_SM)+deltagloop
//    gslpp::complex deltagVl=0.0;
//    gslpp::complex deltagAl=0.0;
//    calculate the Deltas:
//    double DeltaAl=AlSM*(deltagVl.real()/gVlTL
//                         +deltagAl.real()/gAlTL
//                         -2*(gVlTL*deltagVl.real()+gAlTL*deltagAl.real())/(gVlTL*gVlTL+gAlTL*gAlTL));
////    double DeltaPpoltau=0.0;
//    double DeltaAc=AcSM*(deltagVc.real()/gVuTL
//                         +deltagAc.real()/gAuTL
//                         -2*(gVuTL*deltagVc.real()+gAuTL*deltagAc.real())/(gVuTL*gVuTL+gAuTL*gAuTL));
//    double DeltaAb=AbSM*(deltagVb.real()/gVdTL
//                         +deltagAb.real()/gAdTL
//                         -2*(gVdTL*deltagVb.real()+gAdTL*deltagAb.real())/(gVdTL*gVdTL+gAdTL*gAdTL));
////    double DeltaAFBl0=0.0;
////    double DeltaAFBc0=0.0;
////    double DeltaAFBb0=0.0;
////    double DeltaGammaZ=0.0;
////    double DeltaRl0=0.0;
////    double DeltaRc0=0.0;
////    double DeltaRb0=0.0;
////    double DeltaSigmahad=0.0;
////    double DeltaGammaW=0.0;
////    double Deltasinthetaeffl_2=0.0;
    
    
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

AlTHDM::AlTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AlTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAl=0.0;
    double AlSM=myTHDM->A_f(SM.getLeptons(SM.ELECTRON));
    return AlSM+DeltaAl;
}

PpoltauTHDM::PpoltauTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double PpoltauTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaPpoltau=0.0;
    double PpoltauSM=myTHDM->A_f(SM.getLeptons(SM.TAU));
    return PpoltauSM+DeltaPpoltau;
}

AcTHDM::AcTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AcTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAc=0.0;
    double AcSM=myTHDM->A_f(SM.getQuarks(SM.CHARM));
    return AcSM+DeltaAc;
}

AbTHDM::AbTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AbTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAb=0.0;
    double AbSM=myTHDM->A_f(SM.getQuarks(SM.BOTTOM));
    return AbSM+DeltaAb;
}

AFBl0THDM::AFBl0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AFBl0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAFBl0=0.0;
    double AFBl0SM=myTHDM->AFB(SM.getLeptons(SM.ELECTRON));
    return AFBl0SM+DeltaAFBl0;
}

AFBc0THDM::AFBc0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AFBc0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAFBc0=0.0;
    double AFBc0SM=myTHDM->AFB(SM.getQuarks(SM.CHARM));
    double AFBc0THDM = AFBc0SM+DeltaAFBc0;
    return 0.0;
}

AFBb0THDM::AFBb0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double AFBb0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaAFBb0=0.0;
    double AFBb0SM=myTHDM->AFB(SM.getQuarks(SM.BOTTOM));
    double AFBb0THDM = AFBb0SM+DeltaAFBb0;
    return 0.0;
}

GammaZTHDM::GammaZTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double GammaZTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaGammaZ=0.0;
    double GammaZSM=myTHDM->Gamma_Z();
    double GammaZTHDM = GammaZSM+DeltaGammaZ;
    return 0.0;
}

Rl0THDM::Rl0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double Rl0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaRl0=0.0;
    double Rl0SM=myTHDM->R0_f(SM.getLeptons(SM.ELECTRON));
    double Rl0THDM = Rl0SM+DeltaRl0;
    return 0.0;
}

Rc0THDM::Rc0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double Rc0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaRc0=0.0;
    double Rc0SM=myTHDM->R0_f(SM.getQuarks(SM.CHARM));
    double Rc0THDM = Rc0SM+DeltaRc0;
    return 0.0;
}

Rb0THDM::Rb0THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double Rb0THDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaRb0=0.0;
    double Rb0SM=myTHDM->R0_f(SM.getQuarks(SM.BOTTOM));
    double Rb0THDM = Rb0SM+DeltaRb0;
    return 0.0;
}

SigmahadTHDM::SigmahadTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double SigmahadTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaSigmahad=0.0;
    double SigmahadSM=myTHDM->sigma0_had();
    double SigmahadTHDM = SigmahadSM+DeltaSigmahad;
    return 0.0;
}

GammaWTHDM::GammaWTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double GammaWTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaGammaW=0.0;
    double GammaWSM=myTHDM->GammaW();
    double GammaWTHDM = GammaWSM+DeltaGammaW;
    return 0.0;
}

sinthetaeffl_2THDM::sinthetaeffl_2THDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double sinthetaeffl_2THDM::computeThValue()
{
    computeTHDMcouplings();
    double Deltasinthetaeffl_2=0.0;
    double sinthetaeffl_2SM=myTHDM->sin2thetaEff(SM.getLeptons(SM.ELECTRON));
    double sinthetaeffl_2THDM = sinthetaeffl_2SM+Deltasinthetaeffl_2;
    return 0.0;
}

MWTHDM::MWTHDM(const StandardModel& SM_i)
: EWPO(SM_i)
{}

double MWTHDM::computeThValue()
{
    computeTHDMcouplings();
    double DeltaMW=0.0;
    double MWSM=myTHDM->Mw();
    return MWSM+DeltaMW;
}
