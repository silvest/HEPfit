/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMquantities.h"

tanbetaGM::tanbetaGM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double tanbetaGM::computeThValue()
{
    return myGM.gettanb();
}



m1sqGM::m1sqGM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double m1sqGM::computeThValue()
{
    double vPhi = myGM.getMyGMCache()->vPhi;
    double vDelta = myGM.getMyGMCache()->vDelta;
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double Mu1 = myGM.getMu1();
    return -4.0*lambda1*vPhi*vPhi-(6.0*lambda4+3.0*lambda5)*vDelta*vDelta-1.5*Mu1*vDelta;
}



m2sqGM::m2sqGM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double m2sqGM::computeThValue()
{
    double vPhi = myGM.getMyGMCache()->vPhi;
    double vDelta = myGM.getMyGMCache()->vDelta;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double Mu1 = myGM.getMu1();
    double Mu2 = myGM.getMu2();
    return -(12.0*lambda2+4.0*lambda3)*vDelta*vDelta-(2.0*lambda4+lambda5)*vPhi*vPhi-0.25*Mu1*vPhi*vPhi/vDelta-6.0*Mu2*vDelta;
}



lambda1GM::lambda1GM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double lambda1GM::computeThValue()
{
    return myGM.getMyGMCache()->lambda1;
}



lambda2GM::lambda2GM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double lambda2GM::computeThValue()
{
    return myGM.getMyGMCache()->lambda2;
}



lambda3GM::lambda3GM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double lambda3GM::computeThValue()
{
    return myGM.getMyGMCache()->lambda3;
}



lambda4GM::lambda4GM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double lambda4GM::computeThValue()
{
    return myGM.getMyGMCache()->lambda4;
}



lambda5GM::lambda5GM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double lambda5GM::computeThValue()
{
    return myGM.getMyGMCache()->lambda5;
}



vPhiGM::vPhiGM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double vPhiGM::computeThValue()
{
    return myGM.getMyGMCache()->vPhi;
}



vDeltaGM::vDeltaGM(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double vDeltaGM::computeThValue()
{
    return myGM.getMyGMCache()->vDelta;
}



rh_gaga_GM::rh_gaga_GM(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double rh_gaga_GM::computeThValue()
{
    return myGM.getMyGMCache()->rh_gaga;
}



rh_Zga_GM::rh_Zga_GM(const StandardModel& SM_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double rh_Zga_GM::computeThValue()
{
    return myGM.getMyGMCache()->rh_Zga;
}
