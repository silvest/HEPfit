/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMstability.h"

stability_GTHDM::stability_GTHDM(const StandardModel& SM_i)
: myGTHDM(static_cast<const GeneralTHDM&> (SM_i)), vecMinus1(4,-1.), vecStability(4,0.)
{}

stability_GTHDM::~stability_GTHDM() 
{}

gslpp::vector<double> stability_GTHDM::getStability()
{
    double lambda1 = myGTHDM.getMyGTHDMCache()->lambda1_GTHDM;
    double lambda2 = myGTHDM.getMyGTHDMCache()->lambda2_GTHDM;
    double lambda3 = myGTHDM.getMyGTHDMCache()->lambda3_GTHDM;
    double lambda4 = myGTHDM.getMyGTHDMCache()->lambda4_GTHDM;
    double Relambda5 = myGTHDM.getMyGTHDMCache()->Relambda5_GTHDM;
    double Imlambda5 = myGTHDM.getImlambda5();
    double Relambda6 = myGTHDM.getRelambda6();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda6 = myGTHDM.getImlambda6();
    double Imlambda7 = myGTHDM.getImlambda7();
    
    vecStability(0) = lambda3 + sqrt(lambda1*lambda2);
    vecStability(1) = lambda3 + lambda4 - sqrt(Relambda5*Relambda5+Imlambda5*Imlambda5) + sqrt(lambda1*lambda2);
    vecStability(2) = (lambda1+lambda2)/2. + lambda3 + lambda4 + Relambda5 - 2.*fabs(Relambda6 + Relambda7);
    vecStability(3) = (lambda1+lambda2)/2. + lambda3 + lambda4 - Relambda5 - 2.*fabs(Imlambda6 + Imlambda7);
    
    if (lambda1>0 && lambda2>0)
        return vecStability;
    else
        return vecMinus1;
}


stability1_GTHDM::stability1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability1_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(0);
}


stability2_GTHDM::stability2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability2_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(1);
}


stability3_GTHDM::stability3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability3_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(2);
}


stability4_GTHDM::stability4_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability4_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(3);
}