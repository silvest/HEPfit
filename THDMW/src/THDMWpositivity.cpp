/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWpositivity.h"
#include "StandardModel.h"

THDMWpositivity1::THDMWpositivity1(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity1::computeThValue()
{
    double lambda1=myTHDMW.getTHDMW_lambda1();
    double lambda2=myTHDMW.getTHDMW_lambda2();
    double lambda3=myTHDMW.getTHDMW_lambda3();

    if (lambda1>0 && lambda2>0)
    {
        return lambda3+sqrt(lambda1*lambda2);
    }
    else
    {
        return -1.;
    }
}



THDMWpositivity2::THDMWpositivity2(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity2::computeThValue()
{
    double lambda1=myTHDMW.getTHDMW_lambda1();
    double lambda2=myTHDMW.getTHDMW_lambda2();
    double lambda3=myTHDMW.getTHDMW_lambda3();
    double lambda4=myTHDMW.getTHDMW_lambda4();
    double lambda5=myTHDMW.getTHDMW_lambda5();

    if (lambda1>0 && lambda2>0)
    {
        return lambda3+lambda4-fabs(lambda5)+sqrt(lambda1*lambda2);
    }
    else
    {
        return -1.;
    }
}
