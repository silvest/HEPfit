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



THDMWpositivity3::THDMWpositivity3(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity3::computeThValue()
{
    double mu1=myTHDMW.getTHDMW_mu1();
    double mu2=myTHDMW.getTHDMW_mu2();
    double mu3=myTHDMW.getTHDMW_mu3();
    double mu4=myTHDMW.getTHDMW_mu4();
    double mu5=myTHDMW.getTHDMW_mu5();
    double mu6=myTHDMW.getTHDMW_mu6();

    return mu1+2.0*mu3+2.0*mu5+mu2+2.0*mu4+mu6;
}



THDMWpositivity4::THDMWpositivity4(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity4::computeThValue()
{
    double mu1=myTHDMW.getTHDMW_mu1();
    double mu2=myTHDMW.getTHDMW_mu2();
    double mu3=myTHDMW.getTHDMW_mu3();
    double mu4=myTHDMW.getTHDMW_mu4();
    double mu5=myTHDMW.getTHDMW_mu5();
    double mu6=myTHDMW.getTHDMW_mu6();

    return mu1+2.0*mu3+2.0*mu5+mu2+mu6;
}



THDMWpositivity5::THDMWpositivity5(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity5::computeThValue()
{
    double mu1=myTHDMW.getTHDMW_mu1();
    double mu2=myTHDMW.getTHDMW_mu2();
    double mu3=myTHDMW.getTHDMW_mu3();
    double mu4=myTHDMW.getTHDMW_mu4();
    double mu5=myTHDMW.getTHDMW_mu5();
    double mu6=myTHDMW.getTHDMW_mu6();

    return 3.0*(mu1+2.0*mu3+2.0*mu5)-mu2+mu6;
}



THDMWpositivity6::THDMWpositivity6(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity6::computeThValue()
{
    double lambda1=myTHDMW.getTHDMW_lambda1();
    double mu1=myTHDMW.getTHDMW_mu1();
    double mu2=myTHDMW.getTHDMW_mu2();
    double mu3=myTHDMW.getTHDMW_mu3();
    double mu4=myTHDMW.getTHDMW_mu4();
    double mu5=myTHDMW.getTHDMW_mu5();
    double mu6=myTHDMW.getTHDMW_mu6();
    double mu_sum=mu1+2.0*mu3+2.0*mu5+mu2+2.0*mu4+mu6;
    double nu1=myTHDMW.getTHDMW_nu1();
    double nu2=myTHDMW.getTHDMW_nu2();
    double nu3=myTHDMW.getTHDMW_nu3();
    if (lambda1>0 && mu_sum>0)
    {
        return nu1-fabs(nu2)-fabs(nu3)+sqrt(lambda1*mu_sum);
    }
    else
    {
        return -1.;
    }
}



THDMWpositivity7::THDMWpositivity7(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{}

double THDMWpositivity7::computeThValue()
{
    double lambda2=myTHDMW.getTHDMW_lambda2();
    double mu1=myTHDMW.getTHDMW_mu1();
    double mu2=myTHDMW.getTHDMW_mu2();
    double mu3=myTHDMW.getTHDMW_mu3();
    double mu4=myTHDMW.getTHDMW_mu4();
    double mu5=myTHDMW.getTHDMW_mu5();
    double mu6=myTHDMW.getTHDMW_mu6();
    double mu_sum=mu1+2.0*mu3+2.0*mu5+mu2+2.0*mu4+mu6;
    double omega1=myTHDMW.getTHDMW_omega1();
    double omega2=myTHDMW.getTHDMW_omega2();
    double omega3=myTHDMW.getTHDMW_omega3();
    if (lambda2>0 && mu_sum>0)
    {
        return omega1-fabs(omega2)-fabs(omega3)+sqrt(lambda2*mu_sum);
    }
    else
    {
        return -1.;
    }
}
