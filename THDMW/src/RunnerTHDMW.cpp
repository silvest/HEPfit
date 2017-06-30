/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RunnerTHDMW.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

RunnerTHDMW::RunnerTHDMW(const StandardModel& SM_i) : myTHDMW(static_cast<const THDMW*> (&SM_i))
{}

RunnerTHDMW::~RunnerTHDMW()
{};

int RGEsTHDMW(double t, const double y[], double beta[], void *flags)
{
    (void)(t); /* avoid unused parameter warning */
    int ord = *(int *)flags;
//    int type = flag-ord;
//    double Yb1=0;
//    double Yb2=0;
//    double Ytau1=0;
//    double Ytau2=0;
    double la1=y[0];
    double la2=y[1];
    double la3=y[2];
    double la4=y[3];
    double mu1=y[4];
    double mu3=y[5];
    double mu4=y[6];
    double nu1=y[7];
    double om1=y[8];
    double ka1=y[9];
    double nu2=y[10];
    double om2=y[11];
    double ka2=y[12];
    double nu4=y[13];
    double om4=y[14];

    double pi=M_PI;

    //beta_la1
    beta[0] = (12.0*la1*la1 + 4.0*la3*la3 + 4.0*la3*la4 + 4.0*la4*la4 
               + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
    //beta_la2
    beta[1] = (12.0*la2*la2 + 4.0*la3*la3 + 4.0*la3*la4 + 4.0*la4*la4
               + 8.0*om1*om1 + 8.0*om1*om2 + 8.0*om2*om2)/(16.0*pi*pi);
    //beta_la3
    beta[2] = (4.0*la3*la3 + 4.0*la4*la4 + (la1+la2)*(6.0*la3+2.0*la4) 
               + 8.0*ka2*ka2 + 8.0*nu1*om1 + 4.0*nu2*om1 + 4.0*nu1*om2)/(16.0*pi*pi);
    //beta_la4
    beta[3] = (la1*la4 + la2*la4 + 4.0*la3*la4 + 6.0*la4*la4
               + 4.0*ka1*ka1 + 4.0*ka1*ka2 + 2.0*ka2*ka2 + 2.0*nu2*om2)/(8.0*pi*pi);
    //beta_mu1
    beta[4] = (11.0*mu1*mu1 + 3.0*mu1*mu4 + mu1*(2.0*mu1+6.0*mu3+3.0*mu4)
               + 3.0*nu4*nu4 + 3.0*om4*om4)/(16.0*pi*pi);
    //beta_mu3
    beta[5] = (18.0*ka1*ka1 + 18.0*ka1*ka2 + 134.0*mu1*mu1 + 6.0*mu1*(39.0*mu3 + 22.0*mu4)
               + 3.0*(30.0*mu3*mu3 + 39.0*mu3*mu4 + 9.0*mu4*mu4 
                      + 3.0*nu1*nu1 + 3.0*nu1*nu2 - 5.0*nu4*nu4
                      + 3.0*om1*om1 + 3.0*om1*om2 - 5.0*om4*om4))/(72.0*pi*pi);
    //beta_mu4
    beta[6] = (18.0*ka2*ka2 + 4.0*mu1*mu1 + 156.0*mu1*mu4 + 54.0*mu3*mu4 + 144.0*mu4*mu4
               + 9.0*nu2*nu2 + 6.0*nu4*nu4 + 9.0*om2*om2 + 6.0*om4*om4)/(144.0*pi*pi);
    //beta_nu1
    beta[7] = (6.0*ka1*ka1 + 6.0*ka2*ka2 + 18.0*la1*nu1
               + 78.0*mu1*nu1 + 51.0*mu3*nu1 + 39.0*mu4*nu1 + 6.0*nu1*nu1
               + 6.0*la1*nu2 + 32.0*mu1*nu2 + 24.0*mu3*nu2 + 6.0*mu4*nu2
               + 6.0*nu2*nu2 + 10.0*nu4*nu4
               + 12.0*la3*om1 + 6.0*la4*om1 + 6.0*la3*om2)/(48.0*pi*pi);
    //beta_om1
    beta[8] = (6.0*ka1*ka1 + 6.0*ka2*ka2 
               + 12.0*la3*nu1 + 6.0*la4*nu1 + 6.0*la3*nu2
               + 18.0*la2*om1 + 78.0*mu1*om1 + 51.0*mu3*om1 + 39.0*mu4*om1 + 6.0*om1*om1
               + 6.0*la2*om2 + 32.0*mu1*om2 + 24.0*mu3*om2 + 6.0*mu4*om2 + 6.0*om2*om2
               + 10.0*om4*om4)/(48.0*pi*pi);
    //beta_ka1
    beta[9] = (6.0*ka1*(2.0*la3 + 10.0*la4 + 18.0*mu1 + 17.0*mu3 + 13.0*mu4 + 2.0*nu1 + 2.0*om1)
               + ka2*(24.0*la4 + 64.0*mu1 + 48.0*mu3 + 24.0*mu4 + 9.0*nu2 + 9.0*om2)
               + 20.0*nu4*om4)/(96.0*pi*pi);
    //beta_nu2
    beta[10] = (4.0*ka1*ka2 + 6.0*ka2*ka2 + 2.0*la1*nu2 + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*nu2 
                + 4.0*nu1*nu2 + 6.0*nu2*nu2 + (25.0*nu4*nu4)/3.0 + 2.0*la4*om2)/(16.0*pi*pi);
    //beta_om2
    beta[11] = (4.0*ka1*ka2 + 6.0*ka2*ka2 + 2.0*la4*nu2 + 2.0*la2*om2 
                + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*om2 + 4.0*om1*om2 + 6.0*om2*om2 
                + (25.0*om4*om4)/3.0)/(16.0*pi*pi);
    //beta_ka2
    beta[12] = (ka2*(6.0*la3 + 6.0*la4 + 14.0*mu1 + 3.0*mu3 + 27.0*mu4
                     + 6.0*nu1 + 12.0*nu2 + 6.0*om1 + 12.0*om2)
                + 6.0*ka1*(nu2 + om2) + 42.0*nu4*om4)/(48.0*pi*pi);
    //beta_nu4
    beta[13] = (11.0*mu1*nu4 + 3.0*mu3*nu4 + 9.0*mu4*nu4 + 3.0*nu1*nu4 + 9.0*nu2*nu4 
                + 3.0*ka1*om4 + 9.0*ka2*om4)/(16.0*pi*pi);
    //beta_om4
    beta[14] = (3.0*ka1*nu4 + 9.0*ka2*nu4 
                + (11.0*mu1 + 3.0*(mu3 + 3.0*mu4 + om1 + 3.0*om2))*om4)/(16.0*pi*pi);

        if(ord==1){
            //beta_la1
            beta[0] += 0.0;
            //beta_la2
            beta[1] += 0.0;
            //beta_la3
            beta[2] += 0.0;
            //beta_la4
            beta[3] += 0.0;
            //beta_mu1
            beta[4] += 0.0;
            //beta_mu3
            beta[5] += 0.0;
            //beta_mu4
            beta[6] += 0.0;
            //beta_nu1
            beta[7] += 0.0;
            //beta_om1
            beta[8] += 0.0;
            //beta_ka1
            beta[9] += 0.0;
            //beta_nu2
            beta[10] += 0.0;
            //beta_om2
            beta[11] += 0.0;
            //beta_ka2
            beta[12] += 0.0;
            //beta_nu4
            beta[13] += 0.0;
            //beta_om4
            beta[14] += 0.0;
        }

    return 0;
}

int JacobianTHDMW (double t, const double y[], double *dfdy, double dfdt[], void *order)
{
    return 0;
}

int RGEcheckTHDMW(const double InitialValues[], const double t1, const double Rpeps, const double tNLOuni)
{
    int check=0;

    //positivity checks (empty at the moment)

    //NLO unitarity
//    double la1Q = InitialValues[9];
//    double la2Q = InitialValues[10];
//    double la3Q = InitialValues[11];
//    double la4Q = InitialValues[12];
//    double la5Q = InitialValues[13];

//    double betalambda1 = 
//    ...

//    double         uniLO1   = 
//    gslpp::complex uniNLO1  = 
//    ...

//    diagonalization

//    if( (uniNLOev1 - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
//    ...
//    if( fabs(uniLOev1) > Rpeps && (uniNLOev1/uniLOev1-1.0).abs() > 1.0) check=1;
//    ...

    return check;
}

double RunnerTHDMW::RGERunnerTHDMW(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double NLOuniscale)
{
    //Define which stepping function should be used
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;

    //Allocate space for the stepping function
    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc(T, NumberOfRGEs);

    //Define the absolute (A) and relative (R) error on y at each step.
    //The real error will be compared to the following error estimate:
    //  A + R * |y_i|
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new(1e-6, 0.0);

    //Allocate space for the evolutor
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc(NumberOfRGEs);

    //Definition of the RGE system (the Jacobian is not necessary for the RK4 method; it's an empty function here)
    gsl_odeiv2_system RGEsystem = {RGEsTHDMW, JacobianTHDMW, NumberOfRGEs, &order};

    //Set starting and end point as natural logarithmic scales (conversion from decadic log scale)
    double t1 = Q1*log(10.0);
    double t2 = Q2*log(10.0);
    double tNLOuni = NLOuniscale*log(10.0);

    //Set initial step size
    double InitialStepSize = 1e-6;

    //Run!
    while (t1 < t2)
    {
        int status = gsl_odeiv2_evolve_apply (e, c, s, &RGEsystem, &t1, t2, &InitialStepSize, InitialValues);
        if(status != GSL_SUCCESS) break;

        //intermediate checks if appropriate
        if(RGEcheckTHDMW(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}
