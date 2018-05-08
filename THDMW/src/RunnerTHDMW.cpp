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

int RGEsMW(double t, const double y[], double beta[], void *flags)
{
    (void)(t); /* avoid unused parameter warning */
    int ord = *(int *)flags;
//    int type = flag-ord;
//    double Yb1=0;
//    double Yb2=0;
//    double Ytau1=0;
//    double Ytau2=0;
    double la1=y[0];
    double nu1=y[1];
    double nu2=y[2];
    double nu3=y[3];
    double nu4=y[4];
    double nu5=y[5];
    double mu1=y[6];
    double mu2=y[7];
    double mu3=y[8];
    double mu4=y[9];
    double mu5=y[10];
    double mu6=y[11];

    double pi=M_PI;

    //beta_la1
    //I have no idea whether this is correct, it's just taken from the custodial THDMW!
    beta[0] = (12.0*la1*la1 + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
    //beta_nu1
    beta[1] = (2.0*nu1*nu1 + nu2*nu2 + 4.0*nu3*nu3 + 4.0*la1*(3.0*nu1+nu2)
               + (7.0*nu4*nu4 - 4.0*nu4*nu5 + 7.0*nu5*nu5)/3.0
               + nu1*(8.0*mu1 + 8.0*mu2 + 17.0*mu3 + 10.0*mu4 + 3.0*mu5 + 5.0*mu6)
               + nu2*(8.0*mu1 + 8.0*mu2 + 24.0*mu3 + 3.0*mu4
                      + 3.0*mu5 + 8.0*mu6)/3.0)/(16.0*pi*pi);
    //beta_nu2
    beta[2] = (2.0*nu2*nu2 + 4.0*nu1*nu2 + 16.0*nu3*nu3 + 4.0*la1*nu2
               + (4.0*nu4*nu4 + 17.0*nu4*nu5 + 4.0*nu5*nu5)/3.0
               + nu2*(8.0*mu1 + 8.0*mu2 + 3.0*mu3 + 24.0*mu4
                      + 3.0*mu5 - mu6)/3.0)/(16.0*pi*pi);
    //beta_nu3
    beta[3] = (2.0*nu3*(2.0*la1 + 2.0*nu1 + 3.0*nu2)
               + (17.0*nu4*nu4 + 16.0*nu4*nu5 + 17.0*nu5*nu5)/12.0
               + nu3*(-mu1 - mu2 + 3.0*mu3 + 3.0*mu4
                      + 24.0*mu5 + 8.0*mu6)/3.0)/(16.0*pi*pi);
    //beta_nu4
    beta[4] = (8.0*nu3*nu4 + 2.0*nu3*nu5
               + nu5*(2.0*nu2 - mu2 + 2.0*mu4 + 4.0*mu5 + mu6)
               + nu4*(3.0*nu1 + 2.0*nu2 + 6.0*mu1 + 2.0*mu2 + 3.0*mu3
                      + 2.0*mu4 + mu5 + mu6))/(16.0*pi*pi);
    //beta_nu5
    beta[5] = (2.0*nu3*nu4 + 8.0*nu3*nu5
               + nu4*(2.0*nu2 - mu1 + 2.0*mu4 + 4.0*mu5 + mu6)
               + nu5*(3.0*nu1 + 2.0*nu2 + 6.0*mu1 + 2.0*mu2 + 3.0*mu3
                      + 2.0*mu4 + mu5 + mu6))/(16.0*pi*pi);
    //beta_mu1
    beta[6] = (3.0*nu4*nu4 + 7.0*mu1*mu1
               + mu1*(6.0*mu2 + 6.0*mu3 + 4.0*mu4 - mu5 - 2.0*mu6)
               + mu2*(4.0*mu4 - mu5)
               - 2.0*mu4*mu6 + 2.0*mu5*mu6 + mu6*mu6)/(16.0*pi*pi);
    //beta_mu2
    beta[7] = (3.0*nu5*nu5 + 7.0*mu2*mu2
               + mu2*(6.0*mu1 + 6.0*mu3 + 4.0*mu4 - mu5 - 2.0*mu6)
               + mu1*(4.0*mu4 - mu5)
               - 2.0*mu4*mu6 + 2.0*mu5*mu6 + mu6*mu6)/(16.0*pi*pi);
    //beta_mu3
    beta[8] = (20.0*mu3*mu3
               + mu3*(288.0*mu1 + 288.0*mu2 + 360.0*mu4 + 108.0*mu5 + 180.0*mu6)/18.0
               + (36.0*nu1*nu1 + 36.0*nu1*nu2 - 24.0*nu4*nu4 - 12.0*nu4*nu5
                  - 24.0*nu5*nu5 + 62.0*mu1*mu1 + 64.0*mu1*mu2 + 62.0*mu2*mu2
                  + (96.0*mu4 + 18.0*mu5 + 58.0*mu6)*(mu1 + mu2)
                  + 54.0*mu4*mu4 + 36.0*mu4*mu5 + 132.0*mu4*mu6 + 18.0*mu5*mu5
                  + 18.0*mu5*mu6 + 29.0*mu6*mu6)/18.0)/(16.0*pi*pi);
    //beta_mu4
    beta[9] = (nu2*nu2 - (nu4*nu4 - 4.0*nu4*nu5 + nu5*nu5)/3.0 + 10.0*mu4*mu4 /*mu4??*/
               + mu5*(mu1 + mu2 + mu6)
               + mu4*(4.0*(4.0*mu1 + 4.0*mu2 + mu6)/3.0 + 2.0*mu5 + 6.0*mu4)
               + 4.0*mu5*mu5
               + (mu1*mu1 + mu2*mu2 - 4.0*mu6*(mu1+mu2) - 2.0*mu6*mu6)/9.0
               + 26.0/9.0*mu1*mu2)/(16.0*pi*pi);
    //beta_mu5
    beta[10] = (4.0*nu3*nu3 - (nu4*nu4 - 4.0*nu4*nu5 + nu5*nu5)/3.0
               + mu5*((mu1 + mu2 + 19.0*mu6)/3.0 + 8.0*mu4 + 6.0*mu3)
               + 2.0*mu4*mu6 + 8.0*mu5*mu5
               + (mu1*mu1 + mu2*mu2 - 4.0*mu6*(mu1+mu2) + 7.0*mu6*mu6)/9.0
               - 10.0/9.0*mu1*mu2)/(16.0*pi*pi);
    //beta_mu6
    beta[11] = (0.5*mu6*mu6 + 3.0*nu4*nu4 + 3.0*nu5*nu5
               - 2.0*(mu1*mu1 + mu2*mu2) + 6.0*mu5*(mu1 + mu2)
               + 7.0*mu6*(mu1 + mu2 + mu3))/(16.0*pi*pi);

        if(ord==1){
            //beta_la1
            beta[0] += 0.0;
            //beta_nu1
            beta[1] += 0.0;
            //beta_nu2
            beta[2] += 0.0;
            //beta_nu3
            beta[3] += 0.0;
            //beta_nu4
            beta[4] += 0.0;
            //beta_nu5
            beta[5] += 0.0;
            //beta_mu1
            beta[6] += 0.0;
            //beta_mu2
            beta[7] += 0.0;
            //beta_mu3
            beta[8] += 0.0;
            //beta_mu4
            beta[9] += 0.0;
            //beta_mu5
            beta[10] += 0.0;
            //beta_mu6
            beta[11] += 0.0;
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

    //perturbativity of the Yukawa couplings
//    for(int i=3;i<6;i++)
//    {
//        if(fabs(InitialValues[i])>sqrt(4.0*M_PI)) check=1;
//    }

    //perturbativity of the quartic Higgs couplings
    for(int i=0;i<15;i++)
    {
        if(fabs(InitialValues[i])>4.0*M_PI) check=1;
    }

    double la1Q = InitialValues[0];
    double la2Q = InitialValues[1];
    double la3Q = InitialValues[2];
    double la4Q = InitialValues[3];
    double mu1Q = InitialValues[4];
    double mu3Q = InitialValues[5];
    double mu4Q = InitialValues[6];
    double nu1Q = InitialValues[7];
    double om1Q = InitialValues[8];
    double ka1Q = InitialValues[9];
    double nu2Q = InitialValues[10];
    double om2Q = InitialValues[11];
    double ka2Q = InitialValues[12];
    double nu4Q = InitialValues[13];
    double om4Q = InitialValues[14];

    //positivity checks
    double muAtimes2=4.0*mu1Q+2.0*mu3Q+4.0*mu4Q;
    if(la1Q<0.0) check=1;
    if(la2Q<0.0) check=1;
    if(muAtimes2<0.0) check=1;
    if(5.0*mu1Q+3.0*mu3Q+3.0*mu4Q-fabs(mu1Q)<0.0) check=1;
    if(sqrt(la1Q*muAtimes2)+nu1Q<0.0) check=1;
    if(sqrt(la1Q*muAtimes2)+nu1Q+2.0*nu2Q<0.0) check=1;
    if(la1Q+0.25*muAtimes2+nu1Q+2.0*nu2Q-2.0*fabs(nu4Q)/sqrt(3.0)<0.0) check=1;
    if(la3Q+sqrt(la1Q*la2Q)<0.0) check=1;
    if(la3Q+2.0*la4Q+sqrt(la1Q*la2Q)<0.0) check=1;
    if(sqrt(la2Q*muAtimes2)+om1Q<0.0) check=1;
    if(sqrt(la2Q*muAtimes2)+om1Q+2.0*om2Q<0.0) check=1;
    if(la2Q+0.25*muAtimes2+om1Q+2.0*om2Q-2.0*fabs(om4Q)/sqrt(3.0)<0.0) check=1;

    /////////////////
    //NLO unitarity//
    /////////////////

    double pi=M_PI;
    gslpp::matrix<gslpp::complex> Smatrix1(4,4,0.), Smatrix2(4,4,0.);
    gslpp::matrix<gslpp::complex> Seigenvectors1(4,4,0.), Seigenvectors2(4,4,0.);
    gslpp::vector<double> Seigenvalues1(4,0.), Seigenvalues2(4,0.);
    gslpp::vector<gslpp::complex> unitarityeigenvalues(11,0.);

    if(t1>tNLOuni)
    {

    //LO part
    Smatrix1.assign(0,0, 3.0*la1Q/(16.0*pi));
    Smatrix1.assign(0,1, (2.0*la3Q+la4Q)/(16.0*pi));
    Smatrix1.assign(1,0, Smatrix1(0,1));
    Smatrix1.assign(0,3, (2.0*nu1Q+nu2Q)/(8.0*sqrt(2.0)*pi));
    Smatrix1.assign(3,0, Smatrix1(0,3));
    Smatrix1.assign(1,1, 3.0*la2Q/(16.0*pi));
    Smatrix1.assign(1,3, (2.0*om1Q+om2Q)/(8.0*sqrt(2.0)*pi));
    Smatrix1.assign(3,1, Smatrix1(1,3));
    Smatrix1.assign(2,2, (la3Q+5.0*la4Q)/(16.0*pi));
    Smatrix1.assign(2,3, (4.0*ka1Q+2.0*ka2Q)/(16.0*pi));
    Smatrix1.assign(3,2, Smatrix1(2,3));
    Smatrix1.assign(3,3, (26.0*mu1Q+17.0*mu3Q+13.0*mu4Q)/(32.0*pi));

    Smatrix2.assign(0,0, la1Q/(16.0*pi));
    Smatrix2.assign(0,1, la4Q/(16.0*pi));
    Smatrix2.assign(1,0, Smatrix2(0,1));
    Smatrix2.assign(0,3, nu2Q/(8.0*sqrt(2.0)*pi));
    Smatrix2.assign(3,0, Smatrix2(0,3));
    Smatrix2.assign(1,1, la2Q/(16.0*pi));
    Smatrix2.assign(1,3, om2Q/(8.0*sqrt(2.0)*pi));
    Smatrix2.assign(3,1, Smatrix2(1,3));
    Smatrix2.assign(2,2, (la3Q+la4Q)/(16.0*pi));
    Smatrix2.assign(2,3, ka2Q/(8.0*pi));
    Smatrix2.assign(3,2, Smatrix2(2,3));
    Smatrix2.assign(3,3, (14.0*mu1Q+3.0*mu3Q+27.0*mu4Q)/(96.0*pi));

    Smatrix1.eigensystem(Seigenvectors1, Seigenvalues1);
    Smatrix2.eigensystem(Seigenvectors2, Seigenvalues2);

    for (int i=0; i < 4; i++) {
        unitarityeigenvalues.assign(i, Seigenvalues1(i));
        unitarityeigenvalues.assign(4+i, Seigenvalues2(i));
    }
    unitarityeigenvalues.assign(8, (la3Q-la4Q)/(16.0*pi));
    unitarityeigenvalues.assign(9, sqrt(15.0)*nu4Q/(16.0*pi));
    unitarityeigenvalues.assign(10, sqrt(15.0)*om4Q/(16.0*pi));

    //beta_la1*16pi^2
    double betala1 = 12.0*la1Q*la1Q + 4.0*la3Q*la3Q + 4.0*la3Q*la4Q + 4.0*la4Q*la4Q
                     + 8.0*nu1Q*nu1Q + 8.0*nu1Q*nu2Q + 8.0*nu2Q*nu2Q;
    //beta_la2*16pi^2
    double betala2 = 12.0*la2Q*la2Q + 4.0*la3Q*la3Q + 4.0*la3Q*la4Q + 4.0*la4Q*la4Q
                     + 8.0*om1Q*om1Q + 8.0*om1Q*om2Q + 8.0*om2Q*om2Q;
    //beta_la3*16pi^2
    double betala3 = 4.0*la3Q*la3Q + 4.0*la4Q*la4Q + (la1Q+la2Q)*(6.0*la3Q+2.0*la4Q)
                     + 8.0*ka2Q*ka2Q + 8.0*nu1Q*om1Q + 4.0*nu2Q*om1Q + 4.0*nu1Q*om2Q;
    //beta_la4*16pi^2
    double betala4 = (la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 6.0*la4Q*la4Q
                      + 4.0*ka1Q*ka1Q + 4.0*ka1Q*ka2Q + 2.0*ka2Q*ka2Q + 2.0*nu2Q*om2Q)*2.0;
    //beta_mu1*16pi^2
    double betamu1 = 11.0*mu1Q*mu1Q + 3.0*mu1Q*mu4Q + mu1Q*(2.0*mu1Q+6.0*mu3Q+3.0*mu4Q)
                     + 3.0*nu4Q*nu4Q + 3.0*om4Q*om4Q;
    //beta_mu3*16pi^2
    double betamu3 = (18.0*ka1Q*ka1Q + 18.0*ka1Q*ka2Q + 134.0*mu1Q*mu1Q + 6.0*mu1Q*(39.0*mu3Q + 22.0*mu4Q)
                      + 3.0*(30.0*mu3Q*mu3Q + 39.0*mu3Q*mu4Q + 9.0*mu4Q*mu4Q
                             + 3.0*nu1Q*nu1Q + 3.0*nu1Q*nu2Q - 5.0*nu4Q*nu4Q
                             + 3.0*om1Q*om1Q + 3.0*om1Q*om2Q - 5.0*om4Q*om4Q))/4.5;
    //beta_mu4*16pi^2
    double betamu4 = (18.0*ka2Q*ka2Q + 4.0*mu1Q*mu1Q + 156.0*mu1Q*mu4Q + 54.0*mu3Q*mu4Q + 144.0*mu4Q*mu4Q
                      + 9.0*nu2Q*nu2Q + 6.0*nu4Q*nu4Q + 9.0*om2Q*om2Q + 6.0*om4Q*om4Q)/9.0;
    //beta_nu1*16pi^2
    double betanu1 = (6.0*ka1Q*ka1Q + 6.0*ka2Q*ka2Q + 18.0*la1Q*nu1Q
                      + 78.0*mu1Q*nu1Q + 51.0*mu3Q*nu1Q + 39.0*mu4Q*nu1Q + 6.0*nu1Q*nu1Q
                      + 6.0*la1Q*nu2Q + 32.0*mu1Q*nu2Q + 24.0*mu3Q*nu2Q + 6.0*mu4Q*nu2Q
                      + 6.0*nu2Q*nu2Q + 10.0*nu4Q*nu4Q
                      + 12.0*la3Q*om1Q + 6.0*la4Q*om1Q + 6.0*la3Q*om2Q)/3.0;
    //beta_om1*16pi^2
    double betaom1 = (6.0*ka1Q*ka1Q + 6.0*ka2Q*ka2Q
                      + 12.0*la3Q*nu1Q + 6.0*la4Q*nu1Q + 6.0*la3Q*nu2Q
                      + 18.0*la2Q*om1Q + 78.0*mu1Q*om1Q + 51.0*mu3Q*om1Q + 39.0*mu4Q*om1Q + 6.0*om1Q*om1Q
                      + 6.0*la2Q*om2Q + 32.0*mu1Q*om2Q + 24.0*mu3Q*om2Q + 6.0*mu4Q*om2Q + 6.0*om2Q*om2Q
                      + 10.0*om4Q*om4Q)/3.0;
    //beta_ka1*16pi^2
    double betaka1 = (6.0*ka1Q*(2.0*la3Q + 10.0*la4Q + 18.0*mu1Q + 17.0*mu3Q + 13.0*mu4Q + 2.0*nu1Q + 2.0*om1Q)
                      + ka2Q*(24.0*la4Q + 64.0*mu1Q + 48.0*mu3Q + 24.0*mu4Q + 9.0*nu2Q + 9.0*om2Q)
                      + 20.0*nu4Q*om4Q)/6.0;
    //beta_nu2*16pi^2
    double betanu2 = 4.0*ka1Q*ka2Q + 6.0*ka2Q*ka2Q + 2.0*la1Q*nu2Q + ((14.0*mu1Q)/3.0 + mu3Q + 9.0*mu4Q)*nu2Q
                     + 4.0*nu1Q*nu2Q + 6.0*nu2Q*nu2Q + (25.0*nu4Q*nu4Q)/3.0 + 2.0*la4Q*om2Q;
    //beta_om2*16pi^2
    double betaom2 = 4.0*ka1Q*ka2Q + 6.0*ka2Q*ka2Q + 2.0*la4Q*nu2Q + 2.0*la2Q*om2Q
                     + ((14.0*mu1Q)/3.0 + mu3Q + 9.0*mu4Q)*om2Q + 4.0*om1Q*om2Q + 6.0*om2Q*om2Q
                     + (25.0*om4Q*om4Q)/3.0;
    //beta_ka2*16pi^2
    double betaka2 = (ka2Q*(6.0*la3Q + 6.0*la4Q + 14.0*mu1Q + 3.0*mu3Q + 27.0*mu4Q
                           + 6.0*nu1Q + 12.0*nu2Q + 6.0*om1Q + 12.0*om2Q)
                      + 6.0*ka1Q*(nu2Q + om2Q) + 42.0*nu4Q*om4Q)/3.0;
    //beta_nu4*16pi^2
    double betanu4 = 11.0*mu1Q*nu4Q + 3.0*mu3Q*nu4Q + 9.0*mu4Q*nu4Q + 3.0*nu1Q*nu4Q + 9.0*nu2Q*nu4Q
                     + 3.0*ka1Q*om4Q + 9.0*ka2Q*om4Q;
    //beta_om4*16pi^2
    double betaom4 = 3.0*ka1Q*nu4Q + 9.0*ka2Q*nu4Q
                     + (11.0*mu1Q + 3.0*(mu3Q + 3.0*mu4Q + om1Q + 3.0*om2Q))*om4Q;

//    diagonalization
    gslpp::matrix<gslpp::complex> Sbmatrix1(4,4,0.), Sbmatrix2(4,4,0.);
    gslpp::matrix<gslpp::complex> Seigenvectors1T(4,4,0.), Seigenvectors2T(4,4,0.);
    gslpp::vector<gslpp::complex> Sbeigenvalues1(4,0.), Sbeigenvalues2(4,0.);
    gslpp::vector<gslpp::complex> betaeigenvalues(11,0.);
    gslpp::vector<gslpp::complex> NLOunitarityeigenvalues(11,0.);

    Sbmatrix1.assign(0,0, 3.0*betala1/(16.0*pi));
    Sbmatrix1.assign(0,1, (2.0*betala3+betala4)/(16.0*pi));
    Sbmatrix1.assign(1,0, Sbmatrix1(0,1));
    Sbmatrix1.assign(0,3, (2.0*betanu1+betanu2)/(8.0*sqrt(2.0)*pi));
    Sbmatrix1.assign(3,0, Sbmatrix1(0,3));
    Sbmatrix1.assign(1,1, 3.0*betala2/(16.0*pi));
    Sbmatrix1.assign(1,3, (2.0*betaom1+betaom2)/(8.0*sqrt(2.0)*pi));
    Sbmatrix1.assign(3,1, Sbmatrix1(1,3));
    Sbmatrix1.assign(2,2, (betala3+5.0*betala4)/(16.0*pi));
    Sbmatrix1.assign(2,3, (4.0*betaka1+2.0*betaka2)/(16.0*pi));
    Sbmatrix1.assign(3,2, Sbmatrix1(2,3));
    Sbmatrix1.assign(3,3, (26.0*betamu1+17.0*betamu3+13.0*betamu4)/(32.0*pi));

    Sbmatrix2.assign(0,0, betala1/(16.0*pi));
    Sbmatrix2.assign(0,1, betala4/(16.0*pi));
    Sbmatrix2.assign(1,0, Sbmatrix2(0,1));
    Sbmatrix2.assign(0,3, betanu2/(8.0*sqrt(2.0)*pi));
    Sbmatrix2.assign(3,0, Sbmatrix2(0,3));
    Sbmatrix2.assign(1,1, betala2/(16.0*pi));
    Sbmatrix2.assign(1,3, betaom2/(8.0*sqrt(2.0)*pi));
    Sbmatrix2.assign(3,1, Sbmatrix2(1,3));
    Sbmatrix2.assign(2,2, (betala3+betala4)/(16.0*pi));
    Sbmatrix2.assign(2,3, betaka2/(8.0*pi));
    Sbmatrix2.assign(3,2, Sbmatrix2(2,3));
    Sbmatrix2.assign(3,3, (14.0*betamu1+3.0*betamu3+27.0*betamu4)/(96.0*pi));

    Seigenvectors1T=Seigenvectors1.hconjugate();
    Seigenvectors2T=Seigenvectors2.hconjugate();

    for (int i=0; i < 4; i++) {
        for (int k=0; k < 4; k++) {
            for (int l=0; l < 4; l++) {
                Sbeigenvalues1.assign(i, Sbeigenvalues1(i) + Seigenvectors1T(i,k) * Sbmatrix1(k,l) * Seigenvectors1(l,i) );
                Sbeigenvalues2.assign(i, Sbeigenvalues2(i) + Seigenvectors2T(i,k) * Sbmatrix2(k,l) * Seigenvectors2(l,i) );
            }                
        }
        betaeigenvalues.assign(i, -1.5 * Sbeigenvalues1(i));
        betaeigenvalues.assign(i+4, -1.5 * Sbeigenvalues2(i));
    }

    betaeigenvalues.assign(8, -1.5 * (betala3-betala4)/(16.0*pi));
    betaeigenvalues.assign(9, -1.5 * sqrt(15.0)*betanu4/(16.0*pi));
    betaeigenvalues.assign(10, -1.5 * sqrt(15.0)*betaom4/(16.0*pi));

    for (int i=0; i < 11; i++) {
        NLOunitarityeigenvalues.assign(i, -(gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i) + betaeigenvalues(i) );
        if( ( unitarityeigenvalues(i) + NLOunitarityeigenvalues(i).real() ).abs() > 0.5) check=1;
        if( (unitarityeigenvalues(i)).abs() > Rpeps && (NLOunitarityeigenvalues(i)/unitarityeigenvalues(i)).abs() > 1.0) check=1;
    }

    } //end of the if(t1>tNLOuni)
    return check;
}

int RGEcheckMW(const double InitialValues[], const double t1, const double Rpeps, const double tNLOuni)
{
    int check=0;

    //perturbativity of the quartic Higgs couplings
    for(int i=0;i<12;i++)
    {
        if(fabs(InitialValues[i])>4.0*M_PI) check=1;
    }

    double la1Q = InitialValues[0];
    double nu1Q = InitialValues[1];
    double nu2Q = InitialValues[2];
    double nu3Q = InitialValues[3];
    double nu4Q = InitialValues[4];
    double nu5Q = InitialValues[5];
    double mu1Q = InitialValues[6];
    double mu2Q = InitialValues[7];
    double mu3Q = InitialValues[8];
    double mu4Q = InitialValues[9];
    double mu5Q = InitialValues[10];
    double mu6Q = InitialValues[11];

    //positivity checks
    double muAtimes2=mu1Q+mu2Q+mu6Q+2.0*(mu3Q+mu4Q+mu5Q);
    if(la1Q<0.0) check=1;
    if(muAtimes2<0.0) check=1;
    if(mu1Q+mu2Q+mu3Q+mu4Q<0.0) check=1;
    if(14.0*(mu1Q+mu2Q)+5.0*mu6Q+24.0*(mu3Q+mu4Q)-3.0*fabs(2.0*(mu1Q+mu2Q)-mu6Q)<0.0) check=1;
    if(5.0*(mu1Q+mu2Q+mu6Q)+6.0*(2.0*mu3Q+mu4Q+mu5Q)-fabs(mu1Q+mu2Q+mu6Q)<0.0) check=1;
    if(sqrt(la1Q*muAtimes2)+nu1Q<0.0) check=1;
    if(sqrt(la1Q*muAtimes2)+nu1Q+nu2Q-2.0*fabs(nu3Q)<0.0) check=1;
    if(la1Q+0.25*muAtimes2+nu1Q+nu2Q+2.0*nu3Q-fabs(nu4Q+nu5Q)/sqrt(3.0)<0.0) check=1;

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

double RunnerTHDMW::RGERunnerMW(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double NLOuniscale)
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
    gsl_odeiv2_system RGEsystem = {RGEsMW, JacobianTHDMW, NumberOfRGEs, &order};

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
        if(RGEcheckMW(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}
