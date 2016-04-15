/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Runner.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

Runner::Runner(const StandardModel& SM_i) : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{
    mym11_2=new m11_2(SM_i);
    mym22_2=new m22_2(SM_i);
    mylambda1=new lambda1(SM_i);
    mylambda2=new lambda2(SM_i);
    mylambda3=new lambda3(SM_i);
    mylambda4=new lambda4(SM_i);
    mylambda5=new lambda5(SM_i);
}

Runner::~Runner()
{
  delete mym11_2;
  delete mym22_2;
  delete mylambda1;
  delete mylambda2;
  delete mylambda3;
  delete mylambda4;
  delete mylambda5;
};

double Runner::computeThValue()
{
    return 0.;
}

int RGEs(double t, const double y[], double beta[], void *params)
{
    (void)(t); /* avoid unused parameter warning */

    double g1=y[0];
    double g2=y[1];
    double g3=y[2];
    double Yt=y[3];
    double Yb=y[4];
    double Ytau=y[5];
    double m11_2=y[6];
    double m22_2=y[7];
    double m12_2=y[8];
    double la1=y[9];
    double la2=y[10];
    double la3=y[11];
    double la4=y[12];
    double la5=y[13];

    double pi=M_PI;

    //beta_g1
    beta[0] = 7.0*g1*g1*g1/(16.0*pi*pi);
    //beta_g2
    beta[1] = -3.0*g2*g2*g2/(16.0*pi*pi);
    //beta_g3
    beta[2] = -7.0*g3*g3*g3/(16.0*pi*pi);
    //beta_Yt
    beta[3] = (Yt*(-17.0*g1*g1 - 27.0*g2*g2 + 6.0*(-16.0*g3*g3 + Yb*Yb + 9.0*Yt*Yt)))/(192.0*pi*pi);
    //beta_Yb
    beta[4] = (Yb*(-5.0*g1*g1 - 27.0*g2*g2 + 6.0*(-16*g3*g3 + 9.0*Yb*Yb + Yt*Yt + 2.0*Ytau*Ytau)))/(192.0*pi*pi);
    //beta_Ytau
    beta[5] = (Ytau*(-15.0*g1*g1 - 9.0*g2*g2 + 12.0*Yb*Yb + 10.0*Ytau*Ytau))/(64.0*pi*pi);
    //beta_m11_2
    beta[6] = (-3.0*g1*g1*m11_2 - 9.0*g2*g2*m11_2 + 4.0*((2.0*la3 + la4)*m22_2 + m11_2*(3.0*la1 + 3.0*Yb*Yb + Ytau*Ytau)))/(32.0*pi*pi);
    //beta_m22_2
    beta[7] = (-3.0*g1*g1*m22_2 - 9.0*g2*g2*m22_2 + 4.0*(2.0*la3*m11_2 + la4*m11_2 + 3.0*la2*m22_2 + 3.0*m22_2*Yt*Yt))/(32.0*pi*pi);
    //beta_m12_2
    beta[8] = (m12_2*(-3.0*g1*g1 - 9.0*g2*g2 + 2.0*(2.0*la3 + 4.0*la4 + 6.0*la5 + 3.0*Yb*Yb + 3.0*Yt*Yt + Ytau*Ytau)))/(32.0*pi*pi);
    //beta_lambda_1
    beta[9] = (3.0*g1*g1*g1*g1 + 9.0*g2*g2*g2*g2 + 6.0*g1*g1*(g2*g2 - 2.0*la1) - 36.0*g2*g2*la1 + 8.0*(6.0*la1*la1 + 2.0*la3*la3 + 2.0*la3*la4 + la4*la4 + la5*la5 + 6.0*la1*Yb*Yb - 6.0*Yb*Yb*Yb*Yb + 2.0*la1*Ytau*Ytau - 2.0*Ytau*Ytau*Ytau*Ytau))/(64.0*pi*pi);
    //beta_lambda_2
    beta[10] = (3.0*g1*g1*g1*g1 + 9.0*g2*g2*g2*g2 + 6.0*g1*g1*(g2*g2 - 2.0*la2) - 36.0*g2*g2*la2 + 8.0*(6.0*la2*la2 + 2.0*la3*la3 + 2.0*la3*la4 + la4*la4 + la5*la5 + 6.0*la2*Yt*Yt - 6.0*Yt*Yt*Yt*Yt))/(64.0*pi*pi);
    //beta_lambda_3
    beta[11] = (3.0*g1*g1*g1*g1 + 9.0*g2*g2*g2*g2 - 36.0*g2*g2*la3 - 6.0*g1*g1*(g2*g2 + 2.0*la3) + 8.0*(3.0*la1*la3 + 3.0*la2*la3 + 2.0*la3*la3 + la1*la4 + la2*la4 + la4*la4 + la5*la5 + 3.0*la3*Yt*Yt + Yb*Yb*(3.0*la3 - 6.0*Yt*Yt) + la3*Ytau*Ytau))/(64.0*pi*pi);
    //beta_lambda_4
    beta[12] = (3.0*g1*g1*(g2*g2 - la4) - 9.0*g2*g2*la4 + 2.0*la1*la4 + 2.0*la2*la4 + 8.0*la3*la4 + 4.0*la4*la4 + 8.0*la5*la5 + 6.0*la4*Yt*Yt + 6.0*Yb*Yb*(la4 + 2.0*Yt*Yt) + 2.0*la4*Ytau*Ytau)/(16.0*pi*pi);
    //beta_lambda_5
    beta[13] = (la5*(-3.0*g1*g1 - 9.0*g2*g2 + 2.0*(la1 + la2 + 4.0*la3 + 6.0*la4 + 3.0*Yb*Yb + 3.0*Yt*Yt + Ytau*Ytau)))/(16.0*pi*pi);

    return 0;
}

int Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    return 0;
}

int RGEcheck(const double InitialValues[], const double t1)
{
    int check=0;

    //perturbativity of the Yukawa couplings
    for(int i=3;i<6;i++)
    {
        if(fabs(InitialValues[i])>sqrt(4.0*M_PI)) check=1;
    }
    //perturbativity of the quartic Higgs couplings
    for(int i=9;i<14;i++)
    {
        if(fabs(InitialValues[i])>4.0*M_PI) check=1;
    }
    //positivity
    if(InitialValues[9]<0.0) check=1;
    if(InitialValues[10]<0.0) check=1;
    if(InitialValues[11]+sqrt(fabs(InitialValues[9]*InitialValues[10]))<0.0) check=1;
    if(InitialValues[11]+InitialValues[12]-fabs(InitialValues[13])+sqrt(fabs(InitialValues[9]*InitialValues[10]))<0.0) check=1;
    //NLO unitarity
    double YtQ = InitialValues[3];
    double Yb1Q = 0.;
    double Yb2Q = 0.;
    double Ytau1Q = 0.;
    double Ytau2Q = 0.;
//    if( modelflag == "type1" ) {
//        Yb2Q=InitialValues[4];
//        Ytau2Q=InitialValues[5];
//    }
//    else if( modelflag == "type2" ) {
        Yb1Q=InitialValues[4];
        Ytau1Q=InitialValues[5];
//    }
//    else if( modelflag == "typeX" ) {
//        Yb2Q=InitialValues[4];
//        Ytau1Q=InitialValues[5];
//    }
//    else if( modelflag == "typeY" ) {
//        Yb1Q=InitialValues[4];
//        Ytau2Q=InitialValues[5];
//    }
    
//    if(t1>6.908)    //1 TeV
    if(t1>6.62)   //750 GeV
    {
    double la1Q = InitialValues[9];
    double la2Q = InitialValues[10];
    double la3Q = InitialValues[11];
    double la4Q = InitialValues[12];
    double la5Q = InitialValues[13];

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double         uniLO1   = -3.0*la1Q/(16.0*M_PI);
    gslpp::complex uniNLO1  = -3.0*la1Q/(16.0*M_PI) +9.0*betalambda1/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(9.0*la1Q*la1Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    double         uniLO2   = -3.0*la2Q/(16.0*M_PI);
    gslpp::complex uniNLO2  = -3.0*la2Q/(16.0*M_PI) +9.0*betalambda2/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(9.0*la2Q*la2Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    double         uniLO3   = -(2.0*la3Q+la4Q)/(16.0*M_PI);
    gslpp::complex uniNLO3  = -(2.0*la3Q+la4Q)/(16.0*M_PI) +3.0*(2.0*betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI) +3.0*(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*(2.0*la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO4   = -(la3Q+2.0*la4Q)/(16.0*M_PI);
    gslpp::complex uniNLO4  = -(la3Q+2.0*la4Q)/(16.0*M_PI) +3.0*(betalambda3+2.0*betalambda4)/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+4.0*la3Q*la4Q+4.0*la4Q*la4Q+9.0*la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO5   = -3.0*la5Q/(16.0*M_PI);
    gslpp::complex uniNLO5  = -3.0*la5Q/(16.0*M_PI) +9.0*betalambda5/(256.0*M_PI*M_PI*M_PI) +3.0*(gslpp::complex::i()*M_PI-1.0)*(la3Q+2.0*la4Q)*la5Q/(128.0*M_PI*M_PI*M_PI);
    double         uniLO6   = -la1Q/(16.0*M_PI);
    gslpp::complex uniNLO6  = -la1Q/(16.0*M_PI) +3.0*betalambda1/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO7   = -la2Q/(16.0*M_PI);
    gslpp::complex uniNLO7  = -la2Q/(16.0*M_PI) +3.0*betalambda2/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO8   = -la4Q/(16.0*M_PI);
    gslpp::complex uniNLO8  = -la4Q/(16.0*M_PI) +3.0*betalambda4/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    double         uniLO10  = -la3Q/(16.0*M_PI);
    gslpp::complex uniNLO10 = -la3Q/(16.0*M_PI) +3.0*betalambda3/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO11  = -la5Q/(16.0*M_PI);
    gslpp::complex uniNLO11 = -la5Q/(16.0*M_PI) +3.0*betalambda5/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    double         uniLO14  = -(la3Q-la4Q)/(16.0*M_PI);
    gslpp::complex uniNLO14 = -(la3Q-la4Q)/(16.0*M_PI) +3.0*(betalambda3-betalambda4)/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la3Q-la4Q)*(la3Q-la4Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO15  = -la1Q/(16.0*M_PI);
    gslpp::complex uniNLO15 = -la1Q/(16.0*M_PI) +3.0*betalambda1/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO16  = -la2Q/(16.0*M_PI);
    gslpp::complex uniNLO16 = -la2Q/(16.0*M_PI) +3.0*betalambda2/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO17  = -la5Q/(16.0*M_PI);
    gslpp::complex uniNLO17 = -la5Q/(16.0*M_PI) +3.0*betalambda5/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    double         uniLO24  = -(la3Q+la4Q)/(16.0*M_PI);
    gslpp::complex uniNLO24 = -(la3Q+la4Q)/(16.0*M_PI) +3.0*(betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la3Q+la4Q)*(la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    double         uniLO26  = -(la3Q+la4Q)/(16.0*M_PI);
    gslpp::complex uniNLO26 = -(la3Q+la4Q)/(16.0*M_PI) +3.0*(betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI) +(gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);

    double uniLOev1  = 0.5*(uniLO1+uniLO2+sqrt((uniLO1-uniLO2)*(uniLO1-uniLO2)+4.0*uniLO3*uniLO3));
    double uniLOev2  = 0.5*(uniLO1+uniLO2-sqrt((uniLO1-uniLO2)*(uniLO1-uniLO2)+4.0*uniLO3*uniLO3));
    double uniLOev3  = uniLO4+uniLO5;
    double uniLOev4  = uniLO4-uniLO5;
    double uniLOev5  = 0.5*(uniLO6+uniLO7+sqrt((uniLO6-uniLO7)*(uniLO6-uniLO7)+4.0*uniLO8*uniLO8));
    double uniLOev6  = 0.5*(uniLO6+uniLO7-sqrt((uniLO6-uniLO7)*(uniLO6-uniLO7)+4.0*uniLO8*uniLO8));
    double uniLOev9  = uniLO10+uniLO11;
    double uniLOev10 = uniLO10-uniLO11;
    double uniLOev13 = 0.5*(uniLO15+uniLO16+sqrt((uniLO15-uniLO16)*(uniLO15-uniLO16)+4.0*uniLO17*uniLO17));
    double uniLOev14 = 0.5*(uniLO15+uniLO16-sqrt((uniLO15-uniLO16)*(uniLO15-uniLO16)+4.0*uniLO17*uniLO17));
    gslpp::complex uniNLOev1wowfr  = 0.5*(uniNLO1+uniNLO2+sqrt((uniNLO1-uniNLO2)*(uniNLO1-uniNLO2)+4.0*uniNLO3*uniNLO3));
    gslpp::complex uniNLOev2wowfr  = 0.5*(uniNLO1+uniNLO2-sqrt((uniNLO1-uniNLO2)*(uniNLO1-uniNLO2)+4.0*uniNLO3*uniNLO3));
    gslpp::complex uniNLOev3wowfr  = uniNLO4+uniNLO5;
    gslpp::complex uniNLOev4wowfr  = uniNLO4-uniNLO5;
    gslpp::complex uniNLOev5wowfr  = 0.5*(uniNLO6+uniNLO7+sqrt((uniNLO6-uniNLO7)*(uniNLO6-uniNLO7)+4.0*uniNLO8*uniNLO8));
    gslpp::complex uniNLOev6wowfr  = 0.5*(uniNLO6+uniNLO7-sqrt((uniNLO6-uniNLO7)*(uniNLO6-uniNLO7)+4.0*uniNLO8*uniNLO8));
    gslpp::complex uniNLOev9wowfr  = uniNLO10+uniNLO11;
    gslpp::complex uniNLOev10wowfr  = uniNLO10-uniNLO11;
    gslpp::complex uniNLOev13wowfr  = 0.5*(uniNLO15+uniNLO16+sqrt((uniNLO15-uniNLO16)*(uniNLO15-uniNLO16)+4.0*uniNLO17*uniNLO17));
    gslpp::complex uniNLOev14wowfr = 0.5*(uniNLO15+uniNLO16-sqrt((uniNLO15-uniNLO16)*(uniNLO15-uniNLO16)+4.0*uniNLO17*uniNLO17));

    if( (uniNLOev1wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev2wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev3wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev4wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev5wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev6wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev9wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev10wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev13wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLOev14wowfr - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLO14 - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLO24 - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;
    if( (uniNLO26 - 0.5*gslpp::complex::i()).abs() > 0.5) check=1;

    double Rpeps=0.01;
    if( fabs(uniLOev1) > Rpeps && (uniNLOev1wowfr/uniLOev1-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev2) > Rpeps && (uniNLOev2wowfr/uniLOev2-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev3) > Rpeps && (uniNLOev3wowfr/uniLOev3-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev4) > Rpeps && (uniNLOev4wowfr/uniLOev4-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev5) > Rpeps && (uniNLOev5wowfr/uniLOev5-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev6) > Rpeps && (uniNLOev6wowfr/uniLOev6-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev9) > Rpeps && (uniNLOev9wowfr/uniLOev9-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev10) > Rpeps && (uniNLOev10wowfr/uniLOev10-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev13) > Rpeps && (uniNLOev13wowfr/uniLOev13-1.0).abs() > 1.0) check=1;
    if( fabs(uniLOev14) > Rpeps && (uniNLOev14wowfr/uniLOev14-1.0).abs() > 1.0) check=1;
    if( fabs(uniLO14) > Rpeps && (uniNLO14/uniLO14-1.0).abs() > 1.0) check=1;
    if( fabs(uniLO24) > Rpeps && (uniNLO24/uniLO24-1.0).abs() > 1.0) check=1;
    if( fabs(uniLO26) > Rpeps && (uniNLO26/uniLO26-1.0).abs() > 1.0) check=1;
    
    }

    return check;
}

double Runner::RGERunner(/*int RGEs(), */double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2)
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

    //Possibility to define a set of parameters which the RGE's depend on
    double RGEparameters = 0;

    //Definition of the RGE system (the Jacobian is not necessary for the RK4 method; it's an empty function here)
//    gsl_odeiv2_system RGEsystem = {mySM.RGEs, Jacobian, mySM.NumberOfRGEs, &RGEparameters};
    gsl_odeiv2_system RGEsystem = {RGEs, Jacobian, NumberOfRGEs, &RGEparameters};

    //Set starting and end point as natural logarithmic scales (conversion from decadic log scale)
    double t1 = Q1*log(10.0);
    double t2 = Q2*log(10.0);

    //Set initial step size
    double InitialStepSize = 1e-6;

    //Run!
    while (t1 < t2)
    {
        int status = gsl_odeiv2_evolve_apply (e, c, s, &RGEsystem, &t1, t2, &InitialStepSize, InitialValues);
        if(status != GSL_SUCCESS) break;

        //intermediate checks if appropriate
        if(RGEcheck(InitialValues,t1) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

//    for(int i=0;i<14;i++) std::cout<<InitialValues[i]<<std::endl;

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}
