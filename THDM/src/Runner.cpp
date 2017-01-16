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

int RGEs(double t, const double y[], double beta[], void *flags)
{
    (void)(t); /* avoid unused parameter warning */
    int flag = *(int *)flags;
    int ord = flag%3;
    int type = flag-ord;
    double Yb1=0;
    double Yb2=0;
    double Ytau1=0;
    double Ytau2=0;
    double g1=y[0];
    double g2=y[1];
    double g3=y[2];
    double Yt=y[3];
    if( type == 0 ) {//type I
        Yb2=y[4];
        Ytau2=y[5];
    }
    else if( type == 3 ) {//type II
        Yb1=y[4];
        Ytau1=y[5];
    }
    else if( type == 6 ) {//type X
        Yb2=y[4];
        Ytau1=y[5];
    }
    else if( type == 9 ) {//type Y
        Yb1=y[4];
        Ytau2=y[5];
    }
    else {
        throw std::runtime_error("type should only be any of 0, 3, 6, 9");
    }
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
    beta[3] = Yt*(-17.0*g1*g1+3.0*(-9.0*g2*g2-32.0*g3*g3+2.0*Yb1*Yb1+6.0*Yb2*Yb2+18.0*Yt*Yt+4.0*Ytau2*Ytau2))/(192.0*pi*pi);
    //beta_Yb
    beta[4] = ((-5.0*g1*g1-27.0*g2*g2-96.0*g3*g3)*(Yb1+Yb2)+(6.0*Yb1+18.0*Yb2)*Yt*Yt
               +54.0*(Yb1*Yb1*Yb1+Yb2*Yb2*Yb2+Yb1*Yb2*Yb2)+12.0*(Yb1*Ytau1*Ytau1+Yb2*Ytau2*Ytau2))/(192.0*pi*pi);
    //beta_Ytau
    beta[5] = ((-15.0*g1*g1-9.0*g2*g2)*(Ytau1+Ytau2)+Ytau2*(12.0*Yb2*Yb2+12.0*Yt*Yt+10.0*Ytau2*Ytau2)
               +Ytau1*(12.0*Yb1*Yb1+10.0*Ytau1*Ytau1))/(64.0*pi*pi);
    //beta_m11_2
    beta[6] = (-3.0*g1*g1*m11_2-9.0*g2*g2*m11_2
               +4.0*(m11_2*(3.0*Yb1*Yb1+Ytau1*Ytau1+3.0*la1)+m22_2*(2.0*la3+la4)))/(32.0*pi*pi);
    //beta_m22_2
    beta[7] = (-3.0*g1*g1*m22_2-9.0*g2*g2*m22_2
               +4.0*(m22_2*(3.0*Yb2*Yb2+3.0*Yt*Yt+Ytau2*Ytau2+3.0*la2)+m11_2*(2.0*la3+la4)))/(32.0*pi*pi);
    //beta_m12_2
    beta[8] = m12_2*(-3.0*g1*g1-9.0*g2*g2
                     +2.0*(3.0*Yb1*Yb1+3.0*Yb2*Yb2+3.0*Yt*Yt+Ytau1*Ytau1+Ytau2*Ytau2+2.0*la3+4.0*la4+6.0*la5))/(32.0*pi*pi);
    //beta_lambda_1
    beta[9] = (3.0*g1*g1*g1*g1+9.0*g2*g2*g2*g2+6.0*g1*g1*(g2*g2-2.0*la1)-36.0*g2*g2*la1
               +8.0*(6.0*la1*la1+2.0*la3*la3+2.0*la3*la4+la4*la4+la5*la5)
               -16.0*(3.0*Yb1*Yb1*Yb1*Yb1+Ytau1*Ytau1*Ytau1*Ytau1)+16.0*la1*(3.0*Yb1*Yb1+Ytau1*Ytau1))/(64.0*pi*pi);
    //beta_lambda_2
    beta[10] = (3.0*g1*g1*g1*g1+9.0*g2*g2*g2*g2+6.0*g1*g1*(g2*g2-2.0*la2)-36.0*g2*g2*la2
                -8.0*(6.0*Yb2*Yb2*Yb2*Yb2+6.0*Yt*Yt*Yt*Yt+2.0*Ytau2*Ytau2*Ytau2*Ytau2-6.0*Yb2*Yb2*la2
                      -6.0*Yt*Yt*la2-2.0*Ytau2*Ytau2*la2-6.0*la2*la2-2.0*la3*la3-2.0*la3*la4-la4*la4-la5*la5))/(64.0*pi*pi);
    //beta_lambda_3
    beta[11] = (3.0*g1*g1*g1*g1+9.0*g2*g2*g2*g2-36.0*g2*g2*la3-6.0*g1*g1*(g2*g2+2.0*la3)
                +8.0*(3.0*Yb2*Yb2*la3+3.0*Yt*Yt*la3+Ytau1*Ytau1*la3+Ytau2*Ytau2*la3+3.0*la1*la3+3.0*la2*la3+2.0*la3*la3
                      +Yb1*Yb1*(-6.0*Yt*Yt+3.0*la3)+la1*la4+la2*la4+la4*la4+la5*la5))/(64.0*pi*pi);
    //beta_lambda_4
    beta[12] = (3.0*g1*g1*(g2*g2-la4)-9.0*g2*g2*la4+6.0*Yb2*Yb2*la4+6.0*Yt*Yt*la4+2.0*Ytau1*Ytau1*la4
                +2.0*Ytau2*Ytau2*la4+2.0*la1*la4+2.0*la2*la4+8.0*la3*la4+4.0*la4*la4+6.0*Yb1*Yb1*(2.0*Yt*Yt+la4)
                +8.0*la5*la5)/(16.0*pi*pi);
    //beta_lambda_5
    beta[13] = (-3.0*g1*g1-9.0*g2*g2
                +2.0*(3.0*Yb1*Yb1+3.0*Yb2*Yb2+3.0*Yt*Yt+Ytau1*Ytau1+Ytau2*Ytau2+la1+la2+4.0*la3+6.0*la4))*la5/(16.0*pi*pi);

    if(ord==1||ord==2){
        //beta_g1
        beta[0] += (g1*g1*g1*(88.0*g3*g3-17.0*Yt*Yt))/(1536.0*pi*pi*pi*pi);
        //beta_g2
        beta[1] += (3.0*g2*g2*g2*(8.0*g3*g3-Yt*Yt))/(512.0*pi*pi*pi*pi);
        //beta_g3
        beta[2] += -((g3*g3*g3*(13.0*g3*g3+Yt*Yt))/(128.0*pi*pi*pi*pi));
        //beta_Yt
        beta[3] += (Yt*(-216.0*g3*g3*g3*g3+72.0*g3*g3*Yt*Yt-24.0*Yt*Yt*Yt*Yt-12.0*Yt*Yt*la2
                        +3.0*la2*la2+2.0*la3*la3+2.0*la3*la4+2.0*la4*la4+3.0*la5*la5))/(512.0*pi*pi*pi*pi);
        //beta_Yb
        beta[4] += (-1296.0*g3*g3*g3*g3*(Yb1+Yb2)+16.0*g3*g3*(4.0*Yb1+3.0*Yb2)*Yt*Yt
                    -3.0*(2.0*Yb1*(5.0*Yt*Yt*Yt*Yt-3.0*la1*la1-2.0*la3*la3
                                   +4.0*Yt*Yt*(la3-la4)-2.0*la3*la4
                                   -2.0*la4*la4-3.0*la5*la5)
                          +Yb2*(Yt*Yt*Yt*Yt-2.0*(3.0*la2*la2+2.0*la3*la3+2.0*la3*la4
                                                 +2.0*la4*la4+3.0*la5*la5))))
                   /(3072.0*pi*pi*pi*pi);
        //beta_Ytau
        beta[5] += (80.0*g3*g3*Yt*Yt*Ytau2-27.0*Yt*Yt*Yt*Yt*Ytau2+6.0*Ytau1*la1*la1+6.0*Ytau2*la2*la2
                    +(Ytau1+Ytau2)*(4.0*la3*la3+4.0*la3*la4+4.0*la4*la4+6.0*la5*la5))
                   /(1024.0*pi*pi*pi*pi);
        //beta_m11_2
        beta[6] += -(8.0*m22_2*(2.0*la3*la3+2.0*la3*la4+2.0*la4*la4
                                +3.0*Yt*Yt*(2.0*la3+la4)+3.0*la5*la5)
                     +m11_2*(30.0*la1*la1+4.0*la3*la3+4.0*la3*la4
                             +4.0*la4*la4+6.0*la5*la5))/(512.0*pi*pi*pi*pi);
        //beta_m22_2
        beta[7] += -((-80.0*g3*g3*m22_2*Yt*Yt
                      +8.0*m11_2*(2.0*la3*la3+2.0*la3*la4+2.0*la4*la4+3.0*la5*la5)
                      +m22_2*(27.0*Yt*Yt*Yt*Yt+72.0*Yt*Yt*la2+30.0*la2*la2+4.0*la3*la3
                              +4.0*la3*la4+4.0*la4*la4+6.0*la5*la5))/(512.0*pi*pi*pi*pi));
        //beta_m12_2
        beta[8] += (m12_2*(72.0*(la1*la1+la2*la2-4.0*la3*la4-8.0*la3*la5-8.0*la4*la5
                                 +2.0*la5*la5-4.0*(la1+la2)*(la3+la4+la5))
                           +Yt*Yt*(960.0*g3*g3
                                   -36.0*(9.0*Yt*Yt+8.0*(la3+2.0*la4+3.0*la5)))))
                   /(12288.0*pi*pi*pi*pi);
        //beta_lambda_1
        beta[9] += -(39.0*la1*la1*la1
                     +la1*(10.0*la3*la3+10.0*la3*la4+6.0*la4*la4+7.0*la5*la5)
                     +2.0*(4.0*la3*la3*la3+6.0*la3*la3*la4+8.0*la3*la4*la4
                     +3.0*la4*la4*la4+10.0*la3*la5*la5+11.0*la4*la5*la5
                     +3.0*Yt*Yt*(2.0*la3*la3+2.0*la3*la4+la4*la4+la5*la5)))
                   /(128.0*pi*pi*pi*pi);
        //beta_lambda_2
        beta[10] += -(-60.0*pow(Yt,6)+3.0*Yt*Yt*Yt*Yt*la2+72.0*Yt*Yt*la2*la2
                      +16.0*g3*g3*(4.0*Yt*Yt*Yt*Yt-5.0*Yt*Yt*la2)
                      +2.0*(39.0*la2*la2*la2+10.0*la2*la3*la3+8.0*la3*la3*la3
                            +10.0*la2*la3*la4+12.0*la3*la3*la4+6.0*la2*la4*la4
                            +16.0*la3*la4*la4+6.0*la4*la4*la4+7.0*la2*la5*la5
                            +20.0*la3*la5*la5+22.0*la4*la5*la5))/(256.0*pi*pi*pi*pi);
        //beta_lambda_3
        beta[11] += -(-80.0*g3*g3*Yt*Yt*la3+27.0*Yt*Yt*Yt*Yt*la3
                      +12.0*Yt*Yt*(2.0*la3*la3+la4*la4+2.0*la2*(3.0*la3+la4)+la5*la5)
                      +2.0*(la1*la1*(15.0*la3+4.0*la4)+la2*la2*(15.0*la3+4.0*la4)
                            +2.0*(la1+la2)*(18.0*la3*la3+8.0*la3*la4+7.0*la4*la4+9.0*la5*la5)
                            +2.0*(6.0*la3*la3*la3+2.0*la3*la3*la4+8.0*la3*la4*la4
                                  +6.0*la4*la4*la4+9.0*la3*la5*la5+22.0*la4*la5*la5)))
                     /(512.0*pi*pi*pi*pi);
        //beta_lambda_4
        beta[12] += -(-80.0*g3*g3*Yt*Yt*la4+27.0*Yt*Yt*Yt*Yt*la4
                      +24.0*Yt*Yt*(la2*la4+2.0*la3*la4+la4*la4+2.0*la5*la5)
                      +2.0*(7.0*la1*la1*la4+7.0*la2*la2*la4+28.0*la3*la3*la4
                            +28.0*la3*la4*la4+48.0*la3*la5*la5+26.0*la4*la5*la5
                            +4.0*(la1+la2)*(10.0*la3*la4+5.0*la4*la4+6.0*la5*la5))
                     )/(512.0*pi*pi*pi*pi);
        //beta_lambda_5
        beta[13] += (la5*(80.0*g3*g3*Yt*Yt-3.0*Yt*Yt*Yt*Yt-24.0*Yt*Yt*(la2+2.0*la3+3.0*la4)
                          -2.0*(7.0*la1*la1+7.0*la2*la2+40.0*la1*la3+40.0*la2*la3+28.0*la3*la3
                                +44.0*la1*la4+44.0*la2*la4+76.0*la3*la4+32.0*la4*la4-6.0*la5*la5))
                    )/(512 *pi*pi*pi*pi);

        if(ord==2){
            //beta_g1
            beta[0] += (g1*g1*g1*(208.0*g1*g1+108.0*g2*g2-15.0*Yb1*Yb1-15.0*Yb2*Yb2
                                  -45.0*Ytau1*Ytau1-45.0*Ytau2*Ytau2))/(4608.0*pi*pi*pi*pi);
            //beta_g2
            beta[1] += (g2*g2*g2*(4.0*g1*g1+16.0*g2*g2-3.0*Yb1*Yb1-3.0*Yb2*Yb2
                        -Ytau1*Ytau1-Ytau2*Ytau2))/(512.0*pi*pi*pi*pi);
            //beta_g3
            beta[2] += (g3*g3*g3*(11.0*g1*g1+3.0*(9.0*g2*g2-4.0*(Yb1*Yb1+Yb2*Yb2))))/(1536.0*pi*pi*pi*pi);
            //beta_Yt
            beta[3] += (2534.0*g1*g1*g1*g1
                        +g1*g1*(-324.0*g2*g2+912.0*g3*g3-123.0*Yb1*Yb1+63.0*Yb2*Yb2
                                +3537.0*Yt*Yt+1350.0*Ytau2*Ytau2)
                        -9.0*(252.0*g2*g2*g2*g2
                              -9.0*g2*g2*(48.0*g3*g3+11.0*Yb1*Yb1+33.0*Yb2*Yb2+75.0*Yt*Yt+10.0*Ytau2*Ytau2)
                              -4.0*(16.0*g3*g3*(4.0*Yb1*Yb1+3.0*Yb2*Yb2)
                                    -3.0*(10.0*Yb1*Yb1*Yb1*Yb1+Yb2*Yb2*Yb2*Yb2
                                          +Yb2*Yb2*(11.0*Yt*Yt-5.0*Ytau2*Ytau2)
                                          +9.0*Ytau2*Ytau2*(Yt*Yt+Ytau2*Ytau2)
                                          +Yb1*Yb1*(10.0*Yt*Yt+3.0*Ytau1*Ytau1+8.0*la3-8.0*la4)))))
                       *Yt/(110592.0*pi*pi*pi*pi);
            //beta_Yb
            beta[4] += -(226.0*g1*g1*g1*g1*(Yb1+Yb2)
                         +3.0*g1*g1*(-711.0*Yb1*Yb1*Yb1-711.0*Yb2*Yb2*Yb2+324.0*g2*g2*(Yb1+Yb2)
                                     -496.0*g3*g3*(Yb1+Yb2)+53.0*Yb1*Yt*Yt-273.0*Yb2*Yt*Yt
                                     -450.0*Yb1*Ytau1*Ytau1-450.0*Yb2*Ytau2*Ytau2)
                         +27.0*(84.0*g2*g2*g2*g2*(Yb1+Yb2)
                                -3.0*g2*g2*(75.0*Yb1*Yb1*Yb1+75.0*Yb2*Yb2*Yb2+48.0*g3*g3*(Yb1+Yb2)
                                            +11.0*Yb1*Yt*Yt+33.0*Yb2*Yt*Yt+10.0*Yb1*Ytau1*Ytau1
                                            +10.0*Yb2*Ytau2*Ytau2)
                                +4.0*(48.0*pow(Yb1,5)-144.0*g3*g3*(Yb1*Yb1*Yb1+Yb2*Yb2*Yb2)
                                      +3.0*Yb1*(3.0*Ytau1*Ytau1*Ytau1*Ytau1+Yt*Yt*Ytau2*Ytau2)
                                      +Yb1*Yb1*Yb1*(10.0*Yt*Yt+9.0*Ytau1*Ytau1+24.0*la1)
                                      +Yb2*(48.0*Yb2*Yb2*Yb2*Yb2-5.0*Yt*Yt*Ytau2*Ytau2
                                      +9.0*Ytau2*Ytau2*Ytau2*Ytau2
                                      +Yb2*Yb2*(11.0*Yt*Yt+9.0*Ytau2*Ytau2+24.0*la2)))))
                        /(110592.0*pi*pi*pi*pi);
            //beta_Ytau
            beta[5] += (966.0*g1*g1*g1*g1*(Ytau1+Ytau2)
                        +g1*g1*(50.0*Yb1*Yb1*Ytau1+537.0*Ytau1*Ytau1*Ytau1+50.0*Yb2*Yb2*Ytau2
                                +170.0*Yt*Yt*Ytau2+537.0*Ytau2*Ytau2*Ytau2+108.0*g2*g2*(Ytau1+Ytau2))
                        -3.0*(84.0*g2*g2*g2*g2*(Ytau1+Ytau2)
                              -15.0*g2*g2*(6.0*Yb1*Yb1*Ytau1+11.0*Ytau1*Ytau1*Ytau1+6.0*Yb2*Yb2*Ytau2
                                           +6.0*Yt*Yt*Ytau2+11.0*Ytau2*Ytau2*Ytau2)
                              +4.0*(-80.0*g3*g3*(Yb1*Yb1*Ytau1+Yb2*Yb2*Ytau2)
                                    +3.0*(9.0*Yb1*Yb1*Yb1*Yb1*Ytau1+4.0*pow(Ytau1,5)+9.0*Yb2*Yb2*Yb2*Yb2*Ytau2
                                          -2.0*Yb2*Yb2*Yt*Yt*Ytau2+9.0*Yb2*Yb2*Ytau2*Ytau2*Ytau2
                                          +9.0*Yt*Yt*Ytau2*Ytau2*Ytau2+4.0*pow(Ytau2,5)
                                          +3.0*Yb1*Yb1*(3.0*Ytau1*Ytau1*Ytau1+Yt*Yt*(Ytau1+Ytau2))
                                          +8.0*Ytau1*Ytau1*Ytau1*la1+8.0*Ytau2*Ytau2*Ytau2*la2))))
                       /(12288.0*pi*pi*pi*pi);
            //beta_m11_2
            beta[6] += (3.0*g1*g1*g1*g1*(193.0*m11_2+40.0*m22_2)
                        +2.0*g1*g1*(45.0*g2*g2*m11_2
                                    +2.0*(m11_2*(25.0*Yb1*Yb1+75.0*Ytau1*Ytau1+144.0*la1)
                                          +48.0*m22_2*(2.0*la3+la4)))
                        -3.0*(3.0*g2*g2*g2*g2*(41.0*m11_2-40.0*m22_2)
                              -12.0*g2*g2*(m11_2*(15.0*Yb1*Yb1+5.0*Ytau1*Ytau1+48.0*la1)
                                           +16.0*m22_2*(2.0*la3+la4))
                              +8.0*(-80.0*g3*g3*m11_2*Yb1*Yb1
                                    +3.0*m11_2*(9.0*Yb1*Yb1*Yb1*Yb1+3.0*Ytau1*Ytau1*Ytau1*Ytau1
                                                +8.0*Ytau1*Ytau1*la1+3.0*Yb1*Yb1*(Yt*Yt+8.0*la1))
                                    +8.0*m22_2*(3.0*Yb2*Yb2+Ytau2*Ytau2)*(2.0*la3+la4))))
                       /(12288.0*pi*pi*pi*pi);
            //beta_m22_2
            beta[7] += (3.0*g1*g1*g1*g1*(40.0*m11_2+193.0*m22_2)
                        +2.0*g1*g1*(45.0*g2*g2*m22_2
                                    +2.0*(m22_2*(25.0*Yb2*Yb2+85.0*Yt*Yt+75.0*Ytau2*Ytau2+144.0*la2)
                                          +48.0*m11_2*(2.0*la3+la4)))
                        +3.0*(3.0*g2*g2*g2*g2*(40.0*m11_2-41.0*m22_2)
                              +12.0*g2*g2*(m22_2*(15.0*Yb2*Yb2+15.0*Yt*Yt+5.0*Ytau2*Ytau2+48.0*la2)
                                           +16.0*m11_2*(2.0*la3+la4))
                              +8.0*(80.0*g3*g3*m22_2*Yb2*Yb2
                                    -3.0*m22_2*(9.0*Yb2*Yb2*Yb2*Yb2+3.0*Yb1*Yb1*Yt*Yt
                                                +3.0*Ytau2*Ytau2*Ytau2*Ytau2+8.0*Ytau2*Ytau2*la2
                                                +2.0*Yb2*Yb2*(7.0*Yt*Yt+12.0*la2))
                                    -8.0*m11_2*(3.0*Yb1*Yb1+Ytau1*Ytau1)*(2.0*la3+la4))))
                       /(12288.0*pi*pi*pi*pi);
            //beta_m12_2
            beta[8] += (9.0*(51.0*g1*g1*g1*g1+10.0*g1*g1*g2*g2-81.0*g2*g2*g2*g2)
                        -324.0*(Yb1*Yb1*Yb1*Yb1+Yb2*Yb2*Yb2*Yb2)
                        +2.0*(85.0*g1*g1+135.0*g2*g2-432.0*Yb1*Yb1)*Yt*Yt
                        -108.0*(Ytau1*Ytau1*Ytau1*Ytau1+Ytau2*Ytau2*Ytau2*Ytau2)
                        +2.0*(Yb1*Yb1+Yb2*Yb2)*(25.0*g1*g1
                                                +3.0*(45.0*g2*g2+4.0*(40.0*g3*g3+3.0*Yt*Yt-12.0*la3
                                                                      -24.0*la4-36.0*la5)))
                        +192.0*(g1*g1+3.0*g2*g2)*(la3+2.0*la4+3.0*la5)
                        +6.0*(Ytau1*Ytau1+Ytau2*Ytau2)*(25.0*g1*g1+15.0*g2*g2
                                                        -16.0*(la3+2.0*la4+3.0*la5)))
                       *m12_2/(12288.0*pi*pi*pi*pi);
            //beta_lambda_1
            beta[9] += (-393.0*pow(g1,6)
                        +g1*g1*g1*g1*(-573.0*g2*g2+60.0*Yb1*Yb1-300.0*Ytau1*Ytau1
                                      +651.0*la1+120.0*la3+60.0*la4)
                        +3.0*(291.0*pow(g2,6)
                              -3.0*g2*g2*g2*g2*(12.0*Yb1*Yb1+4.0*Ytau1*Ytau1+17.0*la1-40.0*la3-20.0*la4)
                              +12.0*g2*g2*(15.0*Yb1*Yb1*la1+5.0*Ytau1*Ytau1*la1
                                           +4.0*(9.0*la1*la1+(2.0*la3+la4)*(2.0*la3+la4)))
                              -8.0*(-60.0*pow(Yb1,6)-20.0*pow(Ytau1,6)+Ytau1*Ytau1*Ytau1*Ytau1*la1
                                    +24.0*Ytau1*Ytau1*la1*la1+3.0*Yb1*Yb1*Yb1*Yb1*(-4.0*Yt*Yt+la1)
                                    +9.0*Yb1*Yb1*la1*(Yt*Yt+8.0*la1)
                                    +16.0*g3*g3*(4.0*Yb1*Yb1*Yb1*Yb1-5.0*Yb1*Yb1*la1)
                                    +24.0*Yb2*Yb2*la3*la3+8.0*Ytau2*Ytau2*la3*la3+24.0*Yb2*Yb2*la3*la4
                                    +8.0*Ytau2*Ytau2*la3*la4+12.0*Yb2*Yb2*la4*la4+4.0*Ytau2*Ytau2*la4*la4
                                    +12.0*Yb2*Yb2*la5*la5+4.0*Ytau2*Ytau2*la5*la5))
                        +g1*g1*(-303.0*g2*g2*g2*g2
                                +6.0*g2*g2*(36.0*Yb1*Yb1+44.0*Ytau1*Ytau1+39.0*la1+20.0*la4)
                                +4.0*(16.0*Yb1*Yb1*Yb1*Yb1+25.0*Yb1*Yb1*la1
                                      +3.0*(-16.0*Ytau1*Ytau1*Ytau1*Ytau1+25.0*Ytau1*Ytau1*la1+36.0*la1*la1
                                            +16.0*la3*la3+16.0*la3*la4+8.0*la4*la4-4.0*la5*la5))))
                       /(6144.0*pi*pi*pi*pi);
            //beta_lambda_2
            beta[10] += (-393.0*pow(g1,6)
                         +g1*g1*g1*g1*(-573.0*g2*g2+60.0*Yb2*Yb2-228.0*Yt*Yt-300.0*Ytau2*Ytau2
                                       +651.0*la2+120.0*la3+60.0*la4)
                         +g1*g1*(-303.0*g2*g2*g2*g2
                                 +6.0*g2*g2*(36.0*Yb2*Yb2+84.0*Yt*Yt+44.0*Ytau2*Ytau2+39.0*la2+20.0*la4)
                                 +4.0*(16.0*Yb2*Yb2*Yb2*Yb2-32.0*Yt*Yt*Yt*Yt-48.0*Ytau2*Ytau2*Ytau2*Ytau2
                                       +25.0*Yb2*Yb2*la2+85.0*Yt*Yt*la2+75.0*Ytau2*Ytau2*la2
                                       +108.0*la2*la2+48.0*la3*la3+48.0*la3*la4+24.0*la4*la4-12.0*la5*la5))
                         +3.0*(291.0*pow(g2,6)
                               -3.0*g2*g2*g2*g2*(12.0*Yb2*Yb2+12.0*Yt*Yt+4.0*Ytau2*Ytau2
                                                 +17.0*la2-40.0*la3-20.0*la4)
                               +12.0*g2*g2*(15.0*Yb2*Yb2*la2+15.0*Yt*Yt*la2+5.0*Ytau2*Ytau2*la2
                                            +36.0*la2*la2+16.0*la3*la3+16.0*la3*la4+4.0*la4*la4)
                               -8.0*(-60.0*pow(Yb2,6)-12.0*Yb1*Yb1*Yt*Yt*Yt*Yt-20.0*pow(Ytau2,6)
                                     +9.0*Yb1*Yb1*Yt*Yt*la2+Ytau2*Ytau2*Ytau2*Ytau2*la2
                                     +24.0*Ytau2*Ytau2*la2*la2+3.0*Yb2*Yb2*Yb2*Yb2*(4.0*Yt*Yt+la2)
                                     +16.0*g3*g3*(4.0*Yb2*Yb2*Yb2*Yb2-5.0*Yb2*Yb2*la2)
                                     +6.0*Yb2*Yb2*(2.0*Yt*Yt*Yt*Yt+7.0*Yt*Yt*la2+12.0*la2*la2)
                                     +24.0*Yb1*Yb1*la3*la3+8.0*Ytau1*Ytau1*la3*la3+24.0*Yb1*Yb1*la3*la4
                                     +8.0*Ytau1*Ytau1*la3*la4+12.0*Yb1*Yb1*la4*la4+4.0*Ytau1*Ytau1*la4*la4
                                     +12.0*Yb1*Yb1*la5*la5+4.0*Ytau1*Ytau1*la5*la5)))/(6144.0*pi*pi*pi*pi);
            //beta_lambda_3
            beta[11] += (-393.0*pow(g1,6)
                         +3.0*g1*g1*g1*g1*(101.0*g2*g2+10.0*Yb1*Yb1+10.0*Yb2*Yb2-38.0*Yt*Yt-50.0*Ytau1*Ytau1
                                           -50.0*Ytau2*Ytau2+30.0*la1+30.0*la2+197.0*la3+20.0*la4)
                         +g1*g1*(33.0*g2*g2*g2*g2-6.0*g2*g2*(18.0*Yb1*Yb1+18.0*Yb2*Yb2+42.0*Yt*Yt
                                                             +22.0*Ytau1*Ytau1+22.0*Ytau2*Ytau2+10.0*la1
                                                             +10.0*la2-11.0*la3+12.0*la4)
                                 +2.0*(25.0*Yb2*Yb2*la3+85.0*Yt*Yt*la3+75.0*Ytau1*Ytau1*la3
                                       +75.0*Ytau2*Ytau2*la3+144.0*la1*la3+144.0*la2*la3+24*la3*la3
                                       +Yb1*Yb1*(-16.0*Yt*Yt+25.0*la3)
                                       +48.0*la1*la4+48.0*la2*la4-24.0*la4*la4+48.0*la5*la5))
                         +3.0*(291.0*pow(g2,6)-3.0*g2*g2*g2*g2*(6.0*Yb1*Yb1+6.0*Yb2*Yb2+6.0*Yt*Yt
                                                                +2.0*Ytau1*Ytau1+2.0*Ytau2*Ytau2
                                                                -30.0*la1-30.0*la2+37.0*la3-20.0*la4)
                               +6.0*g2*g2*(15.0*Yb1*Yb1*la3+15.0*Yb2*Yb2*la3+15.0*Yt*Yt*la3
                                           +5.0*Ytau1*Ytau1*la3+5.0*Ytau2*Ytau2*la3+48.0*la1*la3
                                           +48.0*la2*la3+8.0*la3*la3+24.0*la1*la4+24.0*la2*la4
                                           -16.0*la3*la4+8.0*la4*la4)
                               -4.0*(-9.0*Yb1*Yb1*Yb1*Yb1*(8.0*Yt*Yt-3.0*la3)+27.0*Yb2*Yb2*Yb2*Yb2*la3
                                     +42.0*Yb2*Yb2*Yt*Yt*la3+9.0*Ytau1*Ytau1*Ytau1*Ytau1*la3
                                     +9.0*Ytau2*Ytau2*Ytau2*Ytau2*la3+24.0*Ytau1*Ytau1*la1*la3
                                     +72.0*Yb2*Yb2*la2*la3+24.0*Ytau2*Ytau2*la2*la3+24.0*Yb2*Yb2*la3*la3
                                     +8.0*Ytau1*Ytau1*la3*la3+8.0*Ytau2*Ytau2*la3*la3
                                     +16.0*g3*g3*(Yb1*Yb1*(8.0*Yt*Yt-5.0*la3)-5.0*Yb2*Yb2*la3)
                                     +48.0*Yb2*Yb2*Yt*Yt*la4+8.0*Ytau1*Ytau1*la1*la4+24.0*Yb2*Yb2*la2*la4
                                     +8.0*Ytau2*Ytau2*la2*la4+12.0*Yb2*Yb2*la4*la4+4.0*Ytau1*Ytau1*la4*la4
                                     +4.0*Ytau2*Ytau2*la4*la4+12.0*Yb2*Yb2*la5*la5+4.0*Ytau1*Ytau1*la5*la5
                                     +4.0*Ytau2*Ytau2*la5*la5
                                     -6.0*Yb1*Yb1*(12.0*Yt*Yt*Yt*Yt+5.0*Yt*Yt*la3
                                                   -2.0*(6.0*la1*la3+2.0*la3*la3+2.0*la1*la4+la4*la4+la5*la5)))))
                        /(6144.0*pi*pi*pi*pi);
            //beta_lambda_4
            beta[12] += (g1*g1*g1*g1*(-876.0*g2*g2+471.0*la4)
                         +2.0*g1*g1*(-168.0*g2*g2*g2*g2+25.0*Yb2*Yb2*la4+85.0*Yt*Yt*la4
                                     +75.0*Ytau1*Ytau1*la4+75.0*Ytau2*Ytau2*la4
                                     +48.0*la1*la4+48.0*la2*la4+48.0*la3*la4+96.0*la4*la4
                                     +Yb1*Yb1*(16.0*Yt*Yt+25.0*la4)
                                     +3.0*g2*g2*(36.0*Yb1*Yb1+36.0*Yb2*Yb2+84.0*Yt*Yt
                                                 +44.0*Ytau1*Ytau1+44.0*Ytau2*Ytau2
                                                 +20.0*la1+20.0*la2+8.0*la3+51.0*la4)
                                     +192.0*la5*la5)
                         -3.0*(231.0*g2*g2*g2*g2*la4-90.0*g2*g2*Yb2*Yb2*la4
                               +108.0*Yb2*Yb2*Yb2*Yb2*la4-90.0*g2*g2*Yt*Yt*la4
                               -216.0*Yb2*Yb2*Yt*Yt*la4-30.0*g2*g2*Ytau1*Ytau1*la4
                               +36.0*Ytau1*Ytau1*Ytau1*Ytau1*la4-30.0*g2*g2*Ytau2*Ytau2*la4
                               +36.0*Ytau2*Ytau2*Ytau2*Ytau2*la4+32.0*Ytau1*Ytau1*la1*la4
                               +96.0*Yb2*Yb2*la2*la4+32.0*Ytau2*Ytau2*la2*la4-288.0*g2*g2*la3*la4
                               +192.0*Yb2*Yb2*la3*la4+64.0*Ytau1*Ytau1*la3*la4+64.0*Ytau2*Ytau2*la3*la4
                               -144.0*g2*g2*la4*la4+96.0*Yb2*Yb2*la4*la4+32.0*Ytau1*Ytau1*la4*la4
                               +32.0*Ytau2*Ytau2*la4*la4+12.0*Yb1*Yb1*Yb1*Yb1*(16.0*Yt*Yt+9.0*la4)
                               -64.0*g3*g3*(5.0*Yb2*Yb2*la4+Yb1*Yb1*(8.0*Yt*Yt+5.0*la4))
                               -432.0*g2*g2*la5*la5+192.0*Yb2*Yb2*la5*la5+64.0*Ytau1*Ytau1*la5*la5
                               +64.0*Ytau2*Ytau2*la5*la5
                               +6.0*Yb1*Yb1*(32.0*Yt*Yt*Yt*Yt-15.0*g2*g2*la4
                                             +4.0*Yt*Yt*(8.0*la3+11.0*la4)
                                             +16.0*(la1*la4+2.0*la3*la4+la4*la4+2.0*la5*la5))))
                        /(6144.0*pi*pi*pi*pi);
            //beta_lambda_5
            beta[13] += (471.0*g1*g1*g1*g1
                         +2.0*g1*g1*(57.0*g2*g2+25.0*Yb1*Yb1+25.0*Yb2*Yb2+85.0*Yt*Yt+75.0*Ytau1*Ytau1
                                     +75.0*Ytau2*Ytau2-24.0*la1-24.0*la2+192.0*la3+288.0*la4)
                         -3.0*(231.0*g2*g2*g2*g2
                               -6.0*g2*g2*(15.0*Yb1*Yb1+15.0*Yb2*Yb2+15.0*Yt*Yt+5.0*Ytau1*Ytau1
                                           +5.0*Ytau2*Ytau2+48.0*la3+96.0*la4)
                         +4.0*(3.0*Yb1*Yb1*Yb1*Yb1+3.0*Yb2*Yb2*Yb2*Yb2-80.0*g3*g3*(Yb1*Yb1+Yb2*Yb2)
                               -6.0*Yb2*Yb2*Yt*Yt+Ytau1*Ytau1*Ytau1*Ytau1+Ytau2*Ytau2*Ytau2*Ytau2
                               +8.0*Ytau1*Ytau1*la1+24.0*Yb2*Yb2*la2+8.0*Ytau2*Ytau2*la2+48.0*Yb2*Yb2*la3
                               +16.0*Ytau1*Ytau1*la3+16.0*Ytau2*Ytau2*la3+72.0*Yb2*Yb2*la4
                               +24.0*Ytau1*Ytau1*la4+24.0*Ytau2*Ytau2*la4
                               +6.0*Yb1*Yb1*(11.0*Yt*Yt+4.0*la1+8.0*la3+12.0*la4))))
                        *la5/(6144.0*pi*pi*pi*pi);
        }
    }

    return 0;
}

int Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *order)
{
    return 0;
}

int RGEcheck(const double InitialValues[], const double t1, const double Rpeps, const double tNLOuni)
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
    if(t1>tNLOuni)
    {
    double la1Q = InitialValues[9];
    double la2Q = InitialValues[10];
    double la3Q = InitialValues[11];
    double la4Q = InitialValues[12];
    double la5Q = InitialValues[13];

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q;
//                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
//                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q + 6.0*la2Q*YtQ*YtQ - 6.0*YtQ*YtQ*YtQ*YtQ;
//                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q 
//                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q;
    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q + 3.0*la3Q*YtQ*YtQ;
//                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q
//                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q + 3.0*la4Q*YtQ*YtQ;
//                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q
//                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q + 3.0*YtQ*YtQ);
//                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q)
//                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

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
    
    }

    return check;
}

double Runner::RGERunner(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double NLOuniscale)
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
    gsl_odeiv2_system RGEsystem = {RGEs, Jacobian, NumberOfRGEs, &order};

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
        if(RGEcheck(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}
