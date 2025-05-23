/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMZ2Runner.h"
#include "GeneralTHDMcache.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

GeneralTHDMZ2Runner::GeneralTHDMZ2Runner(const StandardModel& SM_i)
: myGTHDMZ2(static_cast<const GeneralTHDMZ2*> (&SM_i)), myZ2_at_Q(3, 5, 0.)
{}

GeneralTHDMZ2Runner::~GeneralTHDMZ2Runner()
{};

int RGEsZ2(double t, const double y[], double beta[], void *flags)
{
    (void)(t); /* avoid unused parameter warning */
    int flag = *(int *)flags;
    int ord = flag%3;
    int type = flag-ord;

    double Yb1   = 0;
    double Yb2   = 0;
    double Ytau1 = 0;
    double Ytau2 = 0;
    double g1 = y[0];
    double g2 = y[1];
    double g3 = y[2];
    double Yt = y[3];

    if( type == 0 ) {//type I
        Yb2   = y[4];
        Ytau2 = y[5];
    }
    else if( type == 3 ) {//type II
        Yb1   = y[4];
        Ytau1 = y[5];
    }
    else if( type == 6 ) {//type X
        Yb2   = y[4];
        Ytau1 = y[5];
    }
    else if( type == 9 ) {//type Y
        Yb1   = y[4];
        Ytau2 = y[5];
    }
    else if( type == 12 ) {//inert
        Yt = 0.;
    }
    else {
        throw std::runtime_error("RGEsZ2: type should only be any of 0, 3, 6, 9, 12");
    }

    double m11_2=y[6];
    double m22_2=y[7];
    double m12_2=y[8];
    double la1=y[9];
    double la2=y[10];
    double la3=y[11];
    double la4=y[12];
    double la5=y[13];

    double pi = M_PI;

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

int JacobianZ2 (double t, const double y[], double *dfdy, double dfdt[], void *order)
{
    return 0;
}

int RGEcheckZ2(const double InitialValues[], const double t1, const double Rpeps, const double tNLOuni)
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

double GeneralTHDMZ2Runner::RGEGeneralTHDMZ2Runner(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order)
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
    gsl_odeiv2_system RGEsystem = {RGEsZ2, JacobianZ2, NumberOfRGEs, &order};

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
        if(RGEcheckZ2(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}

void GeneralTHDMZ2Runner::runGeneralTHDMZ2parameters()
{
    std::string modelflag = myGTHDMZ2->getZ2ModelType();
    std::string RGEorder  = myGTHDMZ2->getRGEorderflag();
    int flag;
    //flag will be used to transport information about model and RGEorder to the Runner:
    //flag=3*(0 for type I, 1 for type II, 2 for type X, 3 for type Y, 4 for inert) + (0 for LO, 1 for approxNLO and 2 for NLO)
    if( RGEorder == "LO" )
        flag = 0;
    else if( RGEorder == "approxNLO" )
        flag = 1;
    else if( RGEorder == "NLO" )
        flag = 2;
    else {
        throw std::runtime_error("GeneralTHDMZ2Runner: RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
    }

    double g1_at_MZ = sqrt(4.0*M_PI*Ale/cW2);
    double g2_at_MZ = sqrt(4.0*M_PI*Ale/(1-cW2));
    double g3_at_MZ = sqrt(4.0*M_PI*Als);

    double Ytop_at_MZ     = (sqrt(2.0)*myGTHDMZ2->getQuarks(QCD::TOP).getMass())/(vev*sinb);
    double Ybottom1_at_MZ = 0.0;
    double Ybottom2_at_MZ = 0.0;
    double Ytau1_at_MZ    = 0.0;
    double Ytau2_at_MZ    = 0.0;

    /*link these to the SM values*/
    double Mb_at_MZ   = 2.96;//GeV
    double Mtau_at_MZ = 1.75;//GeV

    if( modelflag == "type1" ) {
        Ybottom2_at_MZ = (sqrt(2.0)*Mb_at_MZ)/(vev*sinb);
        Ytau2_at_MZ    = (sqrt(2.0)*Mtau_at_MZ)/(vev*sinb);
    }
    else if( modelflag == "type2" ) {
        Ybottom1_at_MZ = (sqrt(2.0)*Mb_at_MZ)/(vev*cosb);
        Ytau1_at_MZ    = (sqrt(2.0)*Mtau_at_MZ)/(vev*cosb);

        flag += 3;
    }
    else if( modelflag == "typeX" ) {
        Ybottom2_at_MZ = (sqrt(2.0)*Mb_at_MZ)/(vev*sinb);
        Ytau1_at_MZ    = (sqrt(2.0)*Mtau_at_MZ)/(vev*cosb);

        flag += 6;
    }
    else if( modelflag == "typeY" ) {
        Ybottom1_at_MZ = (sqrt(2.0)*Mb_at_MZ)/(vev*cosb);
        Ytau2_at_MZ    = (sqrt(2.0)*Mtau_at_MZ)/(vev*sinb);
        flag += 9;
    }
    else if( modelflag == "inert" ) {
        Ytop_at_MZ = 0.0;

        flag += 12;
    }
    else {
        throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\", \"typeY\", \"inert\"");
    }

    double lambda1_at_MZ = myGTHDMZ2->getlambda1_Z2();
    double lambda2_at_MZ = myGTHDMZ2->getlambda2_Z2();
    double lambda3_at_MZ = myGTHDMZ2->getlambda3_Z2();
    double lambda4_at_MZ = myGTHDMZ2->getlambda4_Z2();
    double lambda5_at_MZ = myGTHDMZ2->getlambda5_Z2();

    double lambda345 = lambda3_at_MZ + lambda4_at_MZ + lambda5_at_MZ;

    double m12_2_at_MZ = m12_2;
    double m11_2_at_MZ = tanb*m12_2_at_MZ - vev*vev*(cosb*cosb*lambda1_at_MZ - sinb*sinb*lambda345)/2.;
    double m22_2_at_MZ = m12_2_at_MZ/tanb - vev*vev*(sinb*sinb*lambda2_at_MZ - cosb*cosb*lambda345)/2.;

    if(fabs(Q_GTHDM-log10(MZ))<0.005)   //at MZ scale
    {
        Q_cutoff = log10(MZ);

        g1_at_Q = g1_at_MZ;
        g2_at_Q = g2_at_MZ;
        g3_at_Q = g3_at_MZ;
        Ytop_at_Q = Ytop_at_MZ;
        Ybottom1_at_Q = Ybottom1_at_MZ;
        Ybottom2_at_Q = Ybottom2_at_MZ;
        Ytau1_at_Q = Ytau1_at_MZ;
        Ytau2_at_Q = Ytau2_at_MZ;
        m11_2_at_Q = m11_2_at_MZ;
        m22_2_at_Q = m22_2_at_MZ;
        m12_2_at_Q = m12_2_at_MZ;
        lambda1_at_Q = lambda1_at_MZ;
        lambda2_at_Q = lambda2_at_MZ;
        lambda3_at_Q = lambda3_at_MZ;
        lambda4_at_Q = lambda4_at_MZ;
        lambda5_at_Q = lambda5_at_MZ;
    }
    else   //at some other scale
    {
        double InitVals[14];
        InitVals[0]  = g1_at_MZ;
        InitVals[1]  = g2_at_MZ;
        InitVals[2]  = g3_at_MZ;
        InitVals[3]  = Ytop_at_MZ;
        InitVals[4]  = Ybottom1_at_MZ+Ybottom2_at_MZ;
        InitVals[5]  = Ytau1_at_MZ+Ytau2_at_MZ;
        InitVals[6]  = m11_2_at_MZ;
        InitVals[7]  = m22_2_at_MZ;
        InitVals[8]  = m12_2_at_MZ;
        InitVals[9]  = lambda1_at_MZ;
        InitVals[10] = lambda2_at_MZ;
        InitVals[11] = lambda3_at_MZ;
        InitVals[12] = lambda4_at_MZ;
        InitVals[13] = lambda5_at_MZ;

        // Running up to Q_cutoff <= Q_GTHDM
        Q_cutoff = RGEGeneralTHDMZ2Runner(InitVals, 14, log10(MZ), Q_GTHDM, flag);

        g1_at_Q = InitVals[0];
        g2_at_Q = InitVals[1];
        g3_at_Q = InitVals[2];
        Ytop_at_Q = InitVals[3];
        Ybottom1_at_Q = 0.0;
        Ybottom2_at_Q = 0.0;
        Ytau1_at_Q = 0.0;
        Ytau2_at_Q = 0.0;

        if( modelflag == "type1" ) {
            Ybottom2_at_Q = InitVals[4];
            Ytau2_at_Q    = InitVals[5];
        }
        else if( modelflag == "type2" ) {
            Ybottom1_at_Q = InitVals[4];
            Ytau1_at_Q    = InitVals[5];
        }
        else if( modelflag == "typeX" ) {
            Ybottom2_at_Q = InitVals[4];
            Ytau1_at_Q    = InitVals[5];
        }
        else if( modelflag == "typeY" ) {
            Ybottom1_at_Q = InitVals[4];
            Ytau2_at_Q    = InitVals[5];
        }
        else if( modelflag == "inert" ) {
            Ytop_at_Q = 0.0;
        }
        else {
            throw std::runtime_error("modelflag can be only any of \"type1\", \"type2\", \"typeX\" or \"typeY\"");
        }

        m11_2_at_Q = InitVals[6];
        m22_2_at_Q = InitVals[7];
        m12_2_at_Q = InitVals[8];
        lambda1_at_Q = InitVals[9];
        lambda2_at_Q = InitVals[10];
        lambda3_at_Q = InitVals[11];
        lambda4_at_Q = InitVals[12];
        lambda5_at_Q = InitVals[13];
    }
}


// Wavefunction renormalization contributions to unitarity, from Grinstein:2015rtl
void GeneralTHDMZ2Runner::computeWFR_Z2()
{
    double WFRcomb1a = 0.0;
    double WFRcomb1b = 0.0;
    double WFRcomb2a = 0.0;
    double WFRcomb3a = 0.0;
    double WFRcomb3b = 0.0;
    double WFRcomb4a = 0.0;

    double B000mh     = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_0_mHh2(MZ2,mHl2).real();
    double B000mH     = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_0_mHh2(MZ2,mHh2).real();
    double B00mHpmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double B00mHpmH   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double B00mAmh    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_mA2_mHl2(MZ2,mA2,mHl2).real();
    double B00mAmH    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_0_mA2_mHh2(MZ2,mA2,mHh2).real();
    double B0mh00     = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_0_0(MZ2,mHl2).real();
    double B0mh0mHp   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_0_mHp2(MZ2,mHl2,mHp2).real();
    double B0mh0mA    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_0_mA2(MZ2,mHl2,mA2).real();
    double B0mhmhmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_mHl2_mHl2(MZ2,mHl2).real();
    double B0mhmHmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_mHh2_mHl2(MZ2,mHl2,mHh2).real();
    double B0mhmHmH   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_mHh2_mHh2(MZ2,mHl2,mHh2).real();
    double B0mhmHpmHp = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_mHp2_mHp2(MZ2,mHl2,mHp2).real();
    double B0mhmAmA   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHl2_mA2_mA2(MZ2,mHl2,mA2).real();
    double B0mH00     = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_0_0(MZ2,mHh2).real();
    double B0mH0mHp   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_0_mHp2(MZ2,mHh2,mHp2).real();
    double B0mH0mA    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_0_mA2(MZ2,mHh2,mA2).real();
    double B0mHmhmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_mHl2_mHl2(MZ2,mHh2,mHl2).real();
    double B0mHmHmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_mHh2_mHl2(MZ2,mHh2,mHl2).real();
    double B0mHmHmH   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_mHh2_mHh2(MZ2,mHh2).real();
    double B0mHmHpmHp = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_mHp2_mHp2(MZ2,mHh2,mHp2).real();
    double B0mHmAmA   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHh2_mA2_mA2(MZ2,mHh2,mA2).real();
    double B0mHp0mh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHp2_0_mHl2(MZ2,mHp2,mHl2).real();
    double B0mHp0mH   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHp2_0_mHh2(MZ2,mHp2,mHh2).real();
    double B0mHpmHpmh = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHp2_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double B0mHpmHpmH = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mHp2_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double B0mA0mh    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mA2_0_mHl2(MZ2,mA2,mHl2).real();
    double B0mA0mH    = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mA2_0_mHh2(MZ2,mA2,mHh2).real();
    double B0mAmAmh   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mA2_mA2_mHl2(MZ2,mA2,mHl2).real();
    double B0mAmAmH   = myGTHDMZ2->getMyGTHDMCache()->B0_MZ2_mA2_mA2_mHh2(MZ2,mA2,mHh2).real();

    double ddpB000mh     = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_0_mHl2(MZ2,mHl2).real();
    double ddpB000mH     = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_0_mHh2(MZ2,mHh2).real();
    double ddpB00mHpmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB00mHpmH   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB00mHpmA   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_mHp2_mA2(MZ2,mHp2,mA2).real();
    double ddpB00mAmh    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_mA2_mHl2(MZ2,mA2,mHl2).real();
    double ddpB00mAmH    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_0_mA2_mHh2(MZ2,mA2,mHh2).real();
    double ddpB0mh00     = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_0_0(MZ2,mHl2).real();
    double ddpB0mh0mHp   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_0_mHp2(MZ2,mHl2,mHp2).real();
    double ddpB0mh0mA    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_0_mA2(MZ2,mHl2,mA2).real();
    double ddpB0mhmhmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_mHl2_mHl2(MZ2,mHl2).real();
    double ddpB0mhmHmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_mHh2_mHl2(MZ2,mHl2,mHh2).real();
    double ddpB0mhmHmH   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_mHh2_mHh2(MZ2,mHl2,mHh2).real();
    double ddpB0mhmHpmHp = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_mHp2_mHp2(MZ2,mHl2,mHp2).real();
    double ddpB0mhmAmA   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHl2_mA2_mA2(MZ2,mHl2,mA2).real();
    double ddpB0mH00     = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_0_0(MZ2,mHh2).real();
    double ddpB0mH0mHp   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_0_mHp2(MZ2,mHh2,mHp2).real();
    double ddpB0mH0mA    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_0_mA2(MZ2,mHh2,mA2).real();
    double ddpB0mHmhmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_mHl2_mHl2(MZ2,mHh2,mHl2).real();
    double ddpB0mHmHmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_mHh2_mHl2(MZ2,mHh2,mHl2).real();
    double ddpB0mHmHmH   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_mHh2_mHh2(MZ2,mHh2).real();
    double ddpB0mHmHpmHp = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_mHp2_mHp2(MZ2,mHh2,mHp2).real();
    double ddpB0mHmAmA   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHh2_mA2_mA2(MZ2,mHh2,mA2).real();
    double ddpB0mHp0mh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHp2_0_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB0mHp0mH   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHp2_0_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB0mHp0mA   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHp2_0_mA2(MZ2,mHp2,mA2).real();
    double ddpB0mHpmHpmh = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHp2_mHp2_mHl2(MZ2,mHp2,mHl2).real();
    double ddpB0mHpmHpmH = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mHp2_mHp2_mHh2(MZ2,mHp2,mHh2).real();
    double ddpB0mA0mh    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mA2_0_mHl2(MZ2,mA2,mHl2).real();
    double ddpB0mA0mH    = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mA2_0_mHh2(MZ2,mA2,mHh2).real();
    double ddpB0mA0mHp   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mA2_0_mHp2(MZ2,mA2,mHp2).real();
    double ddpB0mAmAmh   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mA2_mA2_mHl2(MZ2,mA2,mHl2).real();
    double ddpB0mAmAmH   = myGTHDMZ2->getMyGTHDMCache()->B0p_MZ2_mA2_mA2_mHh2(MZ2,mA2,mHh2).real();

    WFRcomb1a = 3.0*mHl2*mHl2*cosb*cosb*sin(bma)*sin(bma) * ddpB000mh
        + 3.0*mHh2*mHh2*cos(bma)*cos(bma)*cosb*cosb * ddpB000mH
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mHpmh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*cosb*cosb * ddpB00mHpmH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*cosb*cosb * ddpB00mHpmA
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mAmh
        + (mA2-mHh2)*(mA2-mHh2)*cosb*cosb*sin(bma)*sin(bma) * ddpB00mAmH
        + 1.5*mHl2*mHl2*sin(alpha)*sin(alpha)*sin(bma)*sin(bma) * ddpB0mh00
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mHp
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mA
        + 9.0*sin(alpha)*sin(alpha)*pow(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
                                        +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + sin(alpha)*sin(alpha)*pow((cos(alpha)/sinb+sin(alpha)/cosb)*(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-(3.0*m12_2)/(cosb*sinb))),2) * ddpB0mhmHmh
        + sin(alpha)*sin(alpha)*sin(bma)*sin(bma)*pow(-2.0*m12_2+(2.0*mHh2+mHl2-(3.0*m12_2)/(cosb*sinb))*sin(2.0*alpha),2)/(8.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + (sin(alpha)*sin(alpha)*pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                     +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb),2))/(64.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + (sin(alpha)*sin(alpha)*pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2))/(128.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 1.5*mHh2*mHh2*cos(alpha)*cos(alpha)*cos(bma)*cos(bma) * ddpB0mH00
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + (mA2-mHh2)*(mA2-mHh2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mA
        + cos(alpha)*cos(alpha)*cos(bma)*cos(bma)*pow(m12_2+cos(alpha)*sin(alpha)*(mHh2-3.0*m12_2/(cosb*sinb))+mHl2*sin(2.0*alpha),2)/(2.0*cosb*cosb*sinb*sinb) * ddpB0mHmhmh
        + cos(alpha)*cos(alpha)*sin(bma)*sin(bma)*pow(m12_2*cosb*sinb+0.5*sin(2.0*alpha)*(3.0*m12_2-(2.0*mHh2+mHl2)*cosb*sinb),2)/pow(cosb*sinb,4) * ddpB0mHmHmh
        + 9.0*cos(alpha)*cos(alpha)*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                                        +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + (cos(alpha)*cos(alpha)*pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)
                                     +2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta),2))/(256.0*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + (cos(alpha)*cos(alpha)*pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)
                                     -2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2))/(512.0*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mHp0mh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mHp0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mHp0mA
        + 2.0*pow((m12_2*cos(alpha+beta))/(sinb*sinb*cosb*cosb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)-(mHl2+2.0*mHp2)*sin(bma),2)*sinb*sinb * ddpB0mHpmHpmh
        + 2.0*pow(sinb*((mHh2+2.0*mHp2)*cos(bma)-mHh2*cos(2.0*beta)*sin(bma)/(cosb*sinb)-m12_2*sin(alpha+beta)/(sinb*sinb*cosb*cosb)),2) * ddpB0mHpmHpmH
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mA0mh
        + (mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mA0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2)/(64.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta) 
               + 3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2)/(256.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmH;

    WFRcomb1b = (mHl2*(mA2*(2.0*mHl2-3.0*mHp2)+mHl2*mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/(2.0*mA2*mHp2) * B000mh
        - (mHh2*(mA2*(2.0*mHh2-3.0*mHp2)+mHh2*mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/(2.0*mA2*mHp2) * B000mH
        + ((mHl2-mHp2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B00mHpmh
        + ((mHp2-mHh2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B00mHpmH
        + ((mHl2-mA2)*cos(bma)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B00mAmh
        + ((mHh2-mA2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                                +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B00mAmH
        + (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mh00
        + ((mHl2-mHp2)*(mHp2-mHh2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mh0mHp
        - ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mh0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHl2-mHh2)*pow(sinb*cosb,3)) * B0mhmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)
           *(4.0*m12_2*m12_2+2.0*m12_2*(mHl2-mHh2)*sin(2.0*alpha)
             -(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2-(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
               +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))/(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mhmHmh
        - 3.0*sin(2.0*alpha)*sin(bma)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mhmHmH
        + sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
             -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmHpmHp
        + sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(2.0*mA2+mHh2)*cos(bma)-2.0*mA2*cos(alpha+3.0*beta)
             -3.0*mHh2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmAmA
        - (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mH00
        + ((mHh2-mHp2)*(mHl2-mHp2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mH0mHp
        + ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mH0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)*(-4.0*m12_2*m12_2+2.0*m12_2*(mHh2-mHl2)*sin(2.0*alpha)
                                        +(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2
                                          -(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
                                          +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))
          /(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mHmHmh
        + 3.0*sin(bma)*sin(2.0*alpha)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmHmH
        - sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmHpmHp
        - sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mA2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mA2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mA2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmAmA
        - (mHl2*(mHl2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mh
        + (mHh2*(mHh2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mH
        + ((mHp2-mHl2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B0mHpmHpmh
        + ((mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B0mHpmHpmH
        + (mHl2*(mA2-mHl2)*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mh
        + ((mHh2-mA2)*mHh2*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mH
        + ((mA2-mHl2)*cos(bma)*((-2.0*mA2+mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B0mAmAmh
        +((mA2-mHh2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                               +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B0mAmAmH;

    WFRcomb2a = 1.5*mHl2*mHl2*sin(bma)*sin(bma) * ddpB000mh
        + 1.5*mHh2*mHh2*cos(bma)*cos(bma) * ddpB000mH
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB00mHpmh
        + (mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma) * ddpB00mHpmH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB00mHpmA
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB00mAmh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB00mAmH
        + 0.75*mHl2*mHl2*sin(bma)*sin(bma) * ddpB0mh00
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mh0mHp
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mh0mA
        + 9.0*pow(16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta)
                  -mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta)),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + 0.5*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mhmHmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(16.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(-4.0*m12_2+(3.0*mHl2+2.0*mHp2)*cosb*sinb),2)/(32.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(4.0*m12_2-(2.0*mA2+3.0*mHl2)*cosb*sinb),2)/(64.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 0.75*mHh2*mHh2*cos(bma)*cos(bma) * ddpB0mH00
        + (mHh2 - mHp2)*(mHh2 - mHp2)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mH0mA
        + 0.25*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mHmhmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(8.0*cosb*cosb*sinb*sinb) * ddpB0mHmHmh
        + 9.0*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                  +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)
           -16.0*m12_2*sin(alpha+beta),2)/(512*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
           +16.0*m12_2*sin(alpha+beta),2)/(1024*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mHp0mh
        + (mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma) * ddpB0mHp0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mHp0mA
        + pow((m12_2*cos(alpha+beta))/(cosb*cosb*sinb*sinb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)
              -(mHl2+2.0*mHp2)*sin(bma),2) * ddpB0mHpmHpmh
        + pow((mHh2+2.0*mHp2)*cos(bma)-(mHh2*cos(2.0*beta)*sin(bma))/(cosb*sinb)
              -(m12_2*sin(alpha+beta))/(cosb*cosb*sinb*sinb),2) * ddpB0mHpmHpmH
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mA0mh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mA0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*cosb*sinb
              +cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*cosb*sinb),2)/(128.0*pow(cosb*sinb,4)) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
              +16.0*m12_2*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mAmAmH;

    WFRcomb3a = 0.5*mHl2*mHl2*(3.0-cos(2.0*beta))*sin(bma)*sin(bma) * ddpB000mh
        + 0.5*mHh2*mHh2*(3.0-cos(2.0*beta))*cos(bma)*cos(bma) * ddpB000mH
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sinb*sinb * ddpB00mHpmh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*sinb*sinb * ddpB00mHpmH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB00mHpmA
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*cosb*cosb * ddpB00mAmh
        + (mA2-mHh2)*(mA2-mHh2)*cosb*cosb*sin(bma)*sin(bma) * ddpB00mAmH
        + 1.5*mHl2*mHl2*sin(alpha)*sin(alpha)*sin(bma)*sin(bma) * ddpB0mh00
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mHp
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sin(alpha)*sin(alpha) * ddpB0mh0mA
        + 9.0*sin(alpha)*sin(alpha)*pow(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
                                        +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + sin(alpha)*sin(alpha)*pow((cos(alpha)/sinb+sin(alpha)/cosb)*(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-(3.0*m12_2)/(cosb*sinb))),2) * ddpB0mhmHmh
        + sin(alpha)*sin(alpha)*sin(bma)*sin(bma)*pow(-2.0*m12_2+(2.0*mHh2+mHl2-(3.0*m12_2)/(cosb*sinb))*sin(2.0*alpha),2)/(8.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + (sin(alpha)*sin(alpha)*pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                     +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb),2))/(64.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + (sin(alpha)*sin(alpha)*pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2))/(128.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 1.5*mHh2*mHh2*cos(alpha)*cos(alpha)*cos(bma)*cos(bma) * ddpB0mH00
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + (mA2-mHh2)*(mA2-mHh2)*cos(alpha)*cos(alpha)*sin(bma)*sin(bma) * ddpB0mH0mA
        + cos(alpha)*cos(alpha)*cos(bma)*cos(bma)*pow(m12_2+cos(alpha)*sin(alpha)*(mHh2-3.0*m12_2/(cosb*sinb))+mHl2*sin(2.0*alpha),2)/(2.0*cosb*cosb*sinb*sinb) * ddpB0mHmhmh
        + cos(alpha)*cos(alpha)*sin(bma)*sin(bma)*pow(m12_2*cosb*sinb+0.5*sin(2.0*alpha)*(3.0*m12_2-(2.0*mHh2+mHl2)*cosb*sinb),2)/pow(cosb*sinb,4) * ddpB0mHmHmh
        + 9.0*cos(alpha)*cos(alpha)*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                                        +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + (cos(alpha)*cos(alpha)*pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)
                                     +2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta),2))/(256.0*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + (cos(alpha)*cos(alpha)*pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)
                                     -2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2))/(512.0*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 2.0*(mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma)*cosb*cosb * ddpB0mHp0mh
        + 2.0*(mHh2-mHp2)*(mHh2-mHp2)*sin(bma)*sin(bma)*cosb*cosb * ddpB0mHp0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*cosb*cosb * ddpB0mHp0mA
        + 2.0*cosb*cosb*pow((m12_2*cos(alpha+beta))/(cosb*cosb*sinb*sinb)-(mHl2*cos(bma)*cos(2.0*beta))/(cosb*sinb)-(mHl2+2.0*mHp2)*sin(bma),2) * ddpB0mHpmHpmh
        + 2.0*cosb*cosb*pow((mHh2+2.0*mHp2)*cos(bma)-(mHh2*cos(2.0*beta)*sin(bma))/(cosb*sinb)-(m12_2*sin(alpha+beta))/(cosb*cosb*sinb*sinb),2) * ddpB0mHpmHpmH
        + (mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma)*sinb*sinb * ddpB0mA0mh
        + (mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma)*sinb*sinb * ddpB0mA0mH
        + 2.0*(mA2-mHp2)*(mA2-mHp2)*sinb*sinb * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb),2)/(64.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta) 
               + 3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta),2)/(256.0*pow(cosb,4)*sinb*sinb) * ddpB0mAmAmH;

    WFRcomb3b = ((mHl2*mHl2*mHp2+mA2*(-2.0*mHl2*mHl2+mHl2*mHp2))*sin(2.0*bma)*sinb*cosb)/(mA2*mHp2) * B000mh
        + ((-mHh2*mHh2*mHp2+mA2*(2.0*mHh2*mHh2-mHh2*mHp2))*sin(2.0*bma)*sinb*cosb)/(mA2*mHp2) * B000mH
        + ((mHp2-mHl2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                 +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb)))/(2.0*mHp2*sinb*cosb) * B00mHpmh
        + ((mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
                                 -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta)))/(4.0*mHp2*sinb*cosb) * B00mHpmH
        + ((mHl2-mA2)*cos(bma)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B00mAmh
        + ((mHh2-mA2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                                +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B00mAmH
        + (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mh00
        + ((mHl2-mHp2)*(mHp2-mHh2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mh0mHp
        - ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mh0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHl2-mHh2)*pow(sinb*cosb,3)) * B0mhmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)
           *(4.0*m12_2*m12_2+2.0*m12_2*(mHl2-mHh2)*sin(2.0*alpha)
             -(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2-(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
               +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))/(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mhmHmh
        - 3.0*sin(2.0*alpha)*sin(bma)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mhmHmH
        + sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
             -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmHpmHp
        + sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(2.0*mA2+mHh2)*cos(bma)-2.0*mA2*cos(alpha+3.0*beta)
             -3.0*mHh2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mhmAmA
        - (3.0*mHh2*mHl2*sin(2.0*alpha)*sin(2.0*bma))/(4.0*(mHh2-mHl2)) * B0mH00
        + ((mHh2-mHp2)*(mHl2-mHp2)*sin(2.0*alpha)*sin(2.0*bma))/(mHh2-mHl2) * B0mH0mHp
        + ((mA2-mHh2)*(mA2-mHl2)*sin(2.0*alpha)*sin(2.0*bma))/(2.0*(mHh2-mHl2)) * B0mH0mA
        + 3.0*cos(bma)*sin(2.0*alpha)*(m12_2+cos(alpha)*sin(alpha)*(mHh2-(3.0*m12_2)/(sinb*cosb))+mHl2*sin(2.0*alpha))
          *(-mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta))
            +16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmhmh
        + (sin(2.0*bma)*sin(2.0*alpha)*(-4.0*m12_2*m12_2+2.0*m12_2*(mHh2-mHl2)*sin(2.0*alpha)
                                        +(2.0*mHh2*mHh2+5.0*mHh2*mHl2+2.0*mHl2*mHl2
                                          -(9.0*m12_2*(mHh2+mHl2))/(sinb*cosb)
                                          +(9.0*m12_2*m12_2)/(sinb*sinb*cosb*cosb))*sin(2.0*alpha)*sin(2.0*alpha)))
          /(8.0*(mHh2-mHl2)*sinb*sinb*cosb*cosb) * B0mHmHmh
        + 3.0*sin(bma)*sin(2.0*alpha)*(m12_2-cos(alpha)*sin(alpha)*(2.0*mHh2+mHl2-(3.0*m12_2)/(sinb*cosb)))
          *(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
            +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta))/(32.0*(mHh2-mHl2)*pow(sinb*cosb,3)) * B0mHmHmH
        - sin(2.0*alpha)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))
          *((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(128.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmHpmHp
        - sin(2.0*alpha)*((mHl2-2.0*mA2)*cos(alpha-3.0*beta)*2.0*sinb*cosb+cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mA2)*2.0*sinb*cosb))
          *((mHh2-2.0*mA2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mA2)*cos(bma)-3.0*mHh2*cos(alpha+3.0*beta)
            -2.0*mA2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(256.0*(mHh2-mHl2)*pow(cosb*sinb,4)) * B0mHmAmA
        + (mHl2*(mHl2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mh
        - (mHh2*(mHh2-mHp2)*sin(2.0*bma)*2.0*sinb*cosb)/mHp2 * B0mHp0mH
        + (mHl2-mHp2)*cos(bma)*((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(3.0*mHl2+2.0*mHp2)*2.0*sinb*cosb))/(2.0*mHp2*sinb*cosb) * B0mHpmHpmh
        - (mHh2-mHp2)*sin(bma)*((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)
                                -3.0*mHh2*cos(alpha+3.0*beta)-2.0*mHp2*cos(alpha+3.0*beta)-16.0*m12_2*sin(alpha+beta))/(4.0*mHp2*sinb*cosb) * B0mHpmHpmH
        + (mHl2*(mA2-mHl2)*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mh
        + ((mHh2-mA2)*mHh2*sin(bma)*cos(bma)*2.0*sinb*cosb)/mA2 * B0mA0mH
        + ((mA2-mHl2)*cos(bma)*((-2.0*mA2+mHl2)*cos(alpha-3.0*beta)*2.0*sinb*cosb
                                +cos(alpha+beta)*(-8.0*m12_2+(2.0*mA2+3.0*mHl2)*2.0*sinb*cosb)))/(4.0*mA2*sinb*cosb) * B0mAmAmh
        +((mA2-mHh2)*sin(bma)*((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+2.0*mA2*cos(alpha+3.0*beta)
                               +3.0*mHh2*cos(alpha+3.0*beta)+16.0*m12_2*sin(alpha+beta)))/(8.0*mA2*sinb*cosb) * B0mAmAmH;

    WFRcomb4a = 0.5*mHl2*mHl2*sin(bma)*sin(bma) * ddpB000mh
        + 0.5*mHh2*mHh2*cos(bma)*cos(bma) * ddpB000mH
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB00mAmh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB00mAmH
        + 0.75*mHl2*mHl2*sin(bma)*sin(bma) * ddpB0mh00
        + (mHl2-mHp2)*(mHl2-mHp2)*cos(bma)*cos(bma) * ddpB0mh0mHp
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mh0mA
        + 9.0*pow(16.0*m12_2*cos(bma)*cos(bma)*cos(alpha+beta)
                  -mHl2*(3.0*sin(bma)+sin(3.0*bma)+sin(3.0*alpha+beta)+3.0*sin(alpha+3.0*beta)),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mhmhmh
        + 0.5*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mhmHmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(16.0*cosb*cosb*sinb*sinb) * ddpB0mhmHmH
        + pow((mHl2-2.0*mHp2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(-4.0*m12_2+(3.0*mHl2+2.0*mHp2)*cosb*sinb),2)/(32.0*pow(cosb*sinb,4)) * ddpB0mhmHpmHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*cosb*sinb
              +cos(alpha+beta)*(4.0*m12_2-(2.0*mA2+3.0*mHl2)*cosb*sinb),2)/(64.0*pow(cosb*sinb,4)) * ddpB0mhmAmA
        + 0.75*mHh2*mHh2*cos(bma)*cos(bma) * ddpB0mH00
        + (mHh2 - mHp2)*(mHh2 - mHp2)*sin(bma)*sin(bma) * ddpB0mH0mHp
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mH0mA
        + 0.25*pow(cos(alpha)/sinb + sin(alpha)/cosb,2)
          *pow(m12_2+cos(alpha)*sin(alpha)*(mHh2+2.0*mHl2-3.0*m12_2/(cosb*sinb)),2) * ddpB0mHmhmh
        + (pow(-2.0*m12_2+(2.0*mHh2+mHl2-3.0*m12_2/(cosb*sinb))*sin(2.0*alpha),2)
           *sin(bma)*sin(bma))/(8.0*cosb*cosb*sinb*sinb) * ddpB0mHmHmh
        + 9.0*pow(mHh2*(-3.0*cos(bma)+cos(3.0*bma)-cos(3.0*alpha+beta)+3.0*cos(alpha+3.0*beta))
                  +16.0*m12_2*sin(bma)*sin(bma)*sin(alpha+beta),2)/(1024.0*pow(cosb*sinb,4)) * ddpB0mHmHmH
        + pow((mHh2-2.0*mHp2)*cos(alpha-5.0*beta)+2.0*(mHh2+2.0*mHp2)*cos(bma)-(3.0*mHh2+2.0*mHp2)*cos(alpha+3.0*beta)
           -16.0*m12_2*sin(alpha+beta),2)/(512*pow(cosb*sinb,4)) * ddpB0mHmHpmHp
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
           +16.0*m12_2*sin(alpha+beta),2)/(1024*pow(cosb*sinb,4)) * ddpB0mHmAmA
        + 0.5*(mA2-mHl2)*(mA2-mHl2)*cos(bma)*cos(bma) * ddpB0mA0mh
        + 0.5*(mA2-mHh2)*(mA2-mHh2)*sin(bma)*sin(bma) * ddpB0mA0mH
        + (mA2-mHp2)*(mA2-mHp2) * ddpB0mA0mHp
        + pow((2.0*mA2-mHl2)*cos(alpha-3.0*beta)*2.0*cosb*sinb
              +cos(alpha+beta)*(8.0*m12_2-(2.0*mA2+3.0*mHl2)*2.0*cosb*sinb),2)/(128.0*pow(cosb*sinb,4)) * ddpB0mAmAmh
        + pow((2.0*mA2-mHh2)*cos(alpha-5.0*beta)-2.0*(2.0*mA2+mHh2)*cos(bma)+(2.0*mA2+3.0*mHh2)*cos(alpha+3.0*beta)
              +16.0*m12_2*sin(alpha+beta),2)/(512.0*pow(cosb*sinb,4)) * ddpB0mAmAmH;

    WFRcomb1=-(WFRcomb1a+WFRcomb1b)/(vev*vev);
    WFRcomb2=-WFRcomb2a/(vev*vev);
    WFRcomb3=-(WFRcomb3a+WFRcomb3b)/(vev*vev);
    WFRcomb4=-WFRcomb4a/(vev*vev);
}


gslpp::matrix<double> GeneralTHDMZ2Runner::getGTHDMZ2_at_Q()
{
    vev = myGTHDMZ2->v();
    cW2 = myGTHDMZ2->c02();
    Ale = myGTHDMZ2->getAle();
    Als = myGTHDMZ2->getAlsMz();
    MZ  = myGTHDMZ2->getMz();
    MZ2 = MZ*MZ;

    m12_2 = myGTHDMZ2->getm12sq_Z2();
    mHl2  = myGTHDMZ2->getmH1sq();
    mHh2  = myGTHDMZ2->getmH2sq();
    mA2   = myGTHDMZ2->getmH3sq();
    mHp2  = myGTHDMZ2->getmHp2();

    beta  = myGTHDMZ2->getbeta_Z2();
    tanb  = myGTHDMZ2->gettanb_Z2();
    sinb  = myGTHDMZ2->getsinb_Z2();
    cosb  = myGTHDMZ2->getcosb_Z2();
    bma   = myGTHDMZ2->getbma_Z2();
    alpha = beta-bma;

    Q_GTHDM     = myGTHDMZ2->getQ_GTHDM();
    Rpeps       = myGTHDMZ2->getRpepsGTHDM();
    NLOuniscale = myGTHDMZ2->getNLOuniscaleGTHDM();

    runGeneralTHDMZ2parameters();

    // Wavefunction renormalization contributions to unitarity, from Grinstein:2015rtl
    computeWFR_Z2();

    // 1st row: (lambda1, lambda2, lambda3, lambda4, lambda5)
    // 2nd row: (Ytop, Ybottom1, Ybottom2, Ytau1, Ytau2)
    // 3rd row: (WFRcomb1, WFRcomb2, WFRcomb3, WFRcomb4, 0.0)
    myZ2_at_Q.assign(0, 0, lambda1_at_Q);
    myZ2_at_Q.assign(0, 1, lambda2_at_Q);
    myZ2_at_Q.assign(0, 2, lambda3_at_Q);
    myZ2_at_Q.assign(0, 3, lambda4_at_Q);
    myZ2_at_Q.assign(0, 4, lambda5_at_Q);
    myZ2_at_Q.assign(1, 0, Ytop_at_Q);
    myZ2_at_Q.assign(1, 1, Ybottom1_at_Q);
    myZ2_at_Q.assign(1, 2, Ybottom2_at_Q);
    myZ2_at_Q.assign(1, 3, Ytau1_at_Q);
    myZ2_at_Q.assign(1, 4, Ytau2_at_Q);
    myZ2_at_Q.assign(2, 0, WFRcomb1);
    myZ2_at_Q.assign(2, 1, WFRcomb2);
    myZ2_at_Q.assign(2, 2, WFRcomb3);
    myZ2_at_Q.assign(2, 3, WFRcomb4);

    return (myZ2_at_Q);
}