/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
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
    //do something with *params

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

//    return GSL_SUCCESS;
    return 0;
}

int Jacobian (double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    return 0;
}

int RGEcheck(const double InitialValues[])
{
    int check=0;
    
    for(int i=9;i<14;i++)
    {
        if(fabs(InitialValues[i])>4.0*M_PI) check=1;
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
        if(RGEcheck(InitialValues) != 0) break;

        printf ("%.5e %.5e\n", t1, InitialValues[9]);
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    for(int i=0;i<14;i++) std::cout<<InitialValues[i]<<std::endl;

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}



lambda1atQ::lambda1atQ(const StandardModel& SM_i)
: Runner(SM_i)
{}

double lambda1atQ::computeThValue()
{

//    Ale=myTHDM->getAle();
//    GF=myTHDM->getGF();
//    Als=myTHDM->getAlsMz();
    double InitVals[14];
    InitVals[0]=0.35611;    //g1
//    double g2=sqrt(8. * mySUSY.getGF() / sqrt(2.)) * mySUSY.Mw_tree();
    InitVals[1]=0.64883;    //g2
    InitVals[2]=3.54491*sqrt(myTHDM->getAlsMz());   //g3
    InitVals[3]=sqrt(2.0)*myTHDM->getQuarks(QCD::TOP).getMass()/myTHDM->v()/myTHDM->getsinb(); //msbar mass?    //Yt
    InitVals[4]=sqrt(2.0)*myTHDM->getQuarks(QCD::BOTTOM).getMass()/myTHDM->v()/myTHDM->getcosb(); //msbar mass  //Yb
    InitVals[5]=sqrt(2.0)*myTHDM->getLeptons(StandardModel::TAU).getMass()/myTHDM->v()/myTHDM->getcosb(); //msbar mass  //Ytau
    InitVals[6]=mym11_2->computeThValue();  //m11_2
    InitVals[7]=mym22_2->computeThValue();  //m22_2
    InitVals[8]=myTHDM->getm12_2(); //m12_2
    InitVals[9]=mylambda1->computeThValue();    //lambda1
    InitVals[10]=mylambda2->computeThValue();    //lambda2
    InitVals[11]=mylambda3->computeThValue();    //lambda3
    InitVals[12]=mylambda4->computeThValue();    //lambda4
    InitVals[13]=mylambda5->computeThValue();    //lambda5
    double Q1=1.959;
    double Q2=5;

    double t_cutoff=RGERunner(InitVals, 14, Q1, Q2);

    return InitVals[9];
}

lambda2atQ::lambda2atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda2atQ::computeThValue()
{
    return 0.;
}

lambda3atQ::lambda3atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda3atQ::computeThValue()
{
    return 0.;
}

lambda4atQ::lambda4atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda4atQ::computeThValue()
{
    return 0.;
}

lambda5atQ::lambda5atQ(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda5atQ::computeThValue()
{
    return 0.;
}
