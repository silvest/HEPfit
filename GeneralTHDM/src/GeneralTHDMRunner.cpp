/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMRunner.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

GeneralTHDMRunner::GeneralTHDMRunner(const StandardModel& SM_i) : myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

GeneralTHDMRunner::~GeneralTHDMRunner()
{};

int RGEsGTHDM(double t, const double y[], double beta[], void *flags)
{
    (void)(t); /* avoid unused parameter warning */
    int ord = *(int *)flags;
    double g1=y[0];
    double g2=y[1];
    double g3=y[2];
    double v1=y[3];
    double v2=y[4];
    double etaU1=y[5];
    double etaU2=y[6];
    double etaD1=y[7];
    double etaD2=y[8];
    double etaL1=y[9];
    double etaL2=y[10];
    double m11sq=y[11];
    double m12sq=y[12];
    double m22sq=y[13];
    double la1=y[14];
    double la2=y[15];
    double la3=y[16];
    double la4=y[17];
    double la5=y[18];
    double la6=y[19];
    double la7=y[20];

    double pi=M_PI;

    //beta_g1
    beta[0] = 7.*g1*g1*g1;
    //beta_g2
    beta[1] = -3.*g2*g2*g2;
    //beta_g3
    beta[2] = -7.*g3*g3*g3;
    //beta_v1
    beta[3] = (-3.*etaD1*etaD1-etaL1*etaL1-3.*etaU1*etaU1+0.75*g1*g1+2.25*g2*g2)*v1
             +(-3.*etaD1*etaD2-etaL1*etaL2-3.*etaU1*etaU2)*v2;
    //beta_v2
    beta[4] = (-3.*etaD1*etaD2-etaL1*etaL2-3.*etaU1*etaU2)*v1
             +(-3.*etaD2*etaD2-etaL2*etaL2-3.*etaU2*etaU2+0.75*g1*g1+2.25*g2*g2)*v2;
    //beta_etaU1
    beta[5] = 1.5*etaD1*etaD1*etaU1+0.5*etaD2*etaD2*etaU1+etaL1*etaL1*etaU1+4.5*etaU1*etaU1*etaU1
             +etaD1*etaD2*etaU2+etaL1*etaL2*etaU2+4.5*etaU1*etaU2*etaU2
             -1.41667*etaU1*g1*g1-2.25*etaU1*g2*g2-8.*etaU1*g3*g3;
    //beta_etaU2
    beta[6] = etaD1*etaD2*etaU1+etaL1*etaL2*etaU1+0.5*etaD1*etaD1*etaU2
             +etaU2*(1.5*etaD2*etaD2+etaL2*etaL2+4.5*etaU1*etaU1+4.5*etaU2*etaU2
                     -1.41667*g1*g1-2.25*g2*g2-8.*g3*g3);
    //beta_etaD1
    beta[7] = 4.5*etaD1*etaD1*etaD1+etaD2*(etaL1*etaL2+etaU1*etaU2)
             +etaD1*(4.5*etaD2*etaD2+etaL1*etaL1+1.5*etaU1*etaU1+0.5*etaU2*etaU2
                     -0.416667*g1*g1-2.25*g2*g2-8.*g3*g3);
    //beta_etaD2
    beta[8] = 4.5*etaD1*etaD1*etaD2+etaD1*(etaL1*etaL2+etaU1*etaU2)
             +etaD2*(4.5*etaD2*etaD2+etaL2*etaL2+0.5*etaU1*etaU1+1.5*etaU2*etaU2
                     -0.416667*g1*g1-2.25*g2*g2-8.*g3*g3);
    //beta_etaL1
    beta[9] = 3.*etaD1*etaD1*etaL1+2.5*etaL1*etaL1*etaL1+3.*etaD1*etaD2*etaL2
             +3.*etaL2*etaU1*etaU2+etaL1*(2.5*etaL2*etaL2+3.*etaU1*etaU1-3.75*g1*g1-2.25*g2*g2);
    //beta_etaL2
    beta[10] = 3.*etaD1*etaD2*etaL1+3.*etaD2*etaD2*etaL2+2.5*etaL1*etaL1*etaL2
              +2.5*etaL2*etaL2*etaL2+3.*etaL1*etaU1*etaU2+3.*etaL2*etaU2*etaU2
              -3.75*etaL2*g1*g1-2.25*etaL2*g2*g2;
    //beta_m11sq
    beta[11] = 6.*etaD1*etaD1*m11sq+2.*etaL1*etaL1*m11sq+6.*etaU1*etaU1*m11sq
              -1.5*g1*g1*m11sq-4.5*g2*g2*m11sq+6.*la1*m11sq
              -6.*etaD1*etaD2*m12sq-2.*etaL1*etaL2*m12sq-6.*etaU1*etaU2*m12sq
              -12.*la6*m12sq+4.*la3*m22sq+2.*la4*m22sq;
    //beta_m12sq
    beta[12] = (3.*etaD1*etaD1+3.*etaD2*etaD2+etaL1*etaL1+etaL2*etaL2+3.*etaU1*etaU1
                +3.*etaU2*etaU2+2.*la3+4.*la4)*m12sq
              -3.*etaD1*etaD2*(m11sq+m22sq)-etaL1*etaL2*(m11sq+m22sq)
              -1.5*(4.*la6*m11sq+g1*g1*m12sq+3.*g2*g2*m12sq-4.*la5*m12sq
                    +4.*la7*m22sq+2.*etaU1*etaU2*(m11sq+m22sq));
    //beta_m22sq
    beta[13] = 4.*la3*m11sq+2.*la4*m11sq-6.*etaD1*etaD2*m12sq-2.*etaL1*etaL2*m12sq
              -6.*etaU1*etaU2*m12sq-12.*la7*m12sq+6.*etaD2*etaD2*m22sq+2.*etaL2*etaL2*m22sq
              +6.*etaU2*etaU2*m22sq-1.5*g1*g1*m22sq-4.5*g2*g2*m22sq+6.*la2*m22sq;
    //beta_la1
    beta[14] = -12.*etaD1*etaD1*etaD1*etaD1-4.*etaL1*etaL1*etaL1*etaL1-12.*etaU1*etaU1*etaU1*etaU1
               +0.75*g1*g1*g1*g1+1.5*g1*g1*g2*g2+2.25*g2*g2*g2*g2
               +12.*etaD1*etaD1*la1+4.*etaL1*etaL1*la1+12.*etaU1*etaU1*la1
               -3.*g1*g1*la1-9.*g2*g2*la1
               +12.*la1*la1+4.*la3*la3+4.*la3*la4+2.*la4*la4
               +12.*etaD1*etaD2*la6+4.*etaL1*etaL2*la6+12.*etaU1*etaU2*la6
               +24.*la6*la6+2.*la5*la5;
    //beta_la2
    beta[15] = -12.*etaD2*etaD2*etaD2*etaD2-4.*etaL2*etaL2*etaL2*etaL2-12.*etaU2*etaU2*etaU2*etaU2
               +0.75*g1*g1*g1*g1+1.5*g1*g1*g2*g2+2.25*g2*g2*g2*g2
               +12.*etaD2*etaD2*la2+4.*etaL2*etaL2*la2+12.*etaU2*etaU2*la2
               -3.*g1*g1*la2-9.*g2*g2*la2
               +12.*la2*la2+4.*la3*la3+4.*la3*la4+2.*la4*la4
               +12.*etaD1*etaD2*la7+4.*etaL1*etaL2*la7+12.*etaU1*etaU2*la7
               +24.*la7*la7+2.*la5*la5;
    //beta_la3
    beta[16] = -12.*etaD1*etaD1*etaD2*etaD2-4.*etaL1*etaL1*etaL2*etaL2-12.*etaD2*etaD2*etaU1*etaU1
               +24.*etaD1*etaD2*etaU1*etaU2-12.*etaD1*etaD1*etaU2*etaU2-12.*etaU1*etaU1*etaU2*etaU2
               +0.75*(g1*g1*g1*g1+3.*g2*g2*g2*g2+g1*g1*(-2.*g2*g2-4.*la3)-12.*g2*g2*la3)
               +2.*(la3*(3.*etaD1*etaD1+3.*etaD2*etaD2+etaL1*etaL1+etaL2*etaL2
                    +3.*(etaU1*etaU1+etaU2*etaU2+la1+la2)+2.*la3)
                    +(la1+la2)*la4+la4*la4)
               +(3.*etaD1*etaD2+etaL1*etaL2+3.*etaU1*etaU2)*(la6+la7)
               +la7*(3.*etaD1*etaD2+etaL1*etaL2+3.*etaU1*etaU2+8.*la6+4.*la7)
               +la6*(3.*etaD1*etaD2+etaL1*etaL2+3.*etaU1*etaU2+4.*la6+8.*la7)+2.*la5*la5;
    //beta_la4
    beta[17] = 12.*etaD2*etaD2*etaU1*etaU1-12.*etaU1*etaU1*etaU2*etaU2
              +3.*g1*g1*g2*g2+6.*etaD2*etaD2*la4+2.*etaL2*etaL2*la4+6.*etaU1*etaU1*la4
              +6.*etaU2*etaU2*la4-3.*g1*g1*la4-9.*g2*g2*la4
              +2.*la1*la4+2.*la2*la4+8.*la3*la4+4.*la4*la4
              +etaL1*etaL1*(-4.*etaL2*etaL2+2.*la4)
              +etaD1*etaD1*(-12.*etaD2*etaD2+12.*etaU2*etaU2+6.*la4)
              +8.*la5*la5+6.*etaU1*etaU2*la6+10.*la6*la6+6.*etaU1*etaU2*la7
              +4.*la6*la7+10.*la7*la7+etaL1*etaL2*(2.*la6+2.*la7)
              +etaD1*etaD2*(-24.*etaU1*etaU2+6.*la6+6.*la7);
    //beta_la5
    beta[18] = -12.*etaU1*etaU1*etaU2*etaU2+6.*etaD2*etaD2*la5+2.*etaL2*etaL2*la5
               +6.*etaU1*etaU1*la5+6.*etaU2*etaU2*la5-3.*g1*g1*la5-9.*g2*g2*la5
               +2.*la1*la5+2.*la2*la5+8.*la3*la5+12.*la4*la5
               +etaL1*etaL1*(-4.*etaL2*etaL2+2.*la5)
               +etaD1*etaD1*(-12.*etaD2*etaD2+6.*la5)
               +6.*etaU1*etaU2*la6+10.*la6*la6+6.*etaU1*etaU2*la7
               +4.*la6*la7+10.*la7*la7+etaL1*etaL2*(2.*la6+2.*la7)
               +etaD1*etaD2*(6.*la6+6.*la7);
    //beta_la6
    beta[19] = -12.*etaD1*etaD1*etaD1*etaD2-4.*etaL1*etaL1*etaL1*etaL2-12.*etaU1*etaU1*etaU1*etaU2
               +3.*etaU1*etaU2*la1+3.*etaU1*etaU2*la3+3.*etaU1*etaU2*la4
               +3.*etaU1*etaU2*la5+etaL1*etaL2*(la1+la3+la4+la5)
               +3.*etaD1*etaD2*(la1+la3+la4+la5)+9.*etaD1*etaD1*la6
               +3.*etaD2*etaD2*la6+3.*etaL1*etaL1*la6+etaL2*etaL2*la6
               +9.*etaU1*etaU1*la6+3.*etaU2*etaU2*la6-3.*g1*g1*la6-9.*g2*g2*la6
               +12.*la1*la6+6.*la3*la6+8.*la4*la6+10.*la5*la6
               +6.*la3*la7+4.*la4*la7+2.*la5*la7;
    //beta_la7
    beta[20] = -12.*etaU1*etaU2*etaU2*etaU2+3.*etaU1*etaU2*la2+3.*etaU1*etaU2*la3
               +3.*etaU1*etaU2*la4+3.*etaU1*etaU2*la5
               +etaL1*etaL2*(-4.*etaL2*etaL2+la2+la3+la4+la5)
               +etaD1*etaD2*(-12.*etaD2*etaD2+3.*la2+3.*la3+3.*la4+3.*la5)
               +6.*la3*la6+4.*la4*la6+2.*la5*la6+3.*etaD1*etaD1*la7
               +9.*etaD2*etaD2*la7+etaL1*etaL1*la7+3.*etaL2*etaL2*la7
               +3.*etaU1*etaU1*la7+9.*etaU2*etaU2*la7
               -3.*g1*g1*la7-9.*g2*g2*la7+12.*la2*la7
               +6.*la3*la7+8.*la4*la7+10.*la5*la7;

        if(ord==1){
            //beta_g1
            beta[0] += 0.;
            //beta_g2
            beta[1] += 0.;
            //beta_g3
            beta[2] += 0.;
            //beta_v1
            beta[3] += 0.;
            //beta_v2
            beta[4] += 0.;
            //beta_etaU1
            beta[5] += 0.;
            //beta_etaU2
            beta[6] += 0.;
            //beta_etaD1
            beta[7] += 0.;
            //beta_etaD2
            beta[8] += 0.;
            //beta_etaL1
            beta[9] += 0.;
            //beta_etaL2
            beta[10] += 0.;
            //beta_m11sq
            beta[11] += 0.;
            //beta_m12sq
            beta[12] += 0.;
            //beta_m22sq
            beta[13] += 0.;
            //beta_la1
            beta[14] += 0.;
            //beta_la2
            beta[15] += 0.;
            //beta_la3
            beta[16] += 0.;
            //beta_la4
            beta[17] += 0.;
            //beta_la5
            beta[18] += 0.;
            //beta_la6
            beta[19] += 0.;
            //beta_la7
            beta[20] += 0.;
        }

    return 0;
}

int JacobianGTHDM (double t, const double y[], double *dfdy, double dfdt[], void *order)
{
    return 0;
}

int RGEcheckGTHDM(const double InitialValues[], const double t1, const double Rpeps, const double tNLOuni)
{
    int check=0;

    //perturbativity of the Yukawa couplings
    for(int i=5;i<11;i++)
    {
        if(fabs(InitialValues[i])>sqrt(4.0*M_PI)) check=1;
    }

    //perturbativity of the quartic Higgs couplings
    for(int i=14;i<21;i++)
    {
        if(fabs(InitialValues[i])>4.0*M_PI) check=1;
    }

    double la1Q = InitialValues[14];
    double la2Q = InitialValues[15];
    double la3Q = InitialValues[16];
    double la4Q = InitialValues[17];
    double la5Q = InitialValues[18];
    double la6Q = InitialValues[19];
    double la7Q = InitialValues[20];

    //positivity checks
    if(la1Q<0.0) check=1;
    if(la2Q<0.0) check=1;
//    if(la3Q+sqrt(la1Q*la2Q)<0.0) check=1;
//    if(la3Q+2.0*la4Q+sqrt(la1Q*la2Q)<0.0) check=1;

    /////////////////
    //NLO unitarity//
    /////////////////

//    double pi=M_PI;
//    gslpp::matrix<gslpp::complex> Smatrix1(4,4,0.), Smatrix2(4,4,0.);
//    gslpp::matrix<gslpp::complex> Seigenvectors1(4,4,0.), Seigenvectors2(4,4,0.);
//    gslpp::vector<double> Seigenvalues1(4,0.), Seigenvalues2(4,0.);
//    gslpp::vector<gslpp::complex> unitarityeigenvalues(11,0.);

    if(t1>tNLOuni)
    {
    } //end of the if(t1>tNLOuni)
    return check;
}

double GeneralTHDMRunner::RGEGeneralTHDMRunner(double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double NLOuniscale)
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
    gsl_odeiv2_system RGEsystem = {RGEsGTHDM, JacobianGTHDM, NumberOfRGEs, &order};

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
        if(RGEcheckGTHDM(InitialValues,t1,Rpeps,tNLOuni) != 0) break;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    //Return the decadic log scale at which the evolution stopped
    return t1/log(10.0);
}
