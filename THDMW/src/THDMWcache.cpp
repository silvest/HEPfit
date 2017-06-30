/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWcache.h"

THDMWcache::THDMWcache(const StandardModel& SM_i)
        : unitarityeigenvalues(11, 0.), NLOunitarityeigenvalues(11, 0.),
        myTHDMW(static_cast<const THDMW*> (&SM_i)), betaeigenvalues(11, 0.)
{
    myRunnerTHDMW=new RunnerTHDMW(SM_i);
}

THDMWcache::~THDMWcache()
{
  delete myRunnerTHDMW;
}
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
//int THDMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
//                          const int NumPar, const double params[]) const {
//    bool bCache;
//    for(int i=0; i<CacheSize; i++) {
//        bCache = true;
//        for(int j=0; j<NumPar; j++)
//            bCache &= (params[j] == cache[j][i].real());
//        if (bCache) return i;
//    }
//    return -1;
//}
//
//int THDMcache::CacheCheckReal(const double cache[][CacheSize], 
//                          const int NumPar, const double params[]) const {
//    bool bCache;
//    for(int i=0; i<CacheSize; i++) {
//        bCache = true;
//        for(int j=0; j<NumPar; j++)
//            bCache &= (params[j] == cache[j][i]);
//        if (bCache) return i;
//    }
//    return -1;
//}
//
//void THDMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
//                           const double params[], const gslpp::complex newResult) const {
//    // shift old parameters and result
//    for(int i=CacheSize-1; i>0; i--)
//        for(int j=0; j<NumPar+1; j++)
//            cache[j][i] = cache[j][i-1];
//
//    // store new parameters and result
//    for(int j=0; j<NumPar; j++) {
//        cache[j][0] = gslpp::complex(params[j], 0.0, false);
//        cache[NumPar][0] = newResult;
//    }
//}
//
//void THDMcache::CacheShiftReal(double cache[][CacheSize], const int NumPar,
//                           const double params[], const double newResult) const {
//    // shift old parameters and result
//    for(int i=CacheSize-1; i>0; i--)
//        for(int j=0; j<NumPar+1; j++)
//            cache[j][i] = cache[j][i-1];
//
//    // store new parameters and result
//    for(int j=0; j<NumPar; j++) {
//        cache[j][0] = params[j];
//        cache[NumPar][0] = newResult;
//    }
//}

void THDMWcache::runTHDMWparameters()
{
//    vev=myTHDM->v();
//    double cosb=myTHDM->getcosb();

    std::string RGEorder=myTHDMW->getRGEorderflag();
    //flag will be used to transport information about model and RGEorder to the Runner:
    //flag=0 for LO, 1 for approxNLO (and 2 for NLO - not implemented yet)
    int flag;
    if( RGEorder == "LO" ) flag=0;
    else if( RGEorder == "approxNLO" ) flag=1;
//    else if( RGEorder == "NLO" ) flag=2;
    else {
        throw std::runtime_error("RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
    }

    double lambda1_at_MZ=lambda1;
    double lambda2_at_MZ=lambda2;
    double lambda3_at_MZ=lambda3;
    double lambda4_at_MZ=lambda4;
    double mu1_at_MZ=mu1;
    double mu3_at_MZ=mu3;
    double mu4_at_MZ=mu4;
    double nu1_at_MZ=nu1;
    double omega1_at_MZ=omega1;
    double kappa1_at_MZ=kappa1;
    double nu2_at_MZ=nu2;
    double omega2_at_MZ=omega2;
    double kappa2_at_MZ=kappa2;
    double nu4_at_MZ=nu4;
    double omega4_at_MZ=omega4;
    double NLOuniscale=myTHDMW->getNLOuniscaleTHDMW();

    if(fabs(Q_THDMW-log10(MZ))<0.005)   //at MZ scale
    {
        Q_cutoff=log10(MZ);

        lambda1_at_Q = lambda1_at_MZ;
        lambda2_at_Q = lambda2_at_MZ;
        lambda3_at_Q = lambda3_at_MZ;
        lambda4_at_Q = lambda4_at_MZ;
        mu1_at_Q = mu1_at_MZ;
        mu3_at_Q = mu3_at_MZ;
        mu4_at_Q = mu4_at_MZ;
        nu1_at_Q = nu1_at_MZ;
        omega1_at_Q = omega1_at_MZ;
        kappa1_at_Q = kappa1_at_MZ;
        nu2_at_Q = nu2_at_MZ;
        omega2_at_Q = omega2_at_MZ;
        kappa2_at_Q = kappa2_at_MZ;
        nu4_at_Q = nu4_at_MZ;
        omega4_at_Q = omega4_at_MZ;
    }
    else   //at some other scale
    {
        double InitVals[15];
        InitVals[0]=lambda1_at_MZ;
        InitVals[1]=lambda2_at_MZ;
        InitVals[2]=lambda3_at_MZ;
        InitVals[3]=lambda4_at_MZ;
        InitVals[4]=mu1_at_MZ;
        InitVals[5]=mu3_at_MZ;
        InitVals[6]=mu4_at_MZ;
        InitVals[7]=nu1_at_MZ;
        InitVals[8]=omega1_at_MZ;
        InitVals[9]=kappa1_at_MZ;
        InitVals[10]=nu2_at_MZ;
        InitVals[11]=omega2_at_MZ;
        InitVals[12]=kappa2_at_MZ;
        InitVals[13]=nu4_at_MZ;
        InitVals[14]=omega4_at_MZ;

        Q_cutoff=myRunnerTHDMW->RGERunnerTHDMW(InitVals, 15, log10(MZ), Q_THDMW, flag, RpepsTHDMW, NLOuniscale);  //Running up to Q_cutoff<=Q_THDM

        lambda1_at_Q = InitVals[0];
        lambda2_at_Q = InitVals[1];
        lambda3_at_Q = InitVals[2];
        lambda4_at_Q = InitVals[3];
        mu1_at_Q=InitVals[4];
        mu3_at_Q=InitVals[5];
        mu4_at_Q = InitVals[6];
        nu1_at_Q = InitVals[7];
        omega1_at_Q = InitVals[8];
        kappa1_at_Q = InitVals[9];
        nu2_at_Q = InitVals[10];
        omega2_at_Q = InitVals[11];
        kappa2_at_Q = InitVals[12];
        nu4_at_Q = InitVals[13];
        omega4_at_Q = InitVals[14];
    }

}

void THDMWcache::computeUnitarity()
{
    std::string ModelType=myTHDMW->getModelTypeTHDMWflag();
    if( ModelType != "custodial1" )
    {
        throw std::runtime_error("THDMW unitarity constraints are only implemented for the \"custodial1\" model.");
    }

    double pi=M_PI;
    gslpp::matrix<gslpp::complex> Smatrix1(4,4,0.), Smatrix2(4,4,0.);
    gslpp::matrix<gslpp::complex> Sbmatrix1(4,4,0.), Sbmatrix2(4,4,0.);
    gslpp::matrix<gslpp::complex> Seigenvectors1(4,4,0.), Seigenvectors2(4,4,0.);
    gslpp::matrix<gslpp::complex> Seigenvectors1T(4,4,0.), Seigenvectors2T(4,4,0.);
    gslpp::vector<double> Seigenvalues1(4,0.), Seigenvalues2(4,0.);
    gslpp::vector<gslpp::complex> Sbeigenvalues1(4,0.), Sbeigenvalues2(4,0.);

    /*
    *******   LO part   *************
    */

    // Definition of the blocks of the S-matrix
    Smatrix1.assign(0,0, 3.0*lambda1/(16.0*pi));
    Smatrix1.assign(0,1, (2.0*lambda3+lambda4)/(16.0*pi));
    Smatrix1.assign(1,0, Smatrix1(0,1));
    Smatrix1.assign(0,3, (2.0*nu1+nu2)/(8.0*sqrt(2.0)*pi));
    Smatrix1.assign(3,0, Smatrix1(0,3));
    Smatrix1.assign(1,1, 3.0*lambda2/(16.0*pi));
    Smatrix1.assign(1,3, (2.0*omega1+omega2)/(8.0*sqrt(2.0)*pi));
    Smatrix1.assign(3,1, Smatrix1(1,3));
    Smatrix1.assign(2,2, (lambda3+5.0*lambda4)/(16.0*pi));
    Smatrix1.assign(2,3, (4.0*kappa1+2.0*kappa2)/(16.0*pi));
    Smatrix1.assign(3,2, Smatrix1(2,3));
    Smatrix1.assign(3,3, (26.0*mu1+17.0*mu3+13.0*mu4)/(32.0*pi));

    Smatrix2.assign(0,0, lambda1/(16.0*pi));
    Smatrix2.assign(0,1, lambda4/(16.0*pi));
    Smatrix2.assign(1,0, Smatrix2(0,1));
    Smatrix2.assign(0,3, nu2/(8.0*sqrt(2.0)*pi));
    Smatrix2.assign(3,0, Smatrix2(0,3));
    Smatrix2.assign(1,1, lambda2/(16.0*pi));
    Smatrix2.assign(1,3, omega2/(8.0*sqrt(2.0)*pi));
    Smatrix2.assign(3,1, Smatrix2(1,3));
    Smatrix2.assign(2,2, (lambda3+lambda4)/(16.0*pi));
    Smatrix2.assign(2,3, kappa2/(8.0*pi));
    Smatrix2.assign(3,2, Smatrix2(2,3));
    Smatrix2.assign(3,3, (14.0*mu1+3.0*mu3+27.0*mu4)/(96.0*pi));

    Smatrix1.eigensystem(Seigenvectors1, Seigenvalues1);
    Smatrix2.eigensystem(Seigenvectors2, Seigenvalues2);

    for (int i=0; i < 4; i++) {
        unitarityeigenvalues.assign(i, Seigenvalues1(i));
        unitarityeigenvalues.assign(4+i, Seigenvalues2(i));
    }
    unitarityeigenvalues.assign(8, (lambda3-lambda4)/(16.0*pi));
    unitarityeigenvalues.assign(9, sqrt(15.0)*nu4/(16.0*pi));
    unitarityeigenvalues.assign(10, sqrt(15.0)*omega4/(16.0*pi));

    
    
    /*
    *******   NLO part   *************
    */

    double blambda1=(12.0*lambda1*lambda1 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4 
                     + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
    double blambda2=(12.0*lambda2*lambda2 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4
                     + 8.0*omega1*omega1 + 8.0*omega1*omega2 + 8.0*omega2*omega2)/(16.0*pi*pi);
    double blambda3=(4.0*lambda3*lambda3 + 4.0*lambda4*lambda4 + (lambda1+lambda2)*(6.0*lambda3+2.0*lambda4) 
                     + 8.0*kappa2*kappa2 + 8.0*nu1*omega1 + 4.0*nu2*omega1 + 4.0*nu1*omega2)/(16.0*pi*pi);
    double blambda4=(lambda1*lambda4 + lambda2*lambda4 + 4.0*lambda3*lambda4 + 6.0*lambda4*lambda4
                     + 4.0*kappa1*kappa1 + 4.0*kappa1*kappa2 + 2.0*kappa2*kappa2 + 2.0*nu2*omega2)/(8.0*pi*pi);
    double bmu1=(11.0*mu1*mu1 + 3.0*mu1*mu4 + mu1*(2.0*mu1+6.0*mu3+3.0*mu4)
               + 3.0*nu4*nu4 + 3.0*omega4*omega4)/(16.0*pi*pi);
    double bmu3=(18.0*kappa1*kappa1 + 18.0*kappa1*kappa2 + 134.0*mu1*mu1 + 6.0*mu1*(39.0*mu3 + 22.0*mu4)
               + 3.0*(30.0*mu3*mu3 + 39.0*mu3*mu4 + 9.0*mu4*mu4 
                      + 3.0*nu1*nu1 + 3.0*nu1*nu2 - 5.0*nu4*nu4
                      + 3.0*omega1*omega1 + 3.0*omega1*omega2 - 5.0*omega4*omega4))/(72.0*pi*pi);
    double bmu4=(18.0*kappa2*kappa2 + 4.0*mu1*mu1 + 156.0*mu1*mu4 + 54.0*mu3*mu4 + 144.0*mu4*mu4
               + 9.0*nu2*nu2 + 6.0*nu4*nu4 + 9.0*omega2*omega2 + 6.0*omega4*omega4)/(144.0*pi*pi);
    double bnu1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 + 18.0*lambda1*nu1
               + 78.0*mu1*nu1 + 51.0*mu3*nu1 + 39.0*mu4*nu1 + 6.0*nu1*nu1
               + 6.0*lambda1*nu2 + 32.0*mu1*nu2 + 24.0*mu3*nu2 + 6.0*mu4*nu2
               + 6.0*nu2*nu2 + 10.0*nu4*nu4
               + 12.0*lambda3*omega1 + 6.0*lambda4*omega1 + 6.0*lambda3*omega2)/(48.0*pi*pi);
    double bomega1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 
               + 12.0*lambda3*nu1 + 6.0*lambda4*nu1 + 6.0*lambda3*nu2
               + 18.0*lambda2*omega1 + 78.0*mu1*omega1 + 51.0*mu3*omega1 + 39.0*mu4*omega1 + 6.0*omega1*omega1
               + 6.0*lambda2*omega2 + 32.0*mu1*omega2 + 24.0*mu3*omega2 + 6.0*mu4*omega2 + 6.0*omega2*omega2
               + 10.0*omega4*omega4)/(48.0*pi*pi);
    double bkappa1=(6.0*kappa1*(2.0*lambda3 + 10.0*lambda4 + 18.0*mu1 + 17.0*mu3 + 13.0*mu4 + 2.0*nu1 + 2.0*omega1)
               + kappa2*(24.0*lambda4 + 64.0*mu1 + 48.0*mu3 + 24.0*mu4 + 9.0*nu2 + 9.0*omega2)
               + 20.0*nu4*omega4)/(96.0*pi*pi);
    double bnu2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda1*nu2 + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*nu2 
                + 4.0*nu1*nu2 + 6.0*nu2*nu2 + (25.0*nu4*nu4)/3.0 + 2.0*lambda4*omega2)/(16.0*pi*pi);
    double bomega2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda4*nu2 + 2.0*lambda2*omega2 
                + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*omega2 + 4.0*omega1*omega2 + 6.0*omega2*omega2 
                + (25.0*omega4*omega4)/3.0)/(16.0*pi*pi);
    double bkappa2=(kappa2*(6.0*lambda3 + 6.0*lambda4 + 14.0*mu1 + 3.0*mu3 + 27.0*mu4
                     + 6.0*nu1 + 12.0*nu2 + 6.0*omega1 + 12.0*omega2)
                + 6.0*kappa1*(nu2 + omega2) + 42.0*nu4*omega4)/(48.0*pi*pi);
    double bnu4=(11.0*mu1*nu4 + 3.0*mu3*nu4 + 9.0*mu4*nu4 + 3.0*nu1*nu4 + 9.0*nu2*nu4 
                + 3.0*kappa1*omega4 + 9.0*kappa2*omega4)/(16.0*pi*pi);
    double bomega4=(3.0*kappa1*nu4 + 9.0*kappa2*nu4 
                + (11.0*mu1 + 3.0*(mu3 + 3.0*mu4 + omega1 + 3.0*omega2))*omega4)/(16.0*pi*pi);

    Sbmatrix1.assign(0,0, 3.0*blambda1/(16.0*pi));
    Sbmatrix1.assign(0,1, (2.0*blambda3+blambda4)/(16.0*pi));
    Sbmatrix1.assign(1,0, Smatrix1(0,1));
    Sbmatrix1.assign(0,3, (2.0*bnu1+bnu2)/(8.0*sqrt(2.0)*pi));
    Sbmatrix1.assign(3,0, Smatrix1(0,3));
    Sbmatrix1.assign(1,1, 3.0*blambda2/(16.0*pi));
    Sbmatrix1.assign(1,3, (2.0*bomega1+bomega2)/(8.0*sqrt(2.0)*pi));
    Sbmatrix1.assign(3,1, Smatrix1(1,3));
    Sbmatrix1.assign(2,2, (blambda3+5.0*blambda4)/(16.0*pi));
    Sbmatrix1.assign(2,3, (4.0*bkappa1+2.0*bkappa2)/(16.0*pi));
    Sbmatrix1.assign(3,2, Smatrix1(2,3));
    Sbmatrix1.assign(3,3, (26.0*bmu1+17.0*bmu3+13.0*bmu4)/(32.0*pi));

    Sbmatrix2.assign(0,0, blambda1/(16.0*pi));
    Sbmatrix2.assign(0,1, blambda4/(16.0*pi));
    Sbmatrix2.assign(1,0, Smatrix2(0,1));
    Sbmatrix2.assign(0,3, bnu2/(8.0*sqrt(2.0)*pi));
    Sbmatrix2.assign(3,0, Smatrix2(0,3));
    Sbmatrix2.assign(1,1, blambda2/(16.0*pi));
    Sbmatrix2.assign(1,3, bomega2/(8.0*sqrt(2.0)*pi));
    Sbmatrix2.assign(3,1, Smatrix2(1,3));
    Sbmatrix2.assign(2,2, (blambda3+blambda4)/(16.0*pi));
    Sbmatrix2.assign(2,3, bkappa2/(8.0*pi));
    Sbmatrix2.assign(3,2, Smatrix2(2,3));
    Sbmatrix2.assign(3,3, (14.0*bmu1+3.0*bmu3+27.0*bmu4)/(96.0*pi));

    Seigenvectors1T=Seigenvectors1.hconjugate();
    Seigenvectors2T=Seigenvectors2.hconjugate();

    for (int i=0; i < 4; i++) {
        for (int k=0; k < 4; k++) {
            for (int l=0; l < 4; l++) {
                Sbeigenvalues1.assign(i, Sbeigenvalues1(i) + Seigenvectors1(i,k) * Sbmatrix1(k,l) * Seigenvectors1T(l,i) );
                Sbeigenvalues2.assign(i, Sbeigenvalues2(i) + Seigenvectors2(i,k) * Sbmatrix2(k,l) * Seigenvectors2T(l,i) );
            }                
        }
        betaeigenvalues.assign(i, Sbeigenvalues1(i));
        betaeigenvalues.assign(i+4, Sbeigenvalues2(i));
    }
    
    betaeigenvalues.assign(8, -1.5 * (blambda3-blambda4)/(16.0*pi));
    betaeigenvalues.assign(9, -1.5 * sqrt(15.0)*bnu4/(16.0*pi));
    betaeigenvalues.assign(10, -1.5 * sqrt(15.0)*bomega4/(16.0*pi));

//    std::cout<<"Seigenvectors1 = "<<Seigenvectors1<<std::endl;
//    std::cout<<"Seigenvectors2 = "<<Seigenvectors2<<std::endl;
//    std::cout<<"unitarityeigenvalues = "<<unitarityeigenvalues<<std::endl;
//    std::cout<<"betaeigenvalues = "<<betaeigenvalues<<std::endl;

    for (int i=0; i < 11; i++) {
        NLOunitarityeigenvalues.assign(i, (gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i)+betaeigenvalues(i) );
    }
}

void THDMWcache::updateCache()
{
    Q_THDMW=myTHDMW->getQ_THDMW();
    MZ=myTHDMW->getMz();
    lambda1=myTHDMW->getTHDMW_lambda1();
    lambda2=myTHDMW->getTHDMW_lambda2();
    lambda3=myTHDMW->getTHDMW_lambda3();
    lambda4=myTHDMW->getTHDMW_lambda4();
    mu1=myTHDMW->getTHDMW_mu1();
    mu3=myTHDMW->getTHDMW_mu3();
    mu4=myTHDMW->getTHDMW_mu4();
    nu1=myTHDMW->getTHDMW_nu1();
    omega1=myTHDMW->getTHDMW_omega1();
    kappa1=myTHDMW->getTHDMW_kappa1();
    nu2=myTHDMW->getTHDMW_nu2();
    omega2=myTHDMW->getTHDMW_omega2();
    kappa2=myTHDMW->getTHDMW_kappa2();
    nu4=myTHDMW->getTHDMW_nu4();
    omega4=myTHDMW->getTHDMW_omega4();
    RpepsTHDMW=myTHDMW->getRpepsTHDMW();

    runTHDMWparameters();
    computeUnitarity();
}
