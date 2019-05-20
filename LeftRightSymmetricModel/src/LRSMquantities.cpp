/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LRSMquantities.h"
#include "gslpp.h"
#include "LeftRightSymmetricModel.h"

LRSMquantities::LRSMquantities(const LeftRightSymmetricModel & LRSM_in)
: myLRSM(LRSM_in), Msqneutral(5, 5, 0.)
{
}

//LRSMquantities::LRSMquantities(const StandardModel& SM_i)
//: myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i)), Msqneutral(5,5,0.)
//{}

LRSMquantities::~LRSMquantities()
{
}

bool LRSMquantities::CalcNeutralMasses(gslpp::matrix<gslpp::complex>& U_i, double mH0sq[5])
{
    double Ale = myLRSM.getAle();
    double cW2 = myLRSM.c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM.v();
    double mH2psq = myLRSM.getmH2p_2();
    double xi = myLRSM.getxi_LRSM();
    double xi2 = xi*xi;
    double kappa = sqrt(0.5) * vev * (1.0 - 0.5 * xi * xi);
    double kappasq = kappa*kappa;
    double mWR = myLRSM.getmWR();
    double mWR2 = mWR*mWR;
    double lambda1 = myLRSM.getlambda1_LRSM();
    double lambda2 = myLRSM.getlambda2_LRSM();
    double lambda3 = myLRSM.getlambda3_LRSM();
    double lambda4 = myLRSM.getlambda4_LRSM();
    double rho1 = myLRSM.getrho1_LRSM();
    double alpha1 = myLRSM.getalpha1_LRSM();
    double alpha2 = myLRSM.getalpha2_LRSM();

    Msqneutral.assign(0, 0, (-((2.0 * mH2psq * mWR2 * xi2) / (g2_2 * kappasq + 2.0 * mWR2)) + 4.0 * kappasq * (xi2 - 1.0)*(lambda1 + xi * (2.0 * lambda4 + 2.0 * lambda2 * xi + lambda3 * xi))) / (xi2 - 1.0));
    Msqneutral.assign(0, 2, ((2.0 * mH2psq * mWR2 * xi) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (xi2 - 1.0)*((lambda1 + 2.0 * lambda2 + lambda3) * xi + lambda4 * (1.0 + xi2))) / (xi2 - 1.0));
    Msqneutral.assign(0, 4, (2.0 * kappa * mWR * (alpha1 + 2.0 * alpha2 * xi)) / g2);
    //
    Msqneutral.assign(1, 1, -((xi2 * ((2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (2.0 * lambda2 - lambda3)*(xi2 - 1.0))) / (xi2 - 1.0)));
    Msqneutral.assign(1, 3, -((xi * ((2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (2.0 * lambda2 - lambda3)*(xi2 - 1.0))) / (xi2 - 1.0)));
    //
    Msqneutral.assign(2, 0, ((2.0 * mH2psq * mWR2 * xi) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (xi2 - 1.0)*((lambda1 + 2.0 * lambda2 + lambda3) * xi + lambda4 * (1.0 + xi2))) / (xi2 - 1.0));
    Msqneutral.assign(2, 2, (-((2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2)) + 4.0 * kappasq * (xi2 - 1.0)*(2.0 * lambda2 + lambda3 + 2.0 * lambda4 * xi + lambda1 * xi2)) / (xi2 - 1.0));
    Msqneutral.assign(2, 4, (2.0 * kappa * mWR * (2.0 * alpha2 + ((2.0 * mH2psq) / (kappasq + (2.0 * mWR2) / g2_2) + alpha1) * xi)) / g2);
    //
    Msqneutral.assign(3, 1, -((xi * ((2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (2.0 * lambda2 - lambda3)*(xi2 - 1.0))) / (xi2 - 1.0)));
    Msqneutral.assign(3, 3, -(((2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2) + 4.0 * kappasq * (2.0 * lambda2 - lambda3)*(xi2 - 1.0)) / (xi2 - 1.0)));
    //
    Msqneutral.assign(4, 0, (2.0 * kappa * mWR * (alpha1 + 2.0 * alpha2 * xi)) / g2);
    Msqneutral.assign(4, 2, (2.0 * kappa * mWR * (2.0 * alpha2 + ((2.0 * mH2psq) / (kappasq + (2.0 * mWR2) / g2_2) + alpha1) * xi)) / g2);
    Msqneutral.assign(4, 4, (4.0 * mWR2 * rho1) / g2_2);

    gslpp::vector<gslpp::complex> mH0sq_i(5, 0.);
    Msqneutral.eigensystem(U_i, mH0sq_i);
    for (int i = 0; i < 5; i++) {
        mH0sq[i] = mH0sq_i(i).real();
    }

    int newIndex[5];
    for (int i = 0; i < 5; i++)
        newIndex[i] = i;

    /* sort sfermion masses in increasing order */
    for (int i = 0; i < 4; i++) {
        for (int k = i + 1; k < 5; k++)
            if (mH0sq[i] > mH0sq[k]) {
                std::swap(mH0sq[i], mH0sq[k]);
                std::swap(newIndex[i], newIndex[k]);
            }
    }

    return true;
}

bool LRSMquantities::CalcNeutralMasses_app(double mH0sq_app[4])
{
    double Ale = myLRSM.getAle();
    double cW2 = myLRSM.c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM.v();
    double mH1psq = myLRSM.getmH1p_2();
    double mH2psq = myLRSM.getmH2p_2();
    double xi = myLRSM.getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM.getmWR();
    double mWR2 = mWR*mWR;
    double lambda1 = myLRSM.getlambda1_LRSM();
    double rho1 = myLRSM.getrho1_LRSM();
    double alpha1 = myLRSM.getalpha1_LRSM();

    mH0sq_app[0] = kappasq * (4.0 * lambda1 - (alpha1 * alpha1) / rho1)+(2.0 * mH2psq * mWR2 * xi * xi) / (g2_2 * kappasq + 2.0 * mWR2);
    mH0sq_app[1] = (4.0 * mWR2 * rho1) / g2_2;
    mH0sq_app[2] = (g2_2 * kappasq * (mH1psq - mH2psq) + 2.0 * mH1psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2);
    mH0sq_app[3] = (2.0 * mH2psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2);

    int newIndex[4];
    for (int i = 0; i < 4; i++)
        newIndex[i] = i;

    /* sort sfermion masses in increasing order */
    for (int i = 0; i < 3; i++) {
        for (int k = i + 1; k < 4; k++)
            if (mH0sq_app[i] > mH0sq_app[k]) {
                std::swap(mH0sq_app[i], mH0sq_app[k]);
                std::swap(newIndex[i], newIndex[k]);
            }
    }

    return true;
}

mu1_2_LRSM::mu1_2_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double mu1_2_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;
    double lambda1 = myLRSM->getlambda1_LRSM();
    double alpha1 = myLRSM->getalpha1_LRSM();
    //    bool CPviolation=myLRSM->getCPVflag();

    //    double mu1_2=0.0;
    //    if(CPviolation)
    //    {
    //    }
    //    else
    //    {
    //    }
    return (mWR2 * alpha1) / g2_2 + 2.0 * kappasq*lambda1;
}

mu2_2_LRSM::mu2_2_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double mu2_2_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH2psq = myLRSM->getmH2p_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;
    double lambda4 = myLRSM->getlambda4_LRSM();
    double alpha2 = myLRSM->getalpha2_LRSM();

    return (mWR2 * alpha2) / g2_2 + kappasq * lambda4 + (mH2psq * mWR2 * xi) / (g2_2 * kappasq + 2.0 * mWR2);
}

mu3_2_LRSM::mu3_2_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double mu3_2_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;
    double rho1 = myLRSM->getrho1_LRSM();
    double alpha1 = myLRSM->getalpha1_LRSM();

    return (2.0 * mWR2 * rho1) / g2_2 + kappasq*alpha1;
}

rho2_LRSM::rho2_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double rho2_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH2psq = myLRSM->getmH2p_2();
    double mdeltappRsq = myLRSM->getmdeltappR_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;

    return g2_2 * (mdeltappRsq - (2.0 * kappasq * mH2psq) / (kappasq + (2.0 * mWR2) / g2_2)) / (4.0 * mWR2);
}

rho3_LRSM::rho3_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double rho3_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH1psq = myLRSM->getmH1p_2();
    double mH2psq = myLRSM->getmH2p_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;
    double rho1 = myLRSM->getrho1_LRSM();

    return 2.0 * rho1 + g2_2 / mWR2 * (mH1psq - (g2_2 * kappasq * mH2psq) / (g2_2 * kappasq + 2.0 * mWR2));
}

alpha3_LRSM::alpha3_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double alpha3_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH2psq = myLRSM->getmH2p_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;

    return (2.0 * mH2psq) / ((2.0 * mWR2) / g2_2 + kappasq);
}

MH05_LRSM::MH05_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double MH05_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH1psq = myLRSM->getmH1p_2();
    double mH2psq = myLRSM->getmH2p_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;

    return sqrt((g2_2 * kappasq * (mH1psq - mH2psq) + 2.0 * mH1psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2));
}

MH06_LRSM::MH06_LRSM(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double MH06_LRSM::computeThValue()
{
    double Ale = myLRSM->getAle();
    double cW2 = myLRSM->c02();
    double g2 = sqrt(4.0 * M_PI * Ale / (1 - cW2));
    double g2_2 = g2*g2;
    double vev = myLRSM->v();
    double mH1psq = myLRSM->getmH1p_2();
    double mH2psq = myLRSM->getmH2p_2();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double mWR = myLRSM->getmWR();
    double mWR2 = mWR*mWR;

    return sqrt((g2_2 * kappasq * (mH1psq - mH2psq) + 2.0 * mH1psq * mWR2) / (g2_2 * kappasq + 2.0 * mWR2));
}

MH01_app1::MH01_app1(const StandardModel& SM_i)
: ThObservable(SM_i), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double MH01_app1::computeThValue()
{
    double vev = myLRSM->v();
    double xi = myLRSM->getxi_LRSM();
    double kappasq = 0.5 * vev * vev * (1.0 - xi * xi);
    double lambda1 = myLRSM->getlambda1_LRSM();
    double rho1 = myLRSM->getrho1_LRSM();
    double alpha1 = myLRSM->getalpha1_LRSM();

    return sqrt(kappasq * (4.0 * lambda1 - (alpha1 * alpha1) / rho1));
}

MH0_LRSM::MH0_LRSM(const StandardModel& SM_i, const int ind)
: ThObservable(SM_i), index(ind), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double MH0_LRSM::computeThValue()
{
    switch (index)
    {
        case 0:
            return sqrt(myLRSM->getmH0sq1());
        case 1:
            return sqrt(myLRSM->getmH0sq2());
        case 2:
            return sqrt(myLRSM->getmH0sq3());
        case 3:
            return sqrt(myLRSM->getmH0sq4());
        case 4:
            return sqrt(myLRSM->getmH0sq5());
        default:
            throw std::runtime_error("MH0_LRSM::computeThValue(): undefined index");
    }
}

MH0_app::MH0_app(const StandardModel& SM_i, const int ind)
: ThObservable(SM_i), index(ind), myLRSM(static_cast<const LeftRightSymmetricModel*> (&SM_i))
{
}

double MH0_app::computeThValue()
{
    switch (index)
    {
        case 0:
            return sqrt(myLRSM->getmH0sq1_app());
        case 1:
            return sqrt(myLRSM->getmH0sq2_app());
        case 2:
            return sqrt(myLRSM->getmH0sq3_app());
        case 3:
            return sqrt(myLRSM->getmH0sq4_app());
        default:
            throw std::runtime_error("MH0_app::computeThValue(): undefined index");
    }
}
