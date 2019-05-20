/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMunitarity.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

unitarity_GTHDM::unitarity_GTHDM(const StandardModel& SM_i)
: myGTHDM(static_cast<const GeneralTHDM&> (SM_i)), Smat21(3,3,0.), Smat01(4,4,0.), Smat00(4,4,0.),
        Seigvec21(3,3,0.), Seigvec01(4,4,0.), Seigvec00(4,4,0.), Seigval21(3,0.), Seigval01(4,0.), Seigval00(4,0.)
{}

unitarity_GTHDM::~unitarity_GTHDM() 
{}

bool unitarity_GTHDM::CalcSeigen21(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i)
{
    double lambda1 = myGTHDM.getMyGTHDMCache()->lambda1;
    double lambda2 = myGTHDM.getMyGTHDMCache()->lambda2;
    double lambda3 = myGTHDM.getMyGTHDMCache()->lambda3;
    double lambda4 = myGTHDM.getMyGTHDMCache()->lambda4;
    double Relambda5 = myGTHDM.getRelambda5();
    double Imlambda5 = myGTHDM.getImlambda5();
    double Relambda6 = myGTHDM.getRelambda6();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda6 = myGTHDM.getMyGTHDMCache()->Imlambda6;
    double Imlambda7 = myGTHDM.getMyGTHDMCache()->Imlambda7;
    
    gslpp::complex i = gslpp::complex::i();
    
    Smat21.assign(0,0, lambda1);
    Smat21.assign(0,1, Relambda5 + i*Imlambda5);
    Smat21.assign(0,2, sqrt(2)*(Relambda6 + i*Imlambda6));
    Smat21.assign(1,0, Relambda5 - i*Imlambda5);
    Smat21.assign(1,1, lambda2);
    Smat21.assign(1,2, sqrt(2)*(Relambda7 - i*Imlambda7));
    Smat21.assign(2,0, sqrt(2)*(Relambda6 - i*Imlambda6));
    Smat21.assign(2,1, sqrt(2)*(Relambda7 + i*Imlambda7));
    Smat21.assign(2,2, lambda3 + lambda4);
    
    Smat21.eigensystem(Seigvec_i, Seigval_i);
    
    return true;
}

bool unitarity_GTHDM::CalcSeigen01(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i)
{
    double lambda1 = myGTHDM.getMyGTHDMCache()->lambda1;
    double lambda2 = myGTHDM.getMyGTHDMCache()->lambda2;
    double lambda3 = myGTHDM.getMyGTHDMCache()->lambda3;
    double lambda4 = myGTHDM.getMyGTHDMCache()->lambda4;
    double Relambda5 = myGTHDM.getRelambda5();
    double Imlambda5 = myGTHDM.getImlambda5();
    double Relambda6 = myGTHDM.getRelambda6();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda6 = myGTHDM.getMyGTHDMCache()->Imlambda6;
    double Imlambda7 = myGTHDM.getMyGTHDMCache()->Imlambda7;
    
    gslpp::complex i = gslpp::complex::i();
    
    Smat01.assign(0,0, lambda1);
    Smat01.assign(0,1, lambda4);
    Smat01.assign(0,2, Relambda6 + i*Imlambda6);
    Smat01.assign(0,3, Relambda6 - i*Imlambda6);
    Smat01.assign(1,0, lambda4);
    Smat01.assign(1,1, lambda2);
    Smat01.assign(1,2, Relambda7 + i*Imlambda7);
    Smat01.assign(1,3, Relambda7 - i*Imlambda7);
    Smat01.assign(2,0, Relambda6 - i*Imlambda6);
    Smat01.assign(2,1, Relambda7 - i*Imlambda7);
    Smat01.assign(2,2, lambda3);
    Smat01.assign(2,3, Relambda5 - i*Imlambda5);
    Smat01.assign(3,0, Relambda6 + i*Imlambda6);
    Smat01.assign(3,1, Relambda7 + i*Imlambda7);
    Smat01.assign(3,2, Relambda5 + i*Imlambda5);
    Smat01.assign(3,3, lambda3);
    
    Smat01.eigensystem(Seigvec_i, Seigval_i);

    return true;
}

bool unitarity_GTHDM::CalcSeigen00(gslpp::matrix<gslpp::complex>& Seigvec_i, gslpp::vector<double>& Seigval_i)
{
    double lambda1 = myGTHDM.getMyGTHDMCache()->lambda1;
    double lambda2 = myGTHDM.getMyGTHDMCache()->lambda2;
    double lambda3 = myGTHDM.getMyGTHDMCache()->lambda3;
    double lambda4 = myGTHDM.getMyGTHDMCache()->lambda4;
    double Relambda5 = myGTHDM.getRelambda5();
    double Imlambda5 = myGTHDM.getImlambda5();
    double Relambda6 = myGTHDM.getRelambda6();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda6 = myGTHDM.getMyGTHDMCache()->Imlambda6;
    double Imlambda7 = myGTHDM.getMyGTHDMCache()->Imlambda7;
    
    gslpp::complex i = gslpp::complex::i();
    
    Smat00.assign(0,0, 3.*lambda1);
    Smat00.assign(0,1, 2.*lambda3 + lambda4);
    Smat00.assign(0,2, 3.*(Relambda6 + i*Imlambda6));
    Smat00.assign(0,3, 3.*(Relambda6 - i*Imlambda6));
    Smat00.assign(1,0, 2.*lambda3 + lambda4);
    Smat00.assign(1,1, 3.*lambda2);
    Smat00.assign(1,2, 3.*(Relambda7 + i*Imlambda7));
    Smat00.assign(1,3, 3.*(Relambda7 - i*Imlambda7));
    Smat00.assign(2,0, 3.*(Relambda6 - i*Imlambda6));
    Smat00.assign(2,1, 3.*(Relambda7 - i*Imlambda7));
    Smat00.assign(2,2, lambda3 + 2.*lambda4);
    Smat00.assign(2,3, 3.*(Relambda5 - i*Imlambda5));
    Smat00.assign(3,0, 3.*(Relambda6 + i*Imlambda6));
    Smat00.assign(3,1, 3.*(Relambda7 + i*Imlambda7));
    Smat00.assign(3,2, 3.*(Relambda5 + i*Imlambda5));
    Smat00.assign(3,3, lambda3 + 2.*lambda4);
    
    Smat00.eigensystem(Seigvec_i, Seigval_i);

    return true;
}

gslpp::vector<double> unitarity_GTHDM::getSeigen21()
{
    CalcSeigen21(Seigvec21,Seigval21);
    
    return Seigval21;
}

gslpp::vector<double> unitarity_GTHDM::getSeigen01()
{
    CalcSeigen01(Seigvec01,Seigval01);
    
    return Seigval01;
}

gslpp::vector<double> unitarity_GTHDM::getSeigen00()
{
    CalcSeigen00(Seigvec00,Seigval00);
    
    return Seigval00;
}


unitarity1_GTHDM::unitarity1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity1_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen21())(0);
}


unitarity2_GTHDM::unitarity2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity2_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen21())(1);
}


unitarity3_GTHDM::unitarity3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity3_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen21())(2);
}


unitarity5_GTHDM::unitarity5_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity5_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen01())(0);
}


unitarity6_GTHDM::unitarity6_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity6_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen01())(1);
}


unitarity7_GTHDM::unitarity7_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity7_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen01())(2);
}


unitarity8_GTHDM::unitarity8_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity8_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen01())(3);
}


unitarity9_GTHDM::unitarity9_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity9_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen00())(0);
}


unitarity10_GTHDM::unitarity10_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity10_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen00())(1);
}


unitarity11_GTHDM::unitarity11_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity11_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen00())(2);
}


unitarity12_GTHDM::unitarity12_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myunitarity_GTHDM(SM_i)
{}

double unitarity12_GTHDM::computeThValue()
{
    return (myunitarity_GTHDM.getSeigen00())(3);
}


unitarity4_GTHDM::unitarity4_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double unitarity4_GTHDM::computeThValue()
{
    double lambda3 = myGTHDM.getMyGTHDMCache()->lambda3;
    double lambda4 = myGTHDM.getMyGTHDMCache()->lambda4;
    return (lambda3-lambda4);
}
