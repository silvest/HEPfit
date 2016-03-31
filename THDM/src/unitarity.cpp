/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "unitarity.h"
#include "StandardModel.h"

unitarity::unitarity(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double unitarity::computeThValue()
{
    return 0.0;
}

unitarity1::unitarity1(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity1::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)/sinb/sinb
            +(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
            +sqrt(4.0*pow(mA2 -m12_2/(cosb*sinb),2)
                  +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                       +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                       +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity2::unitarity2(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity2::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)/sinb/sinb
            +(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
            -sqrt(4.0*pow(mA2 -m12_2/(cosb*sinb),2)
                  +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                       +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                       +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity3::unitarity3(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity3::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)/sinb/sinb 
            +(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
            +sqrt(4.0*pow(mA2 -2*mHp2 +m12_2/(cosb*sinb),2)
                  +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                       +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                       +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity4::unitarity4(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity4::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)/sinb/sinb
            +(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
            -sqrt(4.0*pow(mA2 -2*mHp2 +m12_2/(cosb*sinb),2)
                  +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                       +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                       +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity5::unitarity5(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity5::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (3.0*((mHh2/cosb/cosb +mHl*mHl/sinb/sinb)*cosa*cosa
                 +(mHl*mHl/cosb/cosb + mHh2/sinb/sinb)*sina*sina
                 -m12_2*(1.0/tanb/sinb/sinb + tanb/cosb/cosb))
            +sqrt(4.0*pow(mA2 + 2.0*mHp2
                          -(m12_2 +2.0*(mHl*mHl - mHh2)*cosa*sina)/(cosb*sinb),2)
                  +(9.0*pow((mHl*mHl - mHh2)*(cosa*cosa - sina*sina) 
                            +(mHl*mHl + mHh2 -(2.0*m12_2)/(cosb*sinb))
                             *(cosb*cosb - sinb*sinb),2))
                   /(4.0*pow(cosb,4)*pow(sinb,4))))/(2.0*vev*vev);
}

unitarity6::unitarity6(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity6::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (3.0*((mHh2/cosb/cosb +mHl*mHl/sinb/sinb)*cosa*cosa
                 +(mHl*mHl/cosb/cosb + mHh2/sinb/sinb)*sina*sina
                 -m12_2*(1.0/tanb/sinb/sinb + tanb/cosb/cosb))
            -sqrt(4.0*pow(mA2 + 2.0*mHp2
                          -(m12_2 +2.0*(mHl*mHl - mHh2)*cosa*sina)/(cosb*sinb),2)
                  +(9.0*pow((mHl*mHl - mHh2)*(cosa*cosa - sina*sina) 
                            +(mHl*mHl + mHh2 -(2.0*m12_2)/(cosb*sinb))
                             *(cosb*cosb - sinb*sinb),2))
                   /(4.0*pow(cosb,4)*pow(sinb,4))))/(2.0*vev*vev);
}

unitarity7::unitarity7(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity7::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();

    return (mA2 +((mHh2-mHl*mHl)*cosa*sina)/(cosb*sinb))/(vev*vev);
}

unitarity8::unitarity8(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity8::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (-mA2 +4.0*mHp2
            +((mHh2-mHl*mHl)*cosa*sina -2.0*m12_2)/(cosb*sinb))/(vev*vev);
}

unitarity9::unitarity9(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity9::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();

    return (-mA2 +2.0*mHp2
            +((mHh2-mHl*mHl)*cosa*sina)/(cosb*sinb))/(vev*vev);
}

unitarity10::unitarity10(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity10::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (mA2 +2.0*mHp2
            +((mHh2-mHl*mHl)*cosa*sina -2.0*m12_2)/(cosb*sinb))/(vev*vev);
}

unitarity11::unitarity11(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity11::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (-mA2 -2.0*mHp2
            +((mHh2-mHl*mHl)*cosa*sina +4.0*m12_2)/(cosb*sinb))/(vev*vev);
}

unitarity12::unitarity12(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity12::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mA2=myTHDM->getmA2();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (5.0*mA2 -2.0*mHp2
            +((mHh2-mHl*mHl)*cosa*sina -2.0*m12_2)/(cosb*sinb))/(vev*vev);
}

unitarityNLO1::unitarityNLO1(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO1::computeThValue()
{
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;

    double uniNLO1a = -3.0*la1Q/(16.0*M_PI);
    double uniNLO1b = 9.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la1Q*la1Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1d = -3.0*la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    return (uniNLO1a+uniNLO1b+uniNLO1c+uniNLO1d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO2::unitarityNLO2(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO2::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;

    double uniNLO2a = -3.0*la2Q/(16.0*M_PI);
    double uniNLO2b = 9.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la2Q*la2Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2d = -3.0*la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    return (uniNLO2a+uniNLO2b+uniNLO2c+uniNLO2d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO3::unitarityNLO3(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO3::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO3a = -(2.0*la3Q+la4Q)/(16.0*M_PI);
    double uniNLO3b = 3.0*(2.0*betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*(2.0*la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3d = -(2.0*la3Q+la4Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO3a+uniNLO3b+uniNLO3c+uniNLO3d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO4::unitarityNLO4(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO4::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO4a = -(la3Q+2.0*la4Q)/(16.0*M_PI);
    double uniNLO4b = 3.0*(betalambda3+2.0*betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+4.0*la3Q*la4Q+4.0*la4Q*la4Q+9.0*la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4d = -(la3Q+la4Q+la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO4a+uniNLO4b+uniNLO4c+uniNLO4d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO5::unitarityNLO5(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO5::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO5a = -3.0*la5Q/(16.0*M_PI);
    double uniNLO5b = 9.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la3Q+2.0*la4Q)*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5d = -(la4Q+2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO5a+uniNLO5b+uniNLO5c+uniNLO5d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO6::unitarityNLO6(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO6::computeThValue()
{
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;

    double uniNLO6a = -la1Q/(16.0*M_PI);
    double uniNLO6b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    return (uniNLO6a+uniNLO6b+uniNLO6c+uniNLO6d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO7::unitarityNLO7(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO7::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;

    double uniNLO7a = -la2Q/(16.0*M_PI);
    double uniNLO7b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    return (uniNLO7a+uniNLO7b+uniNLO7c+uniNLO7d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO8::unitarityNLO8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO8::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO8a = -la4Q/(16.0*M_PI);
    double uniNLO8b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO8a+uniNLO8b+uniNLO8c+uniNLO8d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO9::unitarityNLO9(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO9::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO9a = -la4Q/(16.0*M_PI);
    double uniNLO9b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9d = -la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO9a+uniNLO9b+uniNLO9c+uniNLO9d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO10::unitarityNLO10(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO10::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO10a = -la3Q/(16.0*M_PI);
    double uniNLO10b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10d = -(la3Q+la4Q-la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO10a+uniNLO10b+uniNLO10c+uniNLO10d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO11::unitarityNLO11(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO11::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO11a = -la5Q/(16.0*M_PI);
    double uniNLO11b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11d = -(la4Q-2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO11a+uniNLO11b+uniNLO11c+uniNLO11d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO12::unitarityNLO12(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO12::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO12a = -la3Q/(16.0*M_PI);
    double uniNLO12b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12d = -la3Q/(512.0*M_PI*M_PI*M_PI) * WFRc3;

    return (uniNLO12a+uniNLO12b+uniNLO12c+uniNLO12d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO13::unitarityNLO13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO13::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO13a = -la5Q/(16.0*M_PI);
    double uniNLO13b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO13a+uniNLO13b+uniNLO13c+uniNLO13d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO14::unitarityNLO14(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO14::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO14a = -(la3Q-la4Q)/(16.0*M_PI);
    double uniNLO14b = 3.0*(betalambda3-betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO14c = (gslpp::complex::i()*M_PI-1.0)*(la3Q-la4Q)*(la3Q-la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO14d = -(la3Q-la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO14a+uniNLO14b+uniNLO14c+uniNLO14d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO15::unitarityNLO15(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO15::computeThValue()
{
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;

    double uniNLO15a = -la1Q/(16.0*M_PI);
    double uniNLO15b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1-2.0*WFRc2+WFRc3+2.0*WFRc4);

    return (uniNLO15a+uniNLO15b+uniNLO15c+uniNLO15d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO16::unitarityNLO16(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO16::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;

    double uniNLO16a = -la2Q/(16.0*M_PI);
    double uniNLO16b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2-WFRc3+2.0*WFRc4);

    return (uniNLO16a+uniNLO16b+uniNLO16c+uniNLO16d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO17::unitarityNLO17(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO17::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO17a = -la5Q/(16.0*M_PI);
    double uniNLO17b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc4;

    return (uniNLO17a+uniNLO17b+uniNLO17c+uniNLO17d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO18::unitarityNLO18(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO18::computeThValue()
{
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;

    double uniNLO18a = -la1Q/(16.0*M_PI);
    double uniNLO18b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    return (uniNLO18a+uniNLO18b+uniNLO18c+uniNLO18d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO19::unitarityNLO19(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO19::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;

    double uniNLO19a = -la2Q/(16.0*M_PI);
    double uniNLO19b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    return (uniNLO19a+uniNLO19b+uniNLO19c+uniNLO19d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO20::unitarityNLO20(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO20::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO20a = -la5Q/(16.0*M_PI);
    double uniNLO20b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO20a+uniNLO20b+uniNLO20c+uniNLO20d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO21::unitarityNLO21(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO21::computeThValue()
{
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;

    double uniNLO21a = -la1Q/(16.0*M_PI);
    double uniNLO21b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1+2.0*WFRc2-WFRc3-2.0*WFRc4);

    return (uniNLO21a+uniNLO21b+uniNLO21c+uniNLO21d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO22::unitarityNLO22(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO22::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;

    double uniNLO22a = -la2Q/(16.0*M_PI);
    double uniNLO22b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2+WFRc3-2.0*WFRc4);

    return (uniNLO22a+uniNLO22b+uniNLO22c+uniNLO22d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO23::unitarityNLO23(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO23::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO23a = -la5Q/(16.0*M_PI);
    double uniNLO23b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * (WFRc2-WFRc4);

    return (uniNLO23a+uniNLO23b+uniNLO23c+uniNLO23d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO24::unitarityNLO24(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO24::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO24a = -(la3Q+la4Q)/(16.0*M_PI);
    double uniNLO24b = 3.0*(betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO24c = (gslpp::complex::i()*M_PI-1.0)*(la3Q+la4Q)*(la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO24d = -2.0*(la3Q+la4Q)/(512.0*M_PI*M_PI*M_PI) * WFRc4;

    return (uniNLO24a+uniNLO24b+uniNLO24c+uniNLO24d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO25::unitarityNLO25(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO25::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO25a = -(la3Q+la4Q)/(16.0*M_PI);
    double uniNLO25b = 3.0*(betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO25c = (gslpp::complex::i()*M_PI-1.0)*(la3Q+la4Q)*(la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO25d = -(la3Q+la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    return (uniNLO25a+uniNLO25b+uniNLO25c+uniNLO25d - gslpp::complex::i()/2.0).abs();
}

unitarityNLO26::unitarityNLO26(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLO26::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO26a = -(la3Q+la4Q)/(16.0*M_PI);
    double uniNLO26b = 3.0*(betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO26c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO26d = -2.0*(la3Q+la4Q)/(512.0*M_PI*M_PI*M_PI) * (WFRc2-WFRc4);

    return (uniNLO26a+uniNLO26b+uniNLO26c+uniNLO26d - gslpp::complex::i()/2.0).abs();
}

// Eigenvalues:

unitarityNLOev1::unitarityNLOev1(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev1::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

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

    double uniNLO1a = -3.0*la1Q/(16.0*M_PI);
    double uniNLO1b = 9.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la1Q*la1Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1d = -3.0*la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO2a = -3.0*la2Q/(16.0*M_PI);
    double uniNLO2b = 9.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la2Q*la2Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2d = -3.0*la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO3a = -(2.0*la3Q+la4Q)/(16.0*M_PI);
    double uniNLO3b = 3.0*(2.0*betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*(2.0*la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3d = -(2.0*la3Q+la4Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO1a+uniNLO1b+uniNLO1c+uniNLO1d;
    gslpp::complex uniB=uniNLO2a+uniNLO2b+uniNLO2c+uniNLO2d;
    gslpp::complex uniC=uniNLO3a+uniNLO3b+uniNLO3c+uniNLO3d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev2::unitarityNLOev2(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev2::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

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

    double uniNLO1a = -3.0*la1Q/(16.0*M_PI);
    double uniNLO1b = 9.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la1Q*la1Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO1d = -3.0*la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO2a = -3.0*la2Q/(16.0*M_PI);
    double uniNLO2b = 9.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2c = (gslpp::complex::i()*M_PI-1.0)*(9.0*la2Q*la2Q+(2.0*la3Q+la4Q)*(2.0*la3Q+la4Q))/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO2d = -3.0*la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO3a = -(2.0*la3Q+la4Q)/(16.0*M_PI);
    double uniNLO3b = 3.0*(2.0*betalambda3+betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*(2.0*la3Q+la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO3d = -(2.0*la3Q+la4Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO1a+uniNLO1b+uniNLO1c+uniNLO1d;
    gslpp::complex uniB=uniNLO2a+uniNLO2b+uniNLO2c+uniNLO2d;
    gslpp::complex uniC=uniNLO3a+uniNLO3b+uniNLO3c+uniNLO3d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev3::unitarityNLOev3(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev3::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO4a = -(la3Q+2.0*la4Q)/(16.0*M_PI);
    double uniNLO4b = 3.0*(betalambda3+2.0*betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+4.0*la3Q*la4Q+4.0*la4Q*la4Q+9.0*la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4d = -(la3Q+la4Q+la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    double uniNLO5a = -3.0*la5Q/(16.0*M_PI);
    double uniNLO5b = 9.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la3Q+2.0*la4Q)*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5d = -(la4Q+2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO4a+uniNLO4b+uniNLO4c+uniNLO4d;
    gslpp::complex uniC=uniNLO5a+uniNLO5b+uniNLO5c+uniNLO5d;

    return (uniA+uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev4::unitarityNLOev4(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev4::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO4a = -(la3Q+2.0*la4Q)/(16.0*M_PI);
    double uniNLO4b = 3.0*(betalambda3+2.0*betalambda4)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+4.0*la3Q*la4Q+4.0*la4Q*la4Q+9.0*la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO4d = -(la3Q+la4Q+la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    double uniNLO5a = -3.0*la5Q/(16.0*M_PI);
    double uniNLO5b = 9.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5c = 3.0*(gslpp::complex::i()*M_PI-1.0)*(la3Q+2.0*la4Q)*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO5d = -(la4Q+2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO4a+uniNLO4b+uniNLO4c+uniNLO4d;
    gslpp::complex uniC=uniNLO5a+uniNLO5b+uniNLO5c+uniNLO5d;

    return (uniA-uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev5::unitarityNLOev5(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev5::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO6a = -la1Q/(16.0*M_PI);
    double uniNLO6b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO7a = -la2Q/(16.0*M_PI);
    double uniNLO7b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO8a = -la4Q/(16.0*M_PI);
    double uniNLO8b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO6a+uniNLO6b+uniNLO6c+uniNLO6d;
    gslpp::complex uniB=uniNLO7a+uniNLO7b+uniNLO7c+uniNLO7d;
    gslpp::complex uniC=uniNLO8a+uniNLO8b+uniNLO8c+uniNLO8d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev6::unitarityNLOev6(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev6::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO6a = -la1Q/(16.0*M_PI);
    double uniNLO6b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO7a = -la2Q/(16.0*M_PI);
    double uniNLO7b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO8a = -la4Q/(16.0*M_PI);
    double uniNLO8b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO8d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO6a+uniNLO6b+uniNLO6c+uniNLO6d;
    gslpp::complex uniB=uniNLO7a+uniNLO7b+uniNLO7c+uniNLO7d;
    gslpp::complex uniC=uniNLO8a+uniNLO8b+uniNLO8c+uniNLO8d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev7::unitarityNLOev7(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev7::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO6a = -la1Q/(16.0*M_PI);
    double uniNLO6b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO7a = -la2Q/(16.0*M_PI);
    double uniNLO7b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO9a = -la4Q/(16.0*M_PI);
    double uniNLO9b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9d = -la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO6a+uniNLO6b+uniNLO6c+uniNLO6d;
    gslpp::complex uniB=uniNLO7a+uniNLO7b+uniNLO7c+uniNLO7d;
    gslpp::complex uniC=uniNLO9a+uniNLO9b+uniNLO9c+uniNLO9d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev8::unitarityNLOev8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev8::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda4 = la1Q*la4Q + la2Q*la4Q + 4.0*la3Q*la4Q + 2.0*la4Q*la4Q + 4.0*la5Q*la5Q
                          + 3.0*la4Q*Yb1Q*Yb1Q + 3.0*la4Q*Yb2Q*Yb2Q + la4Q*Ytau1Q*Ytau1Q + la4Q*Ytau2Q*Ytau2Q + 3.0*la4Q*YtQ*YtQ
                          - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q - YtQ*YtQ);

    double uniNLO6a = -la1Q/(16.0*M_PI);
    double uniNLO6b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO6d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO7a = -la2Q/(16.0*M_PI);
    double uniNLO7b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la4Q*la4Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO7d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO9a = -la4Q/(16.0*M_PI);
    double uniNLO9b = 3.0*betalambda4/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la4Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO9d = -la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO6a+uniNLO6b+uniNLO6c+uniNLO6d;
    gslpp::complex uniB=uniNLO7a+uniNLO7b+uniNLO7c+uniNLO7d;
    gslpp::complex uniC=uniNLO9a+uniNLO9b+uniNLO9c+uniNLO9d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev9::unitarityNLOev9(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev9::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO10a = -la3Q/(16.0*M_PI);
    double uniNLO10b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10d = -(la3Q+la4Q-la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    double uniNLO11a = -la5Q/(16.0*M_PI);
    double uniNLO11b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11d = -(la4Q-2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO10a+uniNLO10b+uniNLO10c+uniNLO10d;
    gslpp::complex uniC=uniNLO11a+uniNLO11b+uniNLO11c+uniNLO11d;

    return (uniA+uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev10::unitarityNLOev10(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev10::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO10a = -la3Q/(16.0*M_PI);
    double uniNLO10b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO10d = -(la3Q+la4Q-la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    double uniNLO11a = -la5Q/(16.0*M_PI);
    double uniNLO11b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO11d = -(la4Q-2.0*la5Q)/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO10a+uniNLO10b+uniNLO10c+uniNLO10d;
    gslpp::complex uniC=uniNLO11a+uniNLO11b+uniNLO11c+uniNLO11d;

    return (uniA-uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev11::unitarityNLOev11(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev11::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO12a = -la3Q/(16.0*M_PI);
    double uniNLO12b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12d = -la3Q/(512.0*M_PI*M_PI*M_PI) * WFRc3;

    double uniNLO13a = -la5Q/(16.0*M_PI);
    double uniNLO13b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO12a+uniNLO12b+uniNLO12c+uniNLO12d;
    gslpp::complex uniC=uniNLO13a+uniNLO13b+uniNLO13c+uniNLO13d;

    return (uniA+uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev12::unitarityNLOev12(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev12::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;

    double betalambda3 = 3.0*la1Q*la3Q + 3.0*la2Q*la3Q + 2.0*la3Q*la3Q + la1Q*la4Q + la2Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 3.0*la3Q*Yb1Q*Yb1Q + 3.0*la3Q*Yb2Q*Yb2Q + la3Q*Ytau1Q*Ytau1Q + la3Q*Ytau2Q*Ytau2Q + 3.0*la3Q*YtQ*YtQ
                          - 6.0*Yb1Q*Yb1Q*(Yb2Q*Yb2Q + YtQ*YtQ) - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO12a = -la3Q/(16.0*M_PI);
    double uniNLO12b = 3.0*betalambda3/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12c = (gslpp::complex::i()*M_PI-1.0)*(la3Q*la3Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO12d = -la3Q/(512.0*M_PI*M_PI*M_PI) * WFRc3;

    double uniNLO13a = -la5Q/(16.0*M_PI);
    double uniNLO13b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13c = (gslpp::complex::i()*M_PI-1.0)*la3Q*la5Q/(128.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO13d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO12a+uniNLO12b+uniNLO12c+uniNLO12d;
    gslpp::complex uniC=uniNLO13a+uniNLO13b+uniNLO13c+uniNLO13d;

    return (uniA-uniC - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev13::unitarityNLOev13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev13::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO15a = -la1Q/(16.0*M_PI);
    double uniNLO15b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1-2.0*WFRc2+WFRc3+2.0*WFRc4);

    double uniNLO16a = -la2Q/(16.0*M_PI);
    double uniNLO16b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2-WFRc3+2.0*WFRc4);

    double uniNLO17a = -la5Q/(16.0*M_PI);
    double uniNLO17b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc4;

    gslpp::complex uniA=uniNLO15a+uniNLO15b+uniNLO15c+uniNLO15d;
    gslpp::complex uniB=uniNLO16a+uniNLO16b+uniNLO16c+uniNLO16d;
    gslpp::complex uniC=uniNLO17a+uniNLO17b+uniNLO17c+uniNLO17d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev14::unitarityNLOev14(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev14::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO15a = -la1Q/(16.0*M_PI);
    double uniNLO15b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO15d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1-2.0*WFRc2+WFRc3+2.0*WFRc4);

    double uniNLO16a = -la2Q/(16.0*M_PI);
    double uniNLO16b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO16d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2-WFRc3+2.0*WFRc4);

    double uniNLO17a = -la5Q/(16.0*M_PI);
    double uniNLO17b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO17d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * WFRc4;

    gslpp::complex uniA=uniNLO15a+uniNLO15b+uniNLO15c+uniNLO15d;
    gslpp::complex uniB=uniNLO16a+uniNLO16b+uniNLO16c+uniNLO16d;
    gslpp::complex uniC=uniNLO17a+uniNLO17b+uniNLO17c+uniNLO17d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev15::unitarityNLOev15(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev15::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO18a = -la1Q/(16.0*M_PI);
    double uniNLO18b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO19a = -la2Q/(16.0*M_PI);
    double uniNLO19b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO20a = -la5Q/(16.0*M_PI);
    double uniNLO20b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO18a+uniNLO18b+uniNLO18c+uniNLO18d;
    gslpp::complex uniB=uniNLO19a+uniNLO19b+uniNLO19c+uniNLO19d;
    gslpp::complex uniC=uniNLO20a+uniNLO20b+uniNLO20c+uniNLO20d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev16::unitarityNLOev16(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev16::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO18a = -la1Q/(16.0*M_PI);
    double uniNLO18b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO18d = -la1Q/(512.0*M_PI*M_PI*M_PI) * WFRc1;

    double uniNLO19a = -la2Q/(16.0*M_PI);
    double uniNLO19b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO19d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2);

    double uniNLO20a = -la5Q/(16.0*M_PI);
    double uniNLO20b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO20d = -la4Q/(512.0*M_PI*M_PI*M_PI) * WFRc2;

    gslpp::complex uniA=uniNLO18a+uniNLO18b+uniNLO18c+uniNLO18d;
    gslpp::complex uniB=uniNLO19a+uniNLO19b+uniNLO19c+uniNLO19d;
    gslpp::complex uniC=uniNLO20a+uniNLO20b+uniNLO20c+uniNLO20d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev17::unitarityNLOev17(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev17::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO21a = -la1Q/(16.0*M_PI);
    double uniNLO21b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1+2.0*WFRc2-WFRc3-2.0*WFRc4);

    double uniNLO22a = -la2Q/(16.0*M_PI);
    double uniNLO22b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2+WFRc3-2.0*WFRc4);

    double uniNLO23a = -la5Q/(16.0*M_PI);
    double uniNLO23b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * (WFRc2-WFRc4);

    gslpp::complex uniA=uniNLO21a+uniNLO21b+uniNLO21c+uniNLO21d;
    gslpp::complex uniB=uniNLO22a+uniNLO22b+uniNLO22c+uniNLO22d;
    gslpp::complex uniC=uniNLO23a+uniNLO23b+uniNLO23c+uniNLO23d;

    return (0.5*(uniA+uniB+sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}

unitarityNLOev18::unitarityNLOev18(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double unitarityNLOev18::computeThValue()
{
    double YtQ = myTHDM.getMyTHDMCache()->Ytop_at_Q;
    double Yb1Q = myTHDM.getMyTHDMCache()->Ybottom1_at_Q;
    double Yb2Q = myTHDM.getMyTHDMCache()->Ybottom2_at_Q;
    double Ytau1Q = myTHDM.getMyTHDMCache()->Ytau1_at_Q;
    double Ytau2Q = myTHDM.getMyTHDMCache()->Ytau2_at_Q;
    double la1Q = myTHDM.getMyTHDMCache()->lambda1_at_Q;
    double la2Q = myTHDM.getMyTHDMCache()->lambda2_at_Q;
    double la3Q = myTHDM.getMyTHDMCache()->lambda3_at_Q;
    double la4Q = myTHDM.getMyTHDMCache()->lambda4_at_Q;
    double la5Q = myTHDM.getMyTHDMCache()->lambda5_at_Q;
    double WFRc1 = myTHDM.getMyTHDMCache()->WFRcomb1;
    double WFRc2 = myTHDM.getMyTHDMCache()->WFRcomb2;
    double WFRc3 = myTHDM.getMyTHDMCache()->WFRcomb3;
    double WFRc4 = myTHDM.getMyTHDMCache()->WFRcomb4;

    double betalambda1 = 6.0*la1Q*la1Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q
                          + 6.0*la1Q*Yb1Q*Yb1Q + 2.0*la1Q*Ytau1Q*Ytau1Q
                          - 6.0*Yb1Q*Yb1Q*Yb1Q*Yb1Q - 2.0*Ytau1Q*Ytau1Q*Ytau1Q*Ytau1Q;
    double betalambda2 = 6.0*la2Q*la2Q + 2.0*la3Q*la3Q + 2.0*la3Q*la4Q + la4Q*la4Q + la5Q*la5Q 
                          + 6.0*la2Q*Yb2Q*Yb2Q + 2.0*la2Q*Ytau2Q*Ytau2Q + 6.0*la2Q*YtQ*YtQ
                          - 6.0*Yb2Q*Yb2Q*Yb2Q*Yb2Q - 2.0*Ytau2Q*Ytau2Q*Ytau2Q*Ytau2Q - 6.0*YtQ*YtQ*YtQ*YtQ;
    double betalambda5 = la5Q*(la1Q + la2Q + 4.0*la3Q + 6.0*la4Q)
                          + 3.0*la5Q*Yb1Q*Yb1Q + la5Q*Ytau1Q*Ytau1Q + la5Q*(3.0*Yb2Q*Yb2Q + Ytau2Q*Ytau2Q + 3.0*YtQ*YtQ)
                          - 6.0*Yb1Q*Yb1Q*Yb2Q*Yb2Q - 2.0*Ytau1Q*Ytau1Q*Ytau2Q*Ytau2Q;

    double uniNLO21a = -la1Q/(16.0*M_PI);
    double uniNLO21b = 3.0*betalambda1/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21c = (gslpp::complex::i()*M_PI-1.0)*(la1Q*la1Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO21d = -la1Q/(512.0*M_PI*M_PI*M_PI) * (WFRc1+2.0*WFRc2-WFRc3-2.0*WFRc4);

    double uniNLO22a = -la2Q/(16.0*M_PI);
    double uniNLO22b = 3.0*betalambda2/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22c = (gslpp::complex::i()*M_PI-1.0)*(la2Q*la2Q+la5Q*la5Q)/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO22d = -la2Q/(512.0*M_PI*M_PI*M_PI) * (-WFRc1+2.0*WFRc2+WFRc3-2.0*WFRc4);

    double uniNLO23a = -la5Q/(16.0*M_PI);
    double uniNLO23b = 3.0*betalambda5/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23c = (gslpp::complex::i()*M_PI-1.0)*(la1Q+la2Q)*la5Q/(256.0*M_PI*M_PI*M_PI);
    gslpp::complex uniNLO23d = -2.0*la5Q/(512.0*M_PI*M_PI*M_PI) * (WFRc2-WFRc4);

    gslpp::complex uniA=uniNLO21a+uniNLO21b+uniNLO21c+uniNLO21d;
    gslpp::complex uniB=uniNLO22a+uniNLO22b+uniNLO22c+uniNLO22d;
    gslpp::complex uniC=uniNLO23a+uniNLO23b+uniNLO23c+uniNLO23d;

    return (0.5*(uniA+uniB-sqrt((uniA-uniB)*(uniA-uniB)+4.0*uniC*uniC)) - 0.5*gslpp::complex::i()).abs();
}
