/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
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
////    double mHl=myTHDM->getMHl();
////    double mA=myTHDM->getMA();
////    double mHh=myTHDM->getMHh();
////    double mHp=myTHDM->getMHp();
////    double vev=myTHDM->v();
////    double sina=myTHDM->computeSina();
////    double cosa=myTHDM->computeCosa();
////    double tanb=myTHDM->getTanb();
////    double sinb=myTHDM->getSinb();
////    double cosb=myTHDM->getCosb();
////    double m12_2=myTHDM->getM12_2();
////
////    double unitarity1=(3.*((mHh*mHh/cosb/cosb +mHl*mHl/sinb/sinb)*cosa*cosa
////                          +(mHl*mHl/cosb/cosb + mHh*mHh/sinb/sinb)*sina*sina
////                          -m12_2*(1/tanb/sinb/sinb + tanb/cosb/cosb))
////                       +sqrt(4.*pow(mA*mA + 2.*mHp*mHp
////                                      -(m12_2 +2.*(mHl*mHl - mHh*mHh)*cosa*sina)/(cosb*sinb),2)
////                             +(9.*pow((mHl*mHl - mHh*mHh)*(cosa*cosa - sina*sina) 
////                                        +(mHl*mHl + mHh*mHh -(2.*m12_2)/(cosb*sinb))
////                                         *(cosb*cosb - sinb*sinb),2))
////                              /(4.*pow(cosb,4)*pow(sinb,4)))
////                       )/vev/vev/2.;
////
////    double unitarity2=(3.*((mHh*mHh/cosb/cosb +mHl*mHl/sinb/sinb)*cosa*cosa
////                          +(mHl*mHl/cosb/cosb + mHh*mHh/sinb/sinb)*sina*sina
////                          -m12_2*(1/tanb/sinb/sinb + tanb/cosb/cosb))
////                       -sqrt(4.*pow(mA*mA + 2.*mHp*mHp
////                                      -(m12_2 +2.*(mHl*mHl - mHh*mHh)*cosa*sina)/(cosb*sinb),2)
////                             +(9.*pow((mHl*mHl - mHh*mHh)*(cosa*cosa - sina*sina) 
////                                        +(mHl*mHl + mHh*mHh -(2.*m12_2)/(cosb*sinb))
////                                         *(cosb*cosb - sinb*sinb),2))
////                              /(4.*pow(cosb,4)*pow(sinb,4)))
////                       )/vev/vev/2.;
////
////    double unitarity3=((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)/sinb/sinb 
////                      +(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
////                      +sqrt(4.*pow(mA*mA -2*mHp*mHp +m12_2/(cosb*sinb),2)
////                            +pow((mHl*mHl/cosb/cosb -mHh*mHh/sinb/sinb)*sina*sina
////                                   +(mHh*mHh/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
////                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/vev/vev/2.;
////
////    double unitarity4=((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)/sinb/sinb
////                      +(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
////                      -sqrt(4.*pow(mA*mA -2*mHp*mHp +m12_2/(cosb*sinb),2)
////                            +pow((mHl*mHl/cosb/cosb -mHh*mHh/sinb/sinb)*sina*sina
////                                   +(mHh*mHh/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
////                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/vev/vev/2.;
////
////    double unitarity5=((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)/sinb/sinb
////                      +(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
////                      +sqrt(4.*pow(mA*mA -m12_2/(cosb*sinb),2)
////                            +pow((mHl*mHl/cosb/cosb -mHh*mHh/sinb/sinb)*sina*sina
////                                   +(mHh*mHh/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
////                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/vev/vev/2.;
////
////    double unitarity6=((mHl*mHl*cosa*cosa +mHh*mHh*sina*sina -m12_2/tanb)/sinb/sinb
////                      +(mHh*mHh*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
////                      -sqrt(4.*pow(mA*mA -m12_2/(cosb*sinb),2)
////                            +pow((mHl*mHl/cosb/cosb -mHh*mHh/sinb/sinb)*sina*sina
////                                   +(mHh*mHh/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
////                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/vev/vev/2.;
////
////    double unitarity7=(5.*mA*mA -2.*mHp*mHp + ((mHh*mHh-mHl*mHl)*cosa*sina -2.*m12_2)/(cosb*sinb))/vev/vev;
////
////    double unitarity8=(mA*mA +2.*mHp*mHp +((mHh*mHh-mHl*mHl)*cosa*sina -2.*m12_2)/(cosb*sinb))/vev/vev;
////
////    double unitarity9=(-mA*mA -2.*mHp*mHp +((mHh*mHh-mHl*mHl)*cosa*sina +4.*m12_2)/(cosb*sinb))/vev/vev;
////
////    double unitarity10=(-mA*mA +2.*mHp*mHp +((mHh*mHh-mHl*mHl)*cosa*sina)/(cosb*sinb))/vev/vev;
////
////    double unitarity11=(mA*mA +((mHh*mHh-mHl*mHl)*cosa*sina)/(cosb*sinb))/vev/vev;
////
////    double unitarity12=(-mA*mA +4.*mHp*mHp +((mHh*mHh-mHl*mHl)*cosa*sina -2.*m12_2)/(cosb*sinb))/vev/vev;
////
////    if (obs ==1) return( unitarity1);
////    if (obs ==2) return( unitarity2);
////    if (obs ==3) return( unitarity3);
////    if (obs ==4) return( unitarity4);
////    if (obs ==5) return( unitarity5);
////    if (obs ==6) return( unitarity6);
////    if (obs ==7) return( unitarity7);
////    if (obs ==8) return( unitarity8);
////    if (obs ==9) return( unitarity9);
////    if (obs ==10) return( unitarity10);
////    if (obs ==11) return( unitarity11);
////    if (obs ==12) return( unitarity12);
////
////    throw std::runtime_error("unitarity::computeThValue(): Observable type not "
////            "defined. Can be only any of (1,2,3,4,5,6,7,8,9,10,11,12)");
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

unitarity2::unitarity2(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity2::computeThValue()
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
                      +sqrt(4.*pow(mA2 -2*mHp2 +m12_2/(cosb*sinb),2)
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
                      -sqrt(4.*pow(mA2 -2*mHp2 +m12_2/(cosb*sinb),2)
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
    double vev=myTHDM->v();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mHl*mHl*cosa*cosa +mHh2*sina*sina -m12_2/tanb)/sinb/sinb
                      +(mHh2*cosa*cosa +mHl*mHl*sina*sina - m12_2*tanb)/cosb/cosb
                      +sqrt(4.*pow(mA2 -m12_2/(cosb*sinb),2)
                            +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                                   +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity6::unitarity6(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity6::computeThValue()
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
                      -sqrt(4.*pow(mA2 -m12_2/(cosb*sinb),2)
                            +pow((mHl*mHl/cosb/cosb -mHh2/sinb/sinb)*sina*sina
                                   +(mHh2/cosb/cosb -mHl*mHl/sinb/sinb)*cosa*cosa
                                   +m12_2/tanb/sinb/sinb -m12_2*tanb/cosb/cosb,2)))/(2.0*vev*vev);
}

unitarity7::unitarity7(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity7::computeThValue()
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

    return (5.*mA2 -2.*mHp2 + ((mHh2-mHl*mHl)*cosa*sina -2.*m12_2)/(cosb*sinb))/(vev*vev);
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

    return (mA2 +2.*mHp2 +((mHh2-mHl*mHl)*cosa*sina -2.*m12_2)/(cosb*sinb))/(vev*vev);
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
    double m12_2=myTHDM->getm12_2();

    return (-mA2 -2.*mHp2 +((mHh2-mHl*mHl)*cosa*sina +4.*m12_2)/(cosb*sinb))/(vev*vev);
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

    return (-mA2 +2.*mHp2 +((mHh2-mHl*mHl)*cosa*sina)/(cosb*sinb))/(vev*vev);
}

unitarity11::unitarity11(const StandardModel& SM_i)
: unitarity(SM_i)
{}

double unitarity11::computeThValue()
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

    return (-mA2 +4.0*mHp2 +((mHh2-mHl*mHl)*cosa*sina -2.0*m12_2)/(cosb*sinb))/(vev*vev);
}
