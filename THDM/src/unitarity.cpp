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
