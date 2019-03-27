/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMquantities.h"
#include "THDM.h"

tanbeta::tanbeta(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double tanbeta::computeThValue()
{
    return myTHDM->gettanb();
}

mass_mHl::mass_mHl(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double mass_mHl::computeThValue()
{
    return sqrt(myTHDM->getmHl2());
}

mass_mHh::mass_mHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double mass_mHh::computeThValue()
{
    return myTHDM->getmHh();
}


mass_mA::mass_mA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double mass_mA::computeThValue()
{
    return myTHDM->getmA();
}


mass_mHp::mass_mHp(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double mass_mHp::computeThValue()
{
    return myTHDM->getmHp();
}


massdifference_mHlmmA::massdifference_mHlmmA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHlmmA::computeThValue()
{
    return sqrt(myTHDM->getmHl2()) - myTHDM->getmA();
}


massdifference_mAmmHl::massdifference_mAmmHl(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mAmmHl::computeThValue()
{
    return myTHDM->getmA() - sqrt(myTHDM->getmHl2());
}


massdifference_mHlmmHp::massdifference_mHlmmHp(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHlmmHp::computeThValue()
{
    return sqrt(myTHDM->getmHl2()) - myTHDM->getmHp();
}


massdifference_mHpmmHl::massdifference_mHpmmHl(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHpmmHl::computeThValue()
{
    return myTHDM->getmHp() - sqrt(myTHDM->getmHl2());
}


massdifference_mHhmmA::massdifference_mHhmmA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHhmmA::computeThValue()
{
    return myTHDM->getmHh() - myTHDM->getmA();
}


massdifference_mAmmHh::massdifference_mAmmHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mAmmHh::computeThValue()
{
    return myTHDM->getmA() - myTHDM->getmHh();
}


massdifference_mHhmmHp::massdifference_mHhmmHp(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHhmmHp::computeThValue()
{
    return myTHDM->getmHh() - myTHDM->getmHp();
}


massdifference_mHpmmHh::massdifference_mHpmmHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHpmmHh::computeThValue()
{
    return myTHDM->getmHp() - myTHDM->getmHh();
}


massdifference_mAmmHp::massdifference_mAmmHp(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mAmmHp::computeThValue()
{
    return myTHDM->getmA() - myTHDM->getmHp();
}


massdifference_mHpmmA::massdifference_mHpmmA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double massdifference_mHpmmA::computeThValue()
{
    return myTHDM->getmHp() - myTHDM->getmA();
}


m11_2::m11_2(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double m11_2::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double tanb=myTHDM->gettanb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (-(cosa*cosa*mHh2) - mHl2*sina*sina + 2.0*m12_2*tanb + cosa*(mHl2 - mHh2)*sina*tanb)/2.;
}


m22_2::m22_2(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double m22_2::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double tanb=myTHDM->gettanb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (2.0*m12_2 + sina*cosa*(mHl2 - mHh2) - (cosa*cosa*mHl2 + sina*sina*mHh2)*tanb)/(2.0*tanb);
}


lambda1::lambda1(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda1::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double cosb=myTHDM->getcosb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (mHh2*cosa*cosa+mHl2*sina*sina-m12_2*tanb)/(vev*vev*cosb*cosb);
}


lambda2::lambda2(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda2::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (mHh2*sina*sina+mHl2*cosa*cosa-m12_2/tanb)/(vev*vev*sinb*sinb);
}


lambda3::lambda3(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda3::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return ((mHh2-mHl2)*cosa*sina+2.0*mHp2*sinb*cosb-m12_2)/(vev*vev*sinb*cosb);
}


lambda4::lambda4(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda4::computeThValue()
{
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return ((mA2-2.0*mHp2)*sinb*cosb+m12_2)/(vev*vev*sinb*cosb);
}


lambda5::lambda5(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda5::computeThValue()
{
    double mA2=myTHDM->getmA2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double m12_2=myTHDM->getm12_2();

    return (m12_2-mA2*sinb*cosb)/(vev*vev*sinb*cosb);
}


lambda345::lambda345(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double lambda345::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (m12_2+(mHh2-mHl2)*cosa*sina)/(vev*vev*sinb*cosb);
}


g_hhh::g_hhh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_hhh::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return 3.0*( -mHl2*(cosa*cosa*cosa/sinb - sina*sina*sina/cosb)
                 +cos(bma)*m12_2*(cosa*cosa/(sinb*sinb) - sina*sina/(cosb*cosb)) )/vev;
}


g_hhHh::g_hhHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_hhHh::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return -cos(bma)*(sinb*cosb*m12_2 
                      +sina*cosa*(sinb*cosb*(2.0*mHl2 + mHh2) - 3.0*m12_2))/(sinb*sinb*cosb*cosb*vev);
}


g_hHhHh::g_hHhHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_hHhHh::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return -sin(bma)*(sinb*cosb*m12_2
                      -sina*cosa*(sinb*cosb*(mHl2 + 2.0*mHh2) - 3.0*m12_2))/(sinb*sinb*cosb*cosb*vev);
}


g_HhHhHh::g_HhHhHh(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_HhHhHh::computeThValue()
{
    double mHh2=myTHDM->getmHh2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return -3.0*( (cosa*cosa*cosa*sinb + sina*sina*sina*cosb)*sinb*cosb*mHh2
                  -(cosa*cosa*sinb*sinb - cosb*cosb*sina*sina)*m12_2*sin(bma) )/(sinb*sinb*cosb*cosb*vev);
}


g_hAA::g_hAA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_hAA::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mA2=myTHDM->getmA2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ( (cosa*cosb-sina*sinb)*(m12_2-sinb*cosb*mHl2)
             +(mHl2-2.0*mA2)*sinb*sinb*cosb*cosb*sin(bma) )/(sinb*sinb*cosb*cosb*vev);
}


g_HhAA::g_HhAA(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_HhAA::computeThValue()
{
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ( (cosa*sinb+cosb*sina)*(m12_2-sinb*cosb*mHh2)
             +(mHh2-2.0*mA2)*sinb*sinb*cosb*cosb*cos(bma) )/(sinb*sinb*cosb*cosb*vev);
}


g_hHpHm::g_hHpHm(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_hHpHm::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ( (cosa*cosb-sina*sinb)*(m12_2-sinb*cosb*mHl2)
             +(mHl2-2.0*mHp2)*sinb*sinb*cosb*cosb*sin(bma) )/(sinb*sinb*cosb*cosb*vev);
}


g_HhHpHm::g_HhHpHm(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double g_HhHpHm::computeThValue()
{
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ( (cosa*sinb+cosb*sina)*(m12_2-sinb*cosb*mHh2)
             +(mHh2-2.0*mHp2)*sinb*sinb*cosb*cosb*cos(bma) )/(sinb*sinb*cosb*cosb*vev);
}


Y1_THDM::Y1_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Y1_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();

    return 0.25*(-mHl2-mHh2+(mHl2-mHh2)*cos(2.0*bma));
}


Y2_THDM::Y2_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Y2_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double m12_2=myTHDM->getm12_2();

    
    return 0.125*((mHl2-mHh2)*(3.0*sin2a-sin2b*cos(2.0*bma)
                                  -sin(2.0*bma)*(cosb*cosb-sinb*sinb))
                  +2.0*(8.0*m12_2-(mHl2+mHh2)*sin2b))/sin2b;
}


Y3_THDM::Y3_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Y3_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();

    return 0.25*(mHh2-mHl2)*sin(2.0*bma);
}


Z1_THDM::Z1_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z1_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();

    return (mHl2+mHh2+(mHh2-mHl2)*cos(2.0*bma))/(2.0*vev*vev);
}


Z2_THDM::Z2_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z2_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return (-128.0*m12_2*cos2b*cos2b
            +4.0*(mHl2+mHh2)*(5.0+3.0*(cos2b*cos2b-sin2b*sin2b))*sin2b
            +2.0*(mHl2-mHh2)*((6.0+10.0*(cos2b*cos2b-sin2b*sin2b))*cos(2.0*bma)
                                 +16.0*sin2b*cos2b*sin(2.0*bma))*sin2b)
           /(16.0*sin2b*sin2b*sin2b*vev*vev);
}


Z3_THDM::Z3_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z3_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double bma=myTHDM->getbma();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double m12_2=myTHDM->getm12_2();

    return ((mHh2-mHl2)*(3.0*sin2a-sin(2.0*bma)*cos2b-cos(2.0*bma)*sin2b)
            +2.0*(-8.0*m12_2+(mHh2+mHl2+4.0*mHp2)*sin2b))/(4.0*sin2b*vev*vev);
}


Z4_THDM::Z4_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z4_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double bma=myTHDM->getbma();

    return (2.0*mA2+mHl2+mHh2-4.0*mHp2+(mHl2-mHh2)*cos(2.0*bma))/(2.0*vev*vev);
}


Z5_THDM::Z5_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z5_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double bma=myTHDM->getbma();

    return (-2.0*mA2+mHl2+mHh2+(mHl2-mHh2)*cos(2.0*bma))/(2.0*vev*vev);
}


Z6_THDM::Z6_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z6_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();

    return ((mHl2-mHh2)*sin(2.0*bma))/(2.0*vev*vev);
}


Z7_THDM::Z7_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double Z7_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double m12_2=myTHDM->getm12_2();

    return (4.0*cos2b*(-4.0*m12_2+(mHl2+mHh2)*sin2b)
            +(mHl2-mHh2)*(3.0*(cosa*cosa-sina*sina)
                             +cos(2.0*bma)*cos2b-sin(2.0*bma)*sin2b)*sin2b)
           /(4.0*sin2b*sin2b*vev*vev);
}


xi0_THDM::xi0_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double xi0_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double m12_2=myTHDM->getm12_2();

    return -((mHh2+mHl2)*sin2b+(mHh2-mHl2)*sin2a-4.0*m12_2)/(4.0*sin2b);
}


xi1_THDM::xi1_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double xi1_THDM::computeThValue()
{
    double m12_2=myTHDM->getm12_2();

    return -m12_2;
}


xi3_THDM::xi3_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double xi3_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ((mHl2-mHh2)*sin(2.0*bma)-4.0*m12_2*cos2b)/(4.0*sin2b);
}


eta00_THDM::eta00_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double eta00_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double cos2a=cosa*cosa-sina*sina;
    double m12_2=myTHDM->getm12_2();

    return ((mHh2+mHl2-(mHh2-mHl2)*(cos2a*cos2b-sin2a*sin2b)
             -4.0*m12_2/sin2b)/(sin2b*sin2b)
            +2.0*mHp2)/(4.0*vev*vev);
}


eta3_THDM::eta3_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double eta3_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double cos2a=cosa*cosa-sina*sina;
    double m12_2=myTHDM->getm12_2();

    return ((mHh2+mHl2)*(0.5*tanb-0.5/tanb)
            +(mHh2-mHl2)*cos2a/sin2b
            +4.0*m12_2*cos2b/(sin2b*sin2b))/(4.0*vev*vev*sin2b);
}


E11_THDM::E11_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double E11_THDM::computeThValue()
{
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double mHp2=myTHDM->getmHp2();
    double m12_2=myTHDM->getm12_2();

    return (m12_2/sin2b-0.5*mHp2)/(vev*vev);
}


E22_THDM::E22_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double E22_THDM::computeThValue()
{
    double vev=myTHDM->v();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();

    return (mA2-mHp2)/(2.0*vev*vev);
}


E33_THDM::E33_THDM(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double E33_THDM::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double mHp2=myTHDM->getmHp2();
    double bma=myTHDM->getbma();
    double m12_2=myTHDM->getm12_2();

    return ((mHh2+mHl2-(mHh2-mHl2)*cos(2.0*bma)
             -4.0*m12_2*cos2b*cos2b/sin2b)/(sin2b*sin2b)
            -2.0*mHp2)/(4.0*vev*vev);
}



HHlambda1::HHlambda1(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda1::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double cos2a=cosa*cosa-sina*sina;
    double m12_2=myTHDM->getm12_2();
    
    return (mHh2+mHl2+4.0*m12_2*cos2b/sin2b
            +(mHh2-mHl2)*(cos2a-sin2a/tanb))/(4.0*vev*vev*cosb*cosb);
}


HHlambda2::HHlambda2(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda2::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double tanb=myTHDM->gettanb();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double cos2b=cosb*cosb-sinb*sinb;
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double cos2a=cosa*cosa-sina*sina;
    double m12_2=myTHDM->getm12_2();
    
    return (mHh2+mHl2-4.0*m12_2*cos2b/sin2b
            -(mHh2-mHl2)*(cos2a+sin2a*tanb))/(4.0*vev*vev*sinb*sinb);
}


HHlambda3::HHlambda3(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda3::computeThValue()
{
    double mHl2=myTHDM->getmHl2();
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double mHh2=myTHDM->getmHh2();
    double sina=myTHDM->getsina();
    double cosa=myTHDM->getcosa();
    double sin2a=2.0*sina*cosa;
    double m12_2=myTHDM->getm12_2();
    
    return ((mHh2-mHl2)*sin2a-2.0*m12_2)/(2.0*vev*vev*sin2b);
}


HHlambda4::HHlambda4(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda4::computeThValue()
{
    double vev=myTHDM->v();
    double mHp2=myTHDM->getmHp2();
    
    return (2.0*mHp2)/(vev*vev);
}


HHlambda5::HHlambda5(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda5::computeThValue()
{
    double vev=myTHDM->v();
    double sinb=myTHDM->getsinb();
    double cosb=myTHDM->getcosb();
    double sin2b=2.0*sinb*cosb;
    double m12_2=myTHDM->getm12_2();
    
    return 4.0*m12_2/(vev*vev*sin2b);
}


HHlambda6::HHlambda6(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
{}

double HHlambda6::computeThValue()
{
    double vev=myTHDM->v();
    double mA2=myTHDM->getmA2();
    
    return 2.0*mA2/(vev*vev);
}
