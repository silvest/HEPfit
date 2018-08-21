/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMquantities.h"

tanbeta_GTHDM::tanbeta_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double tanbeta_GTHDM::computeThValue()
{
    return myGTHDM->gettanb();
}

mH1_2::mH1_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH1_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH1sq;
}

mH2_2::mH2_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH2sq;
}

mH3_2::mH3_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH3sq;
}

mH1_GTHDM::mH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH1_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH1_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH1_2);
        return sqrt(myGTHDM.getMyGTHDMCache()->mH1sq);
}

mH2_GTHDM::mH2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH2_2 < 0. || (myGTHDM.getMyGTHDMCache()->R32_GTHDM)*
//            ((myGTHDM.getMyGTHDMCache()->R13_GTHDM)*(myGTHDM.getMyGTHDMCache()->R22_GTHDM) -
//            (myGTHDM.getMyGTHDMCache()->R12_GTHDM)*(myGTHDM.getMyGTHDMCache()->R23_GTHDM)) < 0.00001)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH2_2);
        return sqrt(myGTHDM.getMyGTHDMCache()->mH2sq);
}

mH3_GTHDM::mH3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH3_2 < 0. || (myGTHDM.getMyGTHDMCache()->R33_GTHDM)*
//            ((myGTHDM.getMyGTHDMCache()->R13_GTHDM)*(myGTHDM.getMyGTHDMCache()->R22_GTHDM) -
//            (myGTHDM.getMyGTHDMCache()->R12_GTHDM)*(myGTHDM.getMyGTHDMCache()->R23_GTHDM)) < 0.00001)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH3_2);
        return sqrt(myGTHDM.getMyGTHDMCache()->mH3sq);
}

mHlight_GTHDM::mHlight_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHlight_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mHlight_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mHlight_2);
        return 0.0;
}

mHmedium_GTHDM::mHmedium_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHmedium_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mHmedium_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mHmedium_2);
        return 0.0;
}

mHheavy_GTHDM::mHheavy_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHheavy_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mHheavy_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mHheavy_2);
        return 0.0;
}

mHp_GTHDM::mHp_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHp_GTHDM::computeThValue()
{
    return sqrt(myGTHDM.getmHp2());
}

mH3mmH2_GTHDM::mH3mmH2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmH2_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH3_2 < 0. || myGTHDM.getMyGTHDMCache()->mH2_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH3_2) - sqrt(myGTHDM.getMyGTHDMCache()->mH2_2);
        return 0.0;
}

mH3mmHp_GTHDM::mH3mmHp_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmHp_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH3_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH3_2) - sqrt(myGTHDM.getMyGTHDMCache()->mHp2_GTHDM);
        return 0.0;
}

mH3mmH1_GTHDM::mH3mmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmH1_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH3_2 < 0. || myGTHDM.getMyGTHDMCache()->mH1_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH3_2) - sqrt(myGTHDM.getMyGTHDMCache()->mH1_2);
        return 0.0;
}

mH2mmHp_GTHDM::mH2mmHp_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2mmHp_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH2_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH2_2) - sqrt(myGTHDM.getMyGTHDMCache()->mHp2_GTHDM);   
        return 0.0;
}

mH2mmH1_GTHDM::mH2mmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2mmH1_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH2_2 < 0. || myGTHDM.getMyGTHDMCache()->mH1_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mH2_2) - sqrt(myGTHDM.getMyGTHDMCache()->mH1_2);
        return 0.0;
}

mHpmmH1_GTHDM::mHpmmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHpmmH1_GTHDM::computeThValue()
{
//    if(myGTHDM.getMyGTHDMCache()->mH1_2 < 0.)
//        return std::numeric_limits<double>::quiet_NaN();
//    else
//        return sqrt(myGTHDM.getMyGTHDMCache()->mHp2_GTHDM) - sqrt(myGTHDM.getMyGTHDMCache()->mH1_2);   
        return 0.0;
}


mH1sq_GTHDM::mH1sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH1sq_GTHDM::computeThValue()
{
//    return myGTHDM.getMyGTHDMCache()->mH1_2;
        return 0.0;
}

mH2sq_GTHDM::mH2sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2sq_GTHDM::computeThValue()
{
//    return myGTHDM.getMyGTHDMCache()->mH2_2;
        return 0.0;
}

mH3sq_GTHDM::mH3sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3sq_GTHDM::computeThValue()
{
//    return myGTHDM.getMyGTHDMCache()->mH3_2;
        return 0.0;
}

Msq11_GTHDM::Msq11_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq11_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M11_2;
}

Msq12_GTHDM::Msq12_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq12_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M12_2;
}

Msq13_GTHDM::Msq13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq13_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M13_2;
}

Msq22_GTHDM::Msq22_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq22_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M22_2;
}

Msq23_GTHDM::Msq23_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq23_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M23_2;
}

Msq33_GTHDM::Msq33_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Msq33_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->M33_2;
}

M2_GTHDM::M2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double M2_GTHDM::computeThValue()
{
//    return myGTHDM.getMyGTHDMCache()->M2_GTHDM;
        return 0.0;
}

m11_2_GTHDM::m11_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m11_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m11sq;
}

m22_2_GTHDM::m22_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m22_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m22sq;
}

Rem12_2_GTHDM::Rem12_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Rem12_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rem12sq;
}

Imm12_2_GTHDM::Imm12_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Imm12_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Imm12sq;
}

lambda1_GTHDM::lambda1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda1;
}

lambda2_GTHDM::lambda2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda2;
}

lambda3_GTHDM::lambda3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda3_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda3;
}

lambda4_GTHDM::lambda4_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda4_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda4;
}


R11_GTHDM::R11_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R11_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R11_GTHDM;
}


R12_GTHDM::R12_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R12_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R12_GTHDM;
}


R13_GTHDM::R13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R13_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R13_GTHDM;
}


R21_GTHDM::R21_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R21_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R21_GTHDM;
}


R22_GTHDM::R22_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R22_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R22_GTHDM;
}


R23_GTHDM::R23_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R23_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R23_GTHDM;
}


R31_GTHDM::R31_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R31_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R31_GTHDM;
}


R32_GTHDM::R32_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R32_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R32_GTHDM;
}


R33_GTHDM::R33_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R33_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R33_GTHDM;
}



v1_GTHDM::v1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double v1_GTHDM::computeThValue()
{
    double v = myGTHDM->v();
    double cosb = myGTHDM->getcosb();

    return (v*cosb);
}


v2_GTHDM::v2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double v2_GTHDM::computeThValue()
{
    double v = myGTHDM->v();
    double sinb = myGTHDM->getsinb();

    return (v*sinb);
}

//
//Resigmau::Resigmau(const StandardModel& SM_i)
//: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double Resigmau::computeThValue()
//{
//    double v = myGTHDM->v();
//    double Ytu_33r = myGTHDM->getYtu_33r();
//    double cosb = myGTHDM->getcosb();
//    double mtop = myGTHDM->getQuarks(QCD::TOP).getMass();
//    double tanb = myGTHDM->gettanb();
//
//    return v*Ytu_33r/(sqrt(2.)*cosb*mtop)-tanb;
//}

massdifference_mH1mmH2::massdifference_mH1mmH2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double massdifference_mH1mmH2::computeThValue()
{
    return  sqrt(myGTHDM->getMyGTHDMCache()->mH1sq) - sqrt(myGTHDM->getMyGTHDMCache()->mH2sq);
}

massdifference_mH1mmH3::massdifference_mH1mmH3(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}
double massdifference_mH1mmH3::computeThValue()
{
    return  sqrt(myGTHDM->getMyGTHDMCache()->mH1sq) - sqrt(myGTHDM->getMyGTHDMCache()->mH3sq);
}

massdifference_mH2mmH3::massdifference_mH2mmH3(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double massdifference_mH2mmH3::computeThValue()
{
    return  sqrt(myGTHDM->getMyGTHDMCache()->mH2sq) - sqrt(myGTHDM->getMyGTHDMCache()->mH3sq);
}

cosalpha1_GTHDM::cosalpha1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}
double cosalpha1_GTHDM::computeThValue()
{
    return   myGTHDM->getcosalpha1();
}