/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMquantities.h"

mH1_GTHDM::mH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH1_GTHDM::computeThValue()
{
    return sqrt(myGTHDM.getMyGTHDMCache()->mH1_2);
}

mH2_GTHDM::mH2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2_GTHDM::computeThValue()
{
    return sqrt(myGTHDM.getMyGTHDMCache()->mH2_2);
}

mH3_GTHDM::mH3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3_GTHDM::computeThValue()
{
    return sqrt(myGTHDM.getMyGTHDMCache()->mH3_2);
}

mH1sq_GTHDM::mH1sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH1sq_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH1_2;
}

mH2sq_GTHDM::mH2sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2sq_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH2_2;
}

mH3sq_GTHDM::mH3sq_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3sq_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->mH3_2;
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
    return myGTHDM.getMyGTHDMCache()->M2_GTHDM;
}

m11_2_GTHDM::m11_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m11_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m11_2_GTHDM;
}

m22_2_GTHDM::m22_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m22_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m22_2_GTHDM;
}

Imm12_2_GTHDM::Imm12_2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Imm12_2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Imm12_2_GTHDM;
}

lambda1_GTHDM::lambda1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda1_GTHDM;
}

lambda2_GTHDM::lambda2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda2_GTHDM;
}

lambda3_GTHDM::lambda3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda3_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda3_GTHDM;
}

lambda4_GTHDM::lambda4_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda4_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda4_GTHDM;
}

Relambda5_GTHDM::Relambda5_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Relambda5_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Relambda5_GTHDM;
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


//GTHDM_Relambda5::GTHDM_Relambda5(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_Relambda5e::computeThValue()
//{   
//    double M_2 = myGeneralTHDM->getM_2();
//    double M33_2 = myGeneralTHDM->getM33_2();
//    double v = myGeneralTHDM->v();
//    double tanb = myGeneralTHDM->gettanb();
//    double Relambda6 = myGeneralTHDM->getRelambda6();
//    double Relambda7 = myGeneralTHDM->getRelambda7();
//    
//    return ((M_2 - M33_2)/v/v - lambda6_Re/tanb/2. - tanb*lambda7_Re/2.);
//}

//
//GTHDM_lambda1::GTHDM_lambda1(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_lambda1::computeThValue()
//{   
//    double M_2 = myGeneralTHDM->getM_2();
//    double M11_2 = myGeneralTHDM->getM11_2();
//    double M12_2 = myGeneralTHDM->getM12_2();
//    double M22_2 = myGeneralTHDM->getM22_2();
//    double v = myGeneralTHDM->v();
//    double tanb = myGeneralTHDM->gettanb();
//    double lambda6_Re = myGeneralTHDM->getlambda6_Re();
//    double lambda7_Re = myGeneralTHDM->getlambda7_Re();
//    
//    return ((M11_2 + tanb*tanb*(M22_2-M_2) - 2.*tanb*M12_2)/v/v + tanb*(tanb*tanb*lambda7_Re - 3.*lambda6_Re)/2.);
//}
//
//
//GTHDM_lambda2::GTHDM_lambda2(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_lambda2::computeThValue()
//{   
//    double M_2 = myGeneralTHDM->getM_2();
//    double M11_2 = myGeneralTHDM->getM11_2();
//    double M12_2 = myGeneralTHDM->getM12_2();
//    double M22_2 = myGeneralTHDM->getM22_2();
//    double v = myGeneralTHDM->v();
//    double tanb = myGeneralTHDM->gettanb();
//    double lambda6_Re = myGeneralTHDM->getlambda6_Re();
//    double lambda7_Re = myGeneralTHDM->getlambda7_Re();
//    
//    return ((M11_2 + (M22_2-M_2)/tanb/tanb + 2.*M12_2/tanb)/v/v + (lambda6_Re/tanb/tanb - 3.*lambda7_Re)/tanb/2.);
//}
//
//
//GTHDM_lambda3::GTHDM_lambda3(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_lambda3::computeThValue()
//{
//    double mHp_2 = myGeneralTHDM->getmHp_2();
//    double M_2 = myGeneralTHDM->getM_2();
//    double M11_2 = myGeneralTHDM->getM11_2();
//    double M12_2 = myGeneralTHDM->getM12_2();
//    double M22_2 = myGeneralTHDM->getM22_2();
//    double v = myGeneralTHDM->v();
//    double tanb = myGeneralTHDM->gettanb();
//    double lambda6_Re = myGeneralTHDM->getlambda6_Re();
//    double lambda7_Re = myGeneralTHDM->getlambda7_Re();
//    
//    return ((M11_2 - M22_2 - M_2 + (1./tanb - tanb)*M12_2 + 2.*mHp_2)/v/v - (lambda6_Re/tanb + tanb*lambda7_Re)/2.);
//}
//
//
//GTHDM_lambda4::GTHDM_lambda4(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_lambda4::computeThValue()
//{
//    double mHp_2 = myGeneralTHDM->getmHp_2();
//    double M_2 = myGeneralTHDM->getM_2();
//    double M33_2 = myGeneralTHDM->getM33_2();
//    double v = myGeneralTHDM->v();
//    double tanb = myGeneralTHDM->gettanb();
//    double lambda6_Re = myGeneralTHDM->getlambda6_Re();
//    double lambda7_Re = myGeneralTHDM->getlambda7_Re();
//    
//    return ((M_2 + M33_2 - 2.*mHp_2)/v/v - (lambda6_Re/tanb + tanb*lambda7_Re)/2.);
//}
//
//
//GTHDM_M_2::GTHDM_M_2(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_M_2::computeThValue()
//{    
//    return (myGeneralTHDM->getM_2());
//}
//
//
//GTHDM_mH2_2::GTHDM_mH2_2(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_mH2_2::computeThValue()
//{    
//    return (myGeneralTHDM->getmH2_2());
//}
//
//
//GTHDM_mH3_2::GTHDM_mH3_2(const StandardModel& SM_i)
//: ThObservable(SM_i), myGeneralTHDM(static_cast<const GeneralTHDM*> (&SM_i))
//{}
//
//double GTHDM_mH3_2::computeThValue()
//{    
//    return (myGeneralTHDM->getmH3_2());
//}