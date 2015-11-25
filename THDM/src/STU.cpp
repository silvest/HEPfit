/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "STU.h"
#include "StandardModel.h"

STU::STU(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))

{
    mycache = new THDMcache();
};

double STU::computeThValue()
{
////    double vev=myTHDM->v();
////    double mHl=myTHDM->getMHl();
////    double mHh=myTHDM->getMHh();
////    double mA=myTHDM->getMA();
////    double mHp=myTHDM->getMHp();
////    double sina=myTHDM->computeSina();
////    double cosa=myTHDM->computeCosa();
////    double tanb=myTHDM->getTanb();
////    double sinb=myTHDM->getSinb();
////    double cosb=myTHDM->getCosb();
////    double sin_ba=myTHDM->getsin_ba();
////    double m12_2=myTHDM->getM12_2(); 
//
////    double MZ=myTHDM->getMz();
//
//    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHh2;
//    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHl2;
//    gslpp::complex B00prime_MZ2_MZ2_mHh2_mA2;
//    gslpp::complex B00prime_MZ2_MZ2_mHp2_mHp2;
//    gslpp::complex B00prime_MZ2_MZ2_mHl2_mA2;
//    gslpp::complex B00prime_MZ2_MZ2_MZ2_mHh2;
//    gslpp::complex B00prime_MZ2_MZ2_MZ2_mHl2;
//
////    double MZ2 = MZ*MZ;
////    double sin2_ba = sin_ba*sin_ba;
////    double cos2_ba = 1. - sin2_ba;
////
////    double MW = myTHDM->Mw();
////    double MW2 = MW*MW;
////    double s_W2 = myTHDM->sW2();
//
//    B00prime_MZ2_MZ2_mHh2_mA2 = - mycache->B00_MZ2_MZ2_mHh2_mA2(MZ,mHh,mA) + mycache->B00_MZ2_0_mHh2_mA2(MZ,mHh,mA);
//    B00prime_MZ2_MZ2_mHp2_mHp2 = - mycache->B00_MZ2_MZ2_mHp2_mHp2(MZ,mHp) + mycache->B00_MZ2_0_mHp2_mHp2(MZ,mHp);
//    B00prime_MZ2_MZ2_mHl2_mA2 = - mycache->B00_MZ2_MZ2_mHl2_mA2(MZ,mHl,mA) + mycache->B00_MZ2_0_mHl2_mA2(MZ,mHl,mA);
//    B00prime_MZ2_MZ2_MZ2_mHh2 = - mycache->B00_MZ2_MZ2_MZ2_mHh2(MZ,mHh) + mycache->B00_MZ2_0_MZ2_mHh2(MZ,mHh);
//    B00prime_MZ2_MZ2_MZ2_mHl2 = - mycache->B00_MZ2_0_MZ2_mHl2(MZ,mHl) + mycache->B00_MZ2_0_MZ2_mHl2(MZ,mHl);
//    B0prime_MZ2_MZ2_MZ2_mHh2 = mycache->B0_MZ2_MZ2_MZ2_mHh2(MZ,mHh) - mycache->B0_MZ2_0_MZ2_mHh2(MZ,mHh);
//    B0prime_MZ2_MZ2_MZ2_mHl2 = mycache->B0_MZ2_MZ2_MZ2_mHl2(MZ,mHl) - mycache->B0_MZ2_0_MZ2_mHl2(MZ,mHl);
//    
//    double DeltaS = 1./MZ2/M_PI*(sin2_ba * B00prime_MZ2_MZ2_mHh2_mA2.real() - B00prime_MZ2_MZ2_mHp2_mHp2.real()
//           + cos2_ba * (B00prime_MZ2_MZ2_mHl2_mA2.real() + B00prime_MZ2_MZ2_MZ2_mHh2.real()
//           - B00prime_MZ2_MZ2_MZ2_mHl2.real() - MZ2 * B0prime_MZ2_MZ2_MZ2_mHh2.real()
//           + MZ2 * B0prime_MZ2_MZ2_MZ2_mHl2.real()));
//
//    gslpp::complex B0_MZ2_0_MZ2_mHh2;
//    gslpp::complex B0_MZ2_0_MZ2_mHl2;
//    gslpp::complex B0_MZ2_0_MW2_mHh2;
//    gslpp::complex B0_MZ2_0_MW2_mHl2;    
//
//    B0_MZ2_0_MW2_mHh2 = mycache->B0_MZ2_0_MW2_mHh2(MZ,MW,mHh);
//    B0_MZ2_0_MZ2_mHh2 = mycache->B0_MZ2_0_MZ2_mHh2(MZ,mHh);
//    B0_MZ2_0_MZ2_mHl2 = mycache->B0_MZ2_0_MW2_mHl2(MZ,MW,mHl);
//    B0_MZ2_0_MW2_mHl2 = mycache->B0_MZ2_0_MW2_mHl2(MZ,MW,mHl);
//    
//    double DeltaT = 1. / 16. / M_PI / MW2 / s_W2 * (F(mHp,mA)
//           + sin2_ba * (F(mHp,mHh) - F(mA,mHh)) + cos2_ba * (F(mHp,mHl)
//           - F(mA,mHl) + F(MW,mHh) - F(MW,mHl) - F(MZ,mHh)
//           + F(MZ,mHl) + 4. * MZ2 * (B0_MZ2_0_MZ2_mHh2.real() - B0_MZ2_0_MZ2_mHl2.real())
//           - 4. * MW2 * (B0_MZ2_0_MW2_mHh2.real() - B0_MZ2_0_MW2_mHl2.real())));
//
//    gslpp::complex B0prime_MZ2_MW2_MW2_mHh2;
//    gslpp::complex B0prime_MZ2_MW2_MW2_mHl2;
//    gslpp::complex B00prime_MZ2_MW2_mA2_mHp2;
//    gslpp::complex B00prime_MZ2_MW2_mHp2_mHp2;
//    gslpp::complex B00prime_MZ2_MW2_mHh2_mHp2;
//    gslpp::complex B00prime_MZ2_MW2_mHl2_mHp2;
//    gslpp::complex B00prime_MZ2_MW2_MW2_mHh2;
//    gslpp::complex B00prime_MZ2_MW2_MW2_mHl2;
//      
//    B00prime_MZ2_MW2_mA2_mHp2 = - mycache->B00_MZ2_MW2_mA2_mHp2(MZ,MW,mA,mHp) + mycache->B00_MZ2_0_mA2_mHp2(MZ,mA,mHp);
//    B00prime_MZ2_MW2_mHp2_mHp2 = - mycache->B00_MZ2_MW2_mHp2_mHp2(MZ,MW,mHp) + mycache->B00_MZ2_0_mHp2_mHp2(MZ,mHp);
//    B00prime_MZ2_MW2_mHl2_mHp2 = - mycache->B00_MZ2_MW2_mHl2_mHp2(MZ,MW,mHl,mHp) + mycache->B00_MZ2_0_mHl2_mHp2(MZ,mHl,mHp);
//    B00prime_MZ2_MW2_MW2_mHh2 = - mycache->B00_MZ2_MW2_MW2_mHh2(MZ,MW,mHh) + mycache->B00_MZ2_0_MW2_mHh2(MZ,MW,mHh);
//    B00prime_MZ2_MW2_MW2_mHl2 = - mycache->B00_MZ2_MW2_MW2_mHl2(MZ,MW,mHl) + mycache->B00_MZ2_0_MW2_mHl2(MZ,MW,mHl);
//    B0prime_MZ2_MW2_MW2_mHh2 = mycache->B0_MZ2_MW2_MW2_mHh2(MZ,MW,mHh) - mycache->B0_MZ2_0_MW2_mHh2(MZ,MW,mHh);
//    B0prime_MZ2_MW2_MW2_mHl2 = mycache->B0_MZ2_MW2_MW2_mHl2(MZ,MW,mHl) - mycache->B0_MZ2_0_MW2_mHl2(MZ,MW,mHl);
//    B00prime_MZ2_MW2_mHh2_mHp2 = - mycache->B00_MZ2_MW2_mHh2_mHp2(MZ,MW,mHh,mHp) + mycache->B00_MZ2_0_mHh2_mHp2(MZ,mHh,mHp);
//    
//    double DeltaU = - DeltaS + 1. / M_PI / MZ2 * (B00prime_MZ2_MW2_mA2_mHp2.real()
//           - 2. * B00prime_MZ2_MW2_mHp2_mHp2.real() + sin2_ba * B00prime_MZ2_MW2_mHh2_mHp2.real()
//           + cos2_ba * (B00prime_MZ2_MW2_mHl2_mHp2.real() + B00prime_MZ2_MW2_MW2_mHh2.real()
//           - B00prime_MZ2_MW2_MW2_mHl2.real() - MW2 * B0prime_MZ2_MW2_MW2_mHh2.real()
//           + MW2 * B0prime_MZ2_MW2_MW2_mHl2.real()));

    return 0.0;

}

/////////////////////////////////////////////////////////////////////////////

double STU::F(const double m02, const double m12) const {
    double F;

    if(m02 == 0. && m12 != 0.) {
        F=0.5 * m12;
    } else if(m02 != 0. && m12 == 0.){
        F=0.5 * m02;
    } else if((m02 == 0. && m12 == 0.) || (fabs(m02-m12) < LEPS)){
        F=0.;
    } else if (m02 != 0 && m12 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    } else
        throw std::runtime_error("Error in THDM::F()");
    return (F);
}

/////////////////////////////////////////////////////////////////////////////

DeltaS::DeltaS(const StandardModel& SM_i)
: STU(SM_i)
{}

double DeltaS::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mHl2=mHl*mHl;
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double sin_ba=myTHDM->getsin_ba();
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
    double MZ=myTHDM->getMz();
    double MZ2 = MZ*MZ;

    gslpp::complex B00prime_MZ2_MZ2_mHh2_mA2;
    gslpp::complex B00prime_MZ2_MZ2_mHp2_mHp2;
    gslpp::complex B00prime_MZ2_MZ2_mHl2_mA2;
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mHh2;
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mHl2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHh2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHl2;

    B00prime_MZ2_MZ2_mHh2_mA2 = - mycache->B00_MZ2_MZ2_mHh2_mA2(MZ2,mHh2,mA2) + mycache->B00_MZ2_0_mHh2_mA2(MZ2,mHh2,mA2);
    B00prime_MZ2_MZ2_mHp2_mHp2 = - mycache->B00_MZ2_MZ2_mHp2_mHp2(MZ2,mHp2) + mycache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);
    B00prime_MZ2_MZ2_mHl2_mA2 = - mycache->B00_MZ2_MZ2_mHl2_mA2(MZ2,mHl2,mA2) + mycache->B00_MZ2_0_mHl2_mA2(MZ2,mHl2,mA2);
    B00prime_MZ2_MZ2_MZ2_mHh2 = - mycache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mHh2) + mycache->B00_MZ2_0_MZ2_mHh2(MZ2,mHh2);
    B00prime_MZ2_MZ2_MZ2_mHl2 = - mycache->B00_MZ2_MZ2_MZ2_mHl2(MZ2,mHl2) + mycache->B00_MZ2_0_MZ2_mHl2(MZ2,mHl2);
    B0prime_MZ2_MZ2_MZ2_mHh2 = mycache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mHh2) - mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mHh2);
    B0prime_MZ2_MZ2_MZ2_mHl2 = mycache->B0_MZ2_MZ2_MZ2_mHl2(MZ2,mHl2) - mycache->B0_MZ2_0_MZ2_mHl2(MZ2,mHl2);

    return 1./MZ2/M_PI*(sin2_ba * B00prime_MZ2_MZ2_mHh2_mA2.real() - B00prime_MZ2_MZ2_mHp2_mHp2.real()
           + cos2_ba * (B00prime_MZ2_MZ2_mHl2_mA2.real() + B00prime_MZ2_MZ2_MZ2_mHh2.real()
           - B00prime_MZ2_MZ2_MZ2_mHl2.real() - MZ2 * B0prime_MZ2_MZ2_MZ2_mHh2.real()
           + MZ2 * B0prime_MZ2_MZ2_MZ2_mHl2.real()));
}

DeltaT::DeltaT(const StandardModel& SM_i)
: STU(SM_i)
{}

double DeltaT::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mHl2=mHl*mHl;
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double sin_ba=myTHDM->getsin_ba();
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
    double MZ=myTHDM->getMz();
    double MZ2 = MZ*MZ;
    double MW = myTHDM->Mw();
    double MW2 = MW*MW;
    double s_W2 = myTHDM->sW2();

    gslpp::complex B0_MZ2_0_MZ2_mHh2;
    gslpp::complex B0_MZ2_0_MZ2_mHl2;
    gslpp::complex B0_MZ2_0_MW2_mHh2;
    gslpp::complex B0_MZ2_0_MW2_mHl2;    

    B0_MZ2_0_MZ2_mHh2 = mycache->B0_MZ2_0_MZ2_mHh2(MZ2,mHh2);
    B0_MZ2_0_MZ2_mHl2 = mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHl2);
    B0_MZ2_0_MW2_mHh2 = mycache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mHh2);
    B0_MZ2_0_MW2_mHl2 = mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHl2);

    return 1. / 16. / M_PI / MW2 / s_W2 * (F(mHp2,mA2)
           + sin2_ba * (F(mHp2,mHh2) - F(mA2,mHh2)) + cos2_ba * (F(mHp2,mHl2)
           - F(mA2,mHl2) + F(MW2,mHh2) - F(MW2,mHl2) - F(MZ2,mHh2)
           + F(MZ2,mHl2) + 4. * MZ2 * (B0_MZ2_0_MZ2_mHh2.real() - B0_MZ2_0_MZ2_mHl2.real())
           - 4. * MW2 * (B0_MZ2_0_MW2_mHh2.real() - B0_MZ2_0_MW2_mHl2.real())));
}

DeltaU::DeltaU(const StandardModel& SM_i)
: STU(SM_i)
{
    myDeltaS = new DeltaS(SM_i);
}

double DeltaU::computeThValue()
{
    double mHl=myTHDM->getMHl();
    double mHl2=mHl*mHl;
    double mHh2=myTHDM->getmHh2();
    double mA2=myTHDM->getmA2();
    double mHp2=myTHDM->getmHp2();
    double sin_ba=myTHDM->getsin_ba();
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
    double MZ=myTHDM->getMz();
    double MZ2 = MZ*MZ;
    double MW = myTHDM->Mw();
    double MW2 = MW*MW;

    gslpp::complex B00prime_MZ2_MW2_mA2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mHp2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mHh2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mHl2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mHh2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mHl2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mHh2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mHl2;

    B00prime_MZ2_MW2_mA2_mHp2 = - mycache->B00_MZ2_MW2_mA2_mHp2(MZ2,MW2,mA2,mHp2) + mycache->B00_MZ2_0_mA2_mHp2(MZ2,mA2,mHp2);
    B00prime_MZ2_MW2_mHp2_mHp2 = - mycache->B00_MZ2_MW2_mHp2_mHp2(MZ2,MW2,mHp2) + mycache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);
    B00prime_MZ2_MW2_mHh2_mHp2 = - mycache->B00_MZ2_MW2_mHh2_mHp2(MZ2,MW2,mHh2,mHp2) + mycache->B00_MZ2_0_mHh2_mHp2(MZ2,mHh2,mHp2);
    B00prime_MZ2_MW2_mHl2_mHp2 = - mycache->B00_MZ2_MW2_mHl2_mHp2(MZ2,MW2,mHl2,mHp2) + mycache->B00_MZ2_0_mHl2_mHp2(MZ2,mHl2,mHp2);
    B00prime_MZ2_MW2_MW2_mHh2 = - mycache->B00_MZ2_MW2_MW2_mHh2(MZ2,MW2,mHh2) + mycache->B00_MZ2_0_MW2_mHh2(MZ2,MW2,mHh2);
    B00prime_MZ2_MW2_MW2_mHl2 = - mycache->B00_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHl2) + mycache->B00_MZ2_0_MW2_mHl2(MZ2,MW2,mHl2);
    B0prime_MZ2_MW2_MW2_mHh2 = mycache->B0_MZ2_MW2_MW2_mHh2(MZ2,MW2,mHh2) - mycache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mHh2);
    B0prime_MZ2_MW2_MW2_mHl2 = mycache->B0_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHl2) - mycache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHl2);

    return - myDeltaS->computeThValue() + 1. / M_PI / MZ2 * (B00prime_MZ2_MW2_mA2_mHp2.real()
           - 2. * B00prime_MZ2_MW2_mHp2_mHp2.real() + sin2_ba * B00prime_MZ2_MW2_mHh2_mHp2.real()
           + cos2_ba * (B00prime_MZ2_MW2_mHl2_mHp2.real() + B00prime_MZ2_MW2_MW2_mHh2.real()
           - B00prime_MZ2_MW2_MW2_mHl2.real() - MW2 * B0prime_MZ2_MW2_MW2_mHh2.real()
           + MW2 * B0prime_MZ2_MW2_MW2_mHl2.real()));
}
