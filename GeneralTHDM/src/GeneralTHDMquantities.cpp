/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMquantities.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"


//This doesn't make sense anymore. We work with the Higgs basis, anyother basis is just arbitrary
/*
tanbeta_GTHDM::tanbeta_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double tanbeta_GTHDM::computeThValue()
{
    return myGTHDM->gettanb();
}
*/




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

mH3mmH2_GTHDM::mH3mmH2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmH2_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH3sq < 0. || myGTHDM.getMyGTHDMCache()->mH2sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
        return sqrt(myGTHDM.getMyGTHDMCache()->mH3sq) - sqrt(myGTHDM.getMyGTHDMCache()->mH2sq);
}

mH3mmHp_GTHDM::mH3mmHp_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmHp_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH3sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
        return sqrt(myGTHDM.getMyGTHDMCache()->mH3sq) - sqrt(myGTHDM.getMyGTHDMCache()->mHp2);
}

mH3mmH1_GTHDM::mH3mmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH3mmH1_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH3sq < 0. || myGTHDM.getMyGTHDMCache()->mH1sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
        return sqrt(myGTHDM.getMyGTHDMCache()->mH3sq) - sqrt(myGTHDM.getMyGTHDMCache()->mH1sq);
}

mH2mmHp_GTHDM::mH2mmHp_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2mmHp_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH2sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
         return sqrt(myGTHDM.getMyGTHDMCache()->mH2sq) - sqrt(myGTHDM.getMyGTHDMCache()->mHp2);   
}

mH2mmH1_GTHDM::mH2mmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mH2mmH1_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH2sq < 0. || myGTHDM.getMyGTHDMCache()->mH1sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
        return sqrt(myGTHDM.getMyGTHDMCache()->mH2sq) - sqrt(myGTHDM.getMyGTHDMCache()->mH1sq);
}

mHpmmH1_GTHDM::mHpmmH1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double mHpmmH1_GTHDM::computeThValue()
{
    if(myGTHDM.getMyGTHDMCache()->mH1sq < 0.)
        return std::numeric_limits<double>::quiet_NaN();
    else
        return sqrt(myGTHDM.getMyGTHDMCache()->mHp2) - sqrt(myGTHDM.getMyGTHDMCache()->mH1sq);   
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
: ThObservable(SM_i)
//,myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
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

//Since we include this parameter in the map there is no need of defining again this observable
//lambda1::lambda1(const StandardModel& SM_i)
//: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
//{}
//
//double lambda1::computeThValue()
//{
//    return myGTHDM.getMyGTHDMCache()->lambda1;
//}



R11_GTHDM::R11_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R11_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(0,0);
}


R12_GTHDM::R12_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R12_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(0,1);
}


R13_GTHDM::R13_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R13_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(0,2);
}


R21_GTHDM::R21_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R21_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(1,0);
}


R22_GTHDM::R22_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R22_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(1,1);
}


R23_GTHDM::R23_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R23_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(1,2);
}


R31_GTHDM::R31_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R31_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(2,0);
}


R32_GTHDM::R32_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R32_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(2,1);
}


R33_GTHDM::R33_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double R33_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Rij_GTHDM(2,2);
}


//We work with the Higgs basis, there is only one vev which is the SM one, doesn't make sense to keep this observable
/*
v1_GTHDM::v1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double v1_GTHDM::computeThValue()
{
    double v = myGTHDM->v();
    double cosb = myGTHDM->getcosb();

    return (v*cosb);
}
*/


//This doesn't make sense in this model!!!!
/*
v2_GTHDM::v2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}

double v2_GTHDM::computeThValue()
{
    double v = myGTHDM->v();
    double sinb = myGTHDM->getsinb();

    return (v*sinb);
}
*/



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


cosalpha1_GTHDM::cosalpha1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM*> (&SM_i))
{}
double cosalpha1_GTHDM::computeThValue()
{
    return   myGTHDM->getcosalpha1();
}

m1_2::m1_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m1_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m1_2;
}


m2_2::m2_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m2_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m2_2;
}

m3_2::m3_2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double m3_2::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->m3_2;
}




//  Quantities at higher scales



Q_stGTHDM::Q_stGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Q_stGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Q_cutoff;
}


DeltaQ_GTHDM::DeltaQ_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double DeltaQ_GTHDM::computeThValue()
{
    return myGTHDM.getQ_GTHDM() - myGTHDM.getMyGTHDMCache()->Q_cutoff;
}


g1atQGTHDM::g1atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double g1atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->g1_at_Q;
}

g2atQGTHDM::g2atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double g2atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->g2_at_Q;
}

g3atQGTHDM::g3atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double g3atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->g3_at_Q;
}

etaU1atQGTHDM::etaU1atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaU1atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaU1_at_Q;
}

etaU2atQGTHDM::etaU2atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaU2atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaU2_at_Q;
}

etaD1atQGTHDM::etaD1atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaD1atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaD1_at_Q;
}

etaD2atQGTHDM::etaD2atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaD2atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaD2_at_Q;
}

etaL1atQGTHDM::etaL1atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaL1atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaL1_at_Q;
}

etaL2atQGTHDM::etaL2atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double etaL2atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->etaL2_at_Q;
}

lambda1atQGTHDM::lambda1atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda1atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda1_at_Q;
}

lambda2atQGTHDM::lambda2atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda2atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda2_at_Q;
}

lambda3atQGTHDM::lambda3atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda3atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda3_at_Q;
}

lambda4atQGTHDM::lambda4atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double lambda4atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->lambda4_at_Q;
}

Relambda5atQGTHDM::Relambda5atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Relambda5atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Relambda5_at_Q;
}

Relambda6atQGTHDM::Relambda6atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Relambda6atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Relambda6_at_Q;
}

Relambda7atQGTHDM::Relambda7atQGTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Relambda7atQGTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Relambda7_at_Q;
}
