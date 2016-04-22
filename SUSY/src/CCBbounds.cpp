/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CCBbounds.h"

CCBu11::CCBu11(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu11::computeThValue()
{
    //Is the SCKM trilinear coupling the same as the Aij?
    //Other question in the SUSYMatching...
    //Are the mass parameters assigned correctly?
    //

    double mUP = mySUSY.getQuarks(QCD::UP).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau11 = mUP / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return (lambdau11*lambdau11*(msQhat2(0,0).abs()+msUhat2(0,0).abs()+mHusq)-TUhat(0,0).abs2()) / fabs(msQhat2(0,0).abs()+msUhat2(0,0).abs()+mHusq);
}

CCBu22::CCBu22(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu22::computeThValue()
{
    double mCHARM = mySUSY.getQuarks(QCD::CHARM).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau22 = mCHARM / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdau22*lambdau22-TUhat(1,1).abs2()/(msQhat2(1,1).abs()+msUhat2(1,1).abs()+mHusq);
}

CCBu33::CCBu33(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu33::computeThValue()
{
    double mTOP = mySUSY.getQuarks(QCD::TOP).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau33 = mTOP / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdau33*lambdau33-TUhat(2,2).abs2()/(msQhat2(2,2).abs()+msUhat2(2,2).abs()+mHusq);
}

CCBu12::CCBu12(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu12::computeThValue()
{
    double mCHARM = mySUSY.getQuarks(QCD::CHARM).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau22 = mCHARM / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdau22*lambdau22-TUhat(0,1).abs2()/(msQhat2(0,0).abs()+msUhat2(1,1).abs()+mHusq);
}

CCBu13::CCBu13(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu13::computeThValue()
{
    double mTOP = mySUSY.getQuarks(QCD::TOP).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau33 = mTOP / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdau33*lambdau33-TUhat(0,2).abs2()/(msQhat2(0,0).abs()+msUhat2(2,2).abs()+mHusq);
}

CCBu23::CCBu23(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBu23::computeThValue()
{
    double mTOP = mySUSY.getQuarks(QCD::TOP).getMass();
    double v2 = mySUSY.v2();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    gslpp::matrix<gslpp::complex> TUhat( mySUSY.getTUhat() );

    double lambdau33 = mTOP / v2 * sqrt(2.);
    double mHusq = mA*mA * cosb*cosb + 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdau33*lambdau33-TUhat(1,2).abs2()/(msQhat2(1,1).abs()+msUhat2(2,2).abs()+mHusq);
}

CCBd11::CCBd11(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd11::computeThValue()
{
    double mDOWN = mySUSY.getQuarks(QCD::DOWN).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad11 = mDOWN / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad11*lambdad11-TDhat(0,0).abs2()/(msQhat2(0,0).abs()+msDhat2(0,0).abs()+mHdsq);
}

CCBd22::CCBd22(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd22::computeThValue()
{
    double mSTRANGE = mySUSY.getQuarks(QCD::STRANGE).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad22 = mSTRANGE / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad22*lambdad22-TDhat(1,1).abs2()/(msQhat2(1,1).abs()+msDhat2(1,1).abs()+mHdsq);
}

CCBd33::CCBd33(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd33::computeThValue()
{
    double mBOTTOM = mySUSY.getQuarks(QCD::BOTTOM).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad33 = mBOTTOM / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad33*lambdad33-TDhat(2,2).abs2()/(msQhat2(2,2).abs()+msDhat2(2,2).abs()+mHdsq);
}

CCBd12::CCBd12(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd12::computeThValue()
{
    double mSTRANGE = mySUSY.getQuarks(QCD::STRANGE).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad22 = mSTRANGE / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad22*lambdad22-TDhat(0,1).abs2()/(msQhat2(0,0).abs()+msDhat2(1,1).abs()+mHdsq);
}

CCBd13::CCBd13(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd13::computeThValue()
{
    double mBOTTOM = mySUSY.getQuarks(QCD::BOTTOM).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad33 = mBOTTOM / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad33*lambdad33-TDhat(0,2).abs2()/(msQhat2(0,0).abs()+msDhat2(2,2).abs()+mHdsq);
}

CCBd23::CCBd23(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBd23::computeThValue()
{
    double mBOTTOM = mySUSY.getQuarks(QCD::BOTTOM).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    gslpp::matrix<gslpp::complex> TDhat( mySUSY.getTDhat() );

    double lambdad33 = mBOTTOM / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdad33*lambdad33-TDhat(1,2).abs2()/(msQhat2(1,1).abs()+msDhat2(2,2).abs()+mHdsq);
}

CCBe11::CCBe11(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe11::computeThValue()
{
    double mE = mySUSY.getLeptons(StandardModel::ELECTRON).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal11 = mE / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal11*lambdal11-TEhat(0,0).abs2()/(msLhat2(0,0).abs()+msEhat2(0,0).abs()+mHdsq);
}

CCBe22::CCBe22(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe22::computeThValue()
{
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal22 = mMU / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal22*lambdal22-TEhat(1,1).abs2()/(msLhat2(1,1).abs()+msEhat2(1,1).abs()+mHdsq);
}

CCBe33::CCBe33(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe33::computeThValue()
{
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal33 = mTAU / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal33*lambdal33-TEhat(2,2).abs2()/(msLhat2(2,2).abs()+msEhat2(2,2).abs()+mHdsq);
}

CCBe12::CCBe12(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe12::computeThValue()
{
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal22 = mMU / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal22*lambdal22-TEhat(0,1).abs2()/(msLhat2(0,0).abs()+msEhat2(1,1).abs()+mHdsq);
}

CCBe13::CCBe13(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe13::computeThValue()
{
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal33 = mTAU / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal33*lambdal33-TEhat(0,2).abs2()/(msLhat2(0,0).abs()+msEhat2(2,2).abs()+mHdsq);
}

CCBe23::CCBe23(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double CCBe23::computeThValue()
{
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    double v1 = mySUSY.v1();
    double mA = mySUSY.getMHptree();
    double cosb = mySUSY.getCosb();
    double MZ = mySUSY.getMz();
    double sinb = mySUSY.getSinb();
    double cos2b = cosb*cosb-sinb*sinb;
    gslpp::complex MuH = mySUSY.getMuH();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    gslpp::matrix<gslpp::complex> TEhat( mySUSY.getTEhat() );

    double lambdal33 = mTAU / v1 * sqrt(2.);
    double mHdsq = mA*mA * sinb*sinb - 0.5 * MZ*MZ * cos2b - MuH.abs2();

    return lambdal33*lambdal33-TEhat(1,2).abs2()/(msLhat2(1,1).abs()+msEhat2(2,2).abs()+mHdsq);
}
