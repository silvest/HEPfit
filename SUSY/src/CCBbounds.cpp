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

    return lambdau11*lambdau11*(msQhat2(0,0).abs()+msUhat2(0,0).abs()+mHusq)-TUhat(0,0).real()/lambdau11;
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

    return lambdad11*lambdad11*(msQhat2(0,0).abs()+msDhat2(0,0).abs()+mHdsq)-TDhat(0,0).real()/lambdad11;
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

    return lambdal11*lambdal11*(msLhat2(0,0).abs()+msEhat2(0,0).abs()+mHdsq)-TEhat(0,0).real()/lambdal11;
}
