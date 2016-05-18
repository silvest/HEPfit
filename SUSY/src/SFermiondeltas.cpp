/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SFermiondeltas.h"

deltaLL1_q::deltaLL1_q(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL1_q::computeThValue()
{
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    double m1sq=msQhat2(0,0).abs();
    double m2sq=msQhat2(1,1).abs();
    return (m2sq-m1sq)/sqrt(m2sq*m1sq);
}

deltaLL2_q::deltaLL2_q(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL2_q::computeThValue()
{
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    double m1sq=msQhat2(0,0).abs();
    double m3sq=msQhat2(2,2).abs();
    return (m3sq-m1sq)/sqrt(m3sq*m1sq);
}

deltaLL3_q::deltaLL3_q(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL3_q::computeThValue()
{
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    double m2sq=msQhat2(1,1).abs();
    double m3sq=msQhat2(2,2).abs();
    return (m3sq-m2sq)/sqrt(m3sq*m2sq);
}

deltaRR1_u::deltaRR1_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR1_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double m1sq=msUhat2(0,0).abs();
    double m2sq=msUhat2(1,1).abs();
    return (m2sq-m1sq)/sqrt(m2sq*m1sq);
}

deltaRR2_u::deltaRR2_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR2_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double m1sq=msUhat2(0,0).abs();
    double m3sq=msUhat2(2,2).abs();
    return (m3sq-m1sq)/sqrt(m3sq*m1sq);
}

deltaRR3_u::deltaRR3_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR3_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double m2sq=msUhat2(1,1).abs();
    double m3sq=msUhat2(2,2).abs();
    return (m3sq-m2sq)/sqrt(m3sq*m2sq);
}

deltaRR1_d::deltaRR1_d(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR1_d::computeThValue()
{
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    double m1sq=msDhat2(0,0).abs();
    double m2sq=msDhat2(1,1).abs();
    return (m2sq-m1sq)/sqrt(m2sq*m1sq);
}

deltaRR2_d::deltaRR2_d(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR2_d::computeThValue()
{
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    double m1sq=msDhat2(0,0).abs();
    double m3sq=msDhat2(2,2).abs();
    return (m3sq-m1sq)/sqrt(m3sq*m1sq);
}

deltaRR3_d::deltaRR3_d(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR3_d::computeThValue()
{
    gslpp::matrix<gslpp::complex> msDhat2( mySUSY.getMsDhat2() );
    double m2sq=msDhat2(1,1).abs();
    double m3sq=msDhat2(2,2).abs();
    return (m3sq-m2sq)/sqrt(m3sq*m2sq);
}

deltaLL1_l::deltaLL1_l(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL1_l::computeThValue()
{
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    double m1sq=msLhat2(0,0).abs();
    double m2sq=msLhat2(1,1).abs();
    return (m2sq-m1sq)/sqrt(m2sq*m1sq);
}

deltaLL2_l::deltaLL2_l(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL2_l::computeThValue()
{
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    double m1sq=msLhat2(0,0).abs();
    double m3sq=msLhat2(2,2).abs();
    return (m3sq-m1sq)/sqrt(m3sq*m1sq);
}

deltaLL3_l::deltaLL3_l(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLL3_l::computeThValue()
{
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    double m2sq=msLhat2(1,1).abs();
    double m3sq=msLhat2(2,2).abs();
    return (m3sq-m2sq)/sqrt(m3sq*m2sq);
}

deltaRR1_e::deltaRR1_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR1_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double m1sq=msEhat2(0,0).abs();
    double m2sq=msEhat2(1,1).abs();
    return (m2sq-m1sq)/sqrt(m2sq*m1sq);
}

deltaRR2_e::deltaRR2_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR2_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double m1sq=msEhat2(0,0).abs();
    double m3sq=msEhat2(2,2).abs();
    return (m3sq-m1sq)/sqrt(m3sq*m1sq);
}

deltaRR3_e::deltaRR3_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRR3_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double m2sq=msEhat2(1,1).abs();
    double m3sq=msEhat2(2,2).abs();
    return (m3sq-m2sq)/sqrt(m3sq*m2sq);
}

deltaRL_12_u::deltaRL_12_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_12_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> TUhat = mySUSY.getTUhat();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double Mu2average_21=sqrt(msQhat2(1,1)*msUhat2(0,0)).abs();
    double v2 = mySUSY.v2();

    return TUhat(1,0).real()*v2/Mu2average_21;  //deltaRL_12_u
}

deltaRL_13_u::deltaRL_13_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_13_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> TUhat = mySUSY.getTUhat();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double Mu2average_31=sqrt(msQhat2(2,2)*msUhat2(0,0)).abs();
    double v2 = mySUSY.v2();

    return TUhat(2,0).real()*v2/Mu2average_31;  //deltaRL_13_u
}

deltaRL_23_u::deltaRL_23_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_23_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> TUhat = mySUSY.getTUhat();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double Mu2average_32=sqrt(msQhat2(2,2)*msUhat2(1,1)).abs();
    double v2 = mySUSY.v2();

    return TUhat(2,1).real()*v2/Mu2average_32;  //deltaRL_23_u
}

deltaRL_12_e::deltaRL_12_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_12_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_21=sqrt(msLhat2(1,1)*msEhat2(0,0)).abs();
    double v1 = mySUSY.v1();

    return TEhat(1,0).real()*v1/Me2average_21;  //deltaRL_12_e
}

deltaRL_21_e::deltaRL_21_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_21_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_12=sqrt(msLhat2(0,0)*msEhat2(1,1)).abs();
    double v1 = mySUSY.v1();

    return TEhat(0,1).real()*v1/Me2average_12;  //deltaRL_12_e
}

deltaRL_13_e::deltaRL_13_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_13_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_31=sqrt(msLhat2(2,2)*msEhat2(0,0)).abs();
    double v1 = mySUSY.v1();

    return TEhat(2,0).real()*v1/Me2average_31;  //deltaRL_13_e
}

deltaRL_31_e::deltaRL_31_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_31_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_13=sqrt(msLhat2(0,0)*msEhat2(2,2)).abs();
    double v1 = mySUSY.v1();

    return TEhat(0,2).real()*v1/Me2average_13;  //deltaRL_13_e
}

deltaRL_23_e::deltaRL_23_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_23_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_32=sqrt(msLhat2(2,2)*msEhat2(1,1)).abs();
    double v1 = mySUSY.v1();

    return TEhat(2,1).real()*v1/Me2average_32;  //deltaRL_23_e
}

deltaRL_32_e::deltaRL_32_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_32_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_23=sqrt(msLhat2(1,1)*msEhat2(2,2)).abs();
    double v1 = mySUSY.v1();

    return TEhat(1,2).real()*v1/Me2average_23;  //deltaRL_32_e
}

logdeltaRL_13_e::logdeltaRL_13_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double logdeltaRL_13_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_31=sqrt(msLhat2(2,2)*msEhat2(0,0)).abs();
    double v1 = mySUSY.v1();

    return log10(TEhat(2,0).abs()*v1/Me2average_31);  //logdeltaRL_13_e
}

logdeltaRL_23_e::logdeltaRL_23_e(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double logdeltaRL_23_e::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_32=sqrt(msLhat2(2,2)*msEhat2(1,1)).abs();
    double v1 = mySUSY.v1();

    return log10(TEhat(2,1).abs()*v1/Me2average_32);  //logdeltaRL_23_e
}

logmslepton::logmslepton(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double logmslepton::computeThValue()
{
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double Me2average_32=sqrt(msLhat2(2,2)*msEhat2(1,1)).abs();

    return log10(sqrt(Me2average_32));  //logmslepton
}

deltaTEhat23::deltaTEhat23(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaTEhat23::computeThValue()
{
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    double denominator = sqrt(TEhat(2,1).abs()*TEhat(1,2).abs());
    
    return (TEhat(2,1).abs()-TEhat(1,2).abs())/denominator;
}

deltaLLRR_l::deltaLLRR_l(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaLLRR_l::computeThValue()
{
    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );
    gslpp::matrix<gslpp::complex> msEhat2( mySUSY.getMsEhat2() );
    double denominator = sqrt(msLhat2(2,2).abs()*msEhat2(2,2).abs());
    
    return (msLhat2(2,2).abs()-msEhat2(2,2).abs())/denominator;
}
