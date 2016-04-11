/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SFermiondeltas.h"

deltaRL_12_u::deltaRL_12_u(const StandardModel& SM_i)
: ThObservable(SM_i), mySUSY(static_cast<const SUSY&> (SM_i))
{}

double deltaRL_12_u::computeThValue()
{
    gslpp::matrix<gslpp::complex> TUhat = mySUSY.getTUhat();
    double Yukawa_21=1.;
    double Au_21=TUhat(1,0).real()/Yukawa_21;
    double v2=mySUSY.v2();
    gslpp::matrix<gslpp::complex> msQhat2( mySUSY.getMsQhat2() );
    gslpp::matrix<gslpp::complex> msUhat2( mySUSY.getMsUhat2() );
    double Mu2average_21=sqrt(msQhat2(0,0)*msUhat2(1,1)).abs();

//    gslpp::matrix<gslpp::complex> msLhat2( mySUSY.getMsLhat2() );

//    gslpp::complex sLmass=msLhat2(0,0);

//                gslpp::matrix<gslpp::complex> msLhat2modified(3,3,0.);
//                    msLhat2modified.assign(0, 0, sLmass);
//                    msLhat2modified.assign(1, 1, sLmass);
//                    msLhat2modified.assign(2, 2, sLmass);
//                    msLhat2modified.assign(0, 1, 0.);
//                    msLhat2modified.assign(1, 0, 0.);
//                    msLhat2modified.assign(0, 2, delta13*sLmass);
//                    msLhat2modified.assign(2, 0, delta13*sLmass);
//                    msLhat2modified.assign(1, 2, delta23*sLmass);
//                    msLhat2modified.assign(2, 1, delta23*sLmass);
//                    gslpp::matrix<gslpp::complex> nLL( msLhat2modified + cos2b * Mz2 /2.0 * Id3 );






    return Au_21*v2/Mu2average_21;  //deltaRL_12_u
}
