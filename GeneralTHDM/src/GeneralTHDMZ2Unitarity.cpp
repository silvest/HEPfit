/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMZ2Unitarity.h"
#include "GeneralTHDMZ2.h"
#include "GeneralTHDMcache.h"

/************************************/
/* Eigenvalues of the even 00 block */
/************************************/

unitarity00eveP_Z2::unitarity00eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity00eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B1 = -3.*la1Q + 9.*beta1/2. + (i*M_PI - 1.)*(9.*la1Q*la1Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI;

    gslpp::complex B2 = -3.*la2Q + 9.*beta2/2. + (i*M_PI - 1.)*(9.*la2Q*la2Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI;

    gslpp::complex B3 = -2.*la3Q - la4Q + 3.*(2.*beta3 + beta4)/2. + 3.*(i*M_PI -
                        1.)*(la1Q + la2Q)*(2.*la3Q + la4Q)/16./M_PI/M_PI;

    return ((B1 + B2 + sqrt((B1 - B2)*(B1 - B2) + 4.*B3*B3))/32./M_PI - i/2.).abs();
}

unitarity00eveM_Z2::unitarity00eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity00eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B1 = -3.*la1Q + 9.*beta1/2. + (i*M_PI - 1.)*(9.*la1Q*la1Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI;

    gslpp::complex B2 = -3.*la2Q + 9.*beta2/2. + (i*M_PI - 1.)*(9.*la2Q*la2Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI;

    gslpp::complex B3 = -2.*la3Q - la4Q + 3.*(2.*beta3 + beta4)/2. + 3.*(i*M_PI -
                        1.)*(la1Q + la2Q)*(2.*la3Q + la4Q)/16./M_PI/M_PI;

    return ((B1 + B2 - sqrt((B1 - B2)*(B1 - B2) + 4.*B3*B3))/32./M_PI - i/2.).abs();
}


/***********************************/
/* Eigenvalues of the odd 00 block */
/***********************************/

unitarity00oddP_Z2::unitarity00oddP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity00oddP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B4 = -la3Q - 2.*la4Q + 3.*(beta3 + 2.*beta4)/2. + (i*M_PI - 1.)*(la3Q*la3Q +
                        4.*la3Q*la4Q + 4.*la4Q*la4Q + 9.*la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B6 = -3.*la5Q + 9.*beta5/2. + 6.*(i*M_PI - 1.)*(la3Q + 2.*la4Q)*la5Q/16./M_PI/M_PI;

    return ((B4 + B6)/16./M_PI - i/2.).abs();
}

unitarity00oddM_Z2::unitarity00oddM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity00oddM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B4 = -la3Q - 2.*la4Q + 3.*(beta3 + 2.*beta4)/2. + (i*M_PI - 1.)*(la3Q*la3Q +
                        4.*la3Q*la4Q + 4.*la4Q*la4Q + 9.*la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B6 = -3.*la5Q + 9.*beta5/2. + 6.*(i*M_PI - 1.)*(la3Q + 2.*la4Q)*la5Q/16./M_PI/M_PI;

    return ((B4 - B6)/16./M_PI - i/2.).abs();
}


/************************************/
/* Eigenvalues of the even 01 block */
/************************************/

unitarity01eveP_Z2::unitarity01eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity01eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B7 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la4Q*la4Q)/16./M_PI/M_PI;

    gslpp::complex B8 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la4Q*la4Q)/16./M_PI/M_PI;

    gslpp::complex B9 = -la4Q + 3.*beta4/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la4Q/16./M_PI/M_PI;

    return ((B7 + B8 + sqrt((B7 - B8)*(B7 - B8) + 4.*B9*B9))/32./M_PI - i/2.).abs();
}

unitarity01eveM_Z2::unitarity01eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity01eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B7 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la4Q*la4Q)/16./M_PI/M_PI;

    gslpp::complex B8 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la4Q*la4Q)/16./M_PI/M_PI;

    gslpp::complex B9 = -la4Q + 3.*beta4/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la4Q/16./M_PI/M_PI;

    return ((B7 + B8 - sqrt((B7 - B8)*(B7 - B8) + 4.*B9*B9))/32./M_PI - i/2.).abs();
}


/***********************************/
/* Eigenvalues of the odd 01 block */
/***********************************/

unitarity01oddP_Z2::unitarity01oddP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity01oddP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B13 = -la3Q + 3.*beta3/2. + (i*M_PI - 1.)*(la3Q*la3Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B15 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*la3Q*la5Q/16./M_PI/M_PI;

    return ((B13 + B15)/16./M_PI - i/2.).abs();
}

unitarity01oddM_Z2::unitarity01oddM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity01oddM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B13 = -la3Q + 3.*beta3/2. + (i*M_PI - 1.)*(la3Q*la3Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B15 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*la3Q*la5Q/16./M_PI/M_PI;

    return ((B13 - B15)/16./M_PI - i/2.).abs();
}


/**********************************/
/* Eigenvalue of the odd 10 block */
/**********************************/

unitarity10odd_Z2::unitarity10odd_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity10odd_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B19 = -la3Q + la4Q + 3.*(beta3 - beta4)/2. + (i*M_PI - 1.)*(la3Q -
                         la4Q)*(la3Q - la4Q)/16./M_PI/M_PI;

    return (B19/16./M_PI - i/2.).abs();
}


/************************************/
/* Eigenvalues of the even 11 block */
/************************************/

unitarity11eveP_Z2::unitarity11eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity11eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B20 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B21 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B22 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la5Q/16./M_PI/M_PI;

    return ((B20 + B21 + sqrt((B20 - B21)*(B20 - B21) + 4.*B22*B22))/32./M_PI - i/2.).abs();
}

unitarity11eveM_Z2::unitarity11eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity11eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta1 = myGTHDM.getMyGTHDMCache()->betalambda1_Z2(la1Q, la3Q, la4Q, la5Q);
    double beta2 = myGTHDM.getMyGTHDMCache()->betalambda2_Z2(la2Q, la3Q, la4Q, la5Q);
    double beta5 = myGTHDM.getMyGTHDMCache()->betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B20 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B21 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la5Q*la5Q)/16./M_PI/M_PI;

    gslpp::complex B22 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la5Q/16./M_PI/M_PI;

    return ((B20 + B21 - sqrt((B20 - B21)*(B20 - B21) + 4.*B22*B22))/32./M_PI - i/2.).abs();
}


/**********************************/
/* Eigenvalue of the odd 11 block */
/**********************************/

unitarity11odd_Z2::unitarity11odd_Z2(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDMZ2&> (SM_i))
{}

double unitarity11odd_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    double la1Q = myGTHDM.getlambda1_Z2();
    double la2Q = myGTHDM.getlambda2_Z2();
    double la3Q = myGTHDM.getlambda3_Z2();
    double la4Q = myGTHDM.getlambda4_Z2();
    double la5Q = myGTHDM.getlambda5_Z2();

    double beta3 = myGTHDM.getMyGTHDMCache()->betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);
    double beta4 = myGTHDM.getMyGTHDMCache()->betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q);

    gslpp::complex B30 = -la3Q - la4Q + 3.*(beta3 + beta4)/2. + (i*M_PI - 1.)*(la3Q +
                         la4Q)*(la3Q + la4Q)/16./M_PI/M_PI;

    return (B30/16./M_PI - i/2.).abs();
}
