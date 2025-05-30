/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMZ2Unitarity.h"
#include "GeneralTHDMZ2Runner.h"

/************************************/
/* Eigenvalues of the even 00 block */
/************************************/

unitarity00eveP_Z2::unitarity00eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity00eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B1 = -3.*la1Q + 9.*beta1/2. + (i*M_PI - 1.)*(9.*la1Q*la1Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI -
                         3.*la1Q*WFRc1/32./M_PI/M_PI;

    gslpp::complex B2 = -3.*la2Q + 9.*beta2/2. + (i*M_PI - 1.)*(9.*la2Q*la2Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI -
                         3.*la2Q*(-WFRc1 + 2.*WFRc2)/32./M_PI/M_PI;

    gslpp::complex B3 = -2.*la3Q - la4Q + 3.*(2.*beta3 + beta4)/2. + 3.*(i*M_PI -
                         1.)*(la1Q + la2Q)*(2.*la3Q + la4Q)/16./M_PI/M_PI -
                        (2.*la3Q + la4Q)*WFRc2/32./M_PI/M_PI;

    return ((B1 + B2 + sqrt((B1 - B2)*(B1 - B2) + 4.*B3*B3))/32./M_PI - i/2.).abs();
}

unitarity00eveM_Z2::unitarity00eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity00eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B1 = -3.*la1Q + 9.*beta1/2. + (i*M_PI - 1.)*(9.*la1Q*la1Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI -
                         3.*la1Q*WFRc1/32./M_PI/M_PI;

    gslpp::complex B2 = -3.*la2Q + 9.*beta2/2. + (i*M_PI - 1.)*(9.*la2Q*la2Q +
                        (2.*la3Q + la4Q)*(2.*la3Q + la4Q))/16./M_PI/M_PI -
                         3.*la2Q*(-WFRc1 + 2.*WFRc2)/32./M_PI/M_PI;

    gslpp::complex B3 = -2.*la3Q - la4Q + 3.*(2.*beta3 + beta4)/2. + 3.*(i*M_PI -
                        1.)*(la1Q + la2Q)*(2.*la3Q + la4Q)/16./M_PI/M_PI -
                        (2.*la3Q + la4Q)*WFRc2/32./M_PI/M_PI;

    return ((B1 + B2 - sqrt((B1 - B2)*(B1 - B2) + 4.*B3*B3))/32./M_PI - i/2.).abs();
}


/***********************************/
/* Eigenvalues of the odd 00 block */
/***********************************/

unitarity00oddP_Z2::unitarity00oddP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity00oddP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc2 = myZ2_at_Q(2, 1);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B4 = -la3Q - 2.*la4Q + 3.*(beta3 + 2.*beta4)/2. + (i*M_PI - 1.)*(la3Q*la3Q +
                        4.*la3Q*la4Q + 4.*la4Q*la4Q + 9.*la5Q*la5Q)/16./M_PI/M_PI -
                        (la3Q + la4Q + la5Q)*WFRc2/32./M_PI/M_PI;

    gslpp::complex B6 = -3.*la5Q + 9.*beta5/2. + 6.*(i*M_PI - 1.)*(la3Q + 2.*la4Q)*la5Q/16./M_PI/M_PI -
                        (la4Q + 2.*la5Q)*WFRc2/32./M_PI/M_PI;

    return ((B4 + B6)/16./M_PI - i/2.).abs();
}

unitarity00oddM_Z2::unitarity00oddM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity00oddM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc2 = myZ2_at_Q(2, 1);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B4 = -la3Q - 2.*la4Q + 3.*(beta3 + 2.*beta4)/2. + (i*M_PI - 1.)*(la3Q*la3Q +
                        4.*la3Q*la4Q + 4.*la4Q*la4Q + 9.*la5Q*la5Q)/16./M_PI/M_PI -
                        (la3Q + la4Q + la5Q)*WFRc2/32./M_PI/M_PI;

    gslpp::complex B6 = -3.*la5Q + 9.*beta5/2. + 6.*(i*M_PI - 1.)*(la3Q + 2.*la4Q)*la5Q/16./M_PI/M_PI -
                        (la4Q + 2.*la5Q)*WFRc2/32./M_PI/M_PI;

    return ((B4 - B6)/16./M_PI - i/2.).abs();
}


/************************************/
/* Eigenvalues of the even 01 block */
/************************************/

unitarity01eveP_Z2::unitarity01eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity01eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B7 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la4Q*la4Q)/16./M_PI/M_PI -
                         la1Q*WFRc1/32.*M_PI/M_PI;

    gslpp::complex B8 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la4Q*la4Q)/16./M_PI/M_PI -
                         la2Q*(-WFRc1 + 2.*WFRc2)/32./M_PI/M_PI;

    gslpp::complex B9 = -la4Q + 3.*beta4/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la4Q/16./M_PI/M_PI -
                         la4Q*WFRc2/32./M_PI/M_PI;

    return ((B7 + B8 + sqrt((B7 - B8)*(B7 - B8) + 4.*B9*B9))/32./M_PI - i/2.).abs();
}

unitarity01eveM_Z2::unitarity01eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity01eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B7 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la4Q*la4Q)/16./M_PI/M_PI -
                         la1Q*WFRc1/32.*M_PI/M_PI;

    gslpp::complex B8 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la4Q*la4Q)/16./M_PI/M_PI -
                         la2Q*(-WFRc1 + 2.*WFRc2)/32./M_PI/M_PI;

    gslpp::complex B9 = -la4Q + 3.*beta4/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la4Q/16./M_PI/M_PI -
                         la4Q*WFRc2/32./M_PI/M_PI;

    return ((B7 + B8 - sqrt((B7 - B8)*(B7 - B8) + 4.*B9*B9))/32./M_PI - i/2.).abs();
}


/***********************************/
/* Eigenvalues of the odd 01 block */
/***********************************/

unitarity01oddP_Z2::unitarity01oddP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity01oddP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc2 = myZ2_at_Q(2, 1);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B13 = -la3Q + 3.*beta3/2. + (i*M_PI - 1.)*(la3Q*la3Q + la5Q*la5Q)/16./M_PI/M_PI -
                         (la3Q + la4Q - la5Q)*WFRc2/32./M_PI/M_PI;

    gslpp::complex B15 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*la3Q*la5Q/16./M_PI/M_PI -
                         (la4Q - 2.*la5Q)*WFRc2/32./M_PI/M_PI;

    return ((B13 + B15)/16./M_PI - i/2.).abs();
}

unitarity01oddM_Z2::unitarity01oddM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity01oddM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc2 = myZ2_at_Q(2, 1);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B13 = -la3Q + 3.*beta3/2. + (i*M_PI - 1.)*(la3Q*la3Q + la5Q*la5Q)/16./M_PI/M_PI -
                         (la3Q + la4Q - la5Q)*WFRc2/32./M_PI/M_PI;

    gslpp::complex B15 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*la3Q*la5Q/16./M_PI/M_PI -
                         (la4Q - 2.*la5Q)*WFRc2/32./M_PI/M_PI;

    return ((B13 - B15)/16./M_PI - i/2.).abs();
}


/**********************************/
/* Eigenvalue of the odd 10 block */
/**********************************/

unitarity10odd_Z2::unitarity10odd_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity10odd_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc2 = myZ2_at_Q(2, 1);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B19 = -la3Q + la4Q + 3.*(beta3 - beta4)/2. + (i*M_PI - 1.)*(la3Q -
                          la4Q)*(la3Q - la4Q)/16./M_PI/M_PI - (la3Q - la5Q)*WFRc2/32./M_PI/M_PI;

    return (B19/16./M_PI - i/2.).abs();
}


/************************************/
/* Eigenvalues of the even 11 block */
/************************************/

unitarity11eveP_Z2::unitarity11eveP_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity11eveP_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);
    double WFRc3 = myZ2_at_Q(2, 2);
    double WFRc4 = myZ2_at_Q(2, 3);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B20 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la5Q*la5Q)/16./M_PI/M_PI -
                          la1Q*(WFRc1 - 2.*WFRc2 + WFRc3 + 2.*WFRc4)/32./M_PI/M_PI;

    gslpp::complex B21 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la5Q*la5Q)/16./M_PI/M_PI -
                          la2Q*(-WFRc1 + 2.*WFRc2 - WFRc3 + 2.*WFRc4)/32./M_PI/M_PI;

    gslpp::complex B22 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la5Q/16./M_PI/M_PI -
                          la5Q*2.*WFRc4/32./M_PI/M_PI;

    return ((B20 + B21 + sqrt((B20 - B21)*(B20 - B21) + 4.*B22*B22))/32./M_PI - i/2.).abs();
}

unitarity11eveM_Z2::unitarity11eveM_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity11eveM_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc1 = myZ2_at_Q(2, 0);
    double WFRc2 = myZ2_at_Q(2, 1);
    double WFRc3 = myZ2_at_Q(2, 2);
    double WFRc4 = myZ2_at_Q(2, 3);

    double beta1 = myGTHDM.betalambda1_Z2(la1Q, la3Q, la4Q, la5Q, Yb1Q, Ytau1Q);
    double beta2 = myGTHDM.betalambda2_Z2(la2Q, la3Q, la4Q, la5Q, YtQ, Yb2Q, Ytau2Q);
    double beta5 = myGTHDM.betalambda5_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B20 = -la1Q + 3.*beta1/2. + (i*M_PI - 1.)*(la1Q*la1Q + la5Q*la5Q)/16./M_PI/M_PI -
                          la1Q*(WFRc1 - 2.*WFRc2 + WFRc3 + 2.*WFRc4)/32./M_PI/M_PI;

    gslpp::complex B21 = -la2Q + 3.*beta2/2. + (i*M_PI - 1.)*(la2Q*la2Q + la5Q*la5Q)/16./M_PI/M_PI -
                          la2Q*(-WFRc1 + 2.*WFRc2 - WFRc3 + 2.*WFRc4)/32./M_PI/M_PI;

    gslpp::complex B22 = -la5Q + 3.*beta5/2. + (i*M_PI - 1.)*(la1Q + la2Q)*la5Q/16./M_PI/M_PI -
                          la5Q*2.*WFRc4/32./M_PI/M_PI;

    return ((B20 + B21 - sqrt((B20 - B21)*(B20 - B21) + 4.*B22*B22))/32./M_PI - i/2.).abs();
}


/**********************************/
/* Eigenvalue of the odd 11 block */
/**********************************/

unitarity11odd_Z2::unitarity11odd_Z2(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(SM_i), myZ2_at_Q(3, 5, 0.)
{}

double unitarity11odd_Z2::computeThValue()
{
    gslpp::complex i = gslpp::complex::i();

    myZ2_at_Q = myGTHDM.getGTHDMZ2_at_Q();

    double la1Q = myZ2_at_Q(0, 0);
    double la2Q = myZ2_at_Q(0, 1);
    double la3Q = myZ2_at_Q(0, 2);
    double la4Q = myZ2_at_Q(0, 3);
    double la5Q = myZ2_at_Q(0, 4);

    double YtQ    = myZ2_at_Q(1, 0);
    double Yb1Q   = myZ2_at_Q(1, 1);
    double Yb2Q   = myZ2_at_Q(1, 2);
    double Ytau1Q = myZ2_at_Q(1, 3);
    double Ytau2Q = myZ2_at_Q(1, 4);

    double WFRc4 = myZ2_at_Q(2, 3);

    double beta3 = myGTHDM.betalambda3_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);
    double beta4 = myGTHDM.betalambda4_Z2(la1Q, la2Q, la3Q, la4Q, la5Q, YtQ, Yb1Q, Yb2Q, Ytau1Q, Ytau2Q);

    gslpp::complex B30 = -la3Q - la4Q + 3.*(beta3 + beta4)/2. + (i*M_PI - 1.)*(la3Q +
                          la4Q)*(la3Q + la4Q)/16./M_PI/M_PI - (la3Q + la4Q)*2.*WFRc4/32./M_PI/M_PI;

    return (B30/16./M_PI - i/2.).abs();
}
