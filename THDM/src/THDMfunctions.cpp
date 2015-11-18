/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMfunctions.h"
#include "StandardModel.h"

double THDMfunctions::computeThValue()
{
    return 0;
}



gslpp::complex THDMfunctions::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}



gslpp::complex THDMfunctions::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1-x)*(log((1+sqrt(1-x))/(1-sqrt(1-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1)*asin(sqrt(1.0/x));
        return gg1;
    }
}



gslpp::complex THDMfunctions::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
    }



gslpp::complex THDMfunctions::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
    }



int THDMfunctions::HSTheta (const double x) const{
    if(x<0)return 0.0;
    else return 1.0;
}



double THDMfunctions::KaellenFunction (const double a, const double b, const double c) const{
    double x = (a*a-b*b-c*c)*(a*a-b*b-c*c)-4.0*b*b*c*c;
    return HSTheta(x)*sqrt(std::abs(x))/(2.0*a);
}
