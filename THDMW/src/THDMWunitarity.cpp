/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWunitarity.h"

THDMWunitarity::THDMWunitarity(const StandardModel& SM_i)
: myTHDMW(static_cast<const THDMW*> (&SM_i)), Smatrix(18,18,0.),
        Seigenvectors(18,18,0.), Seigenvalues(18,0.)
{}

THDMWunitarity::~THDMWunitarity() 
{}

bool THDMWunitarity::CalcSeigen(gslpp::matrix<gslpp::complex>& Seigenvectors_i, gslpp::vector<double>& Seigenvalues_i)
{
    double lambda1 = myTHDMW->getTHDMW_lambda1();
    
    gslpp::complex i = gslpp::complex::i();
    
    Smatrix.assign(0,0, lambda1);
    
    Smatrix.eigensystem(Seigenvectors_i, Seigenvalues_i);
    
    return true;
}

gslpp::vector<double> THDMWunitarity::getSeigenvalues()
{
    CalcSeigen(Seigenvectors,Seigenvalues);
    
    return Seigenvalues;
}


THDMWunitarity1::THDMWunitarity1(const StandardModel& SM_i)
: ThObservable(SM_i), myTHDMWunitarity(SM_i)
{}
//
double THDMWunitarity1::computeThValue()
{
//    return (myTHDMWunitarity.getSeigenvalues())(0);
    return 0.0;
}

