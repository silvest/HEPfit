/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMWunitarity.h"
#include "THDMW.h"
#include "THDMWcache.h"

THDMWunitarityLO::THDMWunitarityLO(const StandardModel& SM_i, unsigned int index_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{
    index = index_i;
}

THDMWunitarityLO::~THDMWunitarityLO() 
{}

double THDMWunitarityLO::computeThValue()
{
    if( index > 10 ) {
        throw std::runtime_error("Index out of range in THDMWunitarityLO");
    }
    return (myTHDMW.getMyTHDMWCache()->unitarityeigenvalues(index)).abs();
}

THDMWunitarityNLO::THDMWunitarityNLO(const StandardModel& SM_i, unsigned int index_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{
    index = index_i;
}

THDMWunitarityNLO::~THDMWunitarityNLO() 
{}

double THDMWunitarityNLO::computeThValue()
{
    if( index > 10 ) {
        throw std::runtime_error("Index out of range in THDMWunitarityNLO");
    }

    gslpp::complex a0 = myTHDMW.getMyTHDMWCache()->unitarityeigenvalues(index);
    gslpp::complex a1 = myTHDMW.getMyTHDMWCache()->NLOunitarityeigenvalues(index);
    
    return ( a0 * a0 + 2.0 * a0 * a1.real() ).abs();
}

THDMWunitarityNLOp::THDMWunitarityNLOp(const StandardModel& SM_i, unsigned int index_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{
    index = index_i;
}

THDMWunitarityNLOp::~THDMWunitarityNLOp() 
{}

double THDMWunitarityNLOp::computeThValue()
{
    if( index > 10 ) {
        throw std::runtime_error("Index out of range in THDMWunitarityNLOp");
    }

    gslpp::complex a0 = myTHDMW.getMyTHDMWCache()->unitarityeigenvalues(index);
    gslpp::complex a1 = myTHDMW.getMyTHDMWCache()->NLOunitarityeigenvalues(index);
    
    return ( a0 + a1.real() ).abs();
}

THDMWunitarityRp::THDMWunitarityRp(const StandardModel& SM_i, unsigned int index_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i))
{
    index = index_i;
}

THDMWunitarityRp::~THDMWunitarityRp() 
{}

double THDMWunitarityRp::computeThValue()
{
    if( index > 10 ) {
        throw std::runtime_error("Index out of range in THDMWunitarityRp");
    }

    gslpp::complex a0 = myTHDMW.getMyTHDMWCache()->unitarityeigenvalues(index);
    gslpp::complex a1 = myTHDMW.getMyTHDMWCache()->NLOunitarityeigenvalues(index);
    double Rpeps = myTHDMW.getMyTHDMWCache()->RpepsTHDMW;

    double Rp = 0.01;
    if(a0.abs()>Rpeps)
    {
        Rp = ( a1 / a0 ).abs();
    }
    
    return Rp;
}
