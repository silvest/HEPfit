/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMunitarity.h"

GMunitarityLO::GMunitarityLO(const StandardModel& SM_i, unsigned int index_i)
: ThObservable(SM_i),myGM(static_cast<const GeorgiMachacek&> (SM_i))
{
    index = index_i;
}

GMunitarityLO::~GMunitarityLO() 
{}

double GMunitarityLO::computeThValue()
{
    if( index > 3 ) {
        throw std::runtime_error("Index out of range in GMunitarityLO");
    }
    return (myGM.getMyGMCache()->unitarityeigenvalues(index)).abs();
}
