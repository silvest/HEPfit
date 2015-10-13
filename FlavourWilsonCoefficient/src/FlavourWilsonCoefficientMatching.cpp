/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficientMatching.h"
#include "FlavourWilsonCoefficient.h"
#include <stdexcept>

FlavourWilsonCoefficientMatching::FlavourWilsonCoefficientMatching(const FlavourWilsonCoefficient & FlavourWilsonCoefficient_i) :

    StandardModelMatching(FlavourWilsonCoefficient_i),
    myCKM(3, 3, 0.),
    myFlavourWilsonCoefficient(FlavourWilsonCoefficient_i)
{}

void FlavourWilsonCoefficientMatching::updateFlavourWilsonCoefficientParameters()
{
    myCKM = myFlavourWilsonCoefficient.getVCKM();
    
    StandardModelMatching::updateSMParameters();
}

