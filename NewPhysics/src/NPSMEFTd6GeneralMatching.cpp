/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6GeneralMatching.h"
#include "NPSMEFTd6General.h"
#include <stdexcept>

NPSMEFTd6GeneralMatching::NPSMEFTd6GeneralMatching(const NPSMEFTd6General & NPSMEFTd6General_i) :

    StandardModelMatching(NPSMEFTd6General_i),
    mySMEFT(NPSMEFTd6General_i)
{}

void NPSMEFTd6GeneralMatching::updateNPSMEFTd6GeneralParameters()
{
        
    //mySMEFT.getSMEFTEvol().EvolveToBasis("Numeric",mySMEFT.getLambda_NP(),mySMEFT.getMuw(),mySMEFT.getSMEFTBasisFlag());
    
    StandardModelMatching::updateSMParameters();
}

NPSMEFTd6GeneralMatching::~NPSMEFTd6GeneralMatching()
{}
