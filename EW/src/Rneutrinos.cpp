/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rneutrinos.h"
#include "StandardModel.h"

double Rneutrinos::computeThValue()
{
    double Rnue = SM.R0_f(SM.getLeptons(SM.NEUTRINO_1));
    double Rnumu = SM.R0_f(SM.getLeptons(SM.NEUTRINO_2));
    double Rnutau = SM.R0_f(SM.getLeptons(SM.NEUTRINO_3));    
    
    return (Rnue + Rnumu + Rnutau);
}


