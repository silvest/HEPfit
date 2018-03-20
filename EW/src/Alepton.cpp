/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Alepton.h"
#include "StandardModel.h"

double Alepton::computeThValue()
{
    double Ae = SM.A_f(SM.getLeptons(SM.ELECTRON));
    double Amu = SM.A_f(SM.getLeptons(SM.MU));
    double Atau = SM.A_f(SM.getLeptons(SM.TAU));  

    return (Ae + Amu + Atau)/3.;
}

double Aelectron::computeThValue()
{
    return SM.A_f(SM.getLeptons(SM.ELECTRON));
}

double Amuon::computeThValue()
{
    return SM.A_f(SM.getLeptons(SM.MU));
}

double Atau::computeThValue()
{
    return SM.A_f(SM.getLeptons(SM.TAU));
}


