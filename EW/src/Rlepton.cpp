/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rlepton.h"
#include "StandardModel.h"

double Rlepton::computeThValue()
{
    double Re = SM.R0_f(SM.getLeptons(SM.ELECTRON));
    double Rmu = SM.R0_f(SM.getLeptons(SM.MU));
    double Rtau = SM.R0_f(SM.getLeptons(SM.TAU));    
    
    return (Re + Rmu + Rtau)/3.;
}

double Relectron::computeThValue()
{
    return SM.R0_f(SM.getLeptons(SM.ELECTRON));
}

double Rmuon::computeThValue()
{
    return SM.R0_f(SM.getLeptons(SM.MU));
}

double Rtau::computeThValue()
{
    return SM.R0_f(SM.getLeptons(SM.TAU));
}


