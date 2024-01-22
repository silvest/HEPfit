/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LE_ParityViolation.h"
#include "StandardModel.h"

double QWe::computeThValue()
{
    double q2=0.026, y=0.6; // Values from SLAC E158 Moller scattering experiment
    
    return SM.Qwemoller(q2,y);
}


double QWp::computeThValue()
{   
    return SM.Qwp();
}

double QWn::computeThValue()
{  
    return SM.Qwn();
}

double QWAPV::computeThValue()
{
    double qwproton,qwneutron;
    
    qwproton = SM.Qwp();
    qwneutron = SM.Qwn();
    
    return (Z_at * qwproton + N_at * qwneutron);
}



// The following should not be part of this class. Placed here temporarily for testing.

double agminus2muon::computeThValue()
{
    
    return SM.amuon();
}
