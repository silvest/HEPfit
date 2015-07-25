/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonP_O_Epsilon.h"

double EpsilonP_O_Epsilon::computeThValue()
{
    return (M_SQRT1_2 * (SM.getReA2_Kd() / SM.getReA0_Kd() / SM.getReA0_Kd()) * ((1 / SM.getReA2_Kd() / SM.getReA0_Kd()) * AmpDS1pp2(NLO).imag() - (1 - SM.getOmega_eta_etap()) * AmpDS1pp0(NLO).imag()));
}
