/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonP_O_Epsilon.h"

double EpsilonP_O_Epsilon::getThValue(){
    
    double om = SM.getReA2_kd() / SM.getReA0_kd();
    
    double a = 0.707106781*(om/SM.getReA0_kd())*((1/om)*AmpDS1pp2(NLO).imag()-
               (1-SM.getOmega_eta_etap())*AmpDS1pp0(NLO).imag());
    
    return a;
}
