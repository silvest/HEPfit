/* 
 * File:   obliqueEpsilon3.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueEpsilon3.h"


obliqueEpsilon3::obliqueEpsilon3(const EW& EW_i) : ThObservable(EW_i) {
    double DeltaRhoPrime, DeltaKappaPrime, s02, c02;
    DeltaRhoPrime = 2.0*(sqrt(EW_i.getRhoZ_l(SM.ELECTRON).real()) - 1.0);
    s02 = 0.5 - sqrt(0.25 - M_PI*EW_i.getAlphaMZ()/sqrt(2.0)
                            /SM.getGF()/pow(SM.getMz(),2.0) );
    c02 = 1.0 - s02;
    DeltaKappaPrime = EW_i.sin2thetaEff(SM.ELECTRON)/s02 - 1.0;
    epsilon_3 = c02*DeltaRhoPrime + (c02-s02)*DeltaKappaPrime;
}

double obliqueEpsilon3::getThValue() {   
    return epsilon_3;
}
 


