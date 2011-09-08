/* 
 * File:   obliqueEpsilon2.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueEpsilon2.h"


obliqueEpsilon2::obliqueEpsilon2(const EW& EW_i) : ThObservable(EW_i) {
    double DeltaRhoPrime, DeltaKappaPrime, DeltaRW, s02, c02;
    DeltaRhoPrime = 2.0*(sqrt(EW_i.getRhoZ_l(SM.ELECTRON).real()) - 1.0);
    DeltaRW = 1.0 - M_PI*EW_i.getAlphaMz()/sqrt(2.0)/SM.getGF()
                    /( 1.0 - pow(EW_i.getMw()/SM.getMz(),2.0) )
                    /pow(EW_i.getMw(),2.0);
    s02 = 0.5 - sqrt(0.25 - M_PI*EW_i.getAlphaMz()/sqrt(2.0)
                            /SM.getGF()/pow(SM.getMz(),2.0) );
    c02 = 1.0 - s02;
    DeltaKappaPrime = EW_i.sin2thetaEff(SM.ELECTRON)/s02 - 1.0;
    epsilon_2 = c02*DeltaRhoPrime + s02*DeltaRW/(c02-s02) - 2.0*s02*DeltaKappaPrime;
}

double obliqueEpsilon2::getThValue() {   
    return epsilon_2;
}
 

