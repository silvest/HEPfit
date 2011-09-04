/* 
 * File:   obliqueU.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueU.h"


obliqueU::obliqueU(const EW& EW_i) : ThObservable(EW_i), epsilon_2(EW_i) {
    double s02;
    s02 = 0.5 - sqrt(0.25 - M_PI*EW_i.getAlphaMZ()/sqrt(2.0)
                            /SM.getGF()/pow(SM.getMz(),2.0) );
    U = - epsilon_2.getThValue()/SM.getAlpha()*4.0*s02;
}

double obliqueU::getThValue() {   
    return U;
}
 


