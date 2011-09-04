/* 
 * File:   obliqueS.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueS.h"


obliqueS::obliqueS(const EW& EW_i) : ThObservable(EW_i), epsilon_3(EW_i) {
    double s02;
    s02 = 0.5 - sqrt(0.25 - M_PI*EW_i.getAlphaMZ()/sqrt(2.0)
                            /SM.getGF()/pow(SM.getMz(),2.0) );
    S = epsilon_3.getThValue()/SM.getAle()*4.0*s02;
}

double obliqueS::getThValue() {   
    return S;
}
 

