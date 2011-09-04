/* 
 * File:   obliqueEpsilon1.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueEpsilon1.h"


obliqueEpsilon1::obliqueEpsilon1(const EW& EW_i) : ThObservable(EW_i) {
    double DeltaRhoPrime;
    DeltaRhoPrime = 2.0*(sqrt(EW_i.getRhoZ_l(SM.ELECTRON).real()) - 1.0);
    epsilon_1 = DeltaRhoPrime;
}

double obliqueEpsilon1::getThValue() {   
    return epsilon_1;
}
        


