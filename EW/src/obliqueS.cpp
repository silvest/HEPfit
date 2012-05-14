/* 
 * File:   obliqueS.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueS.h"


obliqueS::obliqueS(const EW& EW_i) : ThObservable(EW_i) {
    S = EW_i.getSM().obliqueS();
}

double obliqueS::getThValue() {   
    return S;
}
 

