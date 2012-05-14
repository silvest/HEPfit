/* 
 * File:   obliqueU.cpp
 * Author: mishima
 */

#include <cmath>
#include "obliqueU.h"


obliqueU::obliqueU(const EW& EW_i) : ThObservable(EW_i) {
    U = EW_i.getSM().obliqueU();
}

double obliqueU::getThValue() {   
    return U;
}
 


