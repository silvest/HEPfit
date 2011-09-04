/* 
 * File:   obliqueT.cpp
 * Author: mishima
 */

#include "obliqueT.h"


obliqueT::obliqueT(const EW& EW_i) : ThObservable(EW_i), epsilon_1(EW_i) {
    T = epsilon_1.getThValue()/SM.getAlpha();
}

double obliqueT::getThValue() {   
    return T;
}
 


