/* 
 * File:   obliqueT.cpp
 * Author: mishima
 */

#include "obliqueT.h"


obliqueT::obliqueT(const EW& EW_i) : ThObservable(EW_i) {
    T = EW_i.getSM().obliqueT();
}

double obliqueT::getThValue() {   
    return T;
}
 


