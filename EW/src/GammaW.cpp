/* 
 * File:   GammaW.cpp
 * Author: mishima
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "GammaW.h"


GammaW::GammaW(const EW& EW_i) : ThObservable(EW_i) {
    Gamma_W = EW_i.getSM().GammaW();
}

double GammaW::getThValue() {   
    return Gamma_W;
}
        