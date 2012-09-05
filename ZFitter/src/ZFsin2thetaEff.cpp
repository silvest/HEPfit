/* 
 * File:   ZFsin2thetaEff.cpp
 * Author: mishima
 */

#include "ZFsin2thetaEff.h"


double ZFsin2thetaEff::getThValue() {
    return myZF.s2teff_f(1);
}

