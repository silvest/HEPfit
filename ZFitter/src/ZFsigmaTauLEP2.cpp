/* 
 * File:   ZFsigmaTauLEP2.cpp
 * Author: mishima
 */

#include "ZFsigmaTauLEP2.h"


double ZFsigmaTauLEP2::getThValue() { 
    double myZFsigmaTauLEP2, dummy;
    myZF.calcXS_AFB(3, sqrt_s, &myZFsigmaTauLEP2, &dummy);
    myZFsigmaTauLEP2 *= 1000.0;// nb --> pb
    return myZFsigmaTauLEP2;
}


