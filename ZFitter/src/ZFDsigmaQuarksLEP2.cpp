/* 
 * File:   ZFDsigmaQuarksLEP2.cpp
 * Author: mishima
 */

#include "ZFDsigmaQuarksLEP2.h"


double ZFDsigmaQuarksLEP2::getThValue() { 
    double myZFDsigmaQuarksLEP2 = 0.0, tmp;
    for (int i=4; i<=9; i++) {
        if (i!=8) {
            myZF.calcDXS(i, sqrt_s, cos_theta, &tmp);
            myZFDsigmaQuarksLEP2 += tmp;
        }
    }
    myZFDsigmaQuarksLEP2 *= 1000.0;// nb --> pb
    return myZFDsigmaQuarksLEP2;
}
