/* 
 * File:   ZFAFBmuLEP2.cpp
 * Author: mishima
 */

#include "ZFAFBmuLEP2.h"


double ZFAFBmuLEP2::getThValue() { 
    double myZFAFBmuLEP2, dummy;
    myZF.calcXS_AFB(2, sqrt_s, &dummy, &myZFAFBmuLEP2);
    return myZFAFBmuLEP2;
}


