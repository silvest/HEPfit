/* 
 * File:   ZFAFBtauLEP2.cpp
 * Author: mishima
 */

#include "ZFAFBtauLEP2.h"


double ZFAFBtauLEP2::getThValue() { 
    double myZFAFBtauLEP2, dummy;
    myZF.calcXS_AFB(3, sqrt_s, &dummy, &myZFAFBtauLEP2);
    return myZFAFBtauLEP2;
}


