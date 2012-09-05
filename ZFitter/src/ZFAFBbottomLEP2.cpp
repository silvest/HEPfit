/* 
 * File:   ZFAFBbottomLEP2.cpp
 * Author: mishima
 */

#include "ZFAFBbottomLEP2.h"


double ZFAFBbottomLEP2::getThValue() { 
    double myZFAFBbottomLEP2, dummy;
    myZF.calcXS_AFB(9, sqrt_s, &dummy, &myZFAFBbottomLEP2);
    return myZFAFBbottomLEP2;
}


