/* 
 * File:   ZFAFBcharmLEP2.cpp
 * Author: mishima
 */

#include "ZFAFBcharmLEP2.h"


double ZFAFBcharmLEP2::getThValue() { 
    double myZFAFBcharmLEP2, dummy;
    myZF.calcXS_AFB(6, sqrt_s, &dummy, &myZFAFBcharmLEP2);
    return myZFAFBcharmLEP2;
}


