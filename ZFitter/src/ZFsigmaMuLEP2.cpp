/* 
 * File:   ZFsigmaMuLEP2.cpp
 * Author: mishima
 */

#include "ZFsigmaMuLEP2.h"


double ZFsigmaMuLEP2::getThValue() { 
    double myZFsigmaMuLEP2, dummy;
    //myZF.calcXS_AFB(2, sqrt_s, &myZFsigmaMuLEP2, &dummy);
    myZF.calcDXS(2,sqrt_s,1.0,&myZFsigmaMuLEP2);
    myZFsigmaMuLEP2 *= 1000.0;// nb --> pb
    return myZFsigmaMuLEP2;
}


