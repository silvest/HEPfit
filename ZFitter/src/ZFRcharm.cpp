/* 
 * File:   ZFRcharm.cpp
 * Author: mishima
 */

#include "ZFRcharm.h"


double ZFRcharm::getThValue() {
    return myZF.Gamma_f(6)/myZF.Gamma_had();
}


