/* 
 * File:   epsilon1.cpp
 * Author: mishima
 */

#include "epsilon1.h"


double epsilon1::getThValue() {  
    double DeltaRhoPrime = 2.0*( sqrt(SM.rhoZ_l(SM.ELECTRON).abs()) - 1.0 );
    
    return DeltaRhoPrime;
}


