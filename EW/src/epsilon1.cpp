/* 
 * File:   epsilon1.cpp
 * Author: mishima
 */

#include "epsilon1.h"


double epsilon1::getThValue() {  
    double DeltaRhoPrime = 2.0*( sqrt(SM.rhoZ_l(SM.ELECTRON).abs()) - 1.0 );
    double eps1 = DeltaRhoPrime;
    
    if ( myEW.checkModelForSTU() )
        eps1 += myEW.That() - myEW.W() 
                + 2.0*sqrt(SM.s02())/sqrt(SM.c02())*myEW.X()
                - SM.s02()/SM.c02()*myEW.Y();
    
    return eps1;
}


