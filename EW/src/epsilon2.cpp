/* 
 * File:   epsilon2.cpp
 * Author: mishima
 */

#include "epsilon2.h"


double epsilon2::getThValue() {  
    double eps2 = SM.epsilon2();
    
    if ( myEW.checkModelForSTU() )
        eps2 += myEW.Uhat() - myEW.V() - myEW.W() 
                + 2.0*sqrt(SM.s02())/sqrt(SM.c02())*myEW.X();
    
    return eps2; 
}
 

