/* 
 * File:   epsilon3.cpp
 * Author: mishima
 */

#include "epsilon3.h"


double epsilon3::getThValue() {  
    double eps3 = SM.epsilon3();
    
    if ( myEW.checkModelForSTU() )
        eps3 += myEW.Shat() - myEW.W() 
                + myEW.X()/sqrt(SM.s02())/sqrt(SM.c02()) - myEW.Y();
    
    return eps3;
}


