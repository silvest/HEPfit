/* 
 * File:   Alepton.cpp
 * Author: mishima
 */

#include "Alepton.h"


double Alepton::getThValue() { 
    double A_l = myEW.A_l(SM.ELECTRON);
        
    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        double s4 = s2*s2;
        
        A_l -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
               *( myEW.S() - 4.0*c2*s2*myEW.T() );
    }
  
    return A_l;
}
        

