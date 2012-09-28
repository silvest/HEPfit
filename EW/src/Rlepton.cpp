/* 
 * File:   Rlepton.cpp
 * Author: mishima
 */

#include "Rlepton.h"


double Rlepton::getThValue() {  
    double R0_l = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
    
    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c02();
        double s2 = myEW.s02();
        double s4 = s2*s2;

        R0_l += 8.0*alpha*(3.0-2.0*s2)*(1.0-5.0*s2)
                /3.0/pow(1.0-4.0*s2+8.0*s4, 2.0)/(c2-s2)
                *( myEW.S() - 4.0*c2*s2*myEW.T() );                
    }
 
    return R0_l;
}
        

