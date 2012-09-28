/* 
 * File:   GammaZ.cpp
 * Author: mishima
 */

#include "GammaZ.h"


double GammaZ::getThValue() { 
    double Gamma_Z = myEW.Gamma_Z();

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c02();
        double s2 = myEW.s02();
        double s4 = s2*s2;
        
        Gamma_Z += alpha*alpha*myEW.getSM().getMz()/72.0/c2/s2/(c2-s2)
                   *( -10.0*(3.0-8.0*s2)*myEW.S() 
                      + (63.0-126.0*s2-40.0*s4)*myEW.T() );
    }
  
    return Gamma_Z;
}
        
