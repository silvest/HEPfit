/* 
 * File:   GammaW.cpp
 * Author: mishima
 */

#include "GammaW.h"


double GammaW::getThValue() {  
    double Gamma_W = myEW.getSM().GammaW();
        
    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c = sqrt(myEW.c02());
        double c2 = myEW.c02();
        double s2 = myEW.s02();
    
        double Wbar = 0.0;        
        if (myEW.getSM().ModelName()=="NewPhysicsSTUVWXY") {
            Wbar = (myEW.V() - myEW.W())/myEW.getSM().alphaMz();
        }

        Gamma_W -= 3.0*alpha*alpha*c*myEW.getSM().getMz()/8.0/s2/(c2-s2)
                   *( myEW.S() - 2.0*c2*myEW.T() - (c2-s2)*myEW.U()/2.0/s2 
                      - 2.0*(c2 - s2)*Wbar );
    }
 
    return Gamma_W;
}
