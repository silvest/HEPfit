/* 
 * File:   Mw.cpp
 * Author: mishima
 */

#include "Mw.h"


double Mw::getThValue() {   
    double myMw = myEW.getSM().Mw();    

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c = sqrt(myEW.c2());
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        
        myMw -= alpha*c*myEW.getSM().getMz()/4.0/(c2-s2)
                *( myEW.S() - 2.0*c2*myEW.T() - (c2-s2)*myEW.U()/2.0/s2 );
    }
    
    return myMw;
}

