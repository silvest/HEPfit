/* 
 * File:   Mw.cpp
 * Author: mishima
 */

#include "Mw.h"


double Mw::getThValue() {
    double myMw;
    if (bCHMN)  
        myMw = myEW_CHMN.Mw();
    else {
        myMw = SM.Mw();    

        if ( myEW.checkModelForSTU() ) {
            double alpha = SM.alphaMz();
            double c = sqrt(myEW.c02());
            double c2 = myEW.c02();
            double s2 = myEW.s02();
            
            myMw -= alpha*c*SM.getMz()/4.0/(c2-s2)
                    *( myEW.S() - 2.0*c2*myEW.T() - (c2-s2)*myEW.U()/2.0/s2 );
        }
    }
    
    return myMw;
}

