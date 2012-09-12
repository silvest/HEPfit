/* 
 * File:   PtauPol.cpp
 * Author: mishima
 */

#include "PtauPol.h"


double PtauPol::getThValue() {  
    double P_tau_pol = myEW.A_l(SM.TAU);

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        double s4 = s2*s2;
    
        P_tau_pol -= 4.0*alpha*s2/pow(1.0-4.0*s2+8.0*s4, 2.0)
                     *( myEW.S() - 4.0*c2*s2*myEW.T() );
    }
 
    return P_tau_pol;
}
        
