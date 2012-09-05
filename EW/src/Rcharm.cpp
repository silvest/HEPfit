/* 
 * File:   Rcharm.cpp
 * Author: mishima
 */

#include "Rcharm.h"


double Rcharm::getThValue() {   
    double R0_c = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        double s4 = s2*s2;

        R0_c -= 9.0*alpha*(9.0-36.0*s2+16.0*s4)
                /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                 *( myEW.S() - 4.0*c2*s2*myEW.T() );
    }

    return R0_c;
}
        

