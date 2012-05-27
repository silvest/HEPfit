/* 
 * File:   Rcharm.cpp
 * Author: mishima
 */

#include "Rcharm.h"


Rcharm::Rcharm(const EW& EW_i) : ThObservable(EW_i) {
    R0_c = EW_i.Gamma_q(SM.CHARM)/EW_i.Gamma_had();

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().alphaMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;

        R0_c -= 9.0*alpha*(9.0-36.0*s2+16.0*s4)
                /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                 *( EW_i.S() - 4.0*c2*s2*EW_i.T() );
    }
}

double Rcharm::getThValue() {   
    return R0_c;
}
        

