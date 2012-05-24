/* 
 * File:   sigmaHadron.cpp
 * Author: mishima
 */

#include "sigmaHadron.h"


sigmaHadron::sigmaHadron(const EW& EW_i) : ThObservable(EW_i) {
    sigma_had = EW_i.sigma0_had();

    if ( EW_i.checkModelForSTU() ) {
        double alpha = EW_i.getSM().getAle();
        double Mz = EW_i.getSM().getMz();
        double c2 = EW_i.c2();
        double s2 = EW_i.s2();
        double s4 = s2*s2;
        double s6 = s4*s2;        
        double s8 = s6*s2;

        sigma_had -= 72.0*M_PI*alpha
                     *(729.0-4788.0*s2+8352.0*s4-6176.0*s6+640.0*s8)
                     /Mz/Mz/pow(63.0-120.0*s2+160.0*s4, 3.0)/(c2-s2)
                     *( EW_i.S() - 4.0*c2*s2*EW_i.T() );         
    }
}

double sigmaHadron::getThValue() {   
    return sigma_had;
}
        


