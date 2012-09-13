/* 
 * File:   sigmaHadron.cpp
 * Author: mishima
 */

#include "sigmaHadron.h"


double sigmaHadron::getThValue() { 
    double sigma_had = myEW.sigma0_had();

    if ( myEW.checkModelForSTU() ) {
        double alpha = myEW.getSM().alphaMz();
        double Mz = myEW.getSM().getMz();
        double c2 = myEW.c2();
        double s2 = myEW.s2();
        double s4 = s2*s2;
        double s6 = s4*s2;        
        double s8 = s6*s2;

        sigma_had -= 72.0*M_PI*alpha
                     *(729.0-4788.0*s2+8352.0*s4-6176.0*s6+640.0*s8)
                     /Mz/Mz/pow(63.0-120.0*s2+160.0*s4, 3.0)/(c2-s2)
                     *( myEW.S() - 4.0*c2*s2*myEW.T() );         
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


