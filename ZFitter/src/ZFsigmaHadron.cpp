/* 
 * File:   ZFsigmaHadron.cpp
 * Author: mishima
 */

#include "ZFsigmaHadron.h"


double ZFsigmaHadron::getThValue() {   
    return ( 12.0*M_PI/pow(SM.getMz(),2.0)
             *myZF.Gamma_f(1)*myZF.Gamma_had()/myZF.Gamma_Z()/myZF.Gamma_Z()
             *GeVminus2_to_nb );
}
