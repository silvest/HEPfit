/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFsigmaHadron.h"


double ZFsigmaHadron::computeThValue() {   
    return ( 12.0*M_PI/pow(SM.getMz(),2.0)
             *myZF.Gamma_f(1)*myZF.Gamma_had()/myZF.Gamma_Z()/myZF.Gamma_Z()
             *SM.GeVminus2_to_nb );
}
