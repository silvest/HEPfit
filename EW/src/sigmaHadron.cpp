/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"


double sigmaHadron::computeThValue() 
{ 
    double sigma_had = SM.sigma0_had();

    return ( sigma_had*SM.GeVminus2_to_nb );
}
        


