/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"
#include "StandardModel.h"

double sigmaHadron::computeThValue()
{
    return ( SM.sigma0_had() * SM.GeVminus2_to_nb);
}



