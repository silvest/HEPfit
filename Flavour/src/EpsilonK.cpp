/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonK.h"

double EpsilonK::computeThValue()
{
#if SUSYFIT_DEBUG & 2
    std::cout << "amplitude = " <<  AmpDK(FULLNLO).imag() << std::endl;
#endif
    return(1. / SM.getDeltaMK() * AmpDK(FULLNLO).imag() * SM.getKbarEpsK() * 
            sin(SM.getphiEpsK() * M_PI / 180.));
}