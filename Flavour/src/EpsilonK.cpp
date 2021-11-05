/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EpsilonK.h"
#include "StandardModel.h"
#include "std_make_vector.h"

EpsilonK::EpsilonK(const StandardModel& SM_i): ThObservable(SM_i), AmpDK2(SM_i)
{
    setParametersForObservable(make_vector<std::string>() << "DeltaMK" << "KbarEpsK" << "phiEpsK" << "DeltattEpsK");
}

double EpsilonK::computeThValue()
{
#if SUSYFIT_DEBUG & 2
    std::cout << "amplitude = " <<  AmpDK(FULLNLO).imag() << std::endl;
#endif
    return(SM.getCepsK() / SM.getOptionalParameter("DeltaMK") * AmpDK(FULLNLO).imag() * SM.getOptionalParameter("KbarEpsK") * 
            // Tarantino et al 2021
            1.01 *
            sin(SM.getOptionalParameter("phiEpsK") * M_PI / 180.));
}