/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Asl.h"
#include "AmpDB2.h"
#include "StandardModelMatching.h"

//#include "std_make_vector.h"

Asl::Asl(const StandardModel& SM_i, QCD::lepton lep_i)
: ThObservable(SM_i), AmpDB2(SM_i) {
    lep=lep_i;
    //setParametersForObservable(make_vector<std::string>() << "Asl_test" );
}

Asl::~Asl() {
}

double Asl::computeThValue() {
    //std::cout << SM.getOptionalParameter("Asl_test") << "\n";
    return getAsl(LO, lep);
}