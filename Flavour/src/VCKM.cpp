/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "VCKM.h"

VCKM::VCKM(const StandardModel& SM_i, unsigned int obsFlag_1, unsigned int obsFlag_2) : ThObservable(SM_i) {
    if (obsFlag_1 > 0 && obsFlag_1 < 4 && obsFlag_2 > 0 && obsFlag_2 < 4) {
        obs_1 = obsFlag_1;
        obs_2 = obsFlag_2;
    }
    else throw std::runtime_error("obsFlag in CKM(myFlavour, obsFlag_1, obsFlag_1) called from ThFactory::ThFactory() can only be 1 - 3 corresponding to the CKM matrix");
}

VCKM::~VCKM() {
}

double VCKM::computeThValue() { 
    
    if (obs_1 == 1 && obs_2 == 1) return(SM.getCKM().getVud());
    if (obs_1 == 1 && obs_2 == 2) return(SM.getCKM().getVus());
    if (obs_1 == 1 && obs_2 == 3) return(SM.getCKM().getVub());
    if (obs_1 == 2 && obs_2 == 1) return(SM.getCKM().getVcd());
    if (obs_1 == 2 && obs_2 == 2) return(SM.getCKM().getVcs());
    if (obs_1 == 2 && obs_2 == 3) return(SM.getCKM().getVcb());
    if (obs_1 == 3 && obs_2 == 1) return(SM.getCKM().getVtd());
    if (obs_1 == 3 && obs_2 == 2) return(SM.getCKM().getVts());
    if (obs_1 == 3 && obs_2 == 3) return(SM.getCKM().getVtb());
    else throw std::runtime_error("obsFlag in CKM(myFlavour, obsFlag_1, obsFlag_1) called from ThFactory::ThFactory() can only be 1 - 3 corresponding to the CKM matrix");
}