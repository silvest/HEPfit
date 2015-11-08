/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "VCKM.h"

VCKM::VCKM(const StandardModel& SM_i, unsigned int obsFlag_1, unsigned int obsFlag_2) 
: ThObservable(SM_i) 
{
    if (obsFlag_1 > 0 && obsFlag_1 < 4 && obsFlag_2 > 0 && obsFlag_2 < 4) {
        obs_1 = obsFlag_1;
        obs_2 = obsFlag_2;
    }
    else throw std::runtime_error("obsFlag in CKM(myFlavour, obsFlag_1, obsFlag_1) called from ThFactory::ThFactory() can only be 1 - 3 corresponding to the CKM matrix");
}

VCKM::~VCKM() 
{}

double VCKM::computeThValue() 
{ 
    
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

Abslam_t::Abslam_t(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_t::computeThValue()
{   
    return (SM.getCKM().computelamt()).abs();
}

Abslam_c::Abslam_c(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_c::computeThValue()
{   
    return (SM.getCKM().computelamc()).abs();
}

Abslam_u::Abslam_u(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_u::computeThValue()
{   
    return (SM.getCKM().computelamu()).abs();
}

Abslam_td::Abslam_td(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_td::computeThValue()
{   
    return (SM.getCKM().computelamt_d()).abs();
}

Abslam_cd::Abslam_cd(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_cd::computeThValue()
{   
    return (SM.getCKM().computelamc_d()).abs();
}

Abslam_ud::Abslam_ud(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_ud::computeThValue()
{   
    return (SM.getCKM().computelamu_d()).abs();
}

Abslam_ts::Abslam_ts(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_ts::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).abs();
}

Abslam_cs::Abslam_cs(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_cs::computeThValue()
{   
    return (SM.getCKM().computelamc_s()).abs();
}

Abslam_us::Abslam_us(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Abslam_us::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).abs();
}