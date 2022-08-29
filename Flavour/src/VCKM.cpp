/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "VCKM.h"
#include "StandardModel.h"

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
    
    if (obs_1 == 1 && obs_2 == 1) return(SM.getCKM().getV_ud().abs());
    if (obs_1 == 1 && obs_2 == 2) return(SM.getCKM().getV_us().abs());
    if (obs_1 == 1 && obs_2 == 3) return(SM.getCKM().getV_ub().abs());
    if (obs_1 == 2 && obs_2 == 1) return(SM.getCKM().getV_cd().abs());
    if (obs_1 == 2 && obs_2 == 2) return(SM.getCKM().getV_cs().abs());
    if (obs_1 == 2 && obs_2 == 3) return(SM.getCKM().getV_cb().abs());
    if (obs_1 == 3 && obs_2 == 1) return(SM.getCKM().getV_td().abs());
    if (obs_1 == 3 && obs_2 == 2) return(SM.getCKM().getV_ts().abs());
    if (obs_1 == 3 && obs_2 == 3) return(SM.getCKM().getV_tb().abs());
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

Imlam_t::Imlam_t(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_t::computeThValue()
{   
    return (SM.getCKM().computelamt()).imag();
}

Imlam_c::Imlam_c(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_c::computeThValue()
{   
    return (SM.getCKM().computelamc()).imag();
}

Imlam_u::Imlam_u(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_u::computeThValue()
{   
    return (SM.getCKM().computelamu()).imag();
}

Imlam_td::Imlam_td(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_td::computeThValue()
{   
    return (SM.getCKM().computelamt_d()).imag();
}

Imlam_cd::Imlam_cd(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_cd::computeThValue()
{   
    return (SM.getCKM().computelamc_d()).imag();
}

Imlam_ud::Imlam_ud(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_ud::computeThValue()
{   
    return (SM.getCKM().computelamu_d()).imag();
}

Imlam_ts::Imlam_ts(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_ts::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).imag();
}

Imlam_cs::Imlam_cs(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_cs::computeThValue()
{   
    return (SM.getCKM().computelamc_s()).imag();
}

Imlam_us::Imlam_us(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Imlam_us::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).imag();
}

Relam_t::Relam_t(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_t::computeThValue()
{   
    return (SM.getCKM().computelamt()).real();
}

Relam_c::Relam_c(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_c::computeThValue()
{   
    return (SM.getCKM().computelamc()).real();
}

Relam_u::Relam_u(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_u::computeThValue()
{   
    return (SM.getCKM().computelamu()).real();
}

Relam_td::Relam_td(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_td::computeThValue()
{   
    return (SM.getCKM().computelamt_d()).real();
}

Relam_cd::Relam_cd(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_cd::computeThValue()
{   
    return (SM.getCKM().computelamc_d()).real();
}

Relam_ud::Relam_ud(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_ud::computeThValue()
{   
    return (SM.getCKM().computelamu_d()).real();
}

Relam_ts::Relam_ts(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_ts::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).real();
}

Relam_cs::Relam_cs(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_cs::computeThValue()
{   
    return (SM.getCKM().computelamc_s()).real();
}

Relam_us::Relam_us(const StandardModel& SM_i)
: ThObservable(SM_i)
{}

double Relam_us::computeThValue()
{   
    return (SM.getCKM().computelamt_s()).real();
}

CKM_Alpha::CKM_Alpha(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Alpha::computeThValue()
{
    return SM.getCKM().computeAlpha()/M_PI*180.;
}

CKM_Beta::CKM_Beta(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Beta::computeThValue()
{
    return SM.getCKM().computeBeta()/M_PI*180.;
}

CKM_Betas::CKM_Betas(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Betas::computeThValue()
{
    return SM.getCKM().computeBetas()/M_PI*180.;
}

CKM_Gamma::CKM_Gamma(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Gamma::computeThValue()
{
    return SM.getCKM().computeGamma()/M_PI*180.;
}

CKM_2BpG::CKM_2BpG(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_2BpG::computeThValue()
{
    return ((SM.getCKM().computeGamma() + 2. * SM.getCKM().computeBeta())/M_PI*180.);
}

CKM_S2Beta::CKM_S2Beta(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_S2Beta::computeThValue()
{
    return sin(2. * SM.getCKM().computeBeta());
}

CKM_C2Beta::CKM_C2Beta(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_C2Beta::computeThValue()
{
    return cos(2. * SM.getCKM().computeBeta());
}

CKM_SinTheta12::CKM_SinTheta12(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_SinTheta12::computeThValue()
{
    return SM.getCKM().gets12();
}

CKM_SinTheta13::CKM_SinTheta13(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_SinTheta13::computeThValue()
{
    return SM.getCKM().gets13();
}

CKM_SinTheta23::CKM_SinTheta23(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_SinTheta23::computeThValue()
{
    return SM.getCKM().gets23();
}

CKM_Delta::CKM_Delta(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Delta::computeThValue()
{
    return SM.getCKM().getdelta();
}

J_CP::J_CP(const StandardModel& SM_i) : ThObservable(SM_i) {}

double J_CP::computeThValue()
{
    return SM.getCKM().getJcp();
}

CKM_Rt::CKM_Rt(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Rt::computeThValue()
{
    return SM.getCKM().computeRt();
}

CKM_Rt_dms::CKM_Rt_dms(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Rt_dms::computeThValue()
{
    return SM.getCKM().computeRt()*sqrt(1.+SM.getCKM().getLambda()*SM.getCKM().getLambda()*(1.-2.*SM.getCKM().getRhoBar()));
}

CKM_Rts::CKM_Rts(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Rts::computeThValue()
{
    return SM.getCKM().computeRts();
}

CKM_Rb::CKM_Rb(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_Rb::computeThValue()
{
    return SM.getCKM().computeRb();
}

CKM_VtdoVts::CKM_VtdoVts(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_VtdoVts::computeThValue()
{
    return (SM.getCKM().getV_td().abs()/SM.getCKM().getV_ts().abs());
}

CKM_rho::CKM_rho(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_rho::computeThValue()
{
    return (SM.getCKM().getRho());
}

CKM_eta::CKM_eta(const StandardModel& SM_i) : ThObservable(SM_i) {}

double CKM_eta::computeThValue()
{
    return (SM.getCKM().getEta());
}
