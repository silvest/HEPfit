/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "lightHiggs.h"
#include "StandardModel.h"



THDM_BR_h_bb::THDM_BR_h_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double THDM_BR_h_bb::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THDM_BR_h_bb;
}



THDM_BR_h_gaga::THDM_BR_h_gaga(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double THDM_BR_h_gaga::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THDM_BR_h_gaga;
}



THDM_BR_h_tautau::THDM_BR_h_tautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double THDM_BR_h_tautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THDM_BR_h_tautau;
}



ggF_tth_htobb::ggF_tth_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htobb::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoWW::ggF_tth_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoWW::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htotautau::ggF_tth_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htotautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoZZ::ggF_tth_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoZZ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htogaga::ggF_tth_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htogaga::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
}



VBF_Vh_htobb::VBF_Vh_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double VBF_Vh_htobb::computeThValue()
{
    return myTHDM.getMyTHDMCache()->VBF_Vh*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



VBF_Vh_htoWW::VBF_Vh_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double VBF_Vh_htoWW::computeThValue()
{
    return myTHDM.getMyTHDMCache()->VBF_Vh*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



VBF_Vh_htotautau::VBF_Vh_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double VBF_Vh_htotautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->VBF_Vh*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



VBF_Vh_htoZZ::VBF_Vh_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double VBF_Vh_htoZZ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->VBF_Vh*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}

VBF_Vh_htogaga::VBF_Vh_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double VBF_Vh_htogaga::computeThValue()
{
    return myTHDM.getMyTHDMCache()->VBF_Vh*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
}

Gamma_h_THDM::Gamma_h_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Gamma_h_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->Gamma_h;
}

rh_gaga_THDM::rh_gaga_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double rh_gaga_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gaga;
}

rh_gg_THDM::rh_gg_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double rh_gg_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg;
}
