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



ggF_tth_htobb8::ggF_tth_htobb8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htobb8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth8*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoWW8::ggF_tth_htoWW8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoWW8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth8*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htotautau8::ggF_tth_htotautau8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htotautau8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth8*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoZZ8::ggF_tth_htoZZ8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoZZ8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth8*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htogaga8::ggF_tth_htogaga8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htogaga8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth8*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htobb13::ggF_tth_htobb13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htobb13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth13*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoWW13::ggF_tth_htoWW13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoWW13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth13*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htotautau13::ggF_tth_htotautau13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htotautau13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth13*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htoZZ13::ggF_tth_htoZZ13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htoZZ13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth13*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_tth_htogaga13::ggF_tth_htogaga13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_tth_htogaga13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->ggF_tth13*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
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



ggF_htobb::ggF_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_htobb::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_htoWW::ggF_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_htoWW::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_htotautau::ggF_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_htotautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_htoZZ::ggF_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_htoZZ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



ggF_htogaga::ggF_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double ggF_htogaga::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_gg*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
}



tth_htobb::tth_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double tth_htobb::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_QuQu*myTHDM.getMyTHDMCache()->rh_QdQd/myTHDM.getMyTHDMCache()->sumModBRs;
}



tth_htoWW::tth_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double tth_htoWW::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_QuQu*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



tth_htotautau::tth_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double tth_htotautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_QuQu*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
}



tth_htoZZ::tth_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double tth_htoZZ::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_QuQu*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



tth_htogaga::tth_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double tth_htogaga::computeThValue()
{
    return myTHDM.getMyTHDMCache()->rh_QuQu*myTHDM.getMyTHDMCache()->rh_gaga/myTHDM.getMyTHDMCache()->sumModBRs;
}



mu_htoWW::mu_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double mu_htoWW::computeThValue()
{
    return myTHDM.getMyTHDMCache()->pph13*myTHDM.getMyTHDMCache()->rh_VV/myTHDM.getMyTHDMCache()->sumModBRs;
}



mu_htotautau::mu_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double mu_htotautau::computeThValue()
{
    return myTHDM.getMyTHDMCache()->pph13*myTHDM.getMyTHDMCache()->rh_ll/myTHDM.getMyTHDMCache()->sumModBRs;
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
