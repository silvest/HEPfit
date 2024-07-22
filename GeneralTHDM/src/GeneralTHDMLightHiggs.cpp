/* 
 * Copyright (C) 2018 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMLightHiggs.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"



Gamma_h_GTHDM::Gamma_h_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_h_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gamma_h;
}


BR_h_bb_GTHDM::BR_h_bb_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_h_bb_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->GTHDM_BR_h_bb;
}


BR_h_gaga_GTHDM::BR_h_gaga_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_h_gaga_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->GTHDM_BR_h_gaga;
}


BR_h_tautau_GTHDM::BR_h_tautau_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_h_tautau_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->GTHDM_BR_h_tautau;
}



GTHDM_ggF_tth_htobb8::GTHDM_ggF_tth_htobb8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htobb8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth8*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

GTHDM_ggF_tth_htoWW8::GTHDM_ggF_tth_htoWW8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htoWW8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth8*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}


GTHDM_ggF_tth_htotautau8::GTHDM_ggF_tth_htotautau8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htotautau8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth8*(myGTHDM.getMyGTHDMCache()->rh_QlQlE+myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

GTHDM_ggF_tth_htoZZ8::GTHDM_ggF_tth_htoZZ8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htoZZ8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth8*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

GTHDM_ggF_tth_htogaga8::GTHDM_ggF_tth_htogaga8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htogaga8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth8*myGTHDM.getMyGTHDMCache()->rh_gaga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}


GTHDM_ggF_tth_htobb13::GTHDM_ggF_tth_htobb13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htobb13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth13*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}


GTHDM_ggF_tth_htoWW13::GTHDM_ggF_tth_htoWW13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htoWW13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth13*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_tth_htotautau13::GTHDM_ggF_tth_htotautau13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htotautau13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth13*(myGTHDM.getMyGTHDMCache()->rh_QlQlE+myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

GTHDM_ggF_tth_htoZZ13::GTHDM_ggF_tth_htoZZ13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htoZZ13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth13*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_tth_htogaga13::GTHDM_ggF_tth_htogaga13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_tth_htogaga13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->ggF_tth13*myGTHDM.getMyGTHDMCache()->rh_gaga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

GTHDM_VBF_Vh_htobb::GTHDM_VBF_Vh_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htobb::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htoWW::GTHDM_VBF_Vh_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htoWW::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htotautau::GTHDM_VBF_Vh_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htotautau::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*(myGTHDM.getMyGTHDMCache()->rh_QlQlE + myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htoZZ::GTHDM_VBF_Vh_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htoZZ::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htogaga::GTHDM_VBF_Vh_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htogaga::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*myGTHDM.getMyGTHDMCache()->rh_gaga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htogg::GTHDM_VBF_Vh_htogg(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htogg::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*(myGTHDM.getMyGTHDMCache()->rh_gg)/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_VBF_Vh_htocc::GTHDM_VBF_Vh_htocc(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_VBF_Vh_htocc::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->VBF_Vh*(myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_c*myGTHDM.getMyGTHDMCache()->beta_h_t))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}


GTHDM_ggF_htobb::GTHDM_ggF_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_htobb::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_gg)*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_htoWW::GTHDM_ggF_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_htoWW::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_gg)*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_htotautau::GTHDM_ggF_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_htotautau::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_gg)*(myGTHDM.getMyGTHDMCache()->rh_QlQlE + myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_htoZZ::GTHDM_ggF_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_htoZZ::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_gg)*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_ggF_htogaga::GTHDM_ggF_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_ggF_htogaga::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_gg)*myGTHDM.getMyGTHDMCache()->rh_gaga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_tth_htobb::GTHDM_tth_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_tth_htobb::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_t*myGTHDM.getMyGTHDMCache()->beta_h_t))*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_tth_htoWW::GTHDM_tth_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_tth_htoWW::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_t*myGTHDM.getMyGTHDMCache()->beta_h_t))*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_tth_htotautau::GTHDM_tth_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_tth_htotautau::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_t*myGTHDM.getMyGTHDMCache()->beta_h_t))*(myGTHDM.getMyGTHDMCache()->rh_QlQlE + myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_tth_htoZZ::GTHDM_tth_htoZZ(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_tth_htoZZ::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_t*myGTHDM.getMyGTHDMCache()->beta_h_t))*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_tth_htogaga::GTHDM_tth_htogaga(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_tth_htogaga::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->rh_QuQuE + myGTHDM.getMyGTHDMCache()->rh_QuQuO/(myGTHDM.getMyGTHDMCache()->beta_h_t*myGTHDM.getMyGTHDMCache()->beta_h_t))*myGTHDM.getMyGTHDMCache()->rh_gaga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_mu_htobb::GTHDM_mu_htobb(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_mu_htobb::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->pph13*(myGTHDM.getMyGTHDMCache()->rh_QdQdE + myGTHDM.getMyGTHDMCache()->rh_QdQdO/(myGTHDM.getMyGTHDMCache()->beta_h_b*myGTHDM.getMyGTHDMCache()->beta_h_b))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_mu_htoWW::GTHDM_mu_htoWW(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_mu_htoWW::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->pph13*myGTHDM.getMyGTHDMCache()->rh_VV/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_mu_htotautau::GTHDM_mu_htotautau(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_mu_htotautau::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->pph13*(myGTHDM.getMyGTHDMCache()->rh_QlQlE + myGTHDM.getMyGTHDMCache()->rh_QlQlO/(myGTHDM.getMyGTHDMCache()->beta_h_tau*myGTHDM.getMyGTHDMCache()->beta_h_tau))/myGTHDM.getMyGTHDMCache()->sumModBRs;
}



GTHDM_mu_htoZga::GTHDM_mu_htoZga(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double GTHDM_mu_htoZga::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->pph13*myGTHDM.getMyGTHDMCache()->rh_Zga/myGTHDM.getMyGTHDMCache()->sumModBRs;
}

yu1r_GTHDM::yu1r_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double yu1r_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->yu1R;
}

yd1r_GTHDM::yd1r_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double yd1r_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->yd1R;
}


yl1r_GTHDM::yl1r_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double yl1r_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->yl1R;
}


suR_GTHDM::suR_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double suR_GTHDM::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->su).real();
}


sdR_GTHDM::sdR_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double sdR_GTHDM::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->sd).real();
}


slR_GTHDM::slR_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double slR_GTHDM::computeThValue()
{
    return (myGTHDM.getMyGTHDMCache()->sl).real();
}


rh_gg_GTHDM::rh_gg_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rh_gg_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rh_gg;
}


rh_gaga_GTHDM::rh_gaga_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rh_gaga_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rh_gaga;
}


rh_Zga_GTHDM::rh_Zga_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rh_Zga_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rh_Zga;
}

