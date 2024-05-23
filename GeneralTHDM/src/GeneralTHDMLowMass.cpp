/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMLowMass.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

/*************************************/
/* CMS Observables with phi_3, i.e A */
/*************************************/

Hobs_pp_h_phi3phi3_mumutautau_CMS13::Hobs_pp_h_phi3phi3_mumutautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_mumutautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_mumutautau_CMS13;
}

Hobs_pp_h_phi3phi3_bbtautau_CMS13::Hobs_pp_h_phi3phi3_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbtautau_CMS13;
}

Hobs_pp_h_phi3phi3_bbmumu_CMS13::Hobs_pp_h_phi3phi3_bbmumu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbmumu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbmumu_CMS13;
}

Hobs_pp_h_phi3Z_mumull_CMS13::Hobs_pp_h_phi3Z_mumull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3Z_mumull_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3Z_mumull_CMS13;
}

Hobs_pp_h_phi3phi3_mumumumu_CMS13::Hobs_pp_h_phi3phi3_mumumumu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_mumumumu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_mumumumu_CMS13;
}

Hobs_pp_h_phi3phi3_gagagaga_CMS13::Hobs_pp_h_phi3phi3_gagagaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_gagagaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_gagagaga_CMS13;
}

Hobs_pp_h_phi3phi3_tautautautau_CMS13::Hobs_pp_h_phi3phi3_tautautautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_tautautautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_tautautautau_CMS13;
}

Hobs_pp_bbphi3_bbtautau_CMS13::Hobs_pp_bbphi3_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_bbphi3_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_bbphi3_bbtautau_CMS13;
}

Hobs_pp_h_phi3phi3_tautautautau_CMS8::Hobs_pp_h_phi3phi3_tautautautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_tautautautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_tautautautau_CMS8;
}

Hobs_pp_h_phi3phi3_bbmumu_CMS8::Hobs_pp_h_phi3phi3_bbmumu_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbmumu_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbmumu_CMS8;
}

Hobs_pp_h_phi3phi3_mumutautau_CMS8::Hobs_pp_h_phi3phi3_mumutautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_mumutautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_mumutautau_CMS8;
}

Hobs_pp_bbphi3_bbtautau_CMS8::Hobs_pp_bbphi3_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_bbphi3_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_bbphi3_bbtautau_CMS8;
}

Hobs_pp_bbphi3_bbmumu_CMS8::Hobs_pp_bbphi3_bbmumu_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_bbphi3_bbmumu_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_bbphi3_bbmumu_CMS8;
}

/***************************************/
/* ATLAS Observables with phi_3, i.e A */
/***************************************/

Hobs_pp_h_phi3phi3_bbmumu_ATLAS13::Hobs_pp_h_phi3phi3_bbmumu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbmumu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbmumu_ATLAS13;
}

Hobs_gg_h_phi3phi3_mumumumu_ATLAS13::Hobs_gg_h_phi3phi3_mumumumu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_h_phi3phi3_mumumumu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_h_phi3phi3_mumumumu_ATLAS13;
}

Hobs_gg_h_phi3Z_mumull_ATLAS13::Hobs_gg_h_phi3Z_mumull_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_h_phi3Z_mumull_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_h_phi3Z_mumull_ATLAS13;
}

Hobs_Vh_h_phi3phi3_bbbb_ATLAS13::Hobs_Vh_h_phi3phi3_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_Vh_h_phi3phi3_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_Vh_h_phi3phi3_bbbb_ATLAS13;
}

Hobs_Zh_h_phi3phi3_bbbb_ATLAS13::Hobs_Zh_h_phi3phi3_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_Zh_h_phi3phi3_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_Zh_h_phi3phi3_bbbb_ATLAS13;
}

Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old::Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbmumu_ATLAS13_old;
}

Hobs_pp_h_phi3phi3_gagagg_ATLAS13::Hobs_pp_h_phi3phi3_gagagg_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_gagagg_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_gagagg_ATLAS13;
}

Hobs_pp_ttphi3_ttmumu_ATLAS13::Hobs_pp_ttphi3_ttmumu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_ttphi3_ttmumu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_ttphi3_ttmumu_ATLAS13;
}

Hobs_pp_h_phi3phi3_gagagaga_ATLAS8::Hobs_pp_h_phi3phi3_gagagaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_gagagaga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_gagagaga_ATLAS8;
}

Hobs_gg_h_phi3phi3_tautautautau_ATLAS8::Hobs_gg_h_phi3phi3_tautautautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_h_phi3phi3_tautautautau_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_h_phi3phi3_tautautautau_ATLAS8;
}


/*************************************/
/* CMS Observables with phi_2, i.e H */
/*************************************/

Hobs_pp_h_phi2Z_mumull_CMS13::Hobs_pp_h_phi2Z_mumull_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi2Z_mumull_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi2Z_mumull_CMS13;
}

Hobs_pp_h_phi2phi2_mumumumu_CMS13::Hobs_pp_h_phi2phi2_mumumumu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi2phi2_mumumumu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi2phi2_mumumumu_CMS13;
}

Hobs_pp_phi2_gaga_CMS13::Hobs_pp_phi2_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_CMS13;
}

Hobs_pp_phi2_gaga_CMS8::Hobs_pp_phi2_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_gaga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_CMS8;
}

/***************************************/
/* ATLAS Observables with phi_2, i.e H */
/***************************************/

Hobs_gg_h_phi2phi2_mumumumu_ATLAS13::Hobs_gg_h_phi2phi2_mumumumu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_h_phi2phi2_mumumumu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_h_phi2phi2_mumumumu_ATLAS13;
}

Hobs_gg_h_phi2Z_mumull_ATLAS13::Hobs_gg_h_phi2Z_mumull_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_h_phi2Z_mumull_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_h_phi2Z_mumull_ATLAS13;
}

Hobs_Vh_h_phi2phi2_bbbb_ATLAS13::Hobs_Vh_h_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_Vh_h_phi2phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_Vh_h_phi2phi2_bbbb_ATLAS13;
}

Hobs_Zh_h_phi2phi2_bbbb_ATLAS13::Hobs_Zh_h_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_Zh_h_phi2phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_Zh_h_phi2phi2_bbbb_ATLAS13;
}

Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old::Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi2phi2_bbmumu_ATLAS13_old;
}

Hobs_pp_h_phi2phi2_gagagg_ATLAS13::Hobs_pp_h_phi2phi2_gagagg_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi2phi2_gagagg_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi2phi2_gagagg_ATLAS13;
}

Hobs_pp_phi2_gaga_ATLAS13_low::Hobs_pp_phi2_gaga_ATLAS13_low(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_gaga_ATLAS13_low::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_ATLAS13_low;
}


/***************************/
/* CMS observables with Hp */
/***************************/

Hobs_t_Hpb_csb_CMS8::Hobs_t_Hpb_csb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_csb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_csb_CMS8;
}

Hobs_t_Hpb_taunub_CMS8::Hobs_t_Hpb_taunub_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_taunub_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_taunub_CMS8;
}

Hobs_t_Hpb_cbb_CMS8::Hobs_t_Hpb_cbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_cbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_cbb_CMS8;
}

Hobs_t_Hpb_WAb_Wmumub_CMS13::Hobs_t_Hpb_WAb_Wmumub_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_WAb_Wmumub_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_WAb_Wmumub_CMS13;
}

Hobs_t_Hpb_csb_CMS13::Hobs_t_Hpb_csb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_csb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_csb_CMS13;
}

Hobs_t_Hpb_taunub_ATLAS8::Hobs_t_Hpb_taunub_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_taunub_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_taunub_ATLAS8;
}

Hobs_t_Hpb_cbb_ATLAS13::Hobs_t_Hpb_cbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_cbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_cbb_ATLAS13;
}

Hobs_t_Hpb_WAb_Wmumub_ATLAS13::Hobs_t_Hpb_WAb_Wmumub_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_WAb_Wmumub_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_WAb_Wmumub_ATLAS13;
}
