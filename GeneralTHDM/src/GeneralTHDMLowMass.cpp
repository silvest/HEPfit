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
/* CMS observables with phi_3, i.e A */
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
/* ATLAS observables with phi_3, i.e A */
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

Hobs_gg_phi3_tautau_ATLAS13_low::Hobs_gg_phi3_tautau_ATLAS13_low(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_gg_phi3_tautau_ATLAS13_low::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_gg_phi3_tautau_ATLAS13_low;
}

Hobs_pp_h_phi3phi3_gagagaga_ATLAS13::Hobs_pp_h_phi3phi3_gagagaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_gagagaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_gagagaga_ATLAS13;
}

Hobs_pp_h_phi3phi3_bbtautau_ATLAS13::Hobs_pp_h_phi3phi3_bbtautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_h_phi3phi3_bbtautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_h_phi3phi3_bbtautau_ATLAS13;
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
/* CMS observables with phi_2, i.e H */
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
/* ATLAS observables with phi_2, i.e H */
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


/*************************************/
/* LEP observables with phi_2, i.e H */
/*************************************/

Hobs_phi2Z_gagaZ_LEP209::Hobs_phi2Z_gagaZ_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi2Z_gagaZ_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi2Z_gagaZ_LEP209;
}

Hobs_phi2Z_bbZ_LEP209::Hobs_phi2Z_bbZ_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi2Z_bbZ_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi2Z_bbZ_LEP209;
}

Hobs_phi2Z_tautauZ_LEP209::Hobs_phi2Z_tautauZ_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi2Z_tautauZ_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi2Z_tautauZ_LEP209;
}


/**************************************/
/* LEP observables with phi_i + phi_j */
/**************************************/

Hobs_phi2phi3_bbbb_LEP209::Hobs_phi2phi3_bbbb_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi2phi3_bbbb_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi2phi3_bbbb_LEP209;
}

Hobs_phi2phi3_tautautautau_LEP209::Hobs_phi2phi3_tautautautau_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi2phi3_tautautautau_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi2phi3_tautautautau_LEP209;
}

Hobs_phi1phi3_bbbb_LEP209::Hobs_phi1phi3_bbbb_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi1phi3_bbbb_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi1phi3_bbbb_LEP209;
}

Hobs_phi1phi3_tautautautau_LEP209::Hobs_phi1phi3_tautautautau_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_phi1phi3_tautautautau_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_phi1phi3_tautautautau_LEP209;
}


/***************************/
/* LHC observables with Hp */
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

Hobs_t_Hpb_csb_ATLAS13::Hobs_t_Hpb_csb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_t_Hpb_csb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_t_Hpb_csb_ATLAS13;
}

/***************************/
/* LEP observables with Hp */
/***************************/

Hobs_HpHm_taunutaunu_LEP209::Hobs_HpHm_taunutaunu_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_HpHm_taunutaunu_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_HpHm_taunutaunu_LEP209;
}

Hobs_HpHm_qqqq_LEP209::Hobs_HpHm_qqqq_LEP209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_HpHm_qqqq_LEP209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_HpHm_qqqq_LEP209;
}

Hobs_HpHm_qqtaunu_OPAL209::Hobs_HpHm_qqtaunu_OPAL209(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_HpHm_qqtaunu_OPAL209::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_HpHm_qqtaunu_OPAL209;
}

Hobs_HpHm_qqtaunu_OPAL172::Hobs_HpHm_qqtaunu_OPAL172(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_HpHm_qqtaunu_OPAL172::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_HpHm_qqtaunu_OPAL172;
}


/*************************/
/* Invisible decay rates */
/*************************/

BR_h_inv_GTHDM::BR_h_inv_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_h_inv_GTHDM::computeThValue()
{
    double BR_h_inv_theo  = myGTHDM.getMyGTHDMCache()->Gamma_h_inv / myGTHDM.getMyGTHDMCache()->Gamma_h;
    double BR_h_inv_ATLAS = 0.107; // Combined ATLAS (7+8+13 TeV) upper limit at 95% CL, from 2301.10731

    return BR_h_inv_theo / BR_h_inv_ATLAS;
}

Gamma_Z_inv_GTHDM::Gamma_Z_inv_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_Z_inv_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gamma_Z_inv;
}

Gamma_W_inv_GTHDM::Gamma_W_inv_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_W_inv_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gamma_W_inv;
}
