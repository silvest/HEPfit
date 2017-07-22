/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "chargedHiggs.h"
#include "StandardModel.h"



Hobs_pp_Hpm_taunu_ATLAS8::Hobs_pp_Hpm_taunu_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hpm_taunu_ATLAS8;
}

Robs_pp_Hpm_taunu_ATLAS8::Robs_pp_Hpm_taunu_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hpm_taunu_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hpm_taunu_ATLAS8;
}



Hobs_pp_Hp_taunu_CMS8::Hobs_pp_Hp_taunu_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hp_taunu_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hp_taunu_CMS8;
}

Robs_pp_Hp_taunu_CMS8::Robs_pp_Hp_taunu_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hp_taunu_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hp_taunu_CMS8;
}



Hobs_pp_Hpm_tb_ATLAS8::Hobs_pp_Hpm_tb_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hpm_tb_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hpm_tb_ATLAS8;
}

Robs_pp_Hpm_tb_ATLAS8::Robs_pp_Hpm_tb_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hpm_tb_ATLAS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hpm_tb_ATLAS8;
}



Hobs_pp_Hp_tb_CMS8::Hobs_pp_Hp_tb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hp_tb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hp_tb_CMS8;
}

Robs_pp_Hp_tb_CMS8::Robs_pp_Hp_tb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hp_tb_CMS8::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hp_tb_CMS8;
}



Hobs_pp_Hpm_taunu_ATLAS13::Hobs_pp_Hpm_taunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hpm_taunu_ATLAS13;
}

Robs_pp_Hpm_taunu_ATLAS13::Robs_pp_Hpm_taunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hpm_taunu_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hpm_taunu_ATLAS13;
}



Hobs_pp_Hpm_taunu_CMS13::Hobs_pp_Hpm_taunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hpm_taunu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hpm_taunu_CMS13;
}

Robs_pp_Hpm_taunu_CMS13::Robs_pp_Hpm_taunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hpm_taunu_CMS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hpm_taunu_CMS13;
}



Hobs_pp_Hp_tb_ATLAS13::Hobs_pp_Hp_tb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hp_tb_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hp_tb_ATLAS13;
}

Robs_pp_Hp_tb_ATLAS13::Robs_pp_Hp_tb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hp_tb_ATLAS13::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hp_tb_ATLAS13;
}



Hobs_pp_Hp_tb_ATLAS13_1::Hobs_pp_Hp_tb_ATLAS13_1(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hp_tb_ATLAS13_1::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hp_tb_ATLAS13_1;
}

Robs_pp_Hp_tb_ATLAS13_1::Robs_pp_Hp_tb_ATLAS13_1(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hp_tb_ATLAS13_1::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hp_tb_ATLAS13_1;
}



Hobs_pp_Hp_tb_ATLAS13_2::Hobs_pp_Hp_tb_ATLAS13_2(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Hobs_pp_Hp_tb_ATLAS13_2::computeThValue()
{
    return myTHDM.getMyTHDMCache()->THoEX_pp_Hp_tb_ATLAS13_2;
}

Robs_pp_Hp_tb_ATLAS13_2::Robs_pp_Hp_tb_ATLAS13_2(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Robs_pp_Hp_tb_ATLAS13_2::computeThValue()
{
    return myTHDM.getMyTHDMCache()->R_pp_Hp_tb_ATLAS13_2;
}



log10_pp_Hpm_taunu_TH8::log10_pp_Hpm_taunu_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_Hpm_taunu_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_Hpm_taunu_TH8);
}



log10_pp_Hp_tb_TH8::log10_pp_Hp_tb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_Hp_tb_TH8::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_Hp_tb_TH8);
}



log10_pp_Hpm_taunu_TH13::log10_pp_Hpm_taunu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_Hpm_taunu_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_Hpm_taunu_TH13);
}



log10_pp_Hp_tb_TH13::log10_pp_Hp_tb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double log10_pp_Hp_tb_TH13::computeThValue()
{
    return log10(myTHDM.getMyTHDMCache()->pp_Hp_tb_TH13);
}



Gamma_Hp_THDM::Gamma_Hp_THDM(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDM(static_cast<const THDM&> (SM_i))
{}

double Gamma_Hp_THDM::computeThValue()
{
    return myTHDM.getMyTHDMCache()->GammaHptot;
}
