/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMHeavyHiggs.h"
#include "StandardModel.h"

Hobs_ggF_phi3_tautau_ATLAS8::Hobs_ggF_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi3_tautau_ATLAS8::computeThValue() 
{
    return myGTHDM.getMyGTHDMCache() -> THoEX_ggF_phi3_tautau_ATLAS8;
}


Robs_ggF_phi3_tautau_ATLAS8::Robs_ggF_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_ggF_phi3_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_tautau_ATLAS8;
}


Hobs_ggF_phi3_tautau_CMS8::Hobs_ggF_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_tautau_CMS8;
}

Robs_ggF_phi3_tautau_CMS8::Robs_ggF_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_ggF_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_tautau_CMS8;
}

Hobs_bbF_phi3_tautau_ATLAS8::Hobs_bbF_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bbF_phi3_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_tautau_ATLAS8;
}

Robs_bbF_phi3_tautau_ATLAS8::Robs_bbF_phi3_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_bbF_phi3_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_tautau_ATLAS8;
}

Hobs_bbF_phi3_tautau_CMS8::Hobs_bbF_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bbF_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_tautau_CMS8;
}

Robs_bbF_phi3_tautau_CMS8::Robs_bbF_phi3_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_bbF_phi3_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_tautau_CMS8;
}

//tau

Hobs_pp_phi3_gaga_ATLAS8::Hobs_pp_phi3_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_pp_phi3_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_gaga_ATLAS8;
}

Robs_pp_phi3_gaga_ATLAS8::Robs_pp_phi3_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_pp_phi3_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_gaga_ATLAS8;
}

Hobs_ggF_phi3_gaga_CMS8::Hobs_ggF_phi3_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi3_gaga_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_gaga_CMS8;
}

Robs_ggF_phi3_gaga_CMS8::Robs_ggF_phi3_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_gaga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_gaga_CMS8;
}

//gaga

Hobs_pp_phi3_Zga_llga_ATLAS8::Hobs_pp_phi3_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_llga_ATLAS8;
}

Robs_pp_phi3_Zga_llga_ATLAS8::Robs_pp_phi3_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_Zga_llga_ATLAS8;
}

Hobs_pp_phi3_Zga_llga_CMS8::Hobs_pp_phi3_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_llga_CMS8;
}

Robs_pp_phi3_Zga_llga_CMS8::Robs_pp_phi3_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_Zga_llga_CMS8;
}

Hobs_mu_pp_phi3_VV_CMS8::Hobs_mu_pp_phi3_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_mu_pp_phi3_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_mu_pp_phi3_VV_CMS8;
}

Robs_mu_pp_phi3_VV_CMS8::Robs_mu_pp_phi3_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_mu_pp_phi3_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_mu_pp_phi3_VV_CMS8;
}

Hobs_ggF_phi3_ZZ_ATLAS8::Hobs_ggF_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_ATLAS8;
}

Robs_ggF_phi3_ZZ_ATLAS8::Robs_ggF_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_ATLAS8;
}



Hobs_VBF_phi3_ZZ_ATLAS8::Hobs_VBF_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_ATLAS8;
}

Robs_VBF_phi3_ZZ_ATLAS8::Robs_VBF_phi3_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_ATLAS8;
}



Hobs_ggF_phi3_WW_ATLAS8::Hobs_ggF_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_WW_ATLAS8;
}

Robs_ggF_phi3_WW_ATLAS8::Robs_ggF_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_WW_ATLAS8;
}



Hobs_VBF_phi3_WW_ATLAS8::Hobs_VBF_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_WW_ATLAS8;
}

Robs_VBF_phi3_WW_ATLAS8::Robs_VBF_phi3_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_WW_ATLAS8;
}

//ggF, ATLAS 8

Hobs_ggF_phi3_phi1phi1_ATLAS8::Hobs_ggF_phi3_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi1_ATLAS8;
}

Robs_ggF_phi3_phi1phi1_ATLAS8::Robs_ggF_phi3_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi1_ATLAS8;
}

Hobs_ggF_phi3_phi1phi2_ATLAS8::Hobs_ggF_phi3_phi1phi2_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi2_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi2_ATLAS8;
}

Robs_ggF_phi3_phi1phi2_ATLAS8::Robs_ggF_phi3_phi1phi2_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi2_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi2_ATLAS8;
}


Hobs_ggF_phi3_phi2phi2_ATLAS8::Hobs_ggF_phi3_phi2phi2_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi2phi2_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi2phi2_ATLAS8;
}

Robs_ggF_phi3_phi2phi2_ATLAS8::Robs_ggF_phi3_phi2phi2_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi2phi2_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi2phi2_ATLAS8;
}


//pp CMS8

Hobs_pp_phi3_phi1phi1_CMS8::Hobs_pp_phi3_phi1phi1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_CMS8;
}

Robs_pp_phi3_phi1phi1_CMS8::Robs_pp_phi3_phi1phi1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_CMS8;
}

Hobs_pp_phi3_phi1phi2_CMS8::Hobs_pp_phi3_phi1phi2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_CMS8;
}

Robs_pp_phi3_phi1phi2_CMS8::Robs_pp_phi3_phi1phi2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_CMS8;
}


Hobs_pp_phi3_phi2phi2_CMS8::Hobs_pp_phi3_phi2phi2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_CMS8;
}

Robs_pp_phi3_phi2phi2_CMS8::Robs_pp_phi3_phi2phi2_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_CMS8;
}


//ggF bbtautau CMS8


Hobs_ggF_phi3_phi1phi1_bbtautau_CMS8::Hobs_ggF_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi1_bbtautau_CMS8;
}

Robs_ggF_phi3_phi1phi1_bbtautau_CMS8::Robs_ggF_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi1_bbtautau_CMS8;
}

Hobs_ggF_phi3_phi1phi2_bbtautau_CMS8::Hobs_ggF_phi3_phi1phi2_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi2_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi2_bbtautau_CMS8;
}

Robs_ggF_phi3_phi1phi2_bbtautau_CMS8::Robs_ggF_phi3_phi1phi2_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi2_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi2_bbtautau_CMS8;
}


Hobs_ggF_phi3_phi2phi2_bbtautau_CMS8::Hobs_ggF_phi3_phi2phi2_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi2phi2_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi2phi2_bbtautau_CMS8;
}

Robs_ggF_phi3_phi2phi2_bbtautau_CMS8::Robs_ggF_phi3_phi2phi2_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi2phi2_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi2phi2_bbtautau_CMS8;
}

//pp bbbb CMS8

Hobs_pp_phi3_phi1phi1_bbbb_CMS8::Hobs_pp_phi3_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_CMS8;
}

Robs_pp_phi3_phi1phi1_bbbb_CMS8::Robs_pp_phi3_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbbb_CMS8;
}

Hobs_pp_phi3_phi1phi2_bbbb_CMS8::Hobs_pp_phi3_phi1phi2_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbbb_CMS8;
}

Robs_pp_phi3_phi1phi2_bbbb_CMS8::Robs_pp_phi3_phi1phi2_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbbb_CMS8;
}


Hobs_pp_phi3_phi2phi2_bbbb_CMS8::Hobs_pp_phi3_phi2phi2_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbbb_CMS8;
}

Robs_pp_phi3_phi2phi2_bbbb_CMS8::Robs_pp_phi3_phi2phi2_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbbb_CMS8;
}

//pp gagabb CMS8

Hobs_pp_phi3_phi1phi1_gagabb_CMS8::Hobs_pp_phi3_phi1phi1_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_gagabb_CMS8;
}

Robs_pp_phi3_phi1phi1_gagabb_CMS8::Robs_pp_phi3_phi1phi1_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_gagabb_CMS8;
}

Hobs_pp_phi3_phi1phi2_gagabb_CMS8::Hobs_pp_phi3_phi1phi2_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_gagabb_CMS8;
}

Robs_pp_phi3_phi1phi2_gagabb_CMS8::Robs_pp_phi3_phi1phi2_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_gagabb_CMS8;
}


Hobs_pp_phi3_phi2phi2_gagabb_CMS8::Hobs_pp_phi3_phi2phi2_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_gagabb_CMS8;
}

Robs_pp_phi3_phi2phi2_gagabb_CMS8::Robs_pp_phi3_phi2phi2_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_gagabb_CMS8;
}


Hobs_ggF_phi3_tt_ATLAS8::Hobs_ggF_phi3_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_tt_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_tt_ATLAS8;
}

Robs_ggF_phi3_tt_ATLAS8::Robs_ggF_phi3_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_tt_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_tt_ATLAS8;
}



Hobs_bbF_phi3_bb_CMS8::Hobs_bbF_phi3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi3_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_bb_CMS8;
}

Robs_bbF_phi3_bb_CMS8::Robs_bbF_phi3_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi3_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_bb_CMS8;
}

//pp,  AZ_bbll CMS8

Hobs_pp_phi3_phi1Z_bbll_CMS8::Hobs_pp_phi3_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1Z_bbll_CMS8;
}

Robs_pp_phi3_phi1Z_bbll_CMS8::Robs_pp_phi3_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1Z_bbll_CMS8;
}

Hobs_pp_phi3_phi2Z_bbll_CMS8::Hobs_pp_phi3_phi2Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_bbll_CMS8;
}

Robs_pp_phi3_phi2Z_bbll_CMS8::Robs_pp_phi3_phi2Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2Z_bbll_CMS8;
}

// pp tautaull CMS8

Hobs_pp_phi3_phi1Z_tautaull_CMS8::Hobs_pp_phi3_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1Z_tautaull_CMS8;
}

Robs_pp_phi3_phi1Z_tautaull_CMS8::Robs_pp_phi3_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1Z_tautaull_CMS8;
}

Hobs_pp_phi3_phi2Z_tautaull_CMS8::Hobs_pp_phi3_phi2Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2Z_tautaull_CMS8;
}

Robs_pp_phi3_phi2Z_tautaull_CMS8::Robs_pp_phi3_phi2Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2Z_tautaull_CMS8;
}


Hobs_ttF_phi3_tt_ATLAS13::Hobs_ttF_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ttF_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ttF_phi3_tt_ATLAS13;
}

Robs_ttF_phi3_tt_ATLAS13::Robs_ttF_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ttF_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ttF_phi3_tt_ATLAS13;
}



Hobs_bbF_phi3_tt_ATLAS13::Hobs_bbF_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_tt_ATLAS13;
}

Robs_bbF_phi3_tt_ATLAS13::Robs_bbF_phi3_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi3_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_tt_ATLAS13;
}

Hobs_ggF_phi3_tautau_ATLAS13::Hobs_ggF_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_tautau_ATLAS13;
}

Robs_ggF_phi3_tautau_ATLAS13::Robs_ggF_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_tautau_ATLAS13;
}



Hobs_bbF_phi3_tautau_ATLAS13::Hobs_bbF_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_tautau_ATLAS13;
}

Robs_bbF_phi3_tautau_ATLAS13::Robs_bbF_phi3_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi3_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_tautau_ATLAS13;
}

Hobs_ggF_phi3_tautau_CMS13::Hobs_ggF_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_tautau_CMS13;
}

Robs_ggF_phi3_tautau_CMS13::Robs_ggF_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_tautau_CMS13;
}

Hobs_bbF_phi3_tautau_CMS13::Hobs_bbF_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi3_tautau_CMS13;
}

Robs_bbF_phi3_tautau_CMS13::Robs_bbF_phi3_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi3_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi3_tautau_CMS13;
}

Hobs_pp_phi3_gaga_ATLAS13::Hobs_pp_phi3_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_gaga_ATLAS13;
}

Robs_pp_phi3_gaga_ATLAS13::Robs_pp_phi3_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_gaga_ATLAS13;
}



Hobs_ggF_phi3_gaga_CMS13::Hobs_ggF_phi3_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_gaga_CMS13;
}

Robs_ggF_phi3_gaga_CMS13::Robs_ggF_phi3_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_gaga_CMS13;
}

Hobs_pp_phi3_Zga_llga_ATLAS13::Hobs_pp_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_ATLAS13;
}

Robs_pp_phi3_Zga_llga_ATLAS13::Robs_pp_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_Zga_ATLAS13;
}

Hobs_ggF_phi3_Zga_llga_ATLAS13::Hobs_ggF_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_Zga_llga_ATLAS13;
}

Robs_ggF_phi3_Zga_llga_ATLAS13::Robs_ggF_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_Zga_llga_ATLAS13;
}

Hobs_pp_phi3_Zga_llga_CMS13::Hobs_pp_phi3_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_llga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_llga_CMS13;
}

Robs_pp_phi3_Zga_llga_CMS13::Robs_pp_phi3_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_Zga_llga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_Zga_llga_CMS13;
}



Hobs_pp_phi3_Zga_qqga_CMS13::Hobs_pp_phi3_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_Zga_qqga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_Zga_qqga_CMS13;
}

Robs_pp_phi3_Zga_qqga_CMS13::Robs_pp_phi3_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_Zga_qqga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_Zga_qqga_CMS13;
}

Hobs_ggF_phi3_Zga_CMS13::Hobs_ggF_phi3_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_Zga_CMS13;
}

Robs_ggF_phi3_Zga_CMS13::Robs_ggF_phi3_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_Zga_CMS13;
}

Hobs_ggF_phi3_ZZ_llllnunu_ATLAS13::Hobs_ggF_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_llllnunu_ATLAS13;
}

Robs_ggF_phi3_ZZ_llllnunu_ATLAS13::Robs_ggF_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_llllnunu_ATLAS13;
}



Hobs_VBF_phi3_ZZ_llllnunu_ATLAS13::Hobs_VBF_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_llllnunu_ATLAS13;
}

Robs_VBF_phi3_ZZ_llllnunu_ATLAS13::Robs_VBF_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_llllnunu_ATLAS13;
}



Hobs_ggF_phi3_ZZ_llnunu_ATLAS13::Hobs_ggF_phi3_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_llnunu_ATLAS13;
}

Robs_ggF_phi3_ZZ_llnunu_ATLAS13::Robs_ggF_phi3_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_llnunu_ATLAS13;
}



Hobs_pp_phi3_ZZ_llnunu_CMS13::Hobs_pp_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_ZZ_llnunu_CMS13;
}

Robs_pp_phi3_ZZ_llnunu_CMS13::Robs_pp_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_ZZ_llnunu_CMS13;
}

Hobs_ggF_phi3_ZZ_llnunu_CMS13::Hobs_ggF_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_llnunu_CMS13;
}

Robs_ggF_phi3_ZZ_llnunu_CMS13::Robs_ggF_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_llnunu_CMS13;
}



Hobs_VBF_phi3_ZZ_llnunu_CMS13::Hobs_VBF_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_llnunu_CMS13;
}

Robs_VBF_phi3_ZZ_llnunu_CMS13::Robs_VBF_phi3_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_llnunu_CMS13;
}


Hobs_ggF_phi3_ZZ_llll_ATLAS13::Hobs_ggF_phi3_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_llll_ATLAS13;
}

Robs_ggF_phi3_ZZ_llll_ATLAS13::Robs_ggF_phi3_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_llll_ATLAS13;
}



Hobs_VBF_phi3_ZZ_llll_ATLAS13::Hobs_VBF_phi3_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_llll_ATLAS13;
}

Robs_VBF_phi3_ZZ_llll_ATLAS13::Robs_VBF_phi3_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_llll_ATLAS13;
}



Hobs_pp_phi3_ZZ_llll_CMS13::Hobs_pp_phi3_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_ZZ_llll_CMS13;
}

Robs_pp_phi3_ZZ_llll_CMS13::Robs_pp_phi3_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_ZZ_llll_CMS13;
}



Hobs_VBF_VH_phi3_ZZ_llll_CMS13::Hobs_VBF_VH_phi3_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_VH_phi3_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_VH_phi3_ZZ_llll_CMS13;
}

Robs_VBF_VH_phi3_ZZ_llll_CMS13::Robs_VBF_VH_phi3_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_VH_phi3_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_VH_phi3_ZZ_llll_CMS13;
}



Hobs_ggF_phi3_ZZ_qqllnunu_ATLAS13::Hobs_ggF_phi3_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_qqllnunu_ATLAS13;
}

Robs_ggF_phi3_ZZ_qqllnunu_ATLAS13::Robs_ggF_phi3_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_qqllnunu_ATLAS13;
}

Hobs_VBF_phi3_ZZ_qqllnunu_ATLAS13::Hobs_VBF_phi3_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_qqllnunu_ATLAS13;
}

Robs_VBF_phi3_ZZ_qqllnunu_ATLAS13::Robs_VBF_phi3_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_qqllnunu_ATLAS13;
}

Hobs_ggF_phi3_ZZ_llqq_ATLAS13::Hobs_ggF_phi3_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_llqq_ATLAS13;
}

Robs_ggF_phi3_ZZ_llqq_ATLAS13::Robs_ggF_phi3_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_llqq_ATLAS13;
}



Hobs_VBF_phi3_ZZ_llqq_ATLAS13::Hobs_VBF_phi3_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_ZZ_llqq_ATLAS13;
}

Robs_VBF_phi3_ZZ_llqq_ATLAS13::Robs_VBF_phi3_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_ZZ_llqq_ATLAS13;
}

Hobs_ggF_phi3_ZZ_nunuqq_ATLAS13::Hobs_ggF_phi3_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_ZZ_nunuqq_ATLAS13;
}

Robs_ggF_phi3_ZZ_nunuqq_ATLAS13::Robs_ggF_phi3_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_ZZ_nunuqq_ATLAS13;
}


Hobs_pp_phi3_ZZ_llqq_CMS13::Hobs_pp_phi3_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_ZZ_llqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_ZZ_llqq_CMS13;
}

Robs_pp_phi3_ZZ_llqq_CMS13::Robs_pp_phi3_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_ZZ_llqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_ZZ_llqq_CMS13;
}

Hobs_ggF_phi3_WW_lnuqq_ATLAS13::Hobs_ggF_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_WW_lnuqq_ATLAS13;
}

Robs_ggF_phi3_WW_lnuqq_ATLAS13::Robs_ggF_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_WW_lnuqq_ATLAS13;
}


Hobs_VBF_phi3_WW_lnuqq_ATLAS13::Hobs_VBF_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_WW_lnuqq_ATLAS13;
}

Robs_VBF_phi3_WW_lnuqq_ATLAS13::Robs_VBF_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_WW_lnuqq_ATLAS13;
}

Hobs_pp_phi3_VV_qqqq_ATLAS13::Hobs_pp_phi3_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_VV_qqqq_ATLAS13;
}

Robs_pp_phi3_VV_qqqq_ATLAS13::Robs_pp_phi3_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_VV_qqqq_ATLAS13;
}

Hobs_ggF_phi3_WW_enumunu_ATLAS13::Hobs_ggF_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_WW_enumunu_ATLAS13;
}

Robs_ggF_phi3_WW_enumunu_ATLAS13::Robs_ggF_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_WW_enumunu_ATLAS13;
}



Hobs_VBF_phi3_WW_enumunu_ATLAS13::Hobs_VBF_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi3_WW_enumunu_ATLAS13;
}

Robs_VBF_phi3_WW_enumunu_ATLAS13::Robs_VBF_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi3_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi3_WW_enumunu_ATLAS13;
}



Hobs_ggF_VBF_phi3_WW_lnulnu_CMS13::Hobs_ggF_VBF_phi3_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_VBF_phi3_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_VBF_phi3_WW_lnulnu_CMS13;
}

Robs_ggF_VBF_phi3_WW_lnulnu_CMS13::Robs_ggF_VBF_phi3_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_VBF_phi3_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_VBF_phi3_WW_lnulnu_CMS13;
}


Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13::Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbgaga_ATLAS13;
}

Robs_pp_phi3_phi1phi1_bbgaga_ATLAS13::Robs_pp_phi3_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbgaga_ATLAS13;
}


Hobs_pp_phi3_phi1phi1_bbgaga_CMS13::Hobs_pp_phi3_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbgaga_CMS13;
}

Robs_pp_phi3_phi1phi1_bbgaga_CMS13::Robs_pp_phi3_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbgaga_CMS13;
}



Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13::Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_ATLAS13;
}

Robs_pp_phi3_phi1phi1_bbbb_ATLAS13::Robs_pp_phi3_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbbb_ATLAS13;
}



Hobs_pp_phi3_phi1phi1_bbbb_CMS13::Hobs_pp_phi3_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbbb_CMS13;
}

Robs_pp_phi3_phi1phi1_bbbb_CMS13::Robs_pp_phi3_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbbb_CMS13;
}



Hobs_ggF_phi3_phi1phi1_bbbb_CMS13::Hobs_ggF_phi3_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi1_bbbb_CMS13;
}

Robs_ggF_phi3_phi1phi1_bbbb_CMS13::Robs_ggF_phi3_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi1_bbbb_CMS13;
}



Hobs_ggF_phi3_phi1phi1_gagaWW_ATLAS13::Hobs_ggF_phi3_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi1_gagaWW_ATLAS13;
}

Robs_ggF_phi3_phi1phi1_gagaWW_ATLAS13::Robs_ggF_phi3_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi1_gagaWW_ATLAS13;
}



Hobs_pp_phi3_phi1phi1_bbtautau_CMS13::Hobs_pp_phi3_phi1phi1_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau_CMS13;
}

Robs_pp_phi3_phi1phi1_bbtautau_CMS13::Robs_pp_phi3_phi1phi1_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbtautau_CMS13;
}



Hobs_pp_phi3_phi1phi1_bbtautau1_CMS13::Hobs_pp_phi3_phi1phi1_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbtautau1_CMS13;
}

Robs_pp_phi3_phi1phi1_bbtautau1_CMS13::Robs_pp_phi3_phi1phi1_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbtautau1_CMS13;
}



Hobs_pp_phi3_phi1phi1_bblnulnu_CMS13::Hobs_pp_phi3_phi1phi1_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bblnulnu_CMS13;
}

Robs_pp_phi3_phi1phi1_bblnulnu_CMS13::Robs_pp_phi3_phi1phi1_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bblnulnu_CMS13;
}



Hobs_pp_phi3_phi1phi1_bbVV_CMS13::Hobs_pp_phi3_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi1_bbVV_CMS13;
}

Robs_pp_phi3_phi1phi1_bbVV_CMS13::Robs_pp_phi3_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi1_bbVV_CMS13;
}



/////


Hobs_pp_phi3_phi2phi2_bbgaga_ATLAS13::Hobs_pp_phi3_phi2phi2_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbgaga_ATLAS13;
}

Robs_pp_phi3_phi2phi2_bbgaga_ATLAS13::Robs_pp_phi3_phi2phi2_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbgaga_ATLAS13;
}



Hobs_pp_phi3_phi2phi2_bbgaga_CMS13::Hobs_pp_phi3_phi2phi2_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbgaga_CMS13;
}

Robs_pp_phi3_phi2phi2_bbgaga_CMS13::Robs_pp_phi3_phi2phi2_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbgaga_CMS13;
}



Hobs_pp_phi3_phi2phi2_bbbb_ATLAS13::Hobs_pp_phi3_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbbb_ATLAS13;
}

Robs_pp_phi3_phi2phi2_bbbb_ATLAS13::Robs_pp_phi3_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbbb_ATLAS13;
}



Hobs_pp_phi3_phi2phi2_bbbb_CMS13::Hobs_pp_phi3_phi2phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbbb_CMS13;
}

Robs_pp_phi3_phi2phi2_bbbb_CMS13::Robs_pp_phi3_phi2phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbbb_CMS13;
}



Hobs_ggF_phi3_phi2phi2_bbbb_CMS13::Hobs_ggF_phi3_phi2phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi2phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi2phi2_bbbb_CMS13;
}




Hobs_pp_phi3_phi1phi2_bbgaga_ATLAS13::Hobs_pp_phi3_phi1phi2_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbgaga_ATLAS13;
}

Robs_pp_phi3_phi1phi2_bbgaga_ATLAS13::Robs_pp_phi3_phi1phi2_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbgaga_ATLAS13;
}



Hobs_pp_phi3_phi1phi2_bbgaga_CMS13::Hobs_pp_phi3_phi1phi2_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbgaga_CMS13;
}

Robs_pp_phi3_phi1phi2_bbgaga_CMS13::Robs_pp_phi3_phi1phi2_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbgaga_CMS13;
}



Hobs_pp_phi3_phi1phi2_bbbb_ATLAS13::Hobs_pp_phi3_phi1phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbbb_ATLAS13;
}

Robs_pp_phi3_phi1phi2_bbbb_ATLAS13::Robs_pp_phi3_phi1phi2_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbbb_ATLAS13;
}



Hobs_pp_phi3_phi1phi2_bbbb_CMS13::Hobs_pp_phi3_phi1phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbbb_CMS13;
}

Robs_pp_phi3_phi1phi2_bbbb_CMS13::Robs_pp_phi3_phi1phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbbb_CMS13;
}



Hobs_ggF_phi3_phi1phi2_bbbb_CMS13::Hobs_ggF_phi3_phi1phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi2_bbbb_CMS13;
}

Robs_ggF_phi3_phi2phi2_bbbb_CMS13::Robs_ggF_phi3_phi2phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi2phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi2phi2_bbbb_CMS13;
}


Robs_ggF_phi3_phi1phi2_bbbb_CMS13::Robs_ggF_phi3_phi1phi2_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi2_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi2_bbbb_CMS13;
}




Hobs_ggF_phi3_phi2phi2_gagaWW_ATLAS13::Hobs_ggF_phi3_phi2phi2_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi2phi2_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi2phi2_gagaWW_ATLAS13;
}

Robs_ggF_phi3_phi2phi2_gagaWW_ATLAS13::Robs_ggF_phi3_phi2phi2_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi2phi2_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi2phi2_gagaWW_ATLAS13;
}



Hobs_pp_phi3_phi2phi2_bbtautau_CMS13::Hobs_pp_phi3_phi2phi2_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbtautau_CMS13;
}

Robs_pp_phi3_phi2phi2_bbtautau_CMS13::Robs_pp_phi3_phi2phi2_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbtautau_CMS13;
}



Hobs_pp_phi3_phi2phi2_bbtautau1_CMS13::Hobs_pp_phi3_phi2phi2_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbtautau1_CMS13;
}

Robs_pp_phi3_phi2phi2_bbtautau1_CMS13::Robs_pp_phi3_phi2phi2_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbtautau1_CMS13;
}



Hobs_pp_phi3_phi2phi2_bblnulnu_CMS13::Hobs_pp_phi3_phi2phi2_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bblnulnu_CMS13;
}

Robs_pp_phi3_phi2phi2_bblnulnu_CMS13::Robs_pp_phi3_phi2phi2_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bblnulnu_CMS13;
}



Hobs_pp_phi3_phi2phi2_bbVV_CMS13::Hobs_pp_phi3_phi2phi2_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi2phi2_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi2phi2_bbVV_CMS13;
}

Robs_pp_phi3_phi2phi2_bbVV_CMS13::Robs_pp_phi3_phi2phi2_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi2phi2_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi2phi2_bbVV_CMS13;
}




Hobs_ggF_phi3_phi1phi2_gagaWW_ATLAS13::Hobs_ggF_phi3_phi1phi2_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi3_phi1phi2_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi3_phi1phi2_gagaWW_ATLAS13;
}

Robs_ggF_phi3_phi1phi2_gagaWW_ATLAS13::Robs_ggF_phi3_phi1phi2_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi3_phi1phi2_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi3_phi1phi2_gagaWW_ATLAS13;
}



Hobs_pp_phi3_phi1phi2_bbtautau_CMS13::Hobs_pp_phi3_phi1phi2_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbtautau_CMS13;
}

Robs_pp_phi3_phi1phi2_bbtautau_CMS13::Robs_pp_phi3_phi1phi2_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbtautau_CMS13;
}



Hobs_pp_phi3_phi1phi2_bbtautau1_CMS13::Hobs_pp_phi3_phi1phi2_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbtautau1_CMS13;
}

Robs_pp_phi3_phi1phi2_bbtautau1_CMS13::Robs_pp_phi3_phi1phi2_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbtautau1_CMS13;
}



Hobs_pp_phi3_phi1phi2_bblnulnu_CMS13::Hobs_pp_phi3_phi1phi2_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bblnulnu_CMS13;
}

Robs_pp_phi3_phi1phi2_bblnulnu_CMS13::Robs_pp_phi3_phi1phi2_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bblnulnu_CMS13;
}



Hobs_pp_phi3_phi1phi2_bbVV_CMS13::Hobs_pp_phi3_phi1phi2_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_phi1phi2_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_phi1phi2_bbVV_CMS13;
}

Robs_pp_phi3_phi1phi2_bbVV_CMS13::Robs_pp_phi3_phi1phi2_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_phi1phi2_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_phi1phi2_bbVV_CMS13;
}


Hobs_pp_phi3_bb_CMS13::Hobs_pp_phi3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi3_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi3_bb_CMS13;
}

Robs_pp_phi3_bb_CMS13::Robs_pp_phi3_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi3_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi3_bb_CMS13;
}



log10_ggF_phi3_tautau_TH8::log10_ggF_phi3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_tautau_TH8);
}



log10_bbF_phi3_tautau_TH8::log10_bbF_phi3_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi3_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi3_tautau_TH8);
}



log10_pp_phi3_gaga_TH8::log10_pp_phi3_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_gaga_TH8);
}



log10_ggF_phi3_gaga_TH8::log10_ggF_phi3_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_gaga_TH8);
}



log10_pp_phi3_Zga_llga_TH8::log10_pp_phi3_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_Zga_llga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_Zga_llga_TH8);
}



log10_mu_pp_phi3_VV_TH8::log10_mu_pp_phi3_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_mu_pp_phi3_VV_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->mu_pp_phi3_VV_TH8);
}



log10_ggF_phi3_ZZ_TH8::log10_ggF_phi3_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_ZZ_TH8);
}



log10_VBF_phi3_ZZ_TH8::log10_VBF_phi3_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi3_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi3_ZZ_TH8);
}



log10_ggF_phi3_WW_TH8::log10_ggF_phi3_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_WW_TH8);
}



log10_VBF_phi3_WW_TH8::log10_VBF_phi3_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi3_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi3_WW_TH8);
}


//ggF and pp TH8

log10_ggF_phi3_phi1phi1_TH8::log10_ggF_phi3_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi1_TH8);
}



log10_pp_phi3_phi1phi1_TH8::log10_pp_phi3_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_TH8);
}


log10_ggF_phi3_phi2phi2_TH8::log10_ggF_phi3_phi2phi2_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi2phi2_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi2phi2_TH8);
}


log10_ggF_phi3_phi1phi2_TH8::log10_ggF_phi3_phi1phi2_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi2_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi2_TH8);
}



log10_pp_phi3_phi2phi2_TH8::log10_pp_phi3_phi2phi2_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_TH8);
}


log10_pp_phi3_phi1phi2_TH8::log10_pp_phi3_phi1phi2_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_TH8);
}


// ggF  bbtautau

log10_ggF_phi3_phi1phi1_bbtautau_TH8::log10_ggF_phi3_phi1phi1_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi1_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi1_bbtautau_TH8);
}

log10_ggF_phi3_phi2phi2_bbtautau_TH8::log10_ggF_phi3_phi2phi2_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi2phi2_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi2phi2_bbtautau_TH8);
}

log10_ggF_phi3_phi1phi2_bbtautau_TH8::log10_ggF_phi3_phi1phi2_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi2_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi2_bbtautau_TH8);
}


// pp  bbbb

log10_pp_phi3_phi1phi1_bbbb_TH8::log10_pp_phi3_phi1phi1_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbbb_TH8);
}

log10_pp_phi3_phi1phi2_bbbb_TH8::log10_pp_phi3_phi1phi2_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_bbbb_TH8);
}

log10_pp_phi3_phi2phi2_bbbb_TH8::log10_pp_phi3_phi2phi2_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_bbbb_TH8);
}

//pp gagabb


log10_pp_phi3_phi1phi1_gagabb_TH8::log10_pp_phi3_phi1phi1_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_gagabb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_gagabb_TH8);
}

log10_pp_phi3_phi1phi2_gagabb_TH8::log10_pp_phi3_phi1phi2_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_gagabb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_gagabb_TH8);
}


log10_pp_phi3_phi2phi2_gagabb_TH8::log10_pp_phi3_phi2phi2_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_gagabb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_gagabb_TH8);
}


//


log10_ggF_phi3_tt_TH8::log10_ggF_phi3_tt_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_tt_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_tt_TH8);
}



log10_bbF_phi3_bb_TH8::log10_bbF_phi3_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi3_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi3_bb_TH8);
}


//pp AZbbll TH8

log10_pp_phi3_phi1Z_bbll_TH8::log10_pp_phi3_phi1Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1Z_bbll_TH8);
}

log10_pp_phi3_phi2Z_bbll_TH8::log10_pp_phi3_phi2Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2Z_bbll_TH8);
}


//pp AZ tautaull

log10_pp_phi3_phi1Z_tautaull_TH8::log10_pp_phi3_phi1Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1Z_tautaull_TH8);
}

log10_pp_phi3_phi2Z_tautaull_TH8::log10_pp_phi3_phi2Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2Z_tautaull_TH8);
}


//


log10_ggF_phi3_tautau_TH13::log10_ggF_phi3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_tautau_TH13);
}



log10_bbF_phi3_tautau_TH13::log10_bbF_phi3_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi3_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi3_tautau_TH13);
}



log10_pp_phi3_gaga_TH13::log10_pp_phi3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_gaga_TH13);
}



log10_ggF_phi3_gaga_TH13::log10_ggF_phi3_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_gaga_TH13);
}



log10_pp_phi3_Zga_TH13::log10_pp_phi3_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_Zga_TH13);
}



log10_ggF_phi3_Zga_TH13::log10_ggF_phi3_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_Zga_TH13);
}



log10_pp_phi3_ZZ_TH13::log10_pp_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_ZZ_TH13);
}

log10_ggF_phi3_ZZ_TH13::log10_ggF_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_ZZ_TH13);
}



log10_VBF_phi3_ZZ_TH13::log10_VBF_phi3_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi3_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi3_ZZ_TH13);
}



log10_ggF_phi3_ZZ_llll_TH13::log10_ggF_phi3_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_ZZ_llll_TH13);
}



log10_VBF_phi3_ZZ_llll_TH13::log10_VBF_phi3_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi3_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi3_ZZ_llll_TH13);
}



log10_pp_phi3_ZZ_llll_TH13::log10_pp_phi3_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_ZZ_llll_TH13);
}



log10_VBF_VH_phi3_ZZ_llll_TH13::log10_VBF_VH_phi3_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_VH_phi3_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_VH_phi3_ZZ_llll_TH13);
}



log10_ggF_phi3_WW_TH13::log10_ggF_phi3_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_WW_TH13);
}



log10_VBF_phi3_WW_TH13::log10_VBF_phi3_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi3_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi3_WW_TH13);
}



log10_ggF_VBF_phi3_WW_lnulnu_TH13::log10_ggF_VBF_phi3_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_VBF_phi3_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_VBF_phi3_WW_lnulnu_TH13);
}


// ggF, phi1phi1, TH13

log10_ggF_phi3_phi1phi1_TH13::log10_ggF_phi3_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi1_TH13);
}

log10_ggF_phi3_phi2phi2_TH13::log10_ggF_phi3_phi2phi2_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi2phi2_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi2phi2_TH13);
}

log10_ggF_phi3_phi1phi2_TH13::log10_ggF_phi3_phi1phi2_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi2_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi2_TH13);
}

//pp phi1phi1 TH13

log10_pp_phi3_phi1phi1_TH13::log10_pp_phi3_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_TH13);
}

log10_pp_phi3_phi2phi2_TH13::log10_pp_phi3_phi2phi2_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_TH13);
}


log10_pp_phi3_phi1phi2_TH13::log10_pp_phi3_phi1phi2_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_TH13);
}



//pp bbbb

log10_pp_phi3_phi1phi1_bbbb_TH13::log10_pp_phi3_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbbb_TH13);
}

log10_pp_phi3_phi2phi2_bbbb_TH13::log10_pp_phi3_phi2phi2_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_bbbb_TH13);
}

log10_pp_phi3_phi1phi2_bbbb_TH13::log10_pp_phi3_phi1phi2_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_bbbb_TH13);
}



//ggF bbbb

log10_ggF_phi3_phi1phi1_bbbb_TH13::log10_ggF_phi3_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi1_bbbb_TH13);
}


log10_ggF_phi3_phi2phi2_bbbb_TH13::log10_ggF_phi3_phi2phi2_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi2phi2_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi2phi2_bbbb_TH13);
}



log10_ggF_phi3_phi1phi2_bbbb_TH13::log10_ggF_phi3_phi1phi2_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi3_phi1phi2_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi3_phi1phi2_bbbb_TH13);
}


//pp gagabb

log10_pp_phi3_phi1phi1_gagabb_TH13::log10_pp_phi3_phi1phi1_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_gagabb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_gagabb_TH13);
}

log10_pp_phi3_phi2phi2_gagabb_TH13::log10_pp_phi3_phi2phi2_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_gagabb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_gagabb_TH13);
}

log10_pp_phi3_phi1phi2_gagabb_TH13::log10_pp_phi3_phi1phi2_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_gagabb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_gagabb_TH13);
}





//pp bbtautau

log10_pp_phi3_phi1phi1_bbtautau_TH13::log10_pp_phi3_phi1phi1_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbtautau_TH13);
}

log10_pp_phi3_phi2phi2_bbtautau_TH13::log10_pp_phi3_phi2phi2_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_bbtautau_TH13);
}

log10_pp_phi3_phi1phi2_bbtautau_TH13::log10_pp_phi3_phi1phi2_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_bbtautau_TH13);
}

//pp bblnulnu

log10_pp_phi3_phi1phi1_bblnulnu_TH13::log10_pp_phi3_phi1phi1_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bblnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bblnulnu_TH13);
}


log10_pp_phi3_phi2phi2_bblnulnu_TH13::log10_pp_phi3_phi2phi2_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_bblnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_bblnulnu_TH13);
}

log10_pp_phi3_phi1phi2_bblnulnu_TH13::log10_pp_phi3_phi1phi2_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_bblnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_bblnulnu_TH13);
}


//pp bbVV

log10_pp_phi3_phi1phi1_bbVV_TH13::log10_pp_phi3_phi1phi1_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi1_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi1_bbVV_TH13);
}

log10_pp_phi3_phi2phi2_bbVV_TH13::log10_pp_phi3_phi2phi2_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi2phi2_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi2phi2_bbVV_TH13);
}


log10_pp_phi3_phi1phi2_bbVV_TH13::log10_pp_phi3_phi1phi2_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_phi1phi2_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_phi1phi2_bbVV_TH13);
}



//

log10_tt_phi3_tt_TH13::log10_tt_phi3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_tt_phi3_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ttF_phi3_tt_TH13);
}



log10_bb_phi3_tt_TH13::log10_bb_phi3_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi3_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi3_tt_TH13);
}



log10_pp_phi3_bb_TH13::log10_pp_phi3_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi3_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi3_bb_TH13);
}




Gamma_phi3_GTHDM::Gamma_phi3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_phi3_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gammaphi3tot;
}



//rHH_gaga_THDM::rHH_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
//{}
//
//double rHH_gaga_THDM::computeThValue()
//{
//    return myGTHDM.getMyGTHDMCache()->rHH_gaga;
//}



rphi3_ggE_GTHDM::rphi3_ggE_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi3_ggE_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi3_ggE;
}

rphi3_ggO_GTHDM::rphi3_ggO_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi3_ggO_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi3_ggO;
}



// phi3 -> phi1phi1 and phi3 -> phi2phi2

BR_phi3_phi1phi1_GTHDM::BR_phi3_phi1phi1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi1phi1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi1phi1;
}


BR_phi3_phi2phi2_GTHDM::BR_phi3_phi2phi2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi2phi2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi2phi2;
}



//phi3 -> H+H-


BR_phi3_HpHm_GTHDM::BR_phi3_HpHm_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_HpHm_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3toHpHm;
}



BR_phi3_phi1Z_GTHDM::BR_phi3_phi1Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi1Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi1Z;
}

BR_phi3_phi2Z_GTHDM::BR_phi3_phi2Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_phi2Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3tophi2Z;
}



BR_phi3_HpW_GTHDM::BR_phi3_HpW_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi3_HpW_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi3toHpW;
}

Hobs_ggF_phi2_tautau_ATLAS8::Hobs_ggF_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi2_tautau_ATLAS8::computeThValue() 
{
    return myGTHDM.getMyGTHDMCache() -> THoEX_ggF_phi2_tautau_ATLAS8;
}


Robs_ggF_phi2_tautau_ATLAS8::Robs_ggF_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_ggF_phi2_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_tautau_ATLAS8;
}


Hobs_ggF_phi2_tautau_CMS8::Hobs_ggF_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_tautau_CMS8;
}

Robs_ggF_phi2_tautau_CMS8::Robs_ggF_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_ggF_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_tautau_CMS8;
}

Hobs_bbF_phi2_tautau_ATLAS8::Hobs_bbF_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bbF_phi2_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_tautau_ATLAS8;
}

Robs_bbF_phi2_tautau_ATLAS8::Robs_bbF_phi2_tautau_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_bbF_phi2_tautau_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_tautau_ATLAS8;
}

Hobs_bbF_phi2_tautau_CMS8::Hobs_bbF_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_bbF_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_tautau_CMS8;
}

Robs_bbF_phi2_tautau_CMS8::Robs_bbF_phi2_tautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_bbF_phi2_tautau_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_tautau_CMS8;
}

//tau

Hobs_pp_phi2_gaga_ATLAS8::Hobs_pp_phi2_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_pp_phi2_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_ATLAS8;
}

Robs_pp_phi2_gaga_ATLAS8::Robs_pp_phi2_gaga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Robs_pp_phi2_gaga_ATLAS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_gaga_ATLAS8;
}

Hobs_ggF_phi2_gaga_CMS8::Hobs_ggF_phi2_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i), myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) {
}

double Hobs_ggF_phi2_gaga_CMS8::computeThValue() {
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_gaga_CMS8;
}

Robs_ggF_phi2_gaga_CMS8::Robs_ggF_phi2_gaga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_gaga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_gaga_CMS8;
}

//gaga


Hobs_pp_phi2_Zga_llga_ATLAS8::Hobs_pp_phi2_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_llga_ATLAS8;
}

Robs_pp_phi2_Zga_llga_ATLAS8::Robs_pp_phi2_Zga_llga_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_Zga_llga_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_Zga_llga_ATLAS8;
}



Hobs_pp_phi2_Zga_llga_CMS8::Hobs_pp_phi2_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_llga_CMS8;
}

Robs_pp_phi2_Zga_llga_CMS8::Robs_pp_phi2_Zga_llga_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_Zga_llga_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_Zga_llga_CMS8;
}



Hobs_mu_pp_phi2_VV_CMS8::Hobs_mu_pp_phi2_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_mu_pp_phi2_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_mu_pp_phi2_VV_CMS8;
}

Robs_mu_pp_phi2_VV_CMS8::Robs_mu_pp_phi2_VV_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_mu_pp_phi2_VV_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_mu_pp_phi2_VV_CMS8;
}

Hobs_ggF_phi2_ZZ_ATLAS8::Hobs_ggF_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_ATLAS8;
}

Robs_ggF_phi2_ZZ_ATLAS8::Robs_ggF_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_ATLAS8;
}

Hobs_VBF_phi2_ZZ_ATLAS8::Hobs_VBF_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_ATLAS8;
}

Robs_VBF_phi2_ZZ_ATLAS8::Robs_VBF_phi2_ZZ_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_ATLAS8;
}

Hobs_ggF_phi2_WW_ATLAS8::Hobs_ggF_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_WW_ATLAS8;
}

Robs_ggF_phi2_WW_ATLAS8::Robs_ggF_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_WW_ATLAS8;
}



Hobs_VBF_phi2_WW_ATLAS8::Hobs_VBF_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_WW_ATLAS8;
}

Robs_VBF_phi2_WW_ATLAS8::Robs_VBF_phi2_WW_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_WW_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_WW_ATLAS8;
}

//ggF, ATLAS 8

Hobs_ggF_phi2_phi1phi1_ATLAS8::Hobs_ggF_phi2_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_phi1phi1_ATLAS8;
}

Robs_ggF_phi2_phi1phi1_ATLAS8::Robs_ggF_phi2_phi1phi1_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_phi1phi1_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_phi1phi1_ATLAS8;
}

//pp CMS8

Hobs_pp_phi2_phi1phi1_CMS8::Hobs_pp_phi2_phi1phi1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_CMS8;
}

Robs_pp_phi2_phi1phi1_CMS8::Robs_pp_phi2_phi1phi1_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_CMS8;
}
//ggF bbtautau CMS8


Hobs_ggF_phi2_phi1phi1_bbtautau_CMS8::Hobs_ggF_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_phi1phi1_bbtautau_CMS8;
}

Robs_ggF_phi2_phi1phi1_bbtautau_CMS8::Robs_ggF_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_phi1phi1_bbtautau_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_phi1phi1_bbtautau_CMS8;
}

//pp bbbb CMS8

Hobs_pp_phi2_phi1phi1_bbbb_CMS8::Hobs_pp_phi2_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_CMS8;
}

Robs_pp_phi2_phi1phi1_bbbb_CMS8::Robs_pp_phi2_phi1phi1_bbbb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbbb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbbb_CMS8;
}


//pp gagabb CMS8

Hobs_pp_phi2_phi1phi1_gagabb_CMS8::Hobs_pp_phi2_phi1phi1_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_gagabb_CMS8;
}

Robs_pp_phi2_phi1phi1_gagabb_CMS8::Robs_pp_phi2_phi1phi1_gagabb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_gagabb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_gagabb_CMS8;
}

//

Hobs_ggF_phi2_tt_ATLAS8::Hobs_ggF_phi2_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_tt_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_tt_ATLAS8;
}

Robs_ggF_phi2_tt_ATLAS8::Robs_ggF_phi2_tt_ATLAS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_tt_ATLAS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_tt_ATLAS8;
}



Hobs_bbF_phi2_bb_CMS8::Hobs_bbF_phi2_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi2_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_bb_CMS8;
}

Robs_bbF_phi2_bb_CMS8::Robs_bbF_phi2_bb_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi2_bb_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_bb_CMS8;
}

//pp,  AZ_bbll CMS8

Hobs_pp_phi2_phi1Z_bbll_CMS8::Hobs_pp_phi2_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1Z_bbll_CMS8;
}

Robs_pp_phi2_phi1Z_bbll_CMS8::Robs_pp_phi2_phi1Z_bbll_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1Z_bbll_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1Z_bbll_CMS8;
}

// pp tautaull CMS8

Hobs_pp_phi2_phi1Z_tautaull_CMS8::Hobs_pp_phi2_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1Z_tautaull_CMS8;
}

Robs_pp_phi2_phi1Z_tautaull_CMS8::Robs_pp_phi2_phi1Z_tautaull_CMS8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1Z_tautaull_CMS8::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1Z_tautaull_CMS8;
}

//

Hobs_ttF_phi2_tt_ATLAS13::Hobs_ttF_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ttF_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ttF_phi2_tt_ATLAS13;
}

Robs_ttF_phi2_tt_ATLAS13::Robs_ttF_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ttF_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ttF_phi2_tt_ATLAS13;
}



Hobs_bbF_phi2_tt_ATLAS13::Hobs_bbF_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_tt_ATLAS13;
}

Robs_bbF_phi2_tt_ATLAS13::Robs_bbF_phi2_tt_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi2_tt_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_tt_ATLAS13;
}



Hobs_ggF_phi2_tautau_ATLAS13::Hobs_ggF_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_tautau_ATLAS13;
}

Robs_ggF_phi2_tautau_ATLAS13::Robs_ggF_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_tautau_ATLAS13;
}



Hobs_bbF_phi2_tautau_ATLAS13::Hobs_bbF_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_tautau_ATLAS13;
}

Robs_bbF_phi2_tautau_ATLAS13::Robs_bbF_phi2_tautau_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi2_tautau_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_tautau_ATLAS13;
}



Hobs_ggF_phi2_tautau_CMS13::Hobs_ggF_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_tautau_CMS13;
}

Robs_ggF_phi2_tautau_CMS13::Robs_ggF_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_tautau_CMS13;
}



Hobs_bbF_phi2_tautau_CMS13::Hobs_bbF_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_bbF_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_bbF_phi2_tautau_CMS13;
}

Robs_bbF_phi2_tautau_CMS13::Robs_bbF_phi2_tautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_bbF_phi2_tautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_bbF_phi2_tautau_CMS13;
}


Hobs_pp_phi2_gaga_ATLAS13::Hobs_pp_phi2_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_gaga_ATLAS13;
}

Robs_pp_phi2_gaga_ATLAS13::Robs_pp_phi2_gaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_gaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_gaga_ATLAS13;
}

Hobs_ggF_phi2_gaga_CMS13::Hobs_ggF_phi2_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_gaga_CMS13;
}

Robs_ggF_phi2_gaga_CMS13::Robs_ggF_phi2_gaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_gaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_gaga_CMS13;
}

Hobs_pp_phi2_Zga_llga_ATLAS13::Hobs_pp_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_ATLAS13;
}

Robs_pp_phi2_Zga_llga_ATLAS13::Robs_pp_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_Zga_ATLAS13;
}

Hobs_ggF_phi2_Zga_llga_ATLAS13::Hobs_ggF_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_Zga_llga_ATLAS13;
}

Robs_ggF_phi2_Zga_llga_ATLAS13::Robs_ggF_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_Zga_llga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_Zga_llga_ATLAS13;
}

Hobs_pp_phi2_Zga_llga_CMS13::Hobs_pp_phi2_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_llga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_llga_CMS13;
}

Robs_pp_phi2_Zga_llga_CMS13::Robs_pp_phi2_Zga_llga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_Zga_llga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_Zga_llga_CMS13;
}



Hobs_pp_phi2_Zga_qqga_CMS13::Hobs_pp_phi2_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_Zga_qqga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_Zga_qqga_CMS13;
}

Robs_pp_phi2_Zga_qqga_CMS13::Robs_pp_phi2_Zga_qqga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_Zga_qqga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_Zga_qqga_CMS13;
}

Hobs_ggF_phi2_Zga_CMS13::Hobs_ggF_phi2_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_Zga_CMS13;
}

Robs_ggF_phi2_Zga_CMS13::Robs_ggF_phi2_Zga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_Zga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_Zga_CMS13;
}



Hobs_ggF_phi2_ZZ_llllnunu_ATLAS13::Hobs_ggF_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_llllnunu_ATLAS13;
}

Robs_ggF_phi2_ZZ_llllnunu_ATLAS13::Robs_ggF_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_llllnunu_ATLAS13;
}

Hobs_VBF_phi2_ZZ_llllnunu_ATLAS13::Hobs_VBF_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_llllnunu_ATLAS13;
}

Robs_VBF_phi2_ZZ_llllnunu_ATLAS13::Robs_VBF_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_llllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_llllnunu_ATLAS13;
}



Hobs_ggF_phi2_ZZ_llnunu_ATLAS13::Hobs_ggF_phi2_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_llnunu_ATLAS13;
}

Robs_ggF_phi2_ZZ_llnunu_ATLAS13::Robs_ggF_phi2_ZZ_llnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_llnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_llnunu_ATLAS13;
}

Hobs_pp_phi2_ZZ_llnunu_CMS13::Hobs_pp_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_ZZ_llnunu_CMS13;
}

Robs_pp_phi2_ZZ_llnunu_CMS13::Robs_pp_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_ZZ_llnunu_CMS13;
}

Hobs_ggF_phi2_ZZ_llnunu_CMS13::Hobs_ggF_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_llnunu_CMS13;
}

Robs_ggF_phi2_ZZ_llnunu_CMS13::Robs_ggF_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_llnunu_CMS13;
}

Hobs_VBF_phi2_ZZ_llnunu_CMS13::Hobs_VBF_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_llnunu_CMS13;
}

Robs_VBF_phi2_ZZ_llnunu_CMS13::Robs_VBF_phi2_ZZ_llnunu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_llnunu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_llnunu_CMS13;
}


Hobs_ggF_phi2_ZZ_llll_ATLAS13::Hobs_ggF_phi2_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_llll_ATLAS13;
}

Robs_ggF_phi2_ZZ_llll_ATLAS13::Robs_ggF_phi2_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_llll_ATLAS13;
}

Hobs_VBF_phi2_ZZ_llll_ATLAS13::Hobs_VBF_phi2_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_llll_ATLAS13;
}

Robs_VBF_phi2_ZZ_llll_ATLAS13::Robs_VBF_phi2_ZZ_llll_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_llll_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_llll_ATLAS13;
}

Hobs_pp_phi2_ZZ_llll_CMS13::Hobs_pp_phi2_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_ZZ_llll_CMS13;
}

Robs_pp_phi2_ZZ_llll_CMS13::Robs_pp_phi2_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_ZZ_llll_CMS13;
}

Hobs_VBF_VH_phi2_ZZ_llll_CMS13::Hobs_VBF_VH_phi2_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_VH_phi2_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_VH_phi2_ZZ_llll_CMS13;
}

Robs_VBF_VH_phi2_ZZ_llll_CMS13::Robs_VBF_VH_phi2_ZZ_llll_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_VH_phi2_ZZ_llll_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_VH_phi2_ZZ_llll_CMS13;
}

Hobs_ggF_phi2_ZZ_qqllnunu_ATLAS13::Hobs_ggF_phi2_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_qqllnunu_ATLAS13;
}

Robs_ggF_phi2_ZZ_qqllnunu_ATLAS13::Robs_ggF_phi2_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_qqllnunu_ATLAS13;
}



Hobs_VBF_phi2_ZZ_qqllnunu_ATLAS13::Hobs_VBF_phi2_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_qqllnunu_ATLAS13;
}

Robs_VBF_phi2_ZZ_qqllnunu_ATLAS13::Robs_VBF_phi2_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_qqllnunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_qqllnunu_ATLAS13;
}

Hobs_ggF_phi2_ZZ_llqq_ATLAS13::Hobs_ggF_phi2_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_llqq_ATLAS13;
}

Robs_ggF_phi2_ZZ_llqq_ATLAS13::Robs_ggF_phi2_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_llqq_ATLAS13;
}

Hobs_VBF_phi2_ZZ_llqq_ATLAS13::Hobs_VBF_phi2_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_ZZ_llqq_ATLAS13;
}

Robs_VBF_phi2_ZZ_llqq_ATLAS13::Robs_VBF_phi2_ZZ_llqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_ZZ_llqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_ZZ_llqq_ATLAS13;
}

Hobs_ggF_phi2_ZZ_nunuqq_ATLAS13::Hobs_ggF_phi2_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_ZZ_nunuqq_ATLAS13;
}

Robs_ggF_phi2_ZZ_nunuqq_ATLAS13::Robs_ggF_phi2_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_ZZ_nunuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_ZZ_nunuqq_ATLAS13;
}



Hobs_pp_phi2_ZZ_llqq_CMS13::Hobs_pp_phi2_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_ZZ_llqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_ZZ_llqq_CMS13;
}

Robs_pp_phi2_ZZ_llqq_CMS13::Robs_pp_phi2_ZZ_llqq_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_ZZ_llqq_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_ZZ_llqq_CMS13;
}

Hobs_ggF_phi2_WW_lnuqq_ATLAS13::Hobs_ggF_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_WW_lnuqq_ATLAS13;
}

Robs_ggF_phi2_WW_lnuqq_ATLAS13::Robs_ggF_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_WW_lnuqq_ATLAS13;
}

Hobs_VBF_phi2_WW_lnuqq_ATLAS13::Hobs_VBF_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_WW_lnuqq_ATLAS13;
}

Robs_VBF_phi2_WW_lnuqq_ATLAS13::Robs_VBF_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_WW_lnuqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_WW_lnuqq_ATLAS13;
}

Hobs_pp_phi2_VV_qqqq_ATLAS13::Hobs_pp_phi2_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_VV_qqqq_ATLAS13;
}

Robs_pp_phi2_VV_qqqq_ATLAS13::Robs_pp_phi2_VV_qqqq_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_VV_qqqq_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_VV_qqqq_ATLAS13;
}

Hobs_ggF_phi2_WW_enumunu_ATLAS13::Hobs_ggF_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_WW_enumunu_ATLAS13;
}

Robs_ggF_phi2_WW_enumunu_ATLAS13::Robs_ggF_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_WW_enumunu_ATLAS13;
}

Hobs_VBF_phi2_WW_enumunu_ATLAS13::Hobs_VBF_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_VBF_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_VBF_phi2_WW_enumunu_ATLAS13;
}

Robs_VBF_phi2_WW_enumunu_ATLAS13::Robs_VBF_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_VBF_phi2_WW_enumunu_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_VBF_phi2_WW_enumunu_ATLAS13;
}

Hobs_ggF_VBF_phi2_WW_lnulnu_CMS13::Hobs_ggF_VBF_phi2_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_VBF_phi2_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_VBF_phi2_WW_lnulnu_CMS13;
}

Robs_ggF_VBF_phi2_WW_lnulnu_CMS13::Robs_ggF_VBF_phi2_WW_lnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_VBF_phi2_WW_lnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_VBF_phi2_WW_lnulnu_CMS13;
}


Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13::Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbgaga_ATLAS13;
}

Robs_pp_phi2_phi1phi1_bbgaga_ATLAS13::Robs_pp_phi2_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbgaga_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbgaga_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbgaga_CMS13::Hobs_pp_phi2_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbgaga_CMS13;
}

Robs_pp_phi2_phi1phi1_bbgaga_CMS13::Robs_pp_phi2_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbgaga_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbgaga_CMS13;
}


Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13::Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_ATLAS13;
}

Robs_pp_phi2_phi1phi1_bbbb_ATLAS13::Robs_pp_phi2_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbbb_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbbb_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbbb_CMS13::Hobs_pp_phi2_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbbb_CMS13;
}

Robs_pp_phi2_phi1phi1_bbbb_CMS13::Robs_pp_phi2_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbbb_CMS13;
}

Hobs_ggF_phi2_phi1phi1_bbbb_CMS13::Hobs_ggF_phi2_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_phi1phi1_bbbb_CMS13;
}

Robs_ggF_phi2_phi1phi1_bbbb_CMS13::Robs_ggF_phi2_phi1phi1_bbbb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_phi1phi1_bbbb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_phi1phi1_bbbb_CMS13;
}



Hobs_ggF_phi2_phi1phi1_gagaWW_ATLAS13::Hobs_ggF_phi2_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_ggF_phi2_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_ggF_phi2_phi1phi1_gagaWW_ATLAS13;
}

Robs_ggF_phi2_phi1phi1_gagaWW_ATLAS13::Robs_ggF_phi2_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_ggF_phi2_phi1phi1_gagaWW_ATLAS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_ggF_phi2_phi1phi1_gagaWW_ATLAS13;
}

Hobs_pp_phi2_phi1phi1_bbtautau_CMS13::Hobs_pp_phi2_phi1phi1_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau_CMS13;
}

Robs_pp_phi2_phi1phi1_bbtautau_CMS13::Robs_pp_phi2_phi1phi1_bbtautau_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbtautau_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbtautau_CMS13;
}

Hobs_pp_phi2_phi1phi1_bbtautau1_CMS13::Hobs_pp_phi2_phi1phi1_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbtautau1_CMS13;
}

Robs_pp_phi2_phi1phi1_bbtautau1_CMS13::Robs_pp_phi2_phi1phi1_bbtautau1_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbtautau1_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbtautau1_CMS13;
}

Hobs_pp_phi2_phi1phi1_bblnulnu_CMS13::Hobs_pp_phi2_phi1phi1_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bblnulnu_CMS13;
}

Robs_pp_phi2_phi1phi1_bblnulnu_CMS13::Robs_pp_phi2_phi1phi1_bblnulnu_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bblnulnu_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bblnulnu_CMS13;
}



Hobs_pp_phi2_phi1phi1_bbVV_CMS13::Hobs_pp_phi2_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_phi1phi1_bbVV_CMS13;
}

Robs_pp_phi2_phi1phi1_bbVV_CMS13::Robs_pp_phi2_phi1phi1_bbVV_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_phi1phi1_bbVV_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_phi1phi1_bbVV_CMS13;
}

Hobs_pp_phi2_bb_CMS13::Hobs_pp_phi2_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Hobs_pp_phi2_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->THoEX_pp_phi2_bb_CMS13;
}

Robs_pp_phi2_bb_CMS13::Robs_pp_phi2_bb_CMS13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Robs_pp_phi2_bb_CMS13::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->R_pp_phi2_bb_CMS13;
}



log10_ggF_phi2_tautau_TH8::log10_ggF_phi2_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_tautau_TH8);
}


log10_bbF_phi2_tautau_TH8::log10_bbF_phi2_tautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi2_tautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi2_tautau_TH8);
}



log10_pp_phi2_gaga_TH8::log10_pp_phi2_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_gaga_TH8);
}



log10_ggF_phi2_gaga_TH8::log10_ggF_phi2_gaga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_gaga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_gaga_TH8);
}



log10_pp_phi2_Zga_llga_TH8::log10_pp_phi2_Zga_llga_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_Zga_llga_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_Zga_llga_TH8);
}



log10_mu_pp_phi2_VV_TH8::log10_mu_pp_phi2_VV_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_mu_pp_phi2_VV_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->mu_pp_phi2_VV_TH8);
}


log10_ggF_phi2_ZZ_TH8::log10_ggF_phi2_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_ZZ_TH8);
}


log10_VBF_phi2_ZZ_TH8::log10_VBF_phi2_ZZ_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi2_ZZ_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi2_ZZ_TH8);
}



log10_ggF_phi2_WW_TH8::log10_ggF_phi2_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_WW_TH8);
}


log10_VBF_phi2_WW_TH8::log10_VBF_phi2_WW_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi2_WW_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi2_WW_TH8);
}

//ggF and pp TH8

log10_ggF_phi2_phi1phi1_TH8::log10_ggF_phi2_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_phi1phi1_TH8);
}

log10_pp_phi2_phi1phi1_TH8::log10_pp_phi2_phi1phi1_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_TH8);
}

// ggF  bbtautau

log10_ggF_phi2_phi1phi1_bbtautau_TH8::log10_ggF_phi2_phi1phi1_bbtautau_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_phi1phi1_bbtautau_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_phi1phi1_bbtautau_TH8);
}


// pp  bbbb

log10_pp_phi2_phi1phi1_bbbb_TH8::log10_pp_phi2_phi1phi1_bbbb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbbb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbbb_TH8);
}

//pp gagabb


log10_pp_phi2_phi1phi1_gagabb_TH8::log10_pp_phi2_phi1phi1_gagabb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_gagabb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_gagabb_TH8);
}


log10_ggF_phi2_tt_TH8::log10_ggF_phi2_tt_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_tt_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_tt_TH8);
}



log10_bbF_phi2_bb_TH8::log10_bbF_phi2_bb_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi2_bb_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi2_bb_TH8);
}


//pp AZbbll TH8

log10_pp_phi2_phi1Z_bbll_TH8::log10_pp_phi2_phi1Z_bbll_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1Z_bbll_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1Z_bbll_TH8);
}



//pp AZ tautaull

log10_pp_phi2_phi1Z_tautaull_TH8::log10_pp_phi2_phi1Z_tautaull_TH8(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1Z_tautaull_TH8::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1Z_tautaull_TH8);
}


log10_ggF_phi2_tautau_TH13::log10_ggF_phi2_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_tautau_TH13);
}



log10_bbF_phi2_tautau_TH13::log10_bbF_phi2_tautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bbF_phi2_tautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi2_tautau_TH13);
}



log10_pp_phi2_gaga_TH13::log10_pp_phi2_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_gaga_TH13);
}



log10_ggF_phi2_gaga_TH13::log10_ggF_phi2_gaga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_gaga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_gaga_TH13);
}



log10_pp_phi2_Zga_TH13::log10_pp_phi2_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_Zga_TH13);
}



log10_ggF_phi2_Zga_TH13::log10_ggF_phi2_Zga_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_Zga_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_Zga_TH13);
}

log10_pp_phi2_ZZ_TH13::log10_pp_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_ZZ_TH13);
}


log10_ggF_phi2_ZZ_TH13::log10_ggF_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_ZZ_TH13);
}


log10_VBF_phi2_ZZ_TH13::log10_VBF_phi2_ZZ_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi2_ZZ_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi2_ZZ_TH13);
}


log10_ggF_phi2_ZZ_llll_TH13::log10_ggF_phi2_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_ZZ_llll_TH13);
}

log10_VBF_phi2_ZZ_llll_TH13::log10_VBF_phi2_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi2_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi2_ZZ_llll_TH13);
}



log10_pp_phi2_ZZ_llll_TH13::log10_pp_phi2_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_ZZ_llll_TH13);
}

log10_VBF_VH_phi2_ZZ_llll_TH13::log10_VBF_VH_phi2_ZZ_llll_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_VH_phi2_ZZ_llll_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_VH_phi2_ZZ_llll_TH13);
}


log10_ggF_phi2_WW_TH13::log10_ggF_phi2_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_WW_TH13);
}



log10_VBF_phi2_WW_TH13::log10_VBF_phi2_WW_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_VBF_phi2_WW_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->VBF_phi2_WW_TH13);
}


log10_ggF_VBF_phi2_WW_lnulnu_TH13::log10_ggF_VBF_phi2_WW_lnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_VBF_phi2_WW_lnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_VBF_phi2_WW_lnulnu_TH13);
}


// ggF, phi1phi1, TH13

log10_ggF_phi2_phi1phi1_TH13::log10_ggF_phi2_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_phi1phi1_TH13);
}


//pp phi1phi1 TH13

log10_pp_phi2_phi1phi1_TH13::log10_pp_phi2_phi1phi1_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_TH13);
}
//pp bbbb

log10_pp_phi2_phi1phi1_bbbb_TH13::log10_pp_phi2_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbbb_TH13);
}

//ggF bbbb

log10_ggF_phi2_phi1phi1_bbbb_TH13::log10_ggF_phi2_phi1phi1_bbbb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_ggF_phi2_phi1phi1_bbbb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ggF_phi2_phi1phi1_bbbb_TH13);
}
//pp gagabb

log10_pp_phi2_phi1phi1_gagabb_TH13::log10_pp_phi2_phi1phi1_gagabb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_gagabb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_gagabb_TH13);
}

//pp bbtautau

log10_pp_phi2_phi1phi1_bbtautau_TH13::log10_pp_phi2_phi1phi1_bbtautau_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbtautau_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbtautau_TH13);
}

//pp bblnulnu

log10_pp_phi2_phi1phi1_bblnulnu_TH13::log10_pp_phi2_phi1phi1_bblnulnu_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bblnulnu_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bblnulnu_TH13);
}

//pp bbVV

log10_pp_phi2_phi1phi1_bbVV_TH13::log10_pp_phi2_phi1phi1_bbVV_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_phi1phi1_bbVV_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_phi1phi1_bbVV_TH13);
}
//

log10_tt_phi2_tt_TH13::log10_tt_phi2_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_tt_phi2_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->ttF_phi2_tt_TH13);
}



log10_bb_phi2_tt_TH13::log10_bb_phi2_tt_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_bb_phi2_tt_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->bbF_phi2_tt_TH13);
}



log10_pp_phi2_bb_TH13::log10_pp_phi2_bb_TH13(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double log10_pp_phi2_bb_TH13::computeThValue()
{
    return log10(myGTHDM.getMyGTHDMCache()->pp_phi2_bb_TH13);
}


Gamma_phi2_GTHDM::Gamma_phi2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double Gamma_phi2_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Gammaphi2tot;
}



//rHH_gaga_THDM::rHH_gaga_THDM(const StandardModel& SM_i)
//: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
//{}
//
//double rHH_gaga_THDM::computeThValue()
//{
//    return myGTHDM.getMyGTHDMCache()->rHH_gaga;
//}



rphi2_ggE_GTHDM::rphi2_ggE_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi2_ggE_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi2_ggE;
}

rphi2_ggO_GTHDM::rphi2_ggO_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double rphi2_ggO_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->rphi2_ggO;
}

// phi2 -> phi1phi1 and phi2 -> phi2phi2

BR_phi2_phi1phi1_GTHDM::BR_phi2_phi1phi1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_phi1phi1_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2tophi1phi1;
}

//phi2 -> H+H-


BR_phi2_HpHm_GTHDM::BR_phi2_HpHm_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_HpHm_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2toHpHm;
}


BR_phi2_phi1Z_GTHDM::BR_phi2_phi1Z_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_phi1Z_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2tophi1Z;
}

BR_phi2_HpW_GTHDM::BR_phi2_HpW_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i))
{}

double BR_phi2_HpW_GTHDM::computeThValue()
{
    return myGTHDM.getMyGTHDMCache()->Br_phi2toHpW;
}
