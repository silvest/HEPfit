/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "EWObservables.h"
#include "DiBosonThObservables.h"
#include "OtherThObservables.h"
using namespace ThObsConst;

void ThObsFactory::registerDiBosonObservables()
{
    //-----  ee -> WW observables: LEP2 total cross section  -----
    obsThFactory["eeWW_LEP2_161"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_161); };
    obsThFactory["eeWW_LEP2_172"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_172); };
    obsThFactory["eeWW_LEP2_183"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_183); };
    obsThFactory["eeWW_LEP2_189"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_189); };
    obsThFactory["eeWW_LEP2_192"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_192); };
    obsThFactory["eeWW_LEP2_196"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_196); };
    obsThFactory["eeWW_LEP2_200"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_200); };
    obsThFactory["eeWW_LEP2_202"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_202); };
    obsThFactory["eeWW_LEP2_205"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_205); };
    obsThFactory["eeWW_LEP2_206"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_206); };
    obsThFactory["eeWW_LEP2_207"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_207); };
    obsThFactory["eeWW_LEP2_208"] = [=](const StandardModel& SM) { return new xseeWW(SM, sqrt_s_LEP2_208); };
    // Similar observables, defined only for the d 6 SMEFT, from arXiv: 1606.06693 [hep-ph].
    obsThFactory["eeWWlept_LEP2_189"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_189); };
    obsThFactory["eeWWlept_LEP2_192"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_192); };
    obsThFactory["eeWWlept_LEP2_196"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_196); };
    obsThFactory["eeWWlept_LEP2_200"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_200); };
    obsThFactory["eeWWlept_LEP2_202"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_202); };
    obsThFactory["eeWWlept_LEP2_205"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_205); };
    obsThFactory["eeWWlept_LEP2_206"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_206); };
    obsThFactory["eeWWlept_LEP2_207"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_207); };
    obsThFactory["eeWWlept_LEP2_208"] = [=](const StandardModel& SM) { return new xseeWWlept(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["eeWWsemil_LEP2_189"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_189); };
    obsThFactory["eeWWsemil_LEP2_192"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_192); };
    obsThFactory["eeWWsemil_LEP2_196"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_196); };
    obsThFactory["eeWWsemil_LEP2_200"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_200); };
    obsThFactory["eeWWsemil_LEP2_202"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_202); };
    obsThFactory["eeWWsemil_LEP2_205"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_205); };
    obsThFactory["eeWWsemil_LEP2_206"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_206); };
    obsThFactory["eeWWsemil_LEP2_207"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_207); };
    obsThFactory["eeWWsemil_LEP2_208"] = [=](const StandardModel& SM) { return new xseeWWsemil(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["eeWWhad_LEP2_189"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_189); };
    obsThFactory["eeWWhad_LEP2_192"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_192); };
    obsThFactory["eeWWhad_LEP2_196"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_196); };
    obsThFactory["eeWWhad_LEP2_200"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_200); };
    obsThFactory["eeWWhad_LEP2_202"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_202); };
    obsThFactory["eeWWhad_LEP2_205"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_205); };
    obsThFactory["eeWWhad_LEP2_206"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_206); };
    obsThFactory["eeWWhad_LEP2_207"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_207); };
    obsThFactory["eeWWhad_LEP2_208"] = [=](const StandardModel& SM) { return new xseeWWhad(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["eeWWtot_LEP2_189"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_189); };
    obsThFactory["eeWWtot_LEP2_192"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_192); };
    obsThFactory["eeWWtot_LEP2_196"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_196); };
    obsThFactory["eeWWtot_LEP2_200"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_200); };
    obsThFactory["eeWWtot_LEP2_202"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_202); };
    obsThFactory["eeWWtot_LEP2_205"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_205); };
    obsThFactory["eeWWtot_LEP2_206"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_206); };
    obsThFactory["eeWWtot_LEP2_207"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_207); };
    obsThFactory["eeWWtot_LEP2_208"] = [=](const StandardModel& SM) { return new xseeWWtot(SM, sqrt_s_LEP2_208); };
    // The same, but only the new physics contribution
    obsThFactory["deltaeeWWlept_LEP2_189"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_189); };
    obsThFactory["deltaeeWWlept_LEP2_192"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_192); };
    obsThFactory["deltaeeWWlept_LEP2_196"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_196); };
    obsThFactory["deltaeeWWlept_LEP2_200"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_200); };
    obsThFactory["deltaeeWWlept_LEP2_202"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_202); };
    obsThFactory["deltaeeWWlept_LEP2_205"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_205); };
    obsThFactory["deltaeeWWlept_LEP2_206"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_206); };
    obsThFactory["deltaeeWWlept_LEP2_207"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_207); };
    obsThFactory["deltaeeWWlept_LEP2_208"] = [=](const StandardModel& SM) { return new deltaxseeWWlept(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["deltaeeWWsemil_LEP2_189"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_189); };
    obsThFactory["deltaeeWWsemil_LEP2_192"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_192); };
    obsThFactory["deltaeeWWsemil_LEP2_196"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_196); };
    obsThFactory["deltaeeWWsemil_LEP2_200"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_200); };
    obsThFactory["deltaeeWWsemil_LEP2_202"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_202); };
    obsThFactory["deltaeeWWsemil_LEP2_205"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_205); };
    obsThFactory["deltaeeWWsemil_LEP2_206"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_206); };
    obsThFactory["deltaeeWWsemil_LEP2_207"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_207); };
    obsThFactory["deltaeeWWsemil_LEP2_208"] = [=](const StandardModel& SM) { return new deltaxseeWWsemil(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["deltaeeWWhad_LEP2_189"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_189); };
    obsThFactory["deltaeeWWhad_LEP2_192"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_192); };
    obsThFactory["deltaeeWWhad_LEP2_196"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_196); };
    obsThFactory["deltaeeWWhad_LEP2_200"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_200); };
    obsThFactory["deltaeeWWhad_LEP2_202"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_202); };
    obsThFactory["deltaeeWWhad_LEP2_205"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_205); };
    obsThFactory["deltaeeWWhad_LEP2_206"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_206); };
    obsThFactory["deltaeeWWhad_LEP2_207"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_207); };
    obsThFactory["deltaeeWWhad_LEP2_208"] = [=](const StandardModel& SM) { return new deltaxseeWWhad(SM, sqrt_s_LEP2_208); };
    //
    obsThFactory["deltaeeWWtot_LEP2_189"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_189); };
    obsThFactory["deltaeeWWtot_LEP2_192"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_192); };
    obsThFactory["deltaeeWWtot_LEP2_196"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_196); };
    obsThFactory["deltaeeWWtot_LEP2_200"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_200); };
    obsThFactory["deltaeeWWtot_LEP2_202"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_202); };
    obsThFactory["deltaeeWWtot_LEP2_205"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_205); };
    obsThFactory["deltaeeWWtot_LEP2_206"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_206); };
    obsThFactory["deltaeeWWtot_LEP2_207"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_207); };
    obsThFactory["deltaeeWWtot_LEP2_208"] = [=](const StandardModel& SM) { return new deltaxseeWWtot(SM, sqrt_s_LEP2_208); };
    //-----  ee -> WW observables: LEP2 differential cross section  -----
    obsThFactory["deeWWdcos_LEP2_183_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos1_LEP2_WW, cos2_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos2_LEP2_WW, cos3_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos3_LEP2_WW, cos4_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos4_LEP2_WW, cos5_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos5_LEP2_WW, cos6_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos6_LEP2_WW, cos7_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos7_LEP2_WW, cos8_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos8_LEP2_WW, cos9_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos9_LEP2_WW, cos10_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_183_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos1, cos10_LEP2_WW, cos11_LEP2_WW); };
    //
    obsThFactory["deeWWdcos_LEP2_189_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos1_LEP2_WW, cos2_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos2_LEP2_WW, cos3_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos3_LEP2_WW, cos4_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos4_LEP2_WW, cos5_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos5_LEP2_WW, cos6_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos6_LEP2_WW, cos7_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos7_LEP2_WW, cos8_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos8_LEP2_WW, cos9_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos9_LEP2_WW, cos10_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_189_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos2, cos10_LEP2_WW, cos11_LEP2_WW); };
    //
    obsThFactory["deeWWdcos_LEP2_198_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos1_LEP2_WW, cos2_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos2_LEP2_WW, cos3_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos3_LEP2_WW, cos4_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos4_LEP2_WW, cos5_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos5_LEP2_WW, cos6_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos6_LEP2_WW, cos7_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos7_LEP2_WW, cos8_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos8_LEP2_WW, cos9_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos9_LEP2_WW, cos10_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_198_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos3, cos10_LEP2_WW, cos11_LEP2_WW); };
    //
    obsThFactory["deeWWdcos_LEP2_206_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos1_LEP2_WW, cos2_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos2_LEP2_WW, cos3_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos3_LEP2_WW, cos4_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos4_LEP2_WW, cos5_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos5_LEP2_WW, cos6_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos6_LEP2_WW, cos7_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos7_LEP2_WW, cos8_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos8_LEP2_WW, cos9_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos9_LEP2_WW, cos10_LEP2_WW); };
    obsThFactory["deeWWdcos_LEP2_206_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_LEP2_WWcos4, cos10_LEP2_WW, cos11_LEP2_WW); };
    // Similar observables, defined only for the d 6 SMEFT, from arXiv: 1606.06693 [hep-ph].
    obsThFactory["deeWW_LEP2_183_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 1); };
    obsThFactory["deeWW_LEP2_183_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 2); };
    obsThFactory["deeWW_LEP2_183_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 3); };
    obsThFactory["deeWW_LEP2_183_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 4); };
    //
    obsThFactory["deeWW_LEP2_206_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 1); };
    obsThFactory["deeWW_LEP2_206_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 2); };
    obsThFactory["deeWW_LEP2_206_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 3); };
    obsThFactory["deeWW_LEP2_206_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 4); };
    // The same but only the NP contribution
    obsThFactory["deltadeeWW_LEP2_183_Bin1"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 1); };
    obsThFactory["deltadeeWW_LEP2_183_Bin2"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 2); };
    obsThFactory["deltadeeWW_LEP2_183_Bin3"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 3); };
    obsThFactory["deltadeeWW_LEP2_183_Bin4"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_183, 4); };
    //
    obsThFactory["deltadeeWW_LEP2_206_Bin1"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 1); };
    obsThFactory["deltadeeWW_LEP2_206_Bin2"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 2); };
    obsThFactory["deltadeeWW_LEP2_206_Bin3"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 3); };
    obsThFactory["deltadeeWW_LEP2_206_Bin4"] = [=](const StandardModel& SM) { return new deltadxseeWWLEP2Bin(SM, sqrt_s_LEP2_206, 4); };
    //-----  ee -> WW observables: Future colliders differential cross section  -----
    obsThFactory["deeWWdcos_161_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos1_ee_WW, cos2_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos2_ee_WW, cos3_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos3_ee_WW, cos4_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos4_ee_WW, cos5_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos5_ee_WW, cos6_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos6_ee_WW, cos7_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos7_ee_WW, cos8_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos8_ee_WW, cos9_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos9_ee_WW, cos10_ee_WW); };
    obsThFactory["deeWWdcos_161_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_161, cos10_ee_WW, cos11_ee_WW); };
    //
    obsThFactory["deeWWdcos_240_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos1_ee_WW, cos2_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos2_ee_WW, cos3_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos3_ee_WW, cos4_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos4_ee_WW, cos5_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos5_ee_WW, cos6_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos6_ee_WW, cos7_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos7_ee_WW, cos8_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos8_ee_WW, cos9_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos9_ee_WW, cos10_ee_WW); };
    obsThFactory["deeWWdcos_240_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_240, cos10_ee_WW, cos11_ee_WW); };
    //
    obsThFactory["deeWWdcos_250_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos1_ee_WW, cos2_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos2_ee_WW, cos3_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos3_ee_WW, cos4_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos4_ee_WW, cos5_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos5_ee_WW, cos6_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos6_ee_WW, cos7_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos7_ee_WW, cos8_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos8_ee_WW, cos9_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos9_ee_WW, cos10_ee_WW); };
    obsThFactory["deeWWdcos_250_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_250, cos10_ee_WW, cos11_ee_WW); };
    //
    obsThFactory["deeWWdcos_350_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos1_ee_WW, cos2_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos2_ee_WW, cos3_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos3_ee_WW, cos4_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos4_ee_WW, cos5_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos5_ee_WW, cos6_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos6_ee_WW, cos7_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos7_ee_WW, cos8_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos8_ee_WW, cos9_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos9_ee_WW, cos10_ee_WW); };
    obsThFactory["deeWWdcos_350_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_350, cos10_ee_WW, cos11_ee_WW); };
    //
    obsThFactory["deeWWdcos_365_Bin1"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos1_ee_WW, cos2_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin2"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos2_ee_WW, cos3_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin3"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos3_ee_WW, cos4_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin4"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos4_ee_WW, cos5_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin5"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos5_ee_WW, cos6_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin6"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos6_ee_WW, cos7_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin7"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos7_ee_WW, cos8_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin8"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos8_ee_WW, cos9_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin9"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos9_ee_WW, cos10_ee_WW); };
    obsThFactory["deeWWdcos_365_Bin10"] = [=](const StandardModel& SM) { return new dxseeWWdcosBin(SM, sqrt_s_leptcoll_365, cos10_ee_WW, cos11_ee_WW); };
    //-----  ee -> WW observables: total rates (ratio with the SM)  -----
    obsThFactory["eeWW161"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_161, 0., 0.); };
    //
    obsThFactory["eeWW230"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeWW230_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_230, pol_80, -pol_30); };
    obsThFactory["eeWW230_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_230, -pol_80, pol_30); };
    //
    obsThFactory["eeWW240"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeWW240_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_240, pol_80, -pol_30); };
    obsThFactory["eeWW240_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_240, -pol_80, pol_30); };
    //
    obsThFactory["eeWW250"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeWW250_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeWW250_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeWW250_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeWW250_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeWW350"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeWW350_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeWW350_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeWW350_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeWW350_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeWW365"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeWW365_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeWW365_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    //
    obsThFactory["eeWW380_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeWW380_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeWW500"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeWW500_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeWW500_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeWW500_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeWW500_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeWW550"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeWW550_p80_m30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeWW550_m80_p30"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeWW550_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeWW550_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeWW1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeWW1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    //
    obsThFactory["eeWW1500_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeWW1500_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeWW3000_p80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeWW3000_m80_0"] = [=](const StandardModel& SM) { return new mueeWW(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //----- High Energy diboson observables at hadron colliders
    obsThFactory["ppZHprobe14"] = [=](const StandardModel& SM) { return new ppZHprobe(SM, sqrt_s_LHC14); };
    obsThFactory["ppZHprobe27"] = [=](const StandardModel& SM) { return new ppZHprobe(SM, sqrt_s_LHC27); };
    obsThFactory["ppZHprobe100"] = [=](const StandardModel& SM) { return new ppZHprobe(SM, sqrt_s_FCC100); };
    //
    obsThFactory["mupTVppWZ_14_Bin1"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 100., 150.); };
    obsThFactory["mupTVppWZ_14_Bin2"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 150., 220.); };
    obsThFactory["mupTVppWZ_14_Bin3"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 220., 300.); };
    obsThFactory["mupTVppWZ_14_Bin4"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 300., 500.); };
    obsThFactory["mupTVppWZ_14_Bin5"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 500., 750.); };
    obsThFactory["mupTVppWZ_14_Bin6"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC14, 750., 1200.); };
    //
    obsThFactory["mupTVppWZ_27_Bin1"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 150., 220.); };
    obsThFactory["mupTVppWZ_27_Bin2"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 220., 300.); };
    obsThFactory["mupTVppWZ_27_Bin3"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 300., 500.); };
    obsThFactory["mupTVppWZ_27_Bin4"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 500., 750.); };
    obsThFactory["mupTVppWZ_27_Bin5"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 750., 1200.); };
    obsThFactory["mupTVppWZ_27_Bin6"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_LHC27, 1200., 1800.); };
    //
    obsThFactory["mupTVppWZ_100_Bin1"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 220., 300.); };
    obsThFactory["mupTVppWZ_100_Bin2"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 300., 500.); };
    obsThFactory["mupTVppWZ_100_Bin3"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 500., 750.); };
    obsThFactory["mupTVppWZ_100_Bin4"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 750., 1200.); };
    obsThFactory["mupTVppWZ_100_Bin5"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 1200., 1800.); };
    obsThFactory["mupTVppWZ_100_Bin6"] = [=](const StandardModel& SM) { return new mupTVppWZ(SM, sqrt_s_FCC100, 1800., 2400.); };
    //
    //-----  Collider observables: LHC dilepton events  ----------
    //----- p p > e e
    obsThFactory["NEvppee13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 1); };
    obsThFactory["NEvppee13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 2); };
    obsThFactory["NEvppee13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 3); };
    obsThFactory["NEvppee13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 4); };
    obsThFactory["NEvppee13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 5); };
    obsThFactory["NEvppee13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 6); };
    obsThFactory["NEvppee13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 7); };
    obsThFactory["NEvppee13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 8); };
    obsThFactory["NEvppee13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 9); };
    obsThFactory["NEvppee13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 10); };
    obsThFactory["NEvppee13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 11); };
    obsThFactory["NEvppee13_Bin12"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 12); };
    obsThFactory["NEvppee13_Bin13"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 13); };
    obsThFactory["NEvppee13_Bin14"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 14); };
    obsThFactory["NEvppee13_Bin15"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 15); };
    obsThFactory["NEvppee13_Bin16"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 16); };
    obsThFactory["NEvppee13_Bin17"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 17); };
    obsThFactory["NEvppee13_Bin18"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 18); };
    obsThFactory["NEvppee13_Bin19"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 19); };
    obsThFactory["NEvppee13_Bin20"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 20); };
    obsThFactory["NEvppee13_Bin21"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 21); };
    obsThFactory["NEvppee13_Bin22"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 22); };
    obsThFactory["NEvppee13_Bin23"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 23); };
    obsThFactory["NEvppee13_Bin24"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 24); };
    obsThFactory["NEvppee13_Bin25"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 25); };
    obsThFactory["NEvppee13_Bin26"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 26); };
    obsThFactory["NEvppee13_Bin27"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 27); };
    obsThFactory["NEvppee13_Bin28"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 28); };
    obsThFactory["NEvppee13_Bin29"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 29); };
    obsThFactory["NEvppee13_Bin30"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 30); };
    obsThFactory["NEvppee13_Bin31"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 31); };
    obsThFactory["NEvppee13_Bin32"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 32); };
    obsThFactory["NEvppee13_Bin33"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 33); };
    obsThFactory["NEvppee13_Bin34"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 34); };
    obsThFactory["NEvppee13_Bin35"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 35); };
    obsThFactory["NEvppee13_Bin36"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 36); };
    obsThFactory["NEvppee13_Bin37"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 37); };
    obsThFactory["NEvppee13_Bin38"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 38); };
    obsThFactory["NEvppee13_Bin39"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 39); };
    obsThFactory["NEvppee13_Bin40"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 40); };
    obsThFactory["NEvppee13_Bin41"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 41); };
    obsThFactory["NEvppee13_Bin42"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 42); };
    obsThFactory["NEvppee13_Bin43"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 43); };
    obsThFactory["NEvppee13_Bin44"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 44); };
    obsThFactory["NEvppee13_Bin45"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 45); };
    obsThFactory["NEvppee13_Bin46"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 46); };
    obsThFactory["NEvppee13_Bin47"] = [=](const StandardModel& SM) { return new NevLHCee13(SM, 47); };
    //
    //----- p p > mu mu
    obsThFactory["NEvppmumu13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 1); };
    obsThFactory["NEvppmumu13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 2); };
    obsThFactory["NEvppmumu13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 3); };
    obsThFactory["NEvppmumu13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 4); };
    obsThFactory["NEvppmumu13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 5); };
    obsThFactory["NEvppmumu13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 6); };
    obsThFactory["NEvppmumu13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 7); };
    obsThFactory["NEvppmumu13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 8); };
    obsThFactory["NEvppmumu13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 9); };
    obsThFactory["NEvppmumu13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 10); };
    obsThFactory["NEvppmumu13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 11); };
    obsThFactory["NEvppmumu13_Bin12"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 12); };
    obsThFactory["NEvppmumu13_Bin13"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 13); };
    obsThFactory["NEvppmumu13_Bin14"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 14); };
    obsThFactory["NEvppmumu13_Bin15"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 15); };
    obsThFactory["NEvppmumu13_Bin16"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 16); };
    obsThFactory["NEvppmumu13_Bin17"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 17); };
    obsThFactory["NEvppmumu13_Bin18"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 18); };
    obsThFactory["NEvppmumu13_Bin19"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 19); };
    obsThFactory["NEvppmumu13_Bin20"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 20); };
    obsThFactory["NEvppmumu13_Bin21"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 21); };
    obsThFactory["NEvppmumu13_Bin22"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 22); };
    obsThFactory["NEvppmumu13_Bin23"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 23); };
    obsThFactory["NEvppmumu13_Bin24"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 24); };
    obsThFactory["NEvppmumu13_Bin25"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 25); };
    obsThFactory["NEvppmumu13_Bin26"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 26); };
    obsThFactory["NEvppmumu13_Bin27"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 27); };
    obsThFactory["NEvppmumu13_Bin28"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 28); };
    obsThFactory["NEvppmumu13_Bin29"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 29); };
    obsThFactory["NEvppmumu13_Bin30"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 30); };
    obsThFactory["NEvppmumu13_Bin31"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 31); };
    obsThFactory["NEvppmumu13_Bin32"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 32); };
    obsThFactory["NEvppmumu13_Bin33"] = [=](const StandardModel& SM) { return new NevLHCmumu13(SM, 33); };
    //
    //----- p p > tau tau
    obsThFactory["NEvpptautau13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 1); };
    obsThFactory["NEvpptautau13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 2); };
    obsThFactory["NEvpptautau13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 3); };
    obsThFactory["NEvpptautau13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 4); };
    obsThFactory["NEvpptautau13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 5); };
    obsThFactory["NEvpptautau13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 6); };
    obsThFactory["NEvpptautau13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 7); };
    obsThFactory["NEvpptautau13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 8); };
    obsThFactory["NEvpptautau13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 9); };
    obsThFactory["NEvpptautau13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 10); };
    obsThFactory["NEvpptautau13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 11); };
    obsThFactory["NEvpptautau13_Bin12"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 12); };
    obsThFactory["NEvpptautau13_Bin13"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 13); };
    obsThFactory["NEvpptautau13_Bin14"] = [=](const StandardModel& SM) { return new NevLHCtautau13(SM, 14); };
    //
    //-----  Collider observables: LHC mono-lepton events  ----------
    //----- p p > e nu
    obsThFactory["NEvppenu13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 1); };
    obsThFactory["NEvppenu13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 2); };
    obsThFactory["NEvppenu13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 3); };
    obsThFactory["NEvppenu13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 4); };
    obsThFactory["NEvppenu13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 5); };
    obsThFactory["NEvppenu13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 6); };
    obsThFactory["NEvppenu13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 7); };
    obsThFactory["NEvppenu13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 8); };
    obsThFactory["NEvppenu13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 9); };
    obsThFactory["NEvppenu13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 10); };
    obsThFactory["NEvppenu13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 11); };
    obsThFactory["NEvppenu13_Bin12"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 12); };
    obsThFactory["NEvppenu13_Bin13"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 13); };
    obsThFactory["NEvppenu13_Bin14"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 14); };
    obsThFactory["NEvppenu13_Bin15"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 15); };
    obsThFactory["NEvppenu13_Bin16"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 16); };
    obsThFactory["NEvppenu13_Bin17"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 17); };
    obsThFactory["NEvppenu13_Bin18"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 18); };
    obsThFactory["NEvppenu13_Bin19"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 19); };
    obsThFactory["NEvppenu13_Bin20"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 20); };
    obsThFactory["NEvppenu13_Bin21"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 21); };
    obsThFactory["NEvppenu13_Bin22"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 22); };
    obsThFactory["NEvppenu13_Bin23"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 23); };
    obsThFactory["NEvppenu13_Bin24"] = [=](const StandardModel& SM) { return new NevLHCenu13(SM, 24); };
    //
    //----- p p > mu nu
    obsThFactory["NEvppmunu13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 1); };
    obsThFactory["NEvppmunu13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 2); };
    obsThFactory["NEvppmunu13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 3); };
    obsThFactory["NEvppmunu13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 4); };
    obsThFactory["NEvppmunu13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 5); };
    obsThFactory["NEvppmunu13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 6); };
    obsThFactory["NEvppmunu13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 7); };
    obsThFactory["NEvppmunu13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 8); };
    obsThFactory["NEvppmunu13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 9); };
    obsThFactory["NEvppmunu13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 10); };
    obsThFactory["NEvppmunu13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 11); };
    obsThFactory["NEvppmunu13_Bin12"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 12); };
    obsThFactory["NEvppmunu13_Bin13"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 13); };
    obsThFactory["NEvppmunu13_Bin14"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 14); };
    obsThFactory["NEvppmunu13_Bin15"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 15); };
    obsThFactory["NEvppmunu13_Bin16"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 16); };
    obsThFactory["NEvppmunu13_Bin17"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 17); };
    obsThFactory["NEvppmunu13_Bin18"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 18); };
    obsThFactory["NEvppmunu13_Bin19"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 19); };
    obsThFactory["NEvppmunu13_Bin20"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 20); };
    obsThFactory["NEvppmunu13_Bin21"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 21); };
    obsThFactory["NEvppmunu13_Bin22"] = [=](const StandardModel& SM) { return new NevLHCmunu13(SM, 22); };
    //
    //----- p p > tau nu
    obsThFactory["NEvpptaunu13_Bin1"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 1); };
    obsThFactory["NEvpptaunu13_Bin2"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 2); };
    obsThFactory["NEvpptaunu13_Bin3"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 3); };
    obsThFactory["NEvpptaunu13_Bin4"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 4); };
    obsThFactory["NEvpptaunu13_Bin5"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 5); };
    obsThFactory["NEvpptaunu13_Bin6"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 6); };
    obsThFactory["NEvpptaunu13_Bin7"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 7); };
    obsThFactory["NEvpptaunu13_Bin8"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 8); };
    obsThFactory["NEvpptaunu13_Bin9"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 9); };
    obsThFactory["NEvpptaunu13_Bin10"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 10); };
    obsThFactory["NEvpptaunu13_Bin11"] = [=](const StandardModel& SM) { return new NevLHCtaunu13(SM, 11); };
    //
}
