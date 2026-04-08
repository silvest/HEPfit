/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "HiggsThObservables.h"
using namespace ThObsConst;

void ThObsFactory::registerHiggsObservables()
{
    //-----  Higgs observables  ----------

    //-----  Production cross sections (ratios with SM)  ----------
    obsThFactory["ggH"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC8); };
    obsThFactory["VBF"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC8); };
    obsThFactory["WH"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC8); };
    obsThFactory["ZH"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC8); };
    obsThFactory["VH"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC8); };
    obsThFactory["ggH+ttH"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_LHC8); };
    obsThFactory["VBF+VH"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_LHC8); };
    obsThFactory["ttH"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC8); };
    obsThFactory["tHq"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC8); };
    //
    obsThFactory["ggH7"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC7); };
    obsThFactory["VBF7"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC7); };
    obsThFactory["WH7"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC7); };
    obsThFactory["ZH7"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC7); };
    obsThFactory["VH7"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC7); };
    obsThFactory["ttH7"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC7); };
    obsThFactory["tHq7"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC7); };
    //
    obsThFactory["ggH8"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC8); };
    obsThFactory["ggH+ttH8"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_LHC8); };
    obsThFactory["VBF8"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC8); };
    obsThFactory["VBF+VH8"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_LHC8); };
    obsThFactory["VBFgamma8"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_LHC8); };
    obsThFactory["VH8"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC8); };
    obsThFactory["WH8"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC8); };
    obsThFactory["ZH8"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC8); };
    obsThFactory["ttH8"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC8); };
    obsThFactory["tHq8"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC8); };
    //
    obsThFactory["ggH13"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC13); };
    obsThFactory["ggH+ttH13"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_LHC13); };
    obsThFactory["VBF13"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC13); };
    obsThFactory["VBF+VH13"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_LHC13); };
    obsThFactory["VBFgamma13"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_LHC13); };
    obsThFactory["VH13"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC13); };
    obsThFactory["WH13"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC13); };
    obsThFactory["ZH13"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC13); };
    obsThFactory["VHpT25013"] = [=](const StandardModel& SM) { return new muVHpT250(SM, sqrt_s_LHC13); };
    obsThFactory["WHpT25013"] = [=](const StandardModel& SM) { return new muWHpT250(SM, sqrt_s_LHC13); };
    obsThFactory["ZHpT25013"] = [=](const StandardModel& SM) { return new muZHpT250(SM, sqrt_s_LHC13); };
    obsThFactory["ttH13"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC13); };
    obsThFactory["tHq13"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC13); };
    //
    obsThFactory["ggH14"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC14); };
    obsThFactory["ggH+ttH14"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_LHC14); };
    obsThFactory["VBF14"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC14); };
    obsThFactory["VBF+VH14"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_LHC14); };
    obsThFactory["VBFgamma14"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_LHC14); };
    obsThFactory["VH14"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC14); };
    obsThFactory["WH14"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC14); };
    obsThFactory["ZH14"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC14); };
    obsThFactory["ttH14"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC14); };
    obsThFactory["tHq14"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC14); };
    //
    obsThFactory["ggH27"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_LHC27); };
    obsThFactory["ggH+ttH27"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_LHC27); };
    obsThFactory["VBF27"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_LHC27); };
    obsThFactory["VBF+VH27"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_LHC27); };
    obsThFactory["VBFgamma27"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_LHC27); };
    obsThFactory["VH27"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_LHC27); };
    obsThFactory["WH27"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_LHC27); };
    obsThFactory["ZH27"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_LHC27); };
    obsThFactory["ttH27"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_LHC27); };
    obsThFactory["tHq27"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_LHC27); };
    //
    obsThFactory["ggH50"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_FCC50); };
    obsThFactory["ggH+ttH50"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_FCC50); };
    obsThFactory["VBF50"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_FCC50); };
    obsThFactory["VBF+VH50"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_FCC50); };
    obsThFactory["VBFgamma50"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_FCC50); };
    obsThFactory["VH50"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_FCC50); };
    obsThFactory["WH50"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_FCC50); };
    obsThFactory["ZH50"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_FCC50); };
    obsThFactory["ttH50"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_FCC50); };
    obsThFactory["tHq50"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_FCC50); };
    //
    obsThFactory["ggH84"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_FCC84); };
    obsThFactory["ggH+ttH84"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_FCC84); };
    obsThFactory["VBF84"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_FCC84); };
    obsThFactory["VBF+VH84"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_FCC84); };
    obsThFactory["VBFgamma84"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_FCC84); };
    obsThFactory["VH84"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_FCC84); };
    obsThFactory["WH84"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_FCC84); };
    obsThFactory["ZH84"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_FCC84); };
    obsThFactory["ttH84"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_FCC84); };
    obsThFactory["tHq84"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_FCC84); };
    //
    obsThFactory["ggH100"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_FCC100); };
    obsThFactory["ggH+ttH100"] = [=](const StandardModel& SM) { return new muggHpttH(SM, sqrt_s_FCC100); };
    obsThFactory["VBF100"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_FCC100); };
    obsThFactory["VBF+VH100"] = [=](const StandardModel& SM) { return new muVBFpVH(SM, sqrt_s_FCC100); };
    obsThFactory["VBFgamma100"] = [=](const StandardModel& SM) { return new muVBFgamma(SM, sqrt_s_FCC100); };
    obsThFactory["VH100"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_FCC100); };
    obsThFactory["WH100"] = [=](const StandardModel& SM) { return new muWH(SM, sqrt_s_FCC100); };
    obsThFactory["ZH100"] = [=](const StandardModel& SM) { return new muZH(SM, sqrt_s_FCC100); };
    obsThFactory["ttH100"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_FCC100); };
    obsThFactory["tHq100"] = [=](const StandardModel& SM) { return new mutHq(SM, sqrt_s_FCC100); };
    //
    obsThFactory["ggH196"] = [=](const StandardModel& SM) { return new muggH(SM, sqrt_s_TeV); };
    obsThFactory["VBF196"] = [=](const StandardModel& SM) { return new muVBF(SM, sqrt_s_TeV); };
    obsThFactory["VH196"] = [=](const StandardModel& SM) { return new muVH(SM, sqrt_s_TeV); };
    obsThFactory["ttH196"] = [=](const StandardModel& SM) { return new muttH(SM, sqrt_s_TeV); };
    //
    // Parameters for inclusive Higgs observables at e+ e-
    const double sqrts_eetoH[8] = {240., 250., 345., 350., 360., 365., 500., 1000.};

    for (int i = 0; i < 8; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrts_eetoH[i]);

        // Unpolarized
        obsThFactory["eeZH_" + sqrt_s_str] = [=](const StandardModel& SM) { return new mueeZHGen(SM, sqrts_eetoH[i], 0., 0.); };

        // Polarized: Pe-: -80% Pe+: +30%
        obsThFactory["eeZH_m80p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new mueeZHGen(SM, sqrts_eetoH[i], -0.8, 0.3); };

        // Polarized: Pe-: 80% Pe+: -30%
        obsThFactory["eeZH_p80m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new mueeZHGen(SM, sqrts_eetoH[i], 0.8, -0.3); };
    }
    //
    obsThFactory["eeZH230"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeZH240"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeZH250"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeZH350"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeZH365"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeZH380"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeZH500"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeZH550"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeZH1000"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeZH1400"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeZH1500"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeZH3000"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mumuZH3000"] = [=](const StandardModel& SM) { return new mummZH(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumuZH10000"] = [=](const StandardModel& SM) { return new mummZH(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["eeZllH230"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeZllH240"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeZllH250"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeZllH350"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeZllH365"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeZllH380"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeZllH500"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeZllH550"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeZllH1000"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeZllH1400"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeZllH1500"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeZllH3000"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["eeZqqH230"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeZqqH240"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeZqqH250"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeZqqH350"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeZqqH365"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeZqqH380"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeZqqH500"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeZqqH550"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeZqqH1000"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeZqqH1400"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeZqqH1500"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeZqqH3000"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["eeZH250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeZH250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeZH250_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeZH250_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeZH350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeZH350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeZH350_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeZH350_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeZH365_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeZH365_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeZH365_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeZH365_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeZH380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeZH380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeZH380_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeZH380_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeZH500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeZH500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeZH500_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeZH500_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeZH550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeZH550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeZH550_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeZH550_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeZH1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeZH1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeZH1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeZH1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeZH1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeZH1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeZH1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeZH1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeZH1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeZH1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeZH1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeZH1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeZH1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeZH1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeZH3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeZH3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeZH3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeZH3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeZllH250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeZllH250_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeZllH250_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeZllH350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeZllH350_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeZllH350_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH365_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeZllH365_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeZllH365_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeZllH365_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeZllH380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeZllH380_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeZllH380_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeZllH500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeZllH500_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeZllH500_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeZllH550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeZllH550_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeZllH550_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeZllH1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeZllH1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeZllH1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeZllH1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeZllH1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeZllH1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeZllH1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeZllH1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeZllH1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeZllH1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeZllH1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeZllH3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeZllH3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeZllH3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeZllH3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZllH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeZqqH250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeZqqH250_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeZqqH250_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeZqqH350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeZqqH350_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeZqqH350_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH365_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeZqqH365_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeZqqH365_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeZqqH365_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeZqqH380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeZqqH380_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeZqqH380_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeZqqH500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeZqqH500_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeZqqH500_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeZqqH550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeZqqH550_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeZqqH550_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeZqqH1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeZqqH1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeZqqH1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeZqqH1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeZqqH1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeZqqH1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeZqqH1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeZqqH1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeZqqH1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeZqqH1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeZqqH1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeZqqH3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeZqqH3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeZqqH3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeZqqH3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZqqH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["aPsk250_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["aPsk250_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    //
    obsThFactory["aPsk350_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["aPsk350_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    //
    obsThFactory["aPsk500_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["aPsk500_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    //
    obsThFactory["bPsk250_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["bPsk250_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    //
    obsThFactory["bPsk350_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["bPsk350_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    //
    obsThFactory["bPsk500_p80_m30"] = [=](const StandardModel& SM) { return new aPsk(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["bPsk500_m80_p30"] = [=](const StandardModel& SM) { return new bPsk(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    //
    obsThFactory["eeWBF230"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeWBF240"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeWBF250"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeWBF350"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeWBF365"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeWBF380"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeWBF500"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeWBF550"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeWBF1000"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeWBF1400"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeWBF1500"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeWBF3000"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["eeWBF250_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeWBF250_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeWBF250_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeWBF250_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF350_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeWBF350_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeWBF350_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeWBF350_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF365_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeWBF365_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeWBF365_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeWBF365_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF380_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeWBF380_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeWBF380_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeWBF380_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF500_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeWBF500_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeWBF500_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeWBF500_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF550_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeWBF550_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeWBF550_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeWBF550_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeWBF1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeWBF1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeWBF1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeWBF1000_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeWBF1000_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeWBF1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeWBF1400_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeWBF1400_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeWBF1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeWBF1500_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeWBF1500_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeWBF3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeWBF3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeWBF3000_p80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeWBF3000_m80_0"] = [=](const StandardModel& SM) { return new mueeWBF(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv230"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeHvv240"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeHvv250"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeHvv350"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeHvv365"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeHvv380"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeHvv500"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeHvv550"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeHvv1000"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeHvv1400"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeHvv1500"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeHvv3000"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mumuHvv3000"] = [=](const StandardModel& SM) { return new mummHvv(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumuHvv10000"] = [=](const StandardModel& SM) { return new mummHvv(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["mumuHmumu3000"] = [=](const StandardModel& SM) { return new mummHmm(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumuHmumu10000"] = [=](const StandardModel& SM) { return new mummHmm(SM, sqrt_s_leptcoll_10000); };
     //
    obsThFactory["eeHvv250_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeHvv250_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeHvv250_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeHvv250_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv350_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeHvv350_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeHvv350_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeHvv350_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv365_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeHvv365_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeHvv365_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeHvv365_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv380_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeHvv380_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeHvv380_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeHvv380_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeHvv500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeHvv500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeHvv500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };  
    //
    obsThFactory["eeHvv550_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeHvv550_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeHvv550_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeHvv550_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeHvv1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeHvv1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeHvv1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeHvv1000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeHvv1000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeHvv1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeHvv1400_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeHvv1400_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeHvv1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeHvv1500_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeHvv1500_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeHvv3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeHvv3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeHvv3000_p80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeHvv3000_m80_0"] = [=](const StandardModel& SM) { return new mueeHvv(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF230"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_230, 0., 0.); };
    obsThFactory["eeZBF240"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_240, 0., 0.); };
    obsThFactory["eeZBF250"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_250, 0., 0.); };
    obsThFactory["eeZBF350"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_350, 0., 0.); };
    obsThFactory["eeZBF365"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_365, 0., 0.); };
    obsThFactory["eeZBF380"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_380, 0., 0.); };
    obsThFactory["eeZBF500"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eeZBF550"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eeZBF1000"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eeZBF1400"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eeZBF1500"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eeZBF3000"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["eeZBF250_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_250, pol_80, -pol_30); };
    obsThFactory["eeZBF250_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_250, -pol_80, pol_30); };
    obsThFactory["eeZBF250_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_250, pol_80, pol_0); };
    obsThFactory["eeZBF250_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_250, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF350_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_350, pol_80, -pol_30); };
    obsThFactory["eeZBF350_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_350, -pol_80, pol_30); };
    obsThFactory["eeZBF350_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_350, pol_80, pol_0); };
    obsThFactory["eeZBF350_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_350, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF365_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_365, pol_80, -pol_30); };
    obsThFactory["eeZBF365_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_365, -pol_80, pol_30); };
    obsThFactory["eeZBF365_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_365, pol_80, pol_0); };
    obsThFactory["eeZBF365_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_365, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF380_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_380, pol_80, -pol_30); };
    obsThFactory["eeZBF380_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_380, -pol_80, pol_30); };
    obsThFactory["eeZBF380_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_380, pol_80, pol_0); };
    obsThFactory["eeZBF380_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_380, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eeZBF500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eeZBF500_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eeZBF500_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF550_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eeZBF550_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eeZBF550_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eeZBF550_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF1000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eeZBF1000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eeZBF1000_p80_m20"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eeZBF1000_m80_p20"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eeZBF1000_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eeZBF1000_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF1400_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eeZBF1400_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eeZBF1400_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eeZBF1400_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF1500_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eeZBF1500_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eeZBF1500_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eeZBF1500_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eeZBF3000_p80_m30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eeZBF3000_m80_p30"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eeZBF3000_p80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eeZBF3000_m80_0"] = [=](const StandardModel& SM) { return new mueeZBF(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["eettH500"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_500, 0., 0.); };
    obsThFactory["eettH550"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_550, 0., 0.); };
    obsThFactory["eettH1000"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, 0., 0.); };
    obsThFactory["eettH1400"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1400, 0., 0.); };
    obsThFactory["eettH1500"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1500, 0., 0.); };
    obsThFactory["eettH3000"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_3000, 0., 0.); };
    //
    obsThFactory["mumuttH1500"] = [=](const StandardModel& SM) { return new mummttH(SM, sqrt_s_leptcoll_3000); };
    obsThFactory["mumuttH3000"] = [=](const StandardModel& SM) { return new mummttH(SM, sqrt_s_leptcoll_10000); };
    //
    obsThFactory["eettH500_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_500, pol_80, -pol_30); };
    obsThFactory["eettH500_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_500, -pol_80, pol_30); };
    obsThFactory["eettH500_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_500, pol_80, pol_0); };
    obsThFactory["eettH500_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_500, -pol_80, pol_0); };
    //
    obsThFactory["eettH550_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_550, pol_80, -pol_30); };
    obsThFactory["eettH550_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_550, -pol_80, pol_30); };
    obsThFactory["eettH550_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_550, pol_80, pol_0); };
    obsThFactory["eettH550_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_550, -pol_80, pol_0); };
    //
    obsThFactory["eettH1000_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_30); };
    obsThFactory["eettH1000_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_30); };
    obsThFactory["eettH1000_p80_m20"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, pol_80, -pol_20); };
    obsThFactory["eettH1000_m80_p20"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_20); };
    obsThFactory["eettH1000_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, pol_80, pol_0); };
    obsThFactory["eettH1000_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1000, -pol_80, pol_0); };
    //
    obsThFactory["eettH1400_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1400, pol_80, -pol_30); };
    obsThFactory["eettH1400_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_30); };
    obsThFactory["eettH1400_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1400, pol_80, pol_0); };
    obsThFactory["eettH1400_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1400, -pol_80, pol_0); };
    //
    obsThFactory["eettH1500_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1500, pol_80, -pol_30); };
    obsThFactory["eettH1500_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_30); };
    obsThFactory["eettH1500_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1500, pol_80, pol_0); };
    obsThFactory["eettH1500_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_1500, -pol_80, pol_0); };
    //
    obsThFactory["eettH3000_p80_m30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_3000, pol_80, -pol_30); };
    obsThFactory["eettH3000_m80_p30"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_30); };
    obsThFactory["eettH3000_p80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_3000, pol_80, pol_0); };
    obsThFactory["eettH3000_m80_0"] = [=](const StandardModel& SM) { return new mueettH(SM, sqrt_s_leptcoll_3000, -pol_80, pol_0); };
    //
    obsThFactory["mumuH125"] = [=](const StandardModel& SM) { return new mummH(SM, sqrt_s_leptcoll_125); };
    //
    obsThFactory["epWBF1200"] = [=](const StandardModel& SM) { return new muepWBF(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["epWBF1300"] = [=](const StandardModel& SM) { return new muepWBF(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["epWBF1800"] = [=](const StandardModel& SM) { return new muepWBF(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["epWBF3500"] = [=](const StandardModel& SM) { return new muepWBF(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["epWBF5000"] = [=](const StandardModel& SM) { return new muepWBF(SM, sqrt_s_FCCep_5); };
    //
    obsThFactory["epZBF1200"] = [=](const StandardModel& SM) { return new muepZBF(SM, sqrt_s_LHeC_1_2); };
    obsThFactory["epZBF1300"] = [=](const StandardModel& SM) { return new muepZBF(SM, sqrt_s_LHeC_1_3); };
    obsThFactory["epZBF1800"] = [=](const StandardModel& SM) { return new muepZBF(SM, sqrt_s_LHeC_1_8); };
    obsThFactory["epZBF3500"] = [=](const StandardModel& SM) { return new muepZBF(SM, sqrt_s_FCCep_3_5); };
    obsThFactory["epZBF5000"] = [=](const StandardModel& SM) { return new muepZBF(SM, sqrt_s_FCCep_5); };
    //-----  Decay width and Branching ratios (ratios with SM)  ----------
    obsThFactory["GammaHggRatio"] = [](const StandardModel& SM) { return new GammaHtoggRatio(SM); };
    obsThFactory["GammaHWWRatio"] = [](const StandardModel& SM) { return new GammaHtoWWRatio(SM); };
    obsThFactory["GammaHZZRatio"] = [](const StandardModel& SM) { return new GammaHtoZZRatio(SM); };
    obsThFactory["GammaHZgaRatio"] = [](const StandardModel& SM) { return new GammaHtoZgaRatio(SM); };
    obsThFactory["GammaHgagaRatio"] = [](const StandardModel& SM) { return new GammaHtogagaRatio(SM); };
    obsThFactory["GammaHmumuRatio"] = [](const StandardModel& SM) { return new GammaHtomumuRatio(SM); };
    obsThFactory["GammaHtautauRatio"] = [](const StandardModel& SM) { return new GammaHtotautauRatio(SM); };
    obsThFactory["GammaHssRatio"] = [](const StandardModel& SM) { return new GammaHtossRatio(SM); };
    obsThFactory["GammaHccRatio"] = [](const StandardModel& SM) { return new GammaHtoccRatio(SM); };
    obsThFactory["GammaHbbRatio"] = [](const StandardModel& SM) { return new GammaHtobbRatio(SM); };
    obsThFactory["GammaHRatio"] = [](const StandardModel& SM) { return new GammaHRatio(SM); };
    //
    obsThFactory["BrHinvisible"] = [](const StandardModel& SM) { return new BrHinvisible(SM); };
    obsThFactory["BrHinvisibleNP"] = [](const StandardModel& SM) { return new BrHinvisibleNP(SM); };
    obsThFactory["BrHexotic"] = [](const StandardModel& SM) { return new BrHexotic(SM); };
    obsThFactory["BrHvisRatio"] = [](const StandardModel& SM) { return new BrHtovisRatio(SM); };
    obsThFactory["BrHtoinvRatio"] = [](const StandardModel& SM) { return new BrHtoinvRatio(SM); };
    obsThFactory["BrHggRatio"] = [](const StandardModel& SM) { return new BrHtoggRatio(SM); };
    obsThFactory["BrHWWRatio"] = [](const StandardModel& SM) { return new BrHtoWWRatio(SM); };
    obsThFactory["BrHZZRatio"] = [](const StandardModel& SM) { return new BrHtoZZRatio(SM); };
    obsThFactory["BrHZgaRatio"] = [](const StandardModel& SM) { return new BrHtoZgaRatio(SM); };
    obsThFactory["BrHZgallRatio"] = [](const StandardModel& SM) { return new BrHtoZgallRatio(SM); };
    obsThFactory["BrHZgaeeRatio"] = [](const StandardModel& SM) { return new BrHtoZgaeeRatio(SM); };
    obsThFactory["BrHZgamumuRatio"] = [](const StandardModel& SM) { return new BrHtoZgamumuRatio(SM); };
    obsThFactory["BrHgagaRatio"] = [](const StandardModel& SM) { return new BrHtogagaRatio(SM); };
    obsThFactory["BrHmumuRatio"] = [](const StandardModel& SM) { return new BrHtomumuRatio(SM); };
    obsThFactory["BrHtautauRatio"] = [](const StandardModel& SM) { return new BrHtotautauRatio(SM); };
    obsThFactory["BrHccRatio"] = [](const StandardModel& SM) { return new BrHtoccRatio(SM); };
    obsThFactory["BrHbbRatio"] = [](const StandardModel& SM) { return new BrHtobbRatio(SM); };
    // Dedicated 4 lepton decays
    obsThFactory["BrHto2l2vRatio"] = [](const StandardModel& SM) { return new BrHto2l2vRatio(SM); };
    obsThFactory["BrHtoevmuvRatio"] = [](const StandardModel& SM) { return new BrHtoevmuvRatio(SM); };
    obsThFactory["BrHto2e2vRatio"] = [](const StandardModel& SM) { return new BrHto2e2vRatio(SM); };
    obsThFactory["BrHto2mu2vRatio"] = [](const StandardModel& SM) { return new BrHto2mu2vRatio(SM); };
    obsThFactory["BrHto4lRatio"] = [](const StandardModel& SM) { return new BrHto4lRatio(SM); };
    obsThFactory["BrHto4eRatio"] = [](const StandardModel& SM) { return new BrHto4eRatio(SM); };
    obsThFactory["BrHto4muRatio"] = [](const StandardModel& SM) { return new BrHto4muRatio(SM); };
    obsThFactory["BrHto2e2muRatio"] = [](const StandardModel& SM) { return new BrHto2e2muRatio(SM); };
    // Other dedicated (semi-)leptonic 4 fermion decays
    obsThFactory["BrHtolljjRatio"] = [](const StandardModel& SM) { return new BrHtolljjRatio(SM); };
    obsThFactory["BrHtolvjjRatio"] = [](const StandardModel& SM) { return new BrHtolvjjRatio(SM); };
    obsThFactory["BrHtolv_lvorjjRatio"] = [](const StandardModel& SM) { return new BrHtolv_lvorjjRatio(SM); };
    obsThFactory["BrHtoll_vvorjjRatio"] = [](const StandardModel& SM) { return new BrHtoll_vvorjjRatio(SM); };
    //-----  Ratios of BR (ratios with SM)  ----------
    obsThFactory["BrHtogaga_over_mumu_Ratio"] = [](const StandardModel& SM) { return new BrHtogaga_over_mumu_Ratio(SM); };
    obsThFactory["BrHtoZga_over_mumu_Ratio"] = [](const StandardModel& SM) { return new BrHtoZga_over_mumu_Ratio(SM); };
    obsThFactory["BrHtoZmumuga_over_mumu_Ratio"] = [](const StandardModel& SM) { return new BrHtoZmumuga_over_mumu_Ratio(SM); };
    obsThFactory["BrHtoZga_over_4mu_Ratio"] = [](const StandardModel& SM) { return new BrHtoZga_over_4mu_Ratio(SM); };
    obsThFactory["BrHtoZmumuga_over_4mu_Ratio"] = [](const StandardModel& SM) { return new BrHtoZmumuga_over_4mu_Ratio(SM); };    
    obsThFactory["BrHtogaga_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHtogaga_over_4l_Ratio(SM); };
    obsThFactory["BrHtobb_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHtobb_over_4l_Ratio(SM); };
    obsThFactory["BrHto2l2v_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHto2l2v_over_4l_Ratio(SM); };
    obsThFactory["BrHtotautau_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHtotautau_over_4l_Ratio(SM); };
    obsThFactory["BrHtogaga_over_2e2mu_Ratio"] = [](const StandardModel& SM) { return new BrHtogaga_over_2e2mu_Ratio(SM); };
    obsThFactory["BrHtoZga_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHtoZga_over_4l_Ratio(SM); };
    obsThFactory["BrHtomumu_over_4l_Ratio"] = [](const StandardModel& SM) { return new BrHtomumu_over_4l_Ratio(SM); };
    obsThFactory["BrHtomumu_over_4mu_Ratio"] = [](const StandardModel& SM) { return new BrHtomumu_over_4mu_Ratio(SM); };
    obsThFactory["BrHto4l_over_gaga_Ratio"] = [](const StandardModel& SM) { return new BrHto4l_over_gaga_Ratio(SM); };
    obsThFactory["BrHtoZga_over_gaga_Ratio"] = [](const StandardModel& SM) { return new BrHtoZga_over_gaga_Ratio(SM); };
    obsThFactory["BrHtomumu_over_gaga_Ratio"] = [](const StandardModel& SM) { return new BrHtomumu_over_gaga_Ratio(SM); };
    obsThFactory["BrHto2l2v_over_gaga_Ratio"] = [](const StandardModel& SM) { return new BrHto2l2v_over_gaga_Ratio(SM); };
    obsThFactory["BrHtobb_over_cc_Ratio"] = [](const StandardModel& SM) { return new BrHtobb_over_cc_Ratio(SM); };
    obsThFactory["BrHtogaga_over_gg_Ratio"] = [](const StandardModel& SM) { return new BrHtogaga_over_gg_Ratio(SM); };
    obsThFactory["BrHtogg_over_bb_Ratio"] = [](const StandardModel& SM) { return new BrHtogg_over_bb_Ratio(SM); };
    obsThFactory["BrHtogg_over_cc_Ratio"] = [](const StandardModel& SM) { return new BrHtogg_over_cc_Ratio(SM); };
    //-----  Special observables --------
    obsThFactory["muttHZbb_boost84"] = [=](const StandardModel& SM) { return new muttHZbbboost(SM, sqrt_s_FCC84); };    
    obsThFactory["muttHgagaZee_boost84"] = [=](const StandardModel& SM) { return new muttHgagaZeeboost(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muttHZbb_boost100"] = [=](const StandardModel& SM) { return new muttHZbbboost(SM, sqrt_s_FCC100); };    
    obsThFactory["muttHgagaZee_boost100"] = [=](const StandardModel& SM) { return new muttHgagaZeeboost(SM, sqrt_s_FCC100); };
    //
    obsThFactory["ggHgagaInt14"] = [=](const StandardModel& SM) { return new muggHgagaInt(SM, sqrt_s_LHC14); };
    //-----  Full Signal strengths per prod and decay: Hadron colliders  ----------
    obsThFactory["muggHgaga"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHgaga"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muVHgaga"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muttHgaga"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muggHZZ"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHZZ"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muVHZZ"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muttHZZ"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muggHWW"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHWW"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muVHWW"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muttHWW"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muggHtautau"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHtautau"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muVHtautau"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muttHtautau"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muggHbb"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHbb"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muVHbb"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muttHbb"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_LHC13); };
    // Indicating the energy explicit in the observable name
    obsThFactory["muggHgaga13"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHgaga13"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muZHgaga13"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muWHgaga13"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muVHgaga13"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muttHgaga13"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muggHZga13"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHZga13"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muZHZga13"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muWHZga13"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muVHZga13"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muttHZga13"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_LHC13); };
    
    obsThFactory["muggHpttHptHpbbH_HZga13"] = [=](const StandardModel& SM) { return new muggHpttHptHpbbH_HZga(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFpVH_HZga13"]    = [=](const StandardModel& SM) { return new muVBFpVH_HZga(SM, sqrt_s_LHC13); };
    
    
    obsThFactory["muggHZZ13"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHZZ13"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muZHZZ13"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muWHZZ13"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muVHZZ13"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muttHZZ13"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muggHZZ4l13"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHZZ4l13"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_LHC13); };
    obsThFactory["muZHZZ4l13"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_LHC13); };
    obsThFactory["muWHZZ4l13"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_LHC13); };
    obsThFactory["muVHZZ4l13"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_LHC13); };
    obsThFactory["muttHZZ4l13"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muggHWW13"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHWW13"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muZHWW13"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muWHWW13"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muVHWW13"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_LHC13); };
    obsThFactory["muttHWW13"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muggHWW2l2v13"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHWW2l2v13"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_LHC13); };
    obsThFactory["muZHWW2l2v13"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_LHC13); };
    obsThFactory["muWHWW2l2v13"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_LHC13); };
    obsThFactory["muVHWW2l2v13"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_LHC13); };
    obsThFactory["muttHWW2l2v13"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muttHVV13"] = [=](const StandardModel& SM) { return new muttHVV(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muggHmumu13"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHmumu13"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muZHmumu13"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muWHmumu13"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muVHmumu13"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muttHmumu13"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muggHtautau13"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHtautau13"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFpVHtautau13"] = [=](const StandardModel& SM) { return new muVBFpVHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muZHtautau13"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muWHtautau13"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muVHtautau13"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muttHtautau13"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_LHC13); };
    obsThFactory["muggHbb13"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFHbb13"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muZHbb13"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muWHbb13"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muVHbb13"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muttHbb13"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_LHC13); };
    obsThFactory["muVHcc13"] = [=](const StandardModel& SM) { return new muVHcc(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muVBFBRinv13"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_LHC13); };
    obsThFactory["muVHBRinv13"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_LHC13); };
    //
    obsThFactory["muVBFHinv13"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_LHC13); };
    obsThFactory["muVHinv13"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_LHC13); };
    //
    //AG:begin
    obsThFactory["ggHgaga8"]   = [=](const StandardModel& SM) { return new ggHgaga(SM, sqrt_s_LHC8); };
    obsThFactory["ggHZZ8"]     = [=](const StandardModel& SM) { return new ggHZZ(SM, sqrt_s_LHC8); };
    obsThFactory["ggHWW8"]     = [=](const StandardModel& SM) { return new ggHWW(SM, sqrt_s_LHC8); };
    obsThFactory["ggHtautau8"] = [=](const StandardModel& SM) { return new ggHtautau(SM, sqrt_s_LHC8); };
    obsThFactory["VBFHgaga8"]  = [=](const StandardModel& SM) { return new VBFHgaga(SM, sqrt_s_LHC8); };
    obsThFactory["VBFHZZ8"]    = [=](const StandardModel& SM) { return new VBFHZZ(SM, sqrt_s_LHC8); };
    obsThFactory["VBFHWW8"]    = [=](const StandardModel& SM) { return new VBFHWW(SM, sqrt_s_LHC8); };
    obsThFactory["VBFHtautau8"]= [=](const StandardModel& SM) { return new VBFHtautau(SM, sqrt_s_LHC8); };
    obsThFactory["WHgaga8"]    = [=](const StandardModel& SM) { return new WHgaga(SM, sqrt_s_LHC8); };
    obsThFactory["WHWW8"]      = [=](const StandardModel& SM) { return new WHWW(SM, sqrt_s_LHC8); };
    obsThFactory["WHtautau8"]  = [=](const StandardModel& SM) { return new WHtautau(SM, sqrt_s_LHC8); };
    obsThFactory["WHbb8"]      = [=](const StandardModel& SM) { return new WHbb(SM, sqrt_s_LHC8); };
    obsThFactory["ZHgaga8"]    = [=](const StandardModel& SM) { return new ZHgaga(SM, sqrt_s_LHC8); };
    obsThFactory["ZHWW8"]      = [=](const StandardModel& SM) { return new ZHWW(SM, sqrt_s_LHC8); };
    obsThFactory["ZHtautau8"]  = [=](const StandardModel& SM) { return new ZHtautau(SM, sqrt_s_LHC8); };
    obsThFactory["ZHbb8"]      = [=](const StandardModel& SM) { return new ZHbb(SM, sqrt_s_LHC8); };
    obsThFactory["ttHgaga8"]   = [=](const StandardModel& SM) { return new ttHgaga(SM, sqrt_s_LHC8); };
    obsThFactory["ttHWW8"]     = [=](const StandardModel& SM) { return new ttHWW(SM, sqrt_s_LHC8); };
    obsThFactory["ttHtautau8"] = [=](const StandardModel& SM) { return new ttHtautau(SM, sqrt_s_LHC8); };
    obsThFactory["ttHbb8"]     = [=](const StandardModel& SM) { return new ttHbb(SM, sqrt_s_LHC8); };

    obsThFactory["muggHpVBFpbbH_Hbb13"] = [=](const StandardModel& SM) { return new muggHpVBFpbbH_Hbb(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_Hbb13"]      = [=](const StandardModel& SM) { return new muttHptH_Hbb(SM, sqrt_s_LHC13); };
    obsThFactory["muggHpbbH_HWW13"]     = [=](const StandardModel& SM) { return new muggHpbbH_HWW(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_HWW13"]      = [=](const StandardModel& SM) { return new muttHptH_HWW(SM, sqrt_s_LHC13); };
    obsThFactory["muggHpbbH_Htautau13"] = [=](const StandardModel& SM) { return new muggHpbbH_Htautau(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_Htautau13"]  = [=](const StandardModel& SM) { return new muttHptH_Htautau(SM, sqrt_s_LHC13); };
    obsThFactory["muggHpbbH_HZZ13"]     = [=](const StandardModel& SM) { return new muggHpbbH_HZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_HZZ13"]      = [=](const StandardModel& SM) { return new muttHptH_HZZ(SM, sqrt_s_LHC13); };
    obsThFactory["muggHpbbH_Hgaga13"]   = [=](const StandardModel& SM) { return new muggHpbbH_Hgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muggHpttHptHpbbH_Hmumu13"] = [=](const StandardModel& SM) { return new muggHpttHptHpbbH_Hmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muVBFpVH_Hmumu13"]    = [=](const StandardModel& SM) { return new muVBFpVH_Hmumu(SM, sqrt_s_LHC13); };
    obsThFactory["mutHgaga13"] = [=](const StandardModel& SM) { return new mutHgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_Hgaga13"]      = [=](const StandardModel& SM) { return new muttHptH_Hgaga(SM, sqrt_s_LHC13); };
    obsThFactory["muttHptH_Hmumu13"]      = [=](const StandardModel& SM) { return new muttHptH_Hmumu(SM, sqrt_s_LHC13); };
    //AG:end
    //
    obsThFactory["muggHgaga14"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHgaga14"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muZHgaga14"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muWHgaga14"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muVHgaga14"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muttHgaga14"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muggHZga14"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHZga14"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muZHZga14"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muWHZga14"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muVHZga14"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muttHZga14"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muggHZZ14"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHZZ14"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muZHZZ14"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muWHZZ14"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muVHZZ14"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muttHZZ14"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muggHZZ4l14"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHZZ4l14"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muZHZZ4l14"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muWHZZ4l14"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muVHZZ4l14"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muttHZZ4l14"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muggHWW14"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHWW14"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muZHWW14"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muWHWW14"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muVHWW14"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muttHWW14"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muggHWW2l2v14"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHWW2l2v14"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muZHWW2l2v14"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muWHWW2l2v14"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muVHWW2l2v14"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muttHWW2l2v14"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muggHmumu14"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHmumu14"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muZHmumu14"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muWHmumu14"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muVHmumu14"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muttHmumu14"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muggHtautau14"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHtautau14"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muZHtautau14"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muWHtautau14"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muVHtautau14"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muttHtautau14"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muggHbb14"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muVBFHbb14"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muZHbb14"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muWHbb14"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muVHbb14"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muttHbb14"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muVBFBRinv14"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_LHC14); };
    obsThFactory["muVHBRinv14"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muVBFHinv14"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_LHC14); };
    obsThFactory["muVHinv14"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muggHgaga27"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHgaga27"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muZHgaga27"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muWHgaga27"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muVHgaga27"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muttHgaga27"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muggHZga27"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHZga27"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muZHZga27"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muWHZga27"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muVHZga27"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muttHZga27"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muggHZZ27"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHZZ27"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muZHZZ27"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muWHZZ27"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muVHZZ27"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muttHZZ27"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muggHZZ4l27"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHZZ4l27"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muZHZZ4l27"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muWHZZ4l27"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muVHZZ4l27"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muttHZZ4l27"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muggHWW27"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHWW27"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muZHWW27"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muWHWW27"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muVHWW27"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muttHWW27"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muggHWW2l2v27"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHWW2l2v27"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muZHWW2l2v27"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muWHWW2l2v27"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muVHWW2l2v27"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muttHWW2l2v27"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muggHmumu27"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHmumu27"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muZHmumu27"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muWHmumu27"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muVHmumu27"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muttHmumu27"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muggHtautau27"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHtautau27"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muZHtautau27"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muWHtautau27"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muVHtautau27"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muttHtautau27"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muggHbb27"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muVBFHbb27"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muZHbb27"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muWHbb27"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muVHbb27"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muttHbb27"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muVBFBRinv27"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_LHC27); };
    obsThFactory["muVHBRinv27"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muVBFHinv27"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_LHC27); };
    obsThFactory["muVHinv27"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muggHgaga50"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHgaga50"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muZHgaga50"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muWHgaga50"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muVHgaga50"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muttHgaga50"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muggHZga50"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muggHZgamumu50"] = [=](const StandardModel& SM) { return new muggHZgamumu(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHZga50"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muZHZga50"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muWHZga50"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muVHZga50"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muttHZga50"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muggHZZ50"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHZZ50"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muZHZZ50"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muWHZZ50"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muVHZZ50"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muttHZZ50"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muggHZZ4l50"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muggHZZ4mu50"] = [=](const StandardModel& SM) { return new muggHZZ4mu(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHZZ4l50"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muZHZZ4l50"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muWHZZ4l50"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muVHZZ4l50"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muttHZZ4l50"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muggHWW50"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHWW50"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muZHWW50"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muWHWW50"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muVHWW50"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muttHWW50"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muggHWW2l2v50"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHWW2l2v50"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muZHWW2l2v50"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muWHWW2l2v50"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muVHWW2l2v50"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muttHWW2l2v50"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muggHmumu50"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHmumu50"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muZHmumu50"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muWHmumu50"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muVHmumu50"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muttHmumu50"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muggHtautau50"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHtautau50"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muZHtautau50"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muWHtautau50"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muVHtautau50"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muttHtautau50"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muggHbb50"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muVBFHbb50"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muZHbb50"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muWHbb50"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muVHbb50"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muttHbb50"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muVBFBRinv50"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_FCC50); };
    obsThFactory["muVHBRinv50"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muVBFHinv50"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_FCC50); };
    obsThFactory["muVHinv50"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muggHgaga84"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHgaga84"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muZHgaga84"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muWHgaga84"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muVHgaga84"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muttHgaga84"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muggHZga84"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muggHZgamumu84"] = [=](const StandardModel& SM) { return new muggHZgamumu(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHZga84"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muZHZga84"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muWHZga84"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muVHZga84"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muttHZga84"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muggHZZ84"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHZZ84"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muZHZZ84"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muWHZZ84"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muVHZZ84"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muttHZZ84"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muggHZZ4l84"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muggHZZ4mu84"] = [=](const StandardModel& SM) { return new muggHZZ4mu(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHZZ4l84"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muZHZZ4l84"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muWHZZ4l84"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muVHZZ4l84"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muttHZZ4l84"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muggHWW84"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHWW84"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muZHWW84"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muWHWW84"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muVHWW84"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muttHWW84"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muggHWW2l2v84"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHWW2l2v84"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muZHWW2l2v84"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muWHWW2l2v84"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muVHWW2l2v84"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muttHWW2l2v84"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muggHmumu84"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHmumu84"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muZHmumu84"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muWHmumu84"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muVHmumu84"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muttHmumu84"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muggHtautau84"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHtautau84"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muZHtautau84"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muWHtautau84"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muVHtautau84"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muttHtautau84"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muggHbb84"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muVBFHbb84"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muZHbb84"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muWHbb84"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muVHbb84"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muttHbb84"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muVBFBRinv84"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_FCC84); };
    obsThFactory["muVHBRinv84"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muVBFHinv84"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_FCC84); };
    obsThFactory["muVHinv84"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muggHgaga100"] = [=](const StandardModel& SM) { return new muggHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHgaga100"] = [=](const StandardModel& SM) { return new muVBFHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muZHgaga100"] = [=](const StandardModel& SM) { return new muZHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muWHgaga100"] = [=](const StandardModel& SM) { return new muWHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muVHgaga100"] = [=](const StandardModel& SM) { return new muVHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muttHgaga100"] = [=](const StandardModel& SM) { return new muttHgaga(SM, sqrt_s_FCC100); };
    obsThFactory["muggHZga100"] = [=](const StandardModel& SM) { return new muggHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muggHZgamumu100"] = [=](const StandardModel& SM) { return new muggHZgamumu(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHZga100"] = [=](const StandardModel& SM) { return new muVBFHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muZHZga100"] = [=](const StandardModel& SM) { return new muZHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muWHZga100"] = [=](const StandardModel& SM) { return new muWHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muVHZga100"] = [=](const StandardModel& SM) { return new muVHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muttHZga100"] = [=](const StandardModel& SM) { return new muttHZga(SM, sqrt_s_FCC100); };
    obsThFactory["muggHZZ100"] = [=](const StandardModel& SM) { return new muggHZZ(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHZZ100"] = [=](const StandardModel& SM) { return new muVBFHZZ(SM, sqrt_s_FCC100); };
    obsThFactory["muZHZZ100"] = [=](const StandardModel& SM) { return new muZHZZ(SM, sqrt_s_FCC100); };
    obsThFactory["muWHZZ100"] = [=](const StandardModel& SM) { return new muWHZZ(SM, sqrt_s_FCC100); };
    obsThFactory["muVHZZ100"] = [=](const StandardModel& SM) { return new muVHZZ(SM, sqrt_s_FCC100); };
    obsThFactory["muttHZZ100"] = [=](const StandardModel& SM) { return new muttHZZ(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muggHZZ4l100"] = [=](const StandardModel& SM) { return new muggHZZ4l(SM, sqrt_s_FCC100); };
    obsThFactory["muggHZZ4mu100"] = [=](const StandardModel& SM) { return new muggHZZ4mu(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHZZ4l100"] = [=](const StandardModel& SM) { return new muVBFHZZ4l(SM, sqrt_s_FCC100); };
    obsThFactory["muZHZZ4l100"] = [=](const StandardModel& SM) { return new muZHZZ4l(SM, sqrt_s_FCC100); };
    obsThFactory["muWHZZ4l100"] = [=](const StandardModel& SM) { return new muWHZZ4l(SM, sqrt_s_FCC100); };
    obsThFactory["muVHZZ4l100"] = [=](const StandardModel& SM) { return new muVHZZ4l(SM, sqrt_s_FCC100); };
    obsThFactory["muttHZZ4l100"] = [=](const StandardModel& SM) { return new muttHZZ4l(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muggHWW100"] = [=](const StandardModel& SM) { return new muggHWW(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHWW100"] = [=](const StandardModel& SM) { return new muVBFHWW(SM, sqrt_s_FCC100); };
    obsThFactory["muZHWW100"] = [=](const StandardModel& SM) { return new muZHWW(SM, sqrt_s_FCC100); };
    obsThFactory["muWHWW100"] = [=](const StandardModel& SM) { return new muWHWW(SM, sqrt_s_FCC100); };
    obsThFactory["muVHWW100"] = [=](const StandardModel& SM) { return new muVHWW(SM, sqrt_s_FCC100); };
    obsThFactory["muttHWW100"] = [=](const StandardModel& SM) { return new muttHWW(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muggHWW2l2v100"] = [=](const StandardModel& SM) { return new muggHWW2l2v(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHWW2l2v100"] = [=](const StandardModel& SM) { return new muVBFHWW2l2v(SM, sqrt_s_FCC100); };
    obsThFactory["muZHWW2l2v100"] = [=](const StandardModel& SM) { return new muZHWW2l2v(SM, sqrt_s_FCC100); };
    obsThFactory["muWHWW2l2v100"] = [=](const StandardModel& SM) { return new muWHWW2l2v(SM, sqrt_s_FCC100); };
    obsThFactory["muVHWW2l2v100"] = [=](const StandardModel& SM) { return new muVHWW2l2v(SM, sqrt_s_FCC100); };
    obsThFactory["muttHWW2l2v100"] = [=](const StandardModel& SM) { return new muttHWW2l2v(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muggHmumu100"] = [=](const StandardModel& SM) { return new muggHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHmumu100"] = [=](const StandardModel& SM) { return new muVBFHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muZHmumu100"] = [=](const StandardModel& SM) { return new muZHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muWHmumu100"] = [=](const StandardModel& SM) { return new muWHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muVHmumu100"] = [=](const StandardModel& SM) { return new muVHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muttHmumu100"] = [=](const StandardModel& SM) { return new muttHmumu(SM, sqrt_s_FCC100); };
    obsThFactory["muggHtautau100"] = [=](const StandardModel& SM) { return new muggHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHtautau100"] = [=](const StandardModel& SM) { return new muVBFHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muZHtautau100"] = [=](const StandardModel& SM) { return new muZHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muWHtautau100"] = [=](const StandardModel& SM) { return new muWHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muVHtautau100"] = [=](const StandardModel& SM) { return new muVHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muttHtautau100"] = [=](const StandardModel& SM) { return new muttHtautau(SM, sqrt_s_FCC100); };
    obsThFactory["muggHbb100"] = [=](const StandardModel& SM) { return new muggHbb(SM, sqrt_s_FCC100); };
    obsThFactory["muVBFHbb100"] = [=](const StandardModel& SM) { return new muVBFHbb(SM, sqrt_s_FCC100); };
    obsThFactory["muZHbb100"] = [=](const StandardModel& SM) { return new muZHbb(SM, sqrt_s_FCC100); };
    obsThFactory["muWHbb100"] = [=](const StandardModel& SM) { return new muWHbb(SM, sqrt_s_FCC100); };
    obsThFactory["muVHbb100"] = [=](const StandardModel& SM) { return new muVHbb(SM, sqrt_s_FCC100); };
    obsThFactory["muttHbb100"] = [=](const StandardModel& SM) { return new muttHbb(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muVBFBRinv100"] = [=](const StandardModel& SM) { return new muVBFBRinv(SM, sqrt_s_FCC100); };
    obsThFactory["muVHBRinv100"] = [=](const StandardModel& SM) { return new muVHBRinv(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muVBFHinv100"] = [=](const StandardModel& SM) { return new muVBFHinv(SM, sqrt_s_FCC100); };
    obsThFactory["muVHinv100"] = [=](const StandardModel& SM) { return new muVHinv(SM, sqrt_s_FCC100); };
    //
    obsThFactory["muppHmumu8"] = [=](const StandardModel& SM) { return new muppHmumu(SM, sqrt_s_LHC8); };
    obsThFactory["muppHmumu13"] = [=](const StandardModel& SM) { return new muppHmumu(SM, sqrt_s_LHC13); };
    obsThFactory["muppHZga8"] = [=](const StandardModel& SM) { return new muppHZga(SM, sqrt_s_LHC8); };
    obsThFactory["muppHZga13"] = [=](const StandardModel& SM) { return new muppHZga(SM, sqrt_s_LHC13); };
    obsThFactory["muggHH2ga2b14"] = [=](const StandardModel& SM) { return new muggHH2ga2b(SM, sqrt_s_LHC14); };
    obsThFactory["muggHH2ga2b50"] = [=](const StandardModel& SM) { return new muggHH2ga2b(SM, sqrt_s_FCC50); };
    obsThFactory["muggHH2ga2b84"] = [=](const StandardModel& SM) { return new muggHH2ga2b(SM, sqrt_s_FCC84); };
    obsThFactory["muggHH2ga2b100"] = [=](const StandardModel& SM) { return new muggHH2ga2b(SM, sqrt_s_FCC100); };
    //
    // Special version of the H signal strength at Hadron collider with separate theory uncertainty per prod x decay channel
    // Only for 14, 27, 50 and 84 TeV
    //
    obsThFactory["muTHUggHgaga14"] = [=](const StandardModel& SM) { return new muTHUggHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHgaga14"] = [=](const StandardModel& SM) { return new muTHUVBFHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHgaga14"] = [=](const StandardModel& SM) { return new muTHUZHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHgaga14"] = [=](const StandardModel& SM) { return new muTHUWHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHgaga14"] = [=](const StandardModel& SM) { return new muTHUVHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHgaga14"] = [=](const StandardModel& SM) { return new muTHUttHgaga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUggHZga14"] = [=](const StandardModel& SM) { return new muTHUggHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHZga14"] = [=](const StandardModel& SM) { return new muTHUVBFHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHZga14"] = [=](const StandardModel& SM) { return new muTHUZHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHZga14"] = [=](const StandardModel& SM) { return new muTHUWHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHZga14"] = [=](const StandardModel& SM) { return new muTHUVHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHZga14"] = [=](const StandardModel& SM) { return new muTHUttHZga(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUggHZZ14"] = [=](const StandardModel& SM) { return new muTHUggHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHZZ14"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHZZ14"] = [=](const StandardModel& SM) { return new muTHUZHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHZZ14"] = [=](const StandardModel& SM) { return new muTHUWHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHZZ14"] = [=](const StandardModel& SM) { return new muTHUVHZZ(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHZZ14"] = [=](const StandardModel& SM) { return new muTHUttHZZ(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUggHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUggHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUZHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUWHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUVHZZ4l(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHZZ4l14"] = [=](const StandardModel& SM) { return new muTHUttHZZ4l(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUggHWW14"] = [=](const StandardModel& SM) { return new muTHUggHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHWW14"] = [=](const StandardModel& SM) { return new muTHUVBFHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHWW14"] = [=](const StandardModel& SM) { return new muTHUZHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHWW14"] = [=](const StandardModel& SM) { return new muTHUWHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHWW14"] = [=](const StandardModel& SM) { return new muTHUVHWW(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHWW14"] = [=](const StandardModel& SM) { return new muTHUttHWW(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUggHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUggHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUVBFHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUZHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUWHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUVHWW2l2v(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHWW2l2v14"] = [=](const StandardModel& SM) { return new muTHUttHWW2l2v(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUggHmumu14"] = [=](const StandardModel& SM) { return new muTHUggHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHmumu14"] = [=](const StandardModel& SM) { return new muTHUVBFHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHmumu14"] = [=](const StandardModel& SM) { return new muTHUZHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHmumu14"] = [=](const StandardModel& SM) { return new muTHUWHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHmumu14"] = [=](const StandardModel& SM) { return new muTHUVHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHmumu14"] = [=](const StandardModel& SM) { return new muTHUttHmumu(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUggHtautau14"] = [=](const StandardModel& SM) { return new muTHUggHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHtautau14"] = [=](const StandardModel& SM) { return new muTHUVBFHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHtautau14"] = [=](const StandardModel& SM) { return new muTHUZHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHtautau14"] = [=](const StandardModel& SM) { return new muTHUWHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHtautau14"] = [=](const StandardModel& SM) { return new muTHUVHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHtautau14"] = [=](const StandardModel& SM) { return new muTHUttHtautau(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUggHbb14"] = [=](const StandardModel& SM) { return new muTHUggHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVBFHbb14"] = [=](const StandardModel& SM) { return new muTHUVBFHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUZHbb14"] = [=](const StandardModel& SM) { return new muTHUZHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUWHbb14"] = [=](const StandardModel& SM) { return new muTHUWHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHbb14"] = [=](const StandardModel& SM) { return new muTHUVHbb(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUttHbb14"] = [=](const StandardModel& SM) { return new muTHUttHbb(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUVBFBRinv14"] = [=](const StandardModel& SM) { return new muTHUVBFBRinv(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHBRinv14"] = [=](const StandardModel& SM) { return new muTHUVHBRinv(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUVBFHinv14"] = [=](const StandardModel& SM) { return new muTHUVBFHinv(SM, sqrt_s_LHC14); };
    obsThFactory["muTHUVHinv14"] = [=](const StandardModel& SM) { return new muTHUVHinv(SM, sqrt_s_LHC14); };
    //
    obsThFactory["muTHUggHgaga27"] = [=](const StandardModel& SM) { return new muTHUggHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHgaga27"] = [=](const StandardModel& SM) { return new muTHUVBFHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHgaga27"] = [=](const StandardModel& SM) { return new muTHUZHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHgaga27"] = [=](const StandardModel& SM) { return new muTHUWHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHgaga27"] = [=](const StandardModel& SM) { return new muTHUVHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHgaga27"] = [=](const StandardModel& SM) { return new muTHUttHgaga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUggHZga27"] = [=](const StandardModel& SM) { return new muTHUggHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHZga27"] = [=](const StandardModel& SM) { return new muTHUVBFHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHZga27"] = [=](const StandardModel& SM) { return new muTHUZHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHZga27"] = [=](const StandardModel& SM) { return new muTHUWHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHZga27"] = [=](const StandardModel& SM) { return new muTHUVHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHZga27"] = [=](const StandardModel& SM) { return new muTHUttHZga(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUggHZZ27"] = [=](const StandardModel& SM) { return new muTHUggHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHZZ27"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHZZ27"] = [=](const StandardModel& SM) { return new muTHUZHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHZZ27"] = [=](const StandardModel& SM) { return new muTHUWHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHZZ27"] = [=](const StandardModel& SM) { return new muTHUVHZZ(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHZZ27"] = [=](const StandardModel& SM) { return new muTHUttHZZ(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUggHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUggHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUZHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUWHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUVHZZ4l(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHZZ4l27"] = [=](const StandardModel& SM) { return new muTHUttHZZ4l(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUggHWW27"] = [=](const StandardModel& SM) { return new muTHUggHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHWW27"] = [=](const StandardModel& SM) { return new muTHUVBFHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHWW27"] = [=](const StandardModel& SM) { return new muTHUZHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHWW27"] = [=](const StandardModel& SM) { return new muTHUWHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHWW27"] = [=](const StandardModel& SM) { return new muTHUVHWW(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHWW27"] = [=](const StandardModel& SM) { return new muTHUttHWW(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUggHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUggHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUVBFHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUZHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUWHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUVHWW2l2v(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHWW2l2v27"] = [=](const StandardModel& SM) { return new muTHUttHWW2l2v(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUggHmumu27"] = [=](const StandardModel& SM) { return new muTHUggHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHmumu27"] = [=](const StandardModel& SM) { return new muTHUVBFHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHmumu27"] = [=](const StandardModel& SM) { return new muTHUZHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHmumu27"] = [=](const StandardModel& SM) { return new muTHUWHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHmumu27"] = [=](const StandardModel& SM) { return new muTHUVHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHmumu27"] = [=](const StandardModel& SM) { return new muTHUttHmumu(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUggHtautau27"] = [=](const StandardModel& SM) { return new muTHUggHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHtautau27"] = [=](const StandardModel& SM) { return new muTHUVBFHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHtautau27"] = [=](const StandardModel& SM) { return new muTHUZHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHtautau27"] = [=](const StandardModel& SM) { return new muTHUWHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHtautau27"] = [=](const StandardModel& SM) { return new muTHUVHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHtautau27"] = [=](const StandardModel& SM) { return new muTHUttHtautau(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUggHbb27"] = [=](const StandardModel& SM) { return new muTHUggHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVBFHbb27"] = [=](const StandardModel& SM) { return new muTHUVBFHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUZHbb27"] = [=](const StandardModel& SM) { return new muTHUZHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUWHbb27"] = [=](const StandardModel& SM) { return new muTHUWHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHbb27"] = [=](const StandardModel& SM) { return new muTHUVHbb(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUttHbb27"] = [=](const StandardModel& SM) { return new muTHUttHbb(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUVBFBRinv27"] = [=](const StandardModel& SM) { return new muTHUVBFBRinv(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHBRinv27"] = [=](const StandardModel& SM) { return new muTHUVHBRinv(SM, sqrt_s_LHC27); };
    //
    obsThFactory["muTHUVBFHinv27"] = [=](const StandardModel& SM) { return new muTHUVBFHinv(SM, sqrt_s_LHC27); };
    obsThFactory["muTHUVHinv27"] = [=](const StandardModel& SM) { return new muTHUVHinv(SM, sqrt_s_LHC27); };
    //
    //
    obsThFactory["muTHUggHgaga50"] = [=](const StandardModel& SM) { return new muTHUggHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHgaga50"] = [=](const StandardModel& SM) { return new muTHUVBFHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHgaga50"] = [=](const StandardModel& SM) { return new muTHUZHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHgaga50"] = [=](const StandardModel& SM) { return new muTHUWHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHgaga50"] = [=](const StandardModel& SM) { return new muTHUVHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHgaga50"] = [=](const StandardModel& SM) { return new muTHUttHgaga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUggHZga50"] = [=](const StandardModel& SM) { return new muTHUggHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHZga50"] = [=](const StandardModel& SM) { return new muTHUVBFHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHZga50"] = [=](const StandardModel& SM) { return new muTHUZHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHZga50"] = [=](const StandardModel& SM) { return new muTHUWHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHZga50"] = [=](const StandardModel& SM) { return new muTHUVHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHZga50"] = [=](const StandardModel& SM) { return new muTHUttHZga(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUggHZZ50"] = [=](const StandardModel& SM) { return new muTHUggHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHZZ50"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHZZ50"] = [=](const StandardModel& SM) { return new muTHUZHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHZZ50"] = [=](const StandardModel& SM) { return new muTHUWHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHZZ50"] = [=](const StandardModel& SM) { return new muTHUVHZZ(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHZZ50"] = [=](const StandardModel& SM) { return new muTHUttHZZ(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUggHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUggHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUZHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUWHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUVHZZ4l(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHZZ4l50"] = [=](const StandardModel& SM) { return new muTHUttHZZ4l(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUggHWW50"] = [=](const StandardModel& SM) { return new muTHUggHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHWW50"] = [=](const StandardModel& SM) { return new muTHUVBFHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHWW50"] = [=](const StandardModel& SM) { return new muTHUZHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHWW50"] = [=](const StandardModel& SM) { return new muTHUWHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHWW50"] = [=](const StandardModel& SM) { return new muTHUVHWW(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHWW50"] = [=](const StandardModel& SM) { return new muTHUttHWW(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUggHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUggHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUVBFHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUZHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUWHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUVHWW2l2v(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHWW2l2v50"] = [=](const StandardModel& SM) { return new muTHUttHWW2l2v(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUggHmumu50"] = [=](const StandardModel& SM) { return new muTHUggHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHmumu50"] = [=](const StandardModel& SM) { return new muTHUVBFHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHmumu50"] = [=](const StandardModel& SM) { return new muTHUZHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHmumu50"] = [=](const StandardModel& SM) { return new muTHUWHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHmumu50"] = [=](const StandardModel& SM) { return new muTHUVHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHmumu50"] = [=](const StandardModel& SM) { return new muTHUttHmumu(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUggHtautau50"] = [=](const StandardModel& SM) { return new muTHUggHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHtautau50"] = [=](const StandardModel& SM) { return new muTHUVBFHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHtautau50"] = [=](const StandardModel& SM) { return new muTHUZHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHtautau50"] = [=](const StandardModel& SM) { return new muTHUWHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHtautau50"] = [=](const StandardModel& SM) { return new muTHUVHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHtautau50"] = [=](const StandardModel& SM) { return new muTHUttHtautau(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUggHbb50"] = [=](const StandardModel& SM) { return new muTHUggHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVBFHbb50"] = [=](const StandardModel& SM) { return new muTHUVBFHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUZHbb50"] = [=](const StandardModel& SM) { return new muTHUZHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUWHbb50"] = [=](const StandardModel& SM) { return new muTHUWHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHbb50"] = [=](const StandardModel& SM) { return new muTHUVHbb(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUttHbb50"] = [=](const StandardModel& SM) { return new muTHUttHbb(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUVBFBRinv50"] = [=](const StandardModel& SM) { return new muTHUVBFBRinv(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHBRinv50"] = [=](const StandardModel& SM) { return new muTHUVHBRinv(SM, sqrt_s_FCC50); };
    //
    obsThFactory["muTHUVBFHinv50"] = [=](const StandardModel& SM) { return new muTHUVBFHinv(SM, sqrt_s_FCC50); };
    obsThFactory["muTHUVHinv50"] = [=](const StandardModel& SM) { return new muTHUVHinv(SM, sqrt_s_FCC50); };
    //
    //
    obsThFactory["muTHUggHgaga84"] = [=](const StandardModel& SM) { return new muTHUggHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHgaga84"] = [=](const StandardModel& SM) { return new muTHUVBFHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHgaga84"] = [=](const StandardModel& SM) { return new muTHUZHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHgaga84"] = [=](const StandardModel& SM) { return new muTHUWHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHgaga84"] = [=](const StandardModel& SM) { return new muTHUVHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHgaga84"] = [=](const StandardModel& SM) { return new muTHUttHgaga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUggHZga84"] = [=](const StandardModel& SM) { return new muTHUggHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHZga84"] = [=](const StandardModel& SM) { return new muTHUVBFHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHZga84"] = [=](const StandardModel& SM) { return new muTHUZHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHZga84"] = [=](const StandardModel& SM) { return new muTHUWHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHZga84"] = [=](const StandardModel& SM) { return new muTHUVHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHZga84"] = [=](const StandardModel& SM) { return new muTHUttHZga(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUggHZZ84"] = [=](const StandardModel& SM) { return new muTHUggHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHZZ84"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHZZ84"] = [=](const StandardModel& SM) { return new muTHUZHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHZZ84"] = [=](const StandardModel& SM) { return new muTHUWHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHZZ84"] = [=](const StandardModel& SM) { return new muTHUVHZZ(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHZZ84"] = [=](const StandardModel& SM) { return new muTHUttHZZ(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUggHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUggHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUVBFHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUZHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUWHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUVHZZ4l(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHZZ4l84"] = [=](const StandardModel& SM) { return new muTHUttHZZ4l(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUggHWW84"] = [=](const StandardModel& SM) { return new muTHUggHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHWW84"] = [=](const StandardModel& SM) { return new muTHUVBFHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHWW84"] = [=](const StandardModel& SM) { return new muTHUZHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHWW84"] = [=](const StandardModel& SM) { return new muTHUWHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHWW84"] = [=](const StandardModel& SM) { return new muTHUVHWW(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHWW84"] = [=](const StandardModel& SM) { return new muTHUttHWW(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUggHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUggHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUVBFHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUZHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUWHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUVHWW2l2v(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHWW2l2v84"] = [=](const StandardModel& SM) { return new muTHUttHWW2l2v(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUggHmumu84"] = [=](const StandardModel& SM) { return new muTHUggHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHmumu84"] = [=](const StandardModel& SM) { return new muTHUVBFHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHmumu84"] = [=](const StandardModel& SM) { return new muTHUZHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHmumu84"] = [=](const StandardModel& SM) { return new muTHUWHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHmumu84"] = [=](const StandardModel& SM) { return new muTHUVHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHmumu84"] = [=](const StandardModel& SM) { return new muTHUttHmumu(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUggHtautau84"] = [=](const StandardModel& SM) { return new muTHUggHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHtautau84"] = [=](const StandardModel& SM) { return new muTHUVBFHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHtautau84"] = [=](const StandardModel& SM) { return new muTHUZHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHtautau84"] = [=](const StandardModel& SM) { return new muTHUWHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHtautau84"] = [=](const StandardModel& SM) { return new muTHUVHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHtautau84"] = [=](const StandardModel& SM) { return new muTHUttHtautau(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUggHbb84"] = [=](const StandardModel& SM) { return new muTHUggHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVBFHbb84"] = [=](const StandardModel& SM) { return new muTHUVBFHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUZHbb84"] = [=](const StandardModel& SM) { return new muTHUZHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUWHbb84"] = [=](const StandardModel& SM) { return new muTHUWHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHbb84"] = [=](const StandardModel& SM) { return new muTHUVHbb(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUttHbb84"] = [=](const StandardModel& SM) { return new muTHUttHbb(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUVBFBRinv84"] = [=](const StandardModel& SM) { return new muTHUVBFBRinv(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHBRinv84"] = [=](const StandardModel& SM) { return new muTHUVHBRinv(SM, sqrt_s_FCC84); };
    //
    obsThFactory["muTHUVBFHinv84"] = [=](const StandardModel& SM) { return new muTHUVBFHinv(SM, sqrt_s_FCC84); };
    obsThFactory["muTHUVHinv84"] = [=](const StandardModel& SM) { return new muTHUVHinv(SM, sqrt_s_FCC84); };
    //    
    //
    //----- STXS bins at hadron colliders
    // Stage 0
    obsThFactory["STXS_0_qqH_13"] = [=](const StandardModel& SM) { return new STXS_0_qqH(SM, sqrt_s_LHC13); };
    // Stage 1: 4l final state
    obsThFactory["STXSggH_VBFtopo_j3v_4l_13"] = [=](const StandardModel& SM) { return new STXSggH_VBFtopo_j3v_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH_VBFtopo_j3_4l_13"] = [=](const StandardModel& SM) { return new STXSggH_VBFtopo_j3_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH0j4l_13"] = [=](const StandardModel& SM) { return new STXSggH0j4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH1j_pTH_0_60_4l_13"] = [=](const StandardModel& SM) { return new STXSggH1j_pTH_0_60_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH1j_pTH_60_120_4l_13"] = [=](const StandardModel& SM) { return new STXSggH1j_pTH_60_120_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH1j_pTH_120_200_4l_13"] = [=](const StandardModel& SM) { return new STXSggH1j_pTH_120_200_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH1j_pTH_200_4l_13"] = [=](const StandardModel& SM) { return new STXSggH1j_pTH_200_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH2j_pTH_0_200_4l_13"] = [=](const StandardModel& SM) { return new STXSggH2j_pTH_0_200_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH2j_pTH_0_60_4l_13"] = [=](const StandardModel& SM) { return new STXSggH2j_pTH_0_60_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH2j_pTH_60_120_4l_13"] = [=](const StandardModel& SM) { return new STXSggH2j_pTH_60_120_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH2j_pTH_120_200_4l_13"] = [=](const StandardModel& SM) { return new STXSggH2j_pTH_120_200_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSggH2j_pTH_200_4l_13"] = [=](const StandardModel& SM) { return new STXSggH2j_pTH_200_4l(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSqqHqq_VBFtopo_Rest_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_VBFtopo_Rest_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_VBFtopo_j3v_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_VBFtopo_j3v_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_VBFtopo_j3_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_VBFtopo_j3_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_nonVHtopo_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_nonVHtopo_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_VHtopo_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_VHtopo_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_Rest_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_Rest_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHqq_pTj_200_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHqq_pTj_200_4l(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSqqHlv_pTV_0_250_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_0_250_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_0_150_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_0_150_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_150_250_0j_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_150_250_0j_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_150_250_1j_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_150_250_1j_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_250_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_250_4l(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSqqHll_pTV_0_150_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_0_150_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_0j_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_0j_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_1j_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_1j_4l(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_250_4l_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_250_4l(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSttHtH4l_13"] = [=](const StandardModel& SM) { return new STXSttHtH4l(SM, sqrt_s_LHC13); };
    // bb
    obsThFactory["STXSqqHlv_pTV_0_250_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_0_250_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_0_150_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_0_150_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_150_250_0j_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_150_250_0j_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_150_250_1j_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_150_250_1j_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHlv_pTV_250_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHlv_pTV_250_bb(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSqqHll_pTV_0_150_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_0_150_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_0j_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_0j_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_150_250_1j_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_150_250_1j_bb(SM, sqrt_s_LHC13); };
    obsThFactory["STXSqqHll_pTV_250_bb_13"] = [=](const StandardModel& SM) { return new STXSqqHll_pTV_250_bb(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSWHqqHqq_VBFtopo_j3v_2b"] = [=](const StandardModel& SM) { return new STXSWHqqHqq_VBFtopo_j3v_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSWHqqHqq_VBFtopo_j3_2b"] = [=](const StandardModel& SM) { return new STXSWHqqHqq_VBFtopo_j3_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSWHqqHqq_VH2j_2b"] = [=](const StandardModel& SM) { return new STXSWHqqHqq_VH2j_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSWHqqHqq_Rest_2b"] = [=](const StandardModel& SM) { return new STXSWHqqHqq_Rest_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSWHqqHqq_pTj1_200_2b"] = [=](const StandardModel& SM) { return new STXSWHqqHqq_pTj1_200_2b(SM, sqrt_s_LHC13); };
    //
    obsThFactory["STXSZHqqHqq_VBFtopo_j3v_2b"] = [=](const StandardModel& SM) { return new STXSZHqqHqq_VBFtopo_j3v_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSZHqqHqq_VBFtopo_j3_2b"] = [=](const StandardModel& SM) { return new STXSZHqqHqq_VBFtopo_j3_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSZHqqHqq_VH2j_2b"] = [=](const StandardModel& SM) { return new STXSZHqqHqq_VH2j_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSZHqqHqq_Rest_2b"] = [=](const StandardModel& SM) { return new STXSZHqqHqq_Rest_2b(SM, sqrt_s_LHC13); };
    obsThFactory["STXSZHqqHqq_pTj1_200_2b"] = [=](const StandardModel& SM) { return new STXSZHqqHqq_pTj1_200_2b(SM, sqrt_s_LHC13); };
    //
    // Stage 1.2: 4l final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300_Nj01(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450_Nj01(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_650_Nj01(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH650_Inf_Nj01(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_Inf_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH60_120_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggHll_pTV0_75_4l"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV0_75(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggHll_pTV75_150_4l"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV75_150(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ggHll_pTV250_Inf_4l"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV250_Inf(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj0_60_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj120_350_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHlv_pTV0_75_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_75(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHlv_pTV75_150_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV75_150(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHlv_pTV250_Inf_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_Inf(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHll_pTV0_75_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_75(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHll_pTV75_150_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV75_150(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_qqHll_pTV250_Inf_4l"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_Inf(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ttH_pTH0_60_4l"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_60(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ttH_pTH60_120_4l"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH60_120(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ttH_pTH120_200_4l"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ttH_pTH200_300_4l"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_ttH_pTH300_Inf_4l"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_Inf(SM, sqrt_s_LHC13, 1); };
    obsThFactory["STXS12_tH_4l"] = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 1); };
    //
    // Stage 1.2: ga ga final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300_Nj01(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450_Nj01(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_650_Nj01(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH650_Inf_Nj01(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_Inf_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH60_120_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggHll_pTV0_75_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV0_75(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggHll_pTV75_150_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV75_150(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggHll_pTV250_Inf_gaga"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV250_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj0_60_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj120_350_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV0_75_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_75(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV75_150_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV75_150(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV250_Inf_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV0_75_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_75(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV75_150_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV75_150(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV250_Inf_gaga"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH0_60_gaga"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_60(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH60_120_gaga"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH60_120(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH120_200_gaga"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH200_300_gaga"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH300_Inf_gaga"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_tH_gaga"] = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 2); };
    //
    // Stage 1.2: bb final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300_Nj01(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450_Nj01(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_650_Nj01(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH650_Inf_Nj01(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_Inf_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH60_120_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggHll_pTV0_75_bb"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV0_75(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggHll_pTV75_150_bb"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV75_150(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggHll_pTV250_Inf_bb"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV250_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj0_60_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj120_350_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV0_75_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_75(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV75_150_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV75_150(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV250_Inf_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV0_75_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_75(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV75_150_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV75_150(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV250_Inf_bb"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH0_60_bb"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_60(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH60_120_bb"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH60_120(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH120_200_bb"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH200_300_bb"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH300_Inf_bb"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_tH_bb"] = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 3); };
    //
    // Stage 1.2: evmuv final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300_Nj01(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450_Nj01(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_650_Nj01(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH650_Inf_Nj01(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_Inf_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH60_120_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggHll_pTV0_75_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV0_75(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggHll_pTV75_150_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV75_150(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ggHll_pTV250_Inf_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ggHll_pTV250_Inf(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj0_60_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj120_350_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHlv_pTV0_75_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_75(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHlv_pTV75_150_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV75_150(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHlv_pTV250_Inf_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_Inf(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHll_pTV0_75_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_75(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHll_pTV75_150_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV75_150(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj1(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_qqHll_pTV250_Inf_evmuv"] = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_Inf(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ttH_pTH0_60_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_60(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ttH_pTH60_120_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH60_120(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ttH_pTH120_200_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ttH_pTH200_300_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_ttH_pTH300_Inf_evmuv"] = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_Inf(SM, sqrt_s_LHC13, 4); };
    obsThFactory["STXS12_tH_evmuv"] = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 4); };
    //
    // AG:begin
    // Stage 1.2: HEPData, 2104706, Figure 7, crossSection[pb]
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_13"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_13"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_13"]                = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_13"]               = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_120_Nj2_13"]        = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_120_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_13"]      = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ggH_mjj350_Inf_pTH0_200_Nj2_13"]      = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_Nj2_13"]    = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2_13"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2_13"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2_13"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2_13"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2_13"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHlv_pTV0_75_13"]                    = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_75(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHlv_pTV75_150_13"]                  = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV75_150(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_13"]             = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHlv_pTV250_400_13"]                 = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_400(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHlv_pTV400_Inf_13"]                 = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV400_Inf(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHll_pTV0_150_13"]                   = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_150(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_13"]             = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_250_Nj0(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHll_pTV250_400_13"]                 = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_400(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_qqHll_pTV400_Inf_13"]                 = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV400_Inf(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH0_60_13"]                      = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_60(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH60_120_13"]                    = [=](const StandardModel& SM) { return new STXS12_ttH_pTH60_120(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH120_200_13"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH200_300_13"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH300_450_13"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_450(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_ttH_pTH450_Inf_13"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH450_Inf(SM, sqrt_s_LHC13, 0); };
    obsThFactory["STXS12_tH_13"]                               = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 0); };

    // Stage 1.2: https://arxiv.org/pdf/2402.05742
//    obsThFactory["STXS12_ggH_pTH0_200_Nj0_WW_placeholderObs"] = 1.;
    obsThFactory["STXS12_ggH_pTH0_200_Nj0_WW_placeholderObs"]  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_200_Nj0(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_WW"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_WW"]                = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_WW"]               = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH60_200_Nj1_WW"]                = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_200_Nj1(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH0_200_Nj2_WW"]                 = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH200_Inf_WW_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_Inf(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH200_300_WW_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_ggH_pTH300_Inf_WW_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_Inf(SM, sqrt_s_LHC13, 6); };

    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_Nj2_WW"]    = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2_WW"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2_WW"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2_WW"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2_WW"]    = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_WW"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_WW"]              = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 6); };
    
    obsThFactory["STXS12_qqHlv_pTV0_150_WW"]                  = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_150(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHlv_pTV150_Inf_WW"]                = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_Inf(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHll_pTV0_150_WW"]                  = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_150(SM, sqrt_s_LHC13, 6); };
    obsThFactory["STXS12_qqHll_pTV150_Inf_WW"]                = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_Inf(SM, sqrt_s_LHC13, 6); };
    
    obsThFactory["STXS12_ttH_tH_WW"]                          = [=](const StandardModel& SM) { return new STXS12_ttH_tH(SM, sqrt_s_LHC13, 6); };
    
    
    
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_ZZ"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH10_200_Nj0_ZZ_placeholderObs"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_200_Nj0(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_ZZ"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_ZZ"]                = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_ZZ"]               = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH0_200_Nj2_ZZ"]                 = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_200_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_pTH200_Inf_ZZ_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_Inf(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_ZZ"]              = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_VHlep_ZZ"]                            = [=](const StandardModel& SM) { return new STXS12_VHlep(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_VHlep_ZZ"]                            = [=](const StandardModel& SM) { return new STXS12_VHlep(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_VHlep_pTV0_150_ZZ"]                   = [=](const StandardModel& SM) { return new STXS12_VHlep_pTV0_150(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_VHlep_pTV150_Inf_ZZ"]                 = [=](const StandardModel& SM) { return new STXS12_VHlep_pTV150_Inf(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ttH_ZZ"]                              = [=](const StandardModel& SM) { return new STXS12_ttH(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ttH_tH_ZZ"]                           = [=](const StandardModel& SM) { return new STXS12_ttH_tH(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_ZZ"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 7); };
    
    obsThFactory["STXS12_ggH_mjj350_Inf_pTH0_200_Nj2_ZZ"]      = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_qqHqq_VH_veto_Nj01_ZZ_placeholderObs"]= [=](const StandardModel& SM) { return new STXS12_qqHqq_VH_veto_Nj01(SM, sqrt_s_LHC13, 7); };
    
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_ZZ_placeholder"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 7); };
    //obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_ZZ_placeholder"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 7);
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_ZZ_placeholder"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2_ZZ_placeholder"] = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2(SM, sqrt_s_LHC13, 7); };
    
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_ZZ"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_ZZ"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH60_120_Nj2(SM, sqrt_s_LHC13, 7); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_ZZ"] = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 7); };
    

    obsThFactory["STXS12_ggH_pTH0_10_Nj0_gaga"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH10_200_Nj0_gaga_placeholderObs"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_200_Nj0(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_120_Nj2_gaga"]       = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_120_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_mjj350_Inf_pTH0_200_Nj2_gaga"]     = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH200_300_gaga_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH300_450_gaga_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ggH_pTH450_Inf_gaga_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_Inf(SM, sqrt_s_LHC13, 2); };

    obsThFactory["STXS12_qqHqq_VH_veto_Nj01_gaga_placeholderObs"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_VH_veto_Nj01(SM, sqrt_s_LHC13, 2); };
    
    
    
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_Nj2_gaga"]   = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2_gaga"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2_gaga"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2_gaga"]= [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2_gaga"]= [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV0_150_gaga"]                  = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_150(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHlv_pTV150_Inf_gaga"]                = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV0_150_gaga"]                  = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_150(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_pTV150_Inf_gaga"]                = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_Inf(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_qqHll_gaga"]                           = [=](const StandardModel& SM) { return new STXS12_qqHll(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_ttH_pTH300_Inf_add_gaga"]              = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_Inf_add(SM, sqrt_s_LHC13, 2); };
    obsThFactory["STXS12_tH_gaga"]                              = [=](const StandardModel& SM) { return new STXS12_tH(SM, sqrt_s_LHC13, 2); };


    
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_tautau"]                  = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_10_Nj0(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH10_200_Nj0_tautau_placeholderObs"] = [=](const StandardModel& SM) { return new STXS12_ggH_pTH10_200_Nj0(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj1_tautau_placeholderObs"]   = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH0_60_Nj1(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_tautau"]            = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_60_Nj1(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_tautau"]          = [=](const StandardModel& SM) { return new STXS12_ggH_pTH60_120_Nj1(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH0_200_Nj2_tautau"]                 = [=](const StandardModel& SM) { return new STXS12_ggH_pTH0_200_Nj2(SM, sqrt_s_LHC13, 5); };
    
    obsThFactory["STXS12_qqHqq_VH_veto_Nj01_tautau_placeholderObs"]= [=](const StandardModel& SM) { return new STXS12_qqHqq_VH_veto_Nj01(SM, sqrt_s_LHC13, 5); };
    
    
    
    obsThFactory["STXS12_ggH_pTH200_300_tautau_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH300_Inf_tautau_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_Inf(SM, sqrt_s_LHC13, 5); };

    
    
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_Nj2_tautau"]    = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2_tautau"]    = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 5); };
    
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_tautau"]  = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(SM, sqrt_s_LHC13, 5); };
    
    
    obsThFactory["STXS12_qqHlv_pTV0_150_tautau"]                  = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV0_150(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHlv_pTV150_Inf_tautau"]                = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_Inf(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHll_pTV0_150_tautau"]                  = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV0_150(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHll_pTV150_Inf_tautau"]                = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV150_Inf(SM, sqrt_s_LHC13, 5); };
    
    obsThFactory["STXS12_ttH_tH_tautau"]                          = [=](const StandardModel& SM) { return new STXS12_ttH_tH(SM, sqrt_s_LHC13, 5); };
    
    
    
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_tautau"]            = [=](const StandardModel& SM) { return new STXS12_ggH_pTH120_200_Nj1(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_tautau"]   = [=](const StandardModel& SM) { return new STXS12_ggH_mjj0_350_pTH120_200_Nj2(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_mjj350_Inf_pTH0_200_Nj2_tautau"]   = [=](const StandardModel& SM) { return new STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_tautau"]           = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj60_120_Nj2(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_qqHqq_mjj350_Inf_Nj2_tautau"]          = [=](const StandardModel& SM) { return new STXS12_qqHqq_mjj350_Inf_Nj2(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ttH_tautau"]                           = [=](const StandardModel& SM) { return new STXS12_ttH(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH200_300_tautau_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH200_300(SM, sqrt_s_LHC13, 5); };
    obsThFactory["STXS12_ggH_pTH300_Inf_tautau_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_Inf(SM, sqrt_s_LHC13, 5); };



    //obsThFactory["STXS12_qqHqq_bb"]                            = bind(boost::factory<STXS12_qqHqq*>(), _1, sqrt_s_LHC13, 3);
    obsThFactory["STXS12_ggH_pTH300_450_bb_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH300_450(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH450_650_bb_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH450_650(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ggH_pTH650_Inf_bb_placeholderObs"]    = [=](const StandardModel& SM) { return new STXS12_ggH_pTH650_Inf(SM, sqrt_s_LHC13, 3); };

    obsThFactory["STXS12_qqHlv_pTV150_250_bb"]             = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV150_250_Nj0(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV250_400_bb"]                 = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_400(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV250_Inf_bb"]                 = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV250_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHlv_pTV400_Inf_bb"]                 = [=](const StandardModel& SM) { return new STXS12_qqHlv_pTV400_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH0_120_bb"]                     = [=](const StandardModel& SM) { return new STXS12_ttH_pTH0_120(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH120_200_bb"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH120_200(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH200_300_bb"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH200_300(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH300_450_bb"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH300_450(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_ttH_pTH450_Inf_bb"]                   = [=](const StandardModel& SM) { return new STXS12_ttH_pTH450_Inf(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV250_400_bb"]                 = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV250_400(SM, sqrt_s_LHC13, 3); };
    obsThFactory["STXS12_qqHll_pTV400_Inf_bb"]                 = [=](const StandardModel& SM) { return new STXS12_qqHll_pTV400_Inf(SM, sqrt_s_LHC13, 3); };

    //AG:end
    //
}
