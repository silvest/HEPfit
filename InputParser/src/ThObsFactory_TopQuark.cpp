/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "NPSMEFT6dtopquark.h"
#include "TopQuarkObservables.h"
#include "OptimizedObservablesSMEFTd6.h"
using namespace ThObsConst;

void ThObsFactory::registerTopQuarkObservables()
{
    //----- NPSMEFT6dtopquark  -----


    //Now that the map is defined these lines should be useless. Remove them and check that is fine
    /*
    obsThFactory["C_phit"] = [](const StandardModel& SM) { return new C_phit(SM); };
    obsThFactory["C_phiQ3"] = [](const StandardModel& SM) { return new C_phiQ3(SM); };
    obsThFactory["C_phiQ1"] = [](const StandardModel& SM) { return new C_phiQ1(SM); };
    obsThFactory["C_phiQm"] = [](const StandardModel& SM) { return new C_phiQm(SM); };
    obsThFactory["C_tW"] = [](const StandardModel& SM) { return new C_tW(SM); };
    obsThFactory["C_tZ"] = [](const StandardModel& SM) { return new C_tZ(SM); };
    obsThFactory["C_tB"] = [](const StandardModel& SM) { return new C_tB(SM); };
    obsThFactory["C_tphi"] = [](const StandardModel& SM) { return new C_tphi(SM); };
    obsThFactory["C_phib"] = [](const StandardModel& SM) { return new C_phib(SM); };
    obsThFactory["C_bW"] = [](const StandardModel& SM) { return new C_bW(SM); };
    obsThFactory["C_bB"] = [](const StandardModel& SM) { return new C_bB(SM); };
    obsThFactory["C_bZ"] = [](const StandardModel& SM) { return new C_bZ(SM); };


    obsThFactory["C_ed"] = [](const StandardModel& SM) { return new C_ed(SM); };
    obsThFactory["C_eq"] = [](const StandardModel& SM) { return new C_eq(SM); };
    obsThFactory["C_ld"] = [](const StandardModel& SM) { return new C_ld(SM); };
    obsThFactory["C_lqP"] = [](const StandardModel& SM) { return new C_lqP(SM); };
    obsThFactory["C_eu"] = [](const StandardModel& SM) { return new C_eu(SM); };
    obsThFactory["C_lu"] = [](const StandardModel& SM) { return new C_lu(SM); };
    obsThFactory["C_lqM"] = [](const StandardModel& SM) { return new C_lqM(SM); };


    obsThFactory["C_phitb"] = [](const StandardModel& SM) { return new C_phitb(SM); };
    obsThFactory["C_tG"] = [](const StandardModel& SM) { return new C_tG(SM); };
    obsThFactory["C_tu8"] = [](const StandardModel& SM) { return new C_tu8(SM); };
    obsThFactory["C_td8"] = [](const StandardModel& SM) { return new C_td8(SM); };
    obsThFactory["C_Qq18"] = [](const StandardModel& SM) { return new C_Qq18(SM); };
    obsThFactory["C_tq8"] = [](const StandardModel& SM) { return new C_tq8(SM); };
    obsThFactory["C_Qq38"] = [](const StandardModel& SM) { return new C_Qq38(SM); };
    obsThFactory["C_Qu8"] = [](const StandardModel& SM) { return new C_Qu8(SM); };
    obsThFactory["C_Qd8"] = [](const StandardModel& SM) { return new C_Qd8(SM); };
    */
   ///////////////////////////////////////////////////////////////////////////////////

    obsThFactory["FB_asymmetry_Tevatron_tt_diff_mtt_NPSMEFT6dtopquark"] = [](const StandardModel& SM) { return new FB_asymmetry_Tevatron_tt_diff_mtt_NPSMEFT6dtopquark(SM); };

    obsThFactory["Rb_NPSMEFT6dtopquark"] = [](const StandardModel& SM) { return new Rb_NPSMEFT6dtopquark(SM); };
    obsThFactory["AFBLR"] = [](const StandardModel& SM) { return new AFBLR(SM); };
    obsThFactory["SigmattZ"] = [](const StandardModel& SM) { return new sigmattZ(SM); };
    obsThFactory["SigmattA"] = [](const StandardModel& SM) { return new sigmattA(SM); };
    obsThFactory["SigmattH"] = [](const StandardModel& SM) { return new sigmattH(SM); };
    obsThFactory["SigmattW"] = [](const StandardModel& SM) { return new sigmattW(SM); };
    obsThFactory["SigmattWqEM"] = [](const StandardModel& SM) { return new ttWqEM(SM); };
    obsThFactory["SigmattWqSUM"] = [](const StandardModel& SM) { return new ttWqSUM(SM); };
    obsThFactory["Sigmatchannel13"] = [](const StandardModel& SM) { return new sigmatchannel13(SM); };
    obsThFactory["Sigmatchannel8"] = [](const StandardModel& SM) { return new sigmatchannel8(SM); };
    obsThFactory["SigmaschannelTev"] = [](const StandardModel& SM) { return new sigmaschannelTev(SM); };
    obsThFactory["Sigmaschannel8"] = [](const StandardModel& SM) { return new sigmaschannel8(SM); };
    obsThFactory["SigmatW"] = [](const StandardModel& SM) { return new sigmatW(SM); };
     obsThFactory["SigmatW_8TeV"] = [](const StandardModel& SM) { return new sigmatW_8TeV(SM); };
    obsThFactory["SigmatqZ"] = [](const StandardModel& SM) { return new sigmatqZ(SM); };
    obsThFactory["SigmatAq"] = [](const StandardModel& SM) { return new sigmatAq(SM); };
    obsThFactory["tH_Theo_Exp"] = [](const StandardModel& SM) { return new tH_tchan(SM); };
    obsThFactory["ttHSUM"] = [](const StandardModel& SM) { return new ttHSUM(SM); };
    obsThFactory["F0"] = [](const StandardModel& SM) { return new F0(SM); };
    obsThFactory["FL"] = [](const StandardModel& SM) { return new FL(SM); };

    obsThFactory["SigmattbarLHC13"] = [](const StandardModel& SM) { return new sigmattbarLHC13(SM); };
    obsThFactory["SigmattbarLHC8"] = [](const StandardModel& SM) { return new sigmattbarLHC8(SM); };
    obsThFactory["SigmattbarTev"] = [](const StandardModel& SM) { return new sigmattbarTev(SM); };



    obsThFactory["sigma_Z_pole_bb"] = [](const StandardModel& SM) { return new sigma_Z_pole_bb(SM); };
    obsThFactory["a_Z_pole_bb"] = [](const StandardModel& SM) { return new a_Z_pole_bb(SM); };
    obsThFactory["sigma_240_bb"] = [](const StandardModel& SM) { return new sigma_240_bb(SM); };
    obsThFactory["a_240_bb"] = [](const StandardModel& SM) { return new a_240_bb(SM); };
    obsThFactory["sigma_360_bb"] = [](const StandardModel& SM) { return new sigma_360_bb(SM); };
    obsThFactory["a_360_bb"] = [](const StandardModel& SM) { return new a_360_bb(SM); };



    obsThFactory["sigma_250_bb_eLpR"] = [](const StandardModel& SM) { return new sigma_250_bb_eLpR(SM); };
    obsThFactory["a_250_bb_eLpR"] = [](const StandardModel& SM) { return new a_250_bb_eLpR(SM); };
    obsThFactory["sigma_250_bb_eRpL"] = [](const StandardModel& SM) { return new sigma_250_bb_eRpL(SM); };
    obsThFactory["a_250_bb_eRpL"] = [](const StandardModel& SM) { return new a_250_bb_eRpL(SM); };
    obsThFactory["sigma_500_bb_eLpR"] = [](const StandardModel& SM) { return new sigma_500_bb_eLpR(SM); };
    obsThFactory["a_500_bb_eLpR"] = [](const StandardModel& SM) { return new a_500_bb_eLpR(SM); };
    obsThFactory["sigma_500_bb_eRpL"] = [](const StandardModel& SM) { return new sigma_500_bb_eRpL(SM); };
    obsThFactory["a_500_bb_eRpL"] = [](const StandardModel& SM) { return new a_500_bb_eRpL(SM); };
    obsThFactory["sigma_500_tt_eLpR"] = [](const StandardModel& SM) { return new sigma_500_tt_eLpR(SM); };
    obsThFactory["a_500_tt_eLpR"] = [](const StandardModel& SM) { return new a_500_tt_eLpR(SM); };
    obsThFactory["sigma_500_tt_eRpL"] = [](const StandardModel& SM) { return new sigma_500_tt_eRpL(SM); };
    obsThFactory["a_500_tt_eRpL"] = [](const StandardModel& SM) { return new a_500_tt_eRpL(SM); };
    obsThFactory["pt_500_tt_eLpR"] = [](const StandardModel& SM) { return new pt_500_tt_eLpR(SM); };
    obsThFactory["pt_500_tt_eRpL"] = [](const StandardModel& SM) { return new pt_500_tt_eRpL(SM); };

    obsThFactory["sigma_1000_bb_eLpR"] = [](const StandardModel& SM) { return new sigma_1000_bb_eLpR(SM); };
    obsThFactory["a_1000_bb_eLpR"] = [](const StandardModel& SM) { return new a_1000_bb_eLpR(SM); };
    obsThFactory["sigma_1000_bb_eRpL"] = [](const StandardModel& SM) { return new sigma_1000_bb_eRpL(SM); };
    obsThFactory["a_1000_bb_eRpL"] = [](const StandardModel& SM) { return new a_1000_bb_eRpL(SM); };


    //ttZ bins

    obsThFactory["sigma_ttz_diff_NLO_ATLAS_210312603"] = [](const StandardModel& SM) { return new sigma_ttz_diff_NLO_ATLAS_210312603(SM); };
    obsThFactory["sigma_ttz_diff_NLO_ATLAS_231204450"] = [](const StandardModel& SM) { return new sigma_ttz_diff_NLO_ATLAS_231204450(SM); };


    //ttA bins

    obsThFactory["sigma_tta_diff_NLO_ATLAS_emu"] = [](const StandardModel& SM) { return new sigma_tta_diff_NLO_ATLAS_emu_200706946(SM); };
    obsThFactory["sigma_tta_diff_NLO_CMS_dilepton_220107301"] = [](const StandardModel& SM) { return new sigma_tta_diff_NLO_CMS_dilepton_220107301(SM); };
    obsThFactory["sigma_tta_diff_NLO_ATLAS_240309452"] = [](const StandardModel& SM) { return new sigma_tta_diff_NLO_ATLAS_240309452(SM); };


    obsThFactory["sigma_tta_diff_NLO_ATLAS_emu"] = [](const StandardModel& SM) { return new sigma_tta_diff_NLO_ATLAS_emu_200706946(SM); };
    obsThFactory["sigma_tta_diff_NLO_CMS_dilepton_220107301"] = [](const StandardModel& SM) { return new sigma_tta_diff_NLO_CMS_dilepton_220107301(SM); };




    //ttH bins

    obsThFactory["sigma_ttH_diff_NLO_ATLAS_220700092"] = [](const StandardModel& SM) { return new sigma_ttH_diff_NLO_ATLAS_220700092(SM); };

    //ttbar bins

    obsThFactory["sigma_tt_diff_NLO"] = [](const StandardModel& SM) { return new sigma_tt_diff_NLO(SM); };

    //Charged asymmetry bins

    obsThFactory["charge_asymmetry_tt_diff_mtt_NLO"] = [](const StandardModel& SM) { return new charge_asymmetry_tt_diff_mtt_NLO(SM); };


    //ttll LHC observables

    obsThFactory["sigma_ttll_diff_LO"] = [](const StandardModel& SM) { return new sigma_ttll_diff_LO(SM); };

    //tt entanglement
    obsThFactory["entang_D_threshold"] = [](const StandardModel& SM) { return new entang_D_threshold(SM); };
    obsThFactory["entang_Dn_boosted"] = [](const StandardModel& SM) { return new entang_Dn_boosted(SM); };


    //ILC OBSERVABLES

    // ILC Z-Pole
    obsThFactory["sigma_Z_pole_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_Z_pole_bb_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_Z_pole_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_Z_pole_bb_eP_M30_eM_P80(SM); };
    obsThFactory["a_Z_pole_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new a_Z_pole_bb_eP_P30_eM_M80(SM); };
    obsThFactory["a_Z_pole_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new a_Z_pole_bb_eP_M30_eM_P80(SM); };

    //ILC 250 GeV

    obsThFactory["sigma_250_bb_eP_M30_eM_M80"] = [](const StandardModel& SM) { return new sigma_250_bb_eP_M30_eM_M80(SM); };
    obsThFactory["sigma_250_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_250_bb_eP_M30_eM_P80(SM); };
    obsThFactory["sigma_250_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_250_bb_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_250_bb_eP_P30_eM_P80"] = [](const StandardModel& SM) { return new sigma_250_bb_eP_P30_eM_P80(SM); };
    obsThFactory["a_250_bb_eP_M30_eM_M80"] = [](const StandardModel& SM) { return new a_250_bb_eP_M30_eM_M80(SM); };
    obsThFactory["a_250_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new a_250_bb_eP_M30_eM_P80(SM); };
    obsThFactory["a_250_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new a_250_bb_eP_P30_eM_M80(SM); };
    obsThFactory["a_250_bb_eP_P30_eM_P80"] = [](const StandardModel& SM) { return new a_250_bb_eP_P30_eM_P80(SM); };



    //ILC 500 GeV

    obsThFactory["sigma_500_bb_eP_M30_eM_M80"] = [](const StandardModel& SM) { return new sigma_500_bb_eP_M30_eM_M80(SM); };
    obsThFactory["sigma_500_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_500_bb_eP_M30_eM_P80(SM); };
    obsThFactory["sigma_500_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_500_bb_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_500_bb_eP_P30_eM_P80"] = [](const StandardModel& SM) { return new sigma_500_bb_eP_P30_eM_P80(SM); };
    obsThFactory["a_500_bb_eP_M30_eM_M80"] = [](const StandardModel& SM) { return new a_500_bb_eP_M30_eM_M80(SM); };
    obsThFactory["a_500_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new a_500_bb_eP_M30_eM_P80(SM); };
    obsThFactory["a_500_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new a_500_bb_eP_P30_eM_M80(SM); };
    obsThFactory["a_500_bb_eP_P30_eM_P80"] = [](const StandardModel& SM) { return new a_500_bb_eP_P30_eM_P80(SM); };
    obsThFactory["sigma_500_ttH_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_500_ttH_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_500_ttH_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_500_ttH_eP_M30_eM_P80(SM); };



    //ILC 1000 GeV

    obsThFactory["sigma_1000_bb_eP_M20_eM_M80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_M20_eM_M80(SM); };
    obsThFactory["sigma_1000_bb_eP_M20_eM_P80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_M20_eM_P80(SM); };
    obsThFactory["sigma_1000_bb_eP_P20_eM_M80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_P20_eM_M80(SM); };
    obsThFactory["sigma_1000_bb_eP_P20_eM_P80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_P20_eM_P80(SM); };
    obsThFactory["a_1000_bb_eP_M20_eM_M80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_M20_eM_M80(SM); };
    obsThFactory["a_1000_bb_eP_M20_eM_P80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_M20_eM_P80(SM); };
    obsThFactory["a_1000_bb_eP_P20_eM_M80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_P20_eM_M80(SM); };
    obsThFactory["a_1000_bb_eP_P20_eM_P80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_P20_eM_P80(SM); };


    obsThFactory["sigma_1000_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_1000_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_1000_bb_eP_M30_eM_P80(SM); };
    obsThFactory["a_1000_bb_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_P30_eM_M80(SM); };
    obsThFactory["a_1000_bb_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new a_1000_bb_eP_M30_eM_P80(SM); };
    obsThFactory["sigma_1000_ttH_eP_P30_eM_M80"] = [](const StandardModel& SM) { return new sigma_1000_ttH_eP_P30_eM_M80(SM); };
    obsThFactory["sigma_1000_ttH_eP_M30_eM_P80"] = [](const StandardModel& SM) { return new sigma_1000_ttH_eP_M30_eM_P80(SM); };


    //CLIC 380 OBSERVABLES

    //CLIC 380 GeV
    obsThFactory["sigma_380_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new sigma_380_bb_eP_0_eM_P80(SM); };
    obsThFactory["sigma_380_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new sigma_380_bb_eP_0_eM_M80(SM); };
    obsThFactory["a_380_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new a_380_bb_eP_0_eM_P80(SM); };
    obsThFactory["a_380_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new a_380_bb_eP_0_eM_M80(SM); };

    //CLIC 1500 GeV
    obsThFactory["sigma_1500_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new sigma_1500_bb_eP_0_eM_P80(SM); };
    obsThFactory["sigma_1500_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new sigma_1500_bb_eP_0_eM_M80(SM); };
    obsThFactory["a_1500_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new a_1500_bb_eP_0_eM_P80(SM); };
    obsThFactory["a_1500_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new a_1500_bb_eP_0_eM_M80(SM); };
    obsThFactory["sigma_1500_ttH_eP_0_eM_M80"] = [](const StandardModel& SM) { return new sigma_1500_ttH_eP_0_eM_M80(SM); };
    obsThFactory["sigma_1500_ttH_eP_0_eM_P80"] = [](const StandardModel& SM) { return new sigma_1500_ttH_eP_0_eM_P80(SM); };

    //CLIC 3000 GeV
    obsThFactory["sigma_3000_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new sigma_3000_bb_eP_0_eM_P80(SM); };
    obsThFactory["sigma_3000_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new sigma_3000_bb_eP_0_eM_M80(SM); };
    obsThFactory["a_3000_bb_eP_0_eM_P80"] = [](const StandardModel& SM) { return new a_3000_bb_eP_0_eM_P80(SM); };
    obsThFactory["a_3000_bb_eP_0_eM_M80"] = [](const StandardModel& SM) { return new a_3000_bb_eP_0_eM_M80(SM); };



    //Muon Collider
        //ttbar VBF
    obsThFactory["sigma_mumu_VBF_3TeV_tt"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_3TeV_tt(SM); };
    obsThFactory["sigma_mumu_VBF_10TeV_tt"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_10TeV_tt(SM); };
    obsThFactory["sigma_mumu_VBF_30TeV_tt"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_30TeV_tt(SM); };

        //ttH
    obsThFactory["sigma_mumu_3TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_3TeV_ttH(SM); };
    obsThFactory["sigma_mumu_10TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_10TeV_ttH(SM); };
    obsThFactory["sigma_mumu_30TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_30TeV_ttH(SM); };


    obsThFactory["sigma_mumu_sch_3TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_sch_3TeV_ttH(SM); };
    obsThFactory["sigma_mumu_sch_10TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_sch_10TeV_ttH(SM); };
    obsThFactory["sigma_mumu_sch_30TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_sch_30TeV_ttH(SM); };


    obsThFactory["sigma_mumu_VBF_3TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_3TeV_ttH(SM); };
    obsThFactory["sigma_mumu_VBF_10TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_10TeV_ttH(SM); };
    obsThFactory["sigma_mumu_VBF_30TeV_ttH"] = [](const StandardModel& SM) { return new sigma_mumu_VBF_30TeV_ttH(SM); };


        //bb
    obsThFactory["sigma_mumu_3TeV_bb"] = [](const StandardModel& SM) { return new sigma_mumu_3TeV_bb(SM); };
    obsThFactory["a_3TeV_mumu_bb"] = [](const StandardModel& SM) { return new a_3TeV_mumu_bb(SM); };
    obsThFactory["sigma_mumu_10TeV_bb"] = [](const StandardModel& SM) { return new sigma_mumu_10TeV_bb(SM); };
    obsThFactory["a_10TeV_mumu_bb"] = [](const StandardModel& SM) { return new a_10TeV_mumu_bb(SM); };
    obsThFactory["sigma_mumu_30TeV_bb"] = [](const StandardModel& SM) { return new sigma_mumu_30TeV_bb(SM); };
    obsThFactory["a_30TeV_mumu_bb"] = [](const StandardModel& SM) { return new a_30TeV_mumu_bb(SM); };


    //Optimal Observables

    obsThFactory["opt_obs_ilc_500_M30_M80"] = [](const StandardModel& SM) { return new opt_obs_ilc_500_M30_M80(SM); };
    obsThFactory["opt_obs_ilc_500_M30_P80"] = [](const StandardModel& SM) { return new opt_obs_ilc_500_M30_P80(SM); };
    obsThFactory["opt_obs_ilc_500_P30_M80"] = [](const StandardModel& SM) { return new opt_obs_ilc_500_P30_M80(SM); };
    obsThFactory["opt_obs_ilc_500_P30_P80"] = [](const StandardModel& SM) { return new opt_obs_ilc_500_P30_P80(SM); };

    obsThFactory["opt_obs_ilc_1000_M30_P80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_M30_P80(SM); };
    obsThFactory["opt_obs_ilc_1000_P30_M80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_P30_M80(SM); };

    obsThFactory["opt_obs_ilc_1000_M20_M80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_M20_M80(SM); };
    obsThFactory["opt_obs_ilc_1000_M20_P80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_M20_P80(SM); };
    obsThFactory["opt_obs_ilc_1000_P20_M80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_P20_M80(SM); };
    obsThFactory["opt_obs_ilc_1000_P20_P80"] = [](const StandardModel& SM) { return new opt_obs_ilc_1000_P20_P80(SM); };


    obsThFactory["opt_obs_clic_380_0_M80"] = [](const StandardModel& SM) { return new opt_obs_clic_380_0_M80(SM); };
    obsThFactory["opt_obs_clic_380_0_P80"] = [](const StandardModel& SM) { return new opt_obs_clic_380_0_P80(SM); };
    obsThFactory["opt_obs_clic_1500_0_M80"] = [](const StandardModel& SM) { return new opt_obs_clic_1500_0_M80(SM); };
    obsThFactory["opt_obs_clic_1500_0_P80"] = [](const StandardModel& SM) { return new opt_obs_clic_1500_0_P80(SM); };
    obsThFactory["opt_obs_clic_3000_0_M80"] = [](const StandardModel& SM) { return new opt_obs_clic_3000_0_M80(SM); };
    obsThFactory["opt_obs_clic_3000_0_P80"] = [](const StandardModel& SM) { return new opt_obs_clic_3000_0_P80(SM); };


    obsThFactory["opt_obs_fcc_350"] = [](const StandardModel& SM) { return new opt_obs_fcc_350(SM); };
    obsThFactory["opt_obs_fcc_365"] = [](const StandardModel& SM) { return new opt_obs_fcc_365(SM); };


    obsThFactory["opt_obs_cepc_350"] = [](const StandardModel& SM) { return new opt_obs_cepc_350(SM); };
    obsThFactory["opt_obs_cepc_360"] = [](const StandardModel& SM) { return new opt_obs_cepc_360(SM); };


    obsThFactory["opt_obs_muon_3TeV"] = [](const StandardModel& SM) { return new opt_obs_muon_3TeV(SM); };
    obsThFactory["opt_obs_muon_10TeV"] = [](const StandardModel& SM) { return new opt_obs_muon_10TeV(SM); };
    obsThFactory["opt_obs_muon_30TeV"] = [](const StandardModel& SM) { return new opt_obs_muon_30TeV(SM); };


    //Observables for CP-violation at tth and thj (prospects and proposed asymmetries)



    obsThFactory["Asymmetry_Dazi_ord_ttH_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_Dazi_ord_ttH_eta_cut_3(SM); };
    obsThFactory["Asymmetry_Dazi_ord_ttH_eta_cut_3_ee"] = [](const StandardModel& SM) { return new Asymmetry_Dazi_ord_ttH_eta_cut_3_ee(SM); };
    obsThFactory["Asymmetry_trip_prod_pt_pe_pp_ttH_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_trip_prod_pt_pe_pp_ttH_eta_cut_3(SM); };
    obsThFactory["Asymmetry_sign_trip_prod_pe_pp_ttH_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_sign_trip_prod_pe_pp_ttH_eta_cut_3(SM); };

    obsThFactory["sigma_ttH_eta_cut_3_diff_LO_mtth"] = [](const StandardModel& SM) { return new sigma_ttH_eta_cut_3_diff_LO_mtth(SM); };
    obsThFactory["b4_ttH_eta_cut_3_LO"] = [](const StandardModel& SM) { return new b4_ttH_eta_cut_3_LO(SM); };




    obsThFactory["Asymmetry_cos_je_tHj_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_cos_je_tHj_eta_cut_3(SM); };
    obsThFactory["Asymmetry_cos_se_tHj_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_cos_se_tHj_eta_cut_3(SM); };
    obsThFactory["Asymmetry_cos_ye_tHj_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_cos_ye_tHj_eta_cut_3(SM); };
    obsThFactory["Asymmetry_trip_prod_ph_pt_pj_tHj_eta_cut_3"] = [](const StandardModel& SM) { return new Asymmetry_trip_prod_ph_pt_pj_tHj_eta_cut_3(SM); };

    obsThFactory["sigma_tHj_eta_cut_3_diff_LO_Del_R_th"] = [](const StandardModel& SM) { return new sigma_tHj_eta_cut_3_diff_LO_Del_R_th(SM); };
    obsThFactory["sigma_tHj_eta_cut_3_diff_LO_mth"] = [](const StandardModel& SM) { return new sigma_tHj_eta_cut_3_diff_LO_mth(SM); };
    obsThFactory["sigma_tHj_eta_cut_3_diff_LO_trip_prod_z_pt_pj"] = [](const StandardModel& SM) { return new sigma_tHj_eta_cut_3_diff_LO_trip_prod_z_pt_pj(SM); };






    obsThFactory["Asymmetry_Dazi_ord_ttH"] = [](const StandardModel& SM) { return new Asymmetry_Dazi_ord_ttH(SM); };
    obsThFactory["Asymmetry_Dazi_ord_ttH_ee"] = [](const StandardModel& SM) { return new Asymmetry_Dazi_ord_ttH_ee(SM); };
    obsThFactory["Asymmetry_trip_prod_pt_pe_pp_ttH"] = [](const StandardModel& SM) { return new Asymmetry_trip_prod_pt_pe_pp_ttH(SM); };
    obsThFactory["Asymmetry_sign_trip_prod_pe_pp_ttH"] = [](const StandardModel& SM) { return new Asymmetry_sign_trip_prod_pe_pp_ttH(SM); };

    obsThFactory["sigma_ttH_diff_LO_mtth"] = [](const StandardModel& SM) { return new sigma_ttH_diff_LO_mtth(SM); };
    obsThFactory["b4_ttH_LO"] = [](const StandardModel& SM) { return new b4_ttH_LO(SM); };




    obsThFactory["Asymmetry_cos_je_tHj"] = [](const StandardModel& SM) { return new Asymmetry_cos_je_tHj(SM); };
    obsThFactory["Asymmetry_cos_se_tHj"] = [](const StandardModel& SM) { return new Asymmetry_cos_se_tHj(SM); };
    obsThFactory["Asymmetry_cos_ye_tHj"] = [](const StandardModel& SM) { return new Asymmetry_cos_ye_tHj(SM); };
    obsThFactory["Asymmetry_trip_prod_ph_pt_pj_tHj"] = [](const StandardModel& SM) { return new Asymmetry_trip_prod_ph_pt_pj_tHj(SM); };

    obsThFactory["sigma_tHj_diff_LO_Del_R_th"] = [](const StandardModel& SM) { return new sigma_tHj_diff_LO_Del_R_th(SM); };
    obsThFactory["sigma_tHj_diff_LO_mth"] = [](const StandardModel& SM) { return new sigma_tHj_diff_LO_mth(SM); };
    obsThFactory["sigma_tHj_diff_LO_trip_prod_z_pt_pj"] = [](const StandardModel& SM) { return new sigma_tHj_diff_LO_trip_prod_z_pt_pj(SM); };




    //obsThFactory["opt_obs_ilc_500_1000"] = boost::factory<opt_obs_ilc_500_1000*>();
    //obsThFactory["test_cov"] = boost::factory<test_cov*>();


    //OPTIMIZED OBSERVABLES OLD!!!!!
    //I don't really like this implementation, these constraints should be included in the prior,
    //as we do now (from 2022). Remove this also in the code and check everything is fine
    obsThFactory["op1"] = [](const StandardModel& SM) { return new op1(SM); };
    obsThFactory["op2"] = [](const StandardModel& SM) { return new op2(SM); };
    obsThFactory["op3"] = [](const StandardModel& SM) { return new op3(SM); };
    obsThFactory["op4"] = [](const StandardModel& SM) { return new op4(SM); };

    //OPTIMIZED OBSERVABLES 1000 GeV
    //Same as above
    obsThFactory["op_1000_1"] = [](const StandardModel& SM) { return new op_1000_1(SM); };
    obsThFactory["op_1000_2"] = [](const StandardModel& SM) { return new op_1000_2(SM); };
    obsThFactory["op_1000_3"] = [](const StandardModel& SM) { return new op_1000_3(SM); };
    obsThFactory["op_1000_4"] = [](const StandardModel& SM) { return new op_1000_4(SM); };
    obsThFactory["op_1000_5"] = [](const StandardModel& SM) { return new op_1000_5(SM); };
    obsThFactory["op_1000_6"] = [](const StandardModel& SM) { return new op_1000_6(SM); };
    obsThFactory["op_1000_7"] = [](const StandardModel& SM) { return new op_1000_7(SM); };
    obsThFactory["op_1000_8"] = [](const StandardModel& SM) { return new op_1000_8(SM); };

    //Some auxiliary observables defining combinations of WC

    obsThFactory["op_eigen_ttll_1"] = [](const StandardModel& SM) { return new op_eigen_ttll_1(SM); };
    obsThFactory["op_eigen_ttll_2"] = [](const StandardModel& SM) { return new op_eigen_ttll_2(SM); };
    obsThFactory["op_eigen_ttll_3"] = [](const StandardModel& SM) { return new op_eigen_ttll_3(SM); };
    obsThFactory["op_eigen_ttll_4"] = [](const StandardModel& SM) { return new op_eigen_ttll_4(SM); };


   obsThFactory["gLt"] = [](const StandardModel& SM) { return new gLt(SM); };
   obsThFactory["gLb"] = [](const StandardModel& SM) { return new gLb(SM); };
   obsThFactory["gRt"] = [](const StandardModel& SM) { return new gRt(SM); };
   obsThFactory["gRb"] = [](const StandardModel& SM) { return new gRb(SM); };





   //----- TopQuarkObservables begin -----

   obsThFactory["FB_asymmetry_Tevatron_tt_diff_mtt_LO"] = [](const StandardModel& SM) { return new FB_asymmetry_Tevatron_tt_diff_mtt_LO(SM); };

   obsThFactory["charge_asymmetry_tt_diff_mtt_LO"] = [](const StandardModel& SM) { return new charge_asymmetry_tt_diff_mtt_LO(SM); };

   obsThFactory["sigma_tt_diff_mtt_LO_CMS_181106625"] = [](const StandardModel& SM) { return new sigma_tt_diff_mtt_LO_CMS_181106625(SM); };
   obsThFactory["sigma_tt_diff_mtt_CMS_LO"] = [](const StandardModel& SM) { return new sigma_tt_diff_mtt_CMS_LO(SM); };
   obsThFactory["sigma_norm_tt_diff_mtt_ATLAS_LO"] = [](const StandardModel& SM) { return new sigma_norm_tt_diff_mtt_ATLAS_LO(SM); };
   obsThFactory["sigma_tt_13_LO"] = [](const StandardModel& SM) { return new sigma_tt_13_LO(SM); };
   obsThFactory["R_tt_8_o_7_LO"] = [](const StandardModel& SM) { return new R_tt_8_o_7_LO(SM); };
   obsThFactory["R_tt_13_o_8_LO"] = [](const StandardModel& SM) { return new R_tt_13_o_8_LO(SM); };

   obsThFactory["sigma_ttz_diff_LO_CMS_190711270"] = [](const StandardModel& SM) { return new sigma_ttz_diff_LO_CMS_190711270(SM); };
   obsThFactory["sigma_ttz_diff_LO_ATLAS_210312603"] = [](const StandardModel& SM) { return new sigma_ttz_diff_LO_ATLAS_210312603(SM); };
   obsThFactory["sigma_ttz_diff_LO_ATLAS_231204450"] = [](const StandardModel& SM) { return new sigma_ttz_diff_LO_ATLAS_231204450(SM); };

   obsThFactory["sigma_tta_diff_LO_ATLAS_emu"] = [](const StandardModel& SM) { return new sigma_tta_diff_LO_ATLAS_emu_200706946(SM); };
   obsThFactory["sigma_tta_diff_LO_CMS_dilepton"] = [](const StandardModel& SM) { return new sigma_tta_diff_LO_CMS_dilepton_220107301(SM); };
   obsThFactory["sigma_tta_diff_LO_CMS_semileptonic"] = [](const StandardModel& SM) { return new sigma_tta_diff_LO_CMS_semileptonic_210701508(SM); };

   obsThFactory["F0_LO"] = [](const StandardModel& SM) { return new F0_LO(SM); };
   obsThFactory["FL_LO"] = [](const StandardModel& SM) { return new FL_LO(SM); };
   obsThFactory["FR_LO"] = [](const StandardModel& SM) { return new FL_LO(SM); };

   obsThFactory["sigma_tttt_13_LO"] = [](const StandardModel& SM) { return new sigma_tttt_13_LO(SM); };
   obsThFactory["sigma_ttbb_13_LO_dilepton"] = [](const StandardModel& SM) { return new sigma_ttbb_13_LO_dilepton(SM); };
   obsThFactory["sigma_ttbb_13_LO_lepjet"] = [](const StandardModel& SM) { return new sigma_ttbb_13_LO_lepjet(SM); };

   obsThFactory["sigma_taq_LO_CMS"] = [](const StandardModel& SM) { return new sigma_taq_LO_CMS(SM); };
   obsThFactory["sigma_taq_LO_ATLAS"] = [](const StandardModel& SM) { return new sigma_taq_LO_ATLAS(SM); };


   obsThFactory["sigma_tzq_LO"] = [](const StandardModel& SM) { return new sigma_tzq_LO(SM); };
   obsThFactory["sigma_ttw_LO"] = [](const StandardModel& SM) { return new sigma_ttw_LO(SM); };
   obsThFactory["R_ttw_LO"] = [](const StandardModel& SM) { return new R_ttw_LO(SM); };
   obsThFactory["sigma_tw_13_LO"] = [](const StandardModel& SM) { return new sigma_tw_13_LO(SM); };
   obsThFactory["sigma_tb_13_LO"] = [](const StandardModel& SM) { return new sigma_tb_13_LO(SM); };
   obsThFactory["sigma_tq_13_LO"] = [](const StandardModel& SM) { return new sigma_tq_13_LO(SM); };
   obsThFactory["sigma_tq_top_13_LO"] = [](const StandardModel& SM) { return new sigma_tq_top_13_LO(SM); };
   obsThFactory["sigma_tq_antitop_13_LO"] = [](const StandardModel& SM) { return new sigma_tq_antitop_13_LO(SM); };

   obsThFactory["sigma_tw_8_LO"] = [](const StandardModel& SM) { return new sigma_tw_8_LO(SM); };
   obsThFactory["sigma_tb_8_LO"] = [](const StandardModel& SM) { return new sigma_tb_8_LO(SM); };
   obsThFactory["sigma_tq_8_LO"] = [](const StandardModel& SM) { return new sigma_tq_8_LO(SM); };
   obsThFactory["sigma_tw_7_LO"] = [](const StandardModel& SM) { return new sigma_tw_7_LO(SM); };
   obsThFactory["sigma_tq_7_LO"] = [](const StandardModel& SM) { return new sigma_tq_7_LO(SM); };







   //----- TopQuarkObservables end -----
}
