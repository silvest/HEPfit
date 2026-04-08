/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "FlavourObservables.h"
#include "LeptonFlavourObservables.h"
#include "LoopMediators.h"
#include "SUSYObservables.h"
using namespace ThObsConst;

void ThObsFactory::registerFlavourObservables()
{
    obsThFactory["Retau_DS1"] = [](const StandardModel& SM) { return new Retau_DS1(SM); };
    obsThFactory["Imtau_DS1"] = [](const StandardModel& SM) { return new Imtau_DS1(SM); };
    //----- B(s) to mu mu  -----
    obsThFactory["BR_Bdmumu"] = [=](const StandardModel& SM) { return new Mll(SM, 1, StandardModel::B_D, StandardModel::MU); };
    obsThFactory["BRbar_Bdmumu"] = [=](const StandardModel& SM) { return new Mll(SM, 2, StandardModel::B_D, StandardModel::MU); };
    obsThFactory["Amumu_Bd"] = [=](const StandardModel& SM) { return new Mll(SM, 3, StandardModel::B_D, StandardModel::MU); };
    obsThFactory["Smumu_Bd"] = [=](const StandardModel& SM) { return new Mll(SM, 4, StandardModel::B_D, StandardModel::MU); };

    obsThFactory["BR_Bsmumu"] = [=](const StandardModel& SM) { return new Mll(SM, 1, StandardModel::B_S, StandardModel::MU); };
    obsThFactory["BRbar_Bsmumu"] = [=](const StandardModel& SM) { return new Mll(SM, 2, StandardModel::B_S, StandardModel::MU); };
    obsThFactory["Amumu_Bs"] = [=](const StandardModel& SM) { return new Mll(SM, 3, StandardModel::B_S, StandardModel::MU); };
    obsThFactory["Smumu_Bs"] = [=](const StandardModel& SM) { return new Mll(SM, 4, StandardModel::B_S, StandardModel::MU); };
    obsThFactory["BR_Bsee"] = [=](const StandardModel& SM) { return new Mll(SM, 1, StandardModel::B_S, StandardModel::ELECTRON); };
    obsThFactory["BRbar_Bsee"] = [=](const StandardModel& SM) { return new Mll(SM, 2, StandardModel::B_S, StandardModel::ELECTRON); };
    obsThFactory["Aee_Bs"] = [=](const StandardModel& SM) { return new Mll(SM, 3, StandardModel::B_S, StandardModel::ELECTRON); };
    obsThFactory["See_Bs"] = [=](const StandardModel& SM) { return new Mll(SM, 4, StandardModel::B_S, StandardModel::ELECTRON); };

    obsThFactory["BR_BdmumuOBR_Bsmumu"] = [](const StandardModel& SM) { return new BdmumuOBsmumu(SM); };
//----- B(s) to mu mu  OLD -----
    obsThFactory["BR_Bdmumu_old"] = [=](const StandardModel& SM) { return new Bdmumu(SM, 1); };
    obsThFactory["BRbar_Bdmumu_old"] = [=](const StandardModel& SM) { return new Bdmumu(SM, 2); };
    obsThFactory["Amumu_Bd_old"] = [=](const StandardModel& SM) { return new Bdmumu(SM, 3); };
    obsThFactory["Smumu_Bd_old"] = [=](const StandardModel& SM) { return new Bdmumu(SM, 4); };
    obsThFactory["BR_Bsmumu_old"] = [=](const StandardModel& SM) { return new Bsmumu(SM, 1); };
    obsThFactory["BRbar_Bsmumu_old"] = [=](const StandardModel& SM) { return new Bsmumu(SM, 2); };
    obsThFactory["Amumu_Bs_old"] = [=](const StandardModel& SM) { return new Bsmumu(SM, 3); };
    obsThFactory["Smumu_Bs_old"] = [=](const StandardModel& SM) { return new Bsmumu(SM, 4); };
   //----- b to q gamma  -----
    obsThFactory["BR_bsgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, StandardModel::STRANGE, 1); };
    obsThFactory["ACP_bsgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, StandardModel::STRANGE, 2); };
    obsThFactory["BR_bdgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, StandardModel::DOWN, 1); };
    obsThFactory["ACP_bdgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, StandardModel::DOWN, 2); };
    obsThFactory["BR_bqgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, 1); };
    obsThFactory["ACP_bqgamma"] = [=](const StandardModel& SM) { return new Bsgamma(SM, 2); };

    //----- Wilson coefficients  -----
    obsThFactory["WC_real_C7g"] = [=](const StandardModel& SM) { return new WC_C7g(SM, 0); };
    obsThFactory["WC_imag_C7g"] = [=](const StandardModel& SM) { return new WC_C7g(SM, 1); };
    obsThFactory["WC_abs_C7g"] = [=](const StandardModel& SM) { return new WC_C7g(SM, 2); };
    obsThFactory["WC_arg_C7g"] = [=](const StandardModel& SM) { return new WC_C7g(SM, 3); };

    obsThFactory["WC_real_C9_mu"] = [=](const StandardModel& SM) { return new WC_C9(SM, 0, StandardModel::MU); };
    obsThFactory["WC_imag_C9_mu"] = [=](const StandardModel& SM) { return new WC_C9(SM, 1, StandardModel::MU); };
    obsThFactory["WC_abs_C9_mu"] = [=](const StandardModel& SM) { return new WC_C9(SM, 2, StandardModel::MU); };
    obsThFactory["WC_arg_C9_mu"] = [=](const StandardModel& SM) { return new WC_C9(SM, 3, StandardModel::MU); };

    obsThFactory["WC_real_C9_el"] = [=](const StandardModel& SM) { return new WC_C9(SM, 0, StandardModel::ELECTRON); };
    obsThFactory["WC_imag_C9_el"] = [=](const StandardModel& SM) { return new WC_C9(SM, 1, StandardModel::ELECTRON); };
    obsThFactory["WC_abs_C9_el"] = [=](const StandardModel& SM) { return new WC_C9(SM, 2, StandardModel::ELECTRON); };
    obsThFactory["WC_arg_C9_el"] = [=](const StandardModel& SM) { return new WC_C9(SM, 3, StandardModel::ELECTRON); };

    obsThFactory["WC_real_C10_mu"] = [=](const StandardModel& SM) { return new WC_C10(SM, 0, StandardModel::MU); };
    obsThFactory["WC_imag_C10_mu"] = [=](const StandardModel& SM) { return new WC_C10(SM, 1, StandardModel::MU); };
    obsThFactory["WC_abs_C10_mu"] = [=](const StandardModel& SM) { return new WC_C10(SM, 2, StandardModel::MU); };
    obsThFactory["WC_arg_C10_mu"] = [=](const StandardModel& SM) { return new WC_C10(SM, 3, StandardModel::MU); };

    obsThFactory["WC_real_C10_el"] = [=](const StandardModel& SM) { return new WC_C10(SM, 0, StandardModel::ELECTRON); };
    obsThFactory["WC_imag_C10_el"] = [=](const StandardModel& SM) { return new WC_C10(SM, 1, StandardModel::ELECTRON); };
    obsThFactory["WC_abs_C10_el"] = [=](const StandardModel& SM) { return new WC_C10(SM, 2, StandardModel::ELECTRON); };
    obsThFactory["WC_arg_C10_el"] = [=](const StandardModel& SM) { return new WC_C10(SM, 3, StandardModel::ELECTRON); };

    obsThFactory["WC_arg_C10_el"] = [=](const StandardModel& SM) { return new WC_C10(SM, 3, StandardModel::ELECTRON); };

    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = [=](const StandardModel& SM) { return new P_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_1_BdKste"] = [=](const StandardModel& SM) { return new P_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_2_BdKstmu"] = [=](const StandardModel& SM) { return new P_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_2_BdKste"] = [=](const StandardModel& SM) { return new P_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_3_BdKstmu"] = [=](const StandardModel& SM) { return new P_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_3_BdKste"] = [=](const StandardModel& SM) { return new P_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_4p_BdKstmu"] = [=](const StandardModel& SM) { return new P_4Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_4p_BdKste"] = [=](const StandardModel& SM) { return new P_4Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_5p_BdKstmu"] = [=](const StandardModel& SM) { return new P_5Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_5p_BdKste"] = [=](const StandardModel& SM) { return new P_5Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_6p_BdKstmu"] = [=](const StandardModel& SM) { return new P_6Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_6p_BdKste"] = [=](const StandardModel& SM) { return new P_6Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["P_8p_BdKstmu"] = [=](const StandardModel& SM) { return new P_8Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_8p_BdKste"] = [=](const StandardModel& SM) { return new P_8Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["Gammap_BdKstmu"] = [=](const StandardModel& SM) { return new GammaPrime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Gammap_BdKste"] = [=](const StandardModel& SM) { return new GammaPrime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["A_FB_BdKstmu"] = [=](const StandardModel& SM) { return new A_FB(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_FB_BdKste"] = [=](const StandardModel& SM) { return new A_FB(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["BR_BdKstmu"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["BR_BdKste"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["RKst_BdKstll"] = [=](const StandardModel& SM) { return new R_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RKstL_BdKstll"] = [=](const StandardModel& SM) { return new RL_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RKstT_BdKstll"] = [=](const StandardModel& SM) { return new RT_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RKstp_BpKstpll"] = [=](const StandardModel& SM) { return new R_MVll(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["R6_BdKstll"] = [=](const StandardModel& SM) { return new R_6(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["ACP_BdKstmu"] = [=](const StandardModel& SM) { return new ACP_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P3CP_BdKstmu"] = [=](const StandardModel& SM) { return new P3CP(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["F_L_BdKstmu"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["F_L_BdKste"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["M_1p_BdKstmu"] = [=](const StandardModel& SM) { return new M_1Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["M_2p_BdKstmu"] = [=](const StandardModel& SM) { return new M_2Prime(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_1c_BdKstmu"] = [=](const StandardModel& SM) { return new S_1c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_1c_BdKste"] = [=](const StandardModel& SM) { return new S_1c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_1s_BdKstmu"] = [=](const StandardModel& SM) { return new S_1s(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_1s_BdKste"] = [=](const StandardModel& SM) { return new S_1s(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_2c_BdKstmu"] = [=](const StandardModel& SM) { return new S_2c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_2c_BdKste"] = [=](const StandardModel& SM) { return new S_2c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_2s_BdKstmu"] = [=](const StandardModel& SM) { return new S_2s(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_2s_BdKste"] = [=](const StandardModel& SM) { return new S_2s(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_3_BdKstmu"] = [=](const StandardModel& SM) { return new S_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_3_BdKste"] = [=](const StandardModel& SM) { return new S_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_4_BdKstmu"] = [=](const StandardModel& SM) { return new S_4(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_4_BdKste"] = [=](const StandardModel& SM) { return new S_4(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_5_BdKstmu"] = [=](const StandardModel& SM) { return new S_5(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_5_BdKste"] = [=](const StandardModel& SM) { return new S_5(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_6c_BdKstmu"] = [=](const StandardModel& SM) { return new S_6c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_6c_BdKste"] = [=](const StandardModel& SM) { return new S_6c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_7_BdKstmu"] = [=](const StandardModel& SM) { return new S_7(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_7_BdKste"] = [=](const StandardModel& SM) { return new S_7(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_8_BdKstmu"] = [=](const StandardModel& SM) { return new S_8(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_8_BdKste"] = [=](const StandardModel& SM) { return new S_8(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["S_9_BdKstmu"] = [=](const StandardModel& SM) { return new S_9(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_9_BdKste"] = [=](const StandardModel& SM) { return new S_9(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON); };
    obsThFactory["A_5_BdKstmu"] = [=](const StandardModel& SM) { return new A_5(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_6_BdKstmu"] = [=](const StandardModel& SM) { return new A_6(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_6c_BdKstmu"] = [=](const StandardModel& SM) { return new A_6c(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_8_BdKstmu"] = [=](const StandardModel& SM) { return new A_8(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_9_BdKstmu"] = [=](const StandardModel& SM) { return new A_9(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["P_1f_BdKstmu"] = [=](const StandardModel& SM) { return new P_1f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_2f_BdKstmu"] = [=](const StandardModel& SM) { return new P_2f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_3f_BdKstmu"] = [=](const StandardModel& SM) { return new P_3f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_4pf_BdKstmu"] = [=](const StandardModel& SM) { return new P_4Primef(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_5pf_BdKstmu"] = [=](const StandardModel& SM) { return new P_5Primef(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_6pf_BdKstmu"] = [=](const StandardModel& SM) { return new P_6Primef(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_8pf_BdKstmu"] = [=](const StandardModel& SM) { return new P_8Primef(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Gammapf_BdKstmu"] = [=](const StandardModel& SM) { return new GammaPrimef(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["BRf_BdKstmu"] = [=](const StandardModel& SM) { return new BRf_MVll(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["A_FBf_BdKstmu"] = [=](const StandardModel& SM) { return new A_FBf(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["F_Lf_BdKstmu"] = [=](const StandardModel& SM) { return new F_Lf(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_3f_BdKstmu"] = [=](const StandardModel& SM) { return new S_3f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_4f_BdKstmu"] = [=](const StandardModel& SM) { return new S_4f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_5f_BdKstmu"] = [=](const StandardModel& SM) { return new S_5f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_7f_BdKstmu"] = [=](const StandardModel& SM) { return new S_7f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_8f_BdKstmu"] = [=](const StandardModel& SM) { return new S_8f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_9f_BdKstmu"] = [=](const StandardModel& SM) { return new S_9f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_relationf"] = [=](const StandardModel& SM) { return new P_relationf(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["P_relation_exactf"] = [=](const StandardModel& SM) { return new P_relation_exactf(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["F_L_BpKstmu"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_3_BpKstmu"] = [=](const StandardModel& SM) { return new S_3(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_4_BpKstmu"] = [=](const StandardModel& SM) { return new S_4(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_5_BpKstmu"] = [=](const StandardModel& SM) { return new S_5(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["A_FB_BpKstmu"] = [=](const StandardModel& SM) { return new A_FB(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_7_BpKstmu"] = [=](const StandardModel& SM) { return new S_7(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_8_BpKstmu"] = [=](const StandardModel& SM) { return new S_8(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["S_9_BpKstmu"] = [=](const StandardModel& SM) { return new S_9(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };

    obsThFactory["P_1_BpKstmu"] = [=](const StandardModel& SM) { return new P_1(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_2_BpKstmu"] = [=](const StandardModel& SM) { return new P_2(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_3_BpKstmu"] = [=](const StandardModel& SM) { return new P_3(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_4p_BpKstmu"] = [=](const StandardModel& SM) { return new P_4Prime(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_5p_BpKstmu"] = [=](const StandardModel& SM) { return new P_5Prime(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_6p_BpKstmu"] = [=](const StandardModel& SM) { return new P_6Prime(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["P_8p_BpKstmu"] = [=](const StandardModel& SM) { return new P_8Prime(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };

    obsThFactory["V0_BdKstmu"] = [=](const StandardModel& SM) { return new V0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Vp_BdKstmu"] = [=](const StandardModel& SM) { return new Vp(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Vm_BdKstmu"] = [=](const StandardModel& SM) { return new Vm(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["T0_BdKstmu"] = [=](const StandardModel& SM) { return new T0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Tp_BdKstmu"] = [=](const StandardModel& SM) { return new Tp(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Tm_BdKstmu"] = [=](const StandardModel& SM) { return new Tm(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["S_BdKstmu"] = [=](const StandardModel& SM) { return new S(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["QCDfC9_1f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9_1f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["QCDfC9_2f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9_2f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["QCDfC9_3f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9_3f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["QCDfC9p_1f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9p_1f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["QCDfC9p_2f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9p_2f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["QCDfC9p_3f_BdKstmu"] = [=](const StandardModel& SM) { return new QCDfC9p_3f(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["Regtilde_1_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };
    obsThFactory["Regtilde_2_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };
    obsThFactory["Regtilde_3_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };

    obsThFactory["Imgtilde_1_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };
    obsThFactory["Imgtilde_2_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };
    obsThFactory["Imgtilde_3_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };

    obsThFactory["Absgtilde_1_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };
    obsThFactory["Absgtilde_2_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };
    obsThFactory["Absgtilde_3_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };

    obsThFactory["Arggtilde_1_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };
    obsThFactory["Arggtilde_2_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };
    obsThFactory["Arggtilde_3_BdKstmu"] = [=](const StandardModel& SM) { return new gtilde_3(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };

    obsThFactory["Reh_0_BdKstmu"] = [=](const StandardModel& SM) { return new h_0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };
    obsThFactory["Reh_p_BdKstmu"] = [=](const StandardModel& SM) { return new h_p(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };
    obsThFactory["Reh_m_BdKstmu"] = [=](const StandardModel& SM) { return new h_m(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1); };

    obsThFactory["Imh_0_BdKstmu"] = [=](const StandardModel& SM) { return new h_0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };
    obsThFactory["Imh_p_BdKstmu"] = [=](const StandardModel& SM) { return new h_p(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };
    obsThFactory["Imh_m_BdKstmu"] = [=](const StandardModel& SM) { return new h_m(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2); };

    obsThFactory["Absh_0_BdKstmu"] = [=](const StandardModel& SM) { return new h_0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };
    obsThFactory["Absh_p_BdKstmu"] = [=](const StandardModel& SM) { return new h_p(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };
    obsThFactory["Absh_m_BdKstmu"] = [=](const StandardModel& SM) { return new h_m(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3); };

    obsThFactory["Argh_0_BdKstmu"] = [=](const StandardModel& SM) { return new h_0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };
    obsThFactory["Argh_p_BdKstmu"] = [=](const StandardModel& SM) { return new h_p(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };
    obsThFactory["Argh_m_BdKstmu"] = [=](const StandardModel& SM) { return new h_m(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4); };

    obsThFactory["BR_MVpsi"] = [=](const StandardModel& SM) { return new BR_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["BR_MVpsi_ratio"] = [=](const StandardModel& SM) { return new BR_MVpsi_ratio(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Abs2Ampar_MVpsi"] = [=](const StandardModel& SM) { return new Abs2Ampar_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["ArgAmpar_MVpsi"] = [=](const StandardModel& SM) { return new ArgAmpar_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Abs2Amperp_MVpsi"] = [=](const StandardModel& SM) { return new Abs2Amperp_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["ArgAmperp_MVpsi"] = [=](const StandardModel& SM) { return new ArgAmperp_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["Abs2Ampzero_MVpsi"] = [=](const StandardModel& SM) { return new Abs2Ampzero_MVpsi(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    obsThFactory["unitarity_bound_1_BdKstmu"] = [=](const StandardModel& SM) { return new unitarity_bound_f_F1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["unitarity_bound_2_BdKstmu"] = [=](const StandardModel& SM) { return new unitarity_bound_g(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["unitarity_bound_3_BdKstmu"] = [=](const StandardModel& SM) { return new unitarity_bound_F2(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["unitarity_bound_4_BdKstmu"] = [=](const StandardModel& SM) { return new unitarity_bound_T1(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };
    obsThFactory["unitarity_bound_5_BdKstmu"] = [=](const StandardModel& SM) { return new unitarity_bound_T2_T0(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::MU); };

    //----- B+ to K*+ ll  -----
    obsThFactory["A_FB_BpKstmu"] = [=](const StandardModel& SM) { return new A_FB(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["F_L_BpKstmu"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };
    obsThFactory["BR_BpKstmu"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU); };

/* BEGIN: REMOVE FROM THE PACKAGE */
    //----- B to X_q ll -----
    obsThFactory["R_BXsee"] = [=](const StandardModel& SM) { return new R_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
    obsThFactory["HT_BXsee"] = [=](const StandardModel& SM) { return new HT_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
    obsThFactory["HL_BXsee"] = [=](const StandardModel& SM) { return new HL_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
    obsThFactory["HA_BXsee"] = [=](const StandardModel& SM) { return new HA_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
    obsThFactory["BR_BXsee"] = [=](const StandardModel& SM) { return new BR_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
    obsThFactory["AFB_BXsee"] = [=](const StandardModel& SM) { return new AFB_BXqll(SM, StandardModel::STRANGE, StandardModel::ELECTRON); };
/* END: REMOVE FROM THE PACKAGE */

    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ACP_BKstgamma"] = [=](const StandardModel& SM) { return new ACP_MVgamma(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["C_BKstgamma"] = [=](const StandardModel& SM) { return new C_MVgamma(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["S_BKstgamma"] = [=](const StandardModel& SM) { return new S_MVgamma(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ADG_BKstgamma"] = [=](const StandardModel& SM) { return new ADG_MVgamma(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["DC7_1"] = [=](const StandardModel& SM) { return new DC7_1(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["DC7_2"] = [=](const StandardModel& SM) { return new DC7_2(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["AbsDC7_L"] = [=](const StandardModel& SM) { return new AbsDC7_L(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["AbsDC7_R"] = [=](const StandardModel& SM) { return new AbsDC7_R(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ReDC7_L_Bd"] = [=](const StandardModel& SM) { return new ReDC7_L(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ReDC7_R_Bd"] = [=](const StandardModel& SM) { return new ReDC7_R(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ImDC7_L_Bd"] = [=](const StandardModel& SM) { return new ImDC7_L(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ImDC7_R_Bd"] = [=](const StandardModel& SM) { return new ImDC7_R(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["hp0_hm0"] = [=](const StandardModel& SM) { return new hp0_hm0(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["AbsDC7_QCDF_Bd"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["AbsDC7_QCDF_Bd_bar"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF_bar(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ReDC7_QCDF_Bd"] = [=](const StandardModel& SM) { return new ReDC7_QCDF(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ReDC7_QCDF_Bd_bar"] = [=](const StandardModel& SM) { return new ReDC7_QCDF_bar(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ImDC7_QCDF_Bd"] = [=](const StandardModel& SM) { return new ImDC7_QCDF(SM, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["ImDC7_QCDF_Bd_bar"] = [=](const StandardModel& SM) { return new ImDC7_QCDF_bar(SM, StandardModel::B_D, StandardModel::K_star); };


    //----- B+ to K*+ gamma  -----
    obsThFactory["BR_BpKstgamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ACP_BpKstgamma"] = [=](const StandardModel& SM) { return new C_MVgamma(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ReDC7_L_Bp"] = [=](const StandardModel& SM) { return new ReDC7_L(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ReDC7_R_Bp"] = [=](const StandardModel& SM) { return new ReDC7_R(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ImDC7_L_Bp"] = [=](const StandardModel& SM) { return new ImDC7_L(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ImDC7_R_Bp"] = [=](const StandardModel& SM) { return new ImDC7_R(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["AbsDC7_QCDF_Bp"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["AbsDC7_QCDF_Bp_bar"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF_bar(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ReDC7_QCDF_Bp"] = [=](const StandardModel& SM) { return new ReDC7_QCDF(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ReDC7_QCDF_Bp_bar"] = [=](const StandardModel& SM) { return new ReDC7_QCDF_bar(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ImDC7_QCDF_Bp"] = [=](const StandardModel& SM) { return new ImDC7_QCDF(SM, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["ImDC7_QCDF_Bp_bar"] = [=](const StandardModel& SM) { return new ImDC7_QCDF_bar(SM, StandardModel::B_P, StandardModel::K_star_P); };

    //----- B to PHI gamma  -----
    obsThFactory["BR_Bsphigamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["C_Bsphigamma"] = [=](const StandardModel& SM) { return new C_MVgamma(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["S_Bsphigamma"] = [=](const StandardModel& SM) { return new S_MVgamma(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ADG_Bsphigamma"] = [=](const StandardModel& SM) { return new ADG_MVgamma(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ReDC7_L_Bs"] = [=](const StandardModel& SM) { return new ReDC7_L(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ReDC7_R_Bs"] = [=](const StandardModel& SM) { return new ReDC7_R(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ImDC7_L_Bs"] = [=](const StandardModel& SM) { return new ImDC7_L(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ImDC7_R_Bs"] = [=](const StandardModel& SM) { return new ImDC7_R(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["AbsDC7_QCDF_Bs"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["AbsDC7_QCDF_Bs_bar"] = [=](const StandardModel& SM) { return new AbsDC7_QCDF_bar(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ReDC7_QCDF_Bs"] = [=](const StandardModel& SM) { return new ReDC7_QCDF(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ReDC7_QCDF_Bs_bar"] = [=](const StandardModel& SM) { return new ReDC7_QCDF_bar(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ImDC7_QCDF_Bs"] = [=](const StandardModel& SM) { return new ImDC7_QCDF(SM, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["ImDC7_QCDF_Bs_bar"] = [=](const StandardModel& SM) { return new ImDC7_QCDF_bar(SM, StandardModel::B_S, StandardModel::PHI); };

    //----- B to RHO gamma  -----
    obsThFactory["BR_Brhogamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_D, StandardModel::RHO); };
    obsThFactory["ACP_Brhogamma"] = [=](const StandardModel& SM) { return new ACP_MVgamma(SM, StandardModel::B_D, StandardModel::RHO); };
    obsThFactory["S_Brhogamma"] = [=](const StandardModel& SM) { return new S_MVgamma(SM, StandardModel::B_D, StandardModel::RHO); };

    //----- B+ to RHO+ gamma  -----
    obsThFactory["BR_Bprhogamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_P, StandardModel::RHO_P); };
    obsThFactory["ACP_Bprhogamma"] = [=](const StandardModel& SM) { return new ACP_MVgamma(SM, StandardModel::B_P, StandardModel::RHO_P); };

    //----- B to OMEGA gamma  -----
    obsThFactory["BR_Bomegagamma"] = [=](const StandardModel& SM) { return new BR_MVgamma(SM, StandardModel::B_D, StandardModel::OMEGA); };
    obsThFactory["ACP_Bomegagamma"] = [=](const StandardModel& SM) { return new ACP_MVgamma(SM, StandardModel::B_D, StandardModel::OMEGA); };

    //----- B to V gamma  -----
    obsThFactory["R_BVgamma"] = [=](const StandardModel& SM) { return new R_MVgamma(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::B_S, StandardModel::PHI); };
    obsThFactory["D0p_BKstgamma"] = [=](const StandardModel& SM) { return new D0p_MVgamma(SM, StandardModel::B_D, StandardModel::K_star, StandardModel::B_P, StandardModel::K_star_P); };
    obsThFactory["DACP_BKstgamma"] = [=](const StandardModel& SM) { return new DACP_MVgamma(SM, StandardModel::B_P, StandardModel::K_star_P, StandardModel::B_D, StandardModel::K_star); };
    obsThFactory["D0p_Brhogamma"] = [=](const StandardModel& SM) { return new D0p_MVgamma(SM, StandardModel::B_D, StandardModel::RHO, StandardModel::B_P, StandardModel::RHO_P); };
    obsThFactory["DACP_Brhogamma"] = [=](const StandardModel& SM) { return new DACP_MVgamma(SM, StandardModel::B_P, StandardModel::RHO_P, StandardModel::B_D, StandardModel::RHO); };

    //----- B to phi ll  -----
    obsThFactory["P_1_Bsphimu"] = [=](const StandardModel& SM) { return new P_1(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_2_Bsphimu"] = [=](const StandardModel& SM) { return new P_2(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_3_Bsphimu"] = [=](const StandardModel& SM) { return new P_3(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_4p_Bsphimu"] = [=](const StandardModel& SM) { return new P_4Prime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_5p_Bsphimu"] = [=](const StandardModel& SM) { return new P_5Prime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_6p_Bsphimu"] = [=](const StandardModel& SM) { return new P_6Prime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_8p_Bsphimu"] = [=](const StandardModel& SM) { return new P_8Prime(SM, StandardModel::B_D, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["Gammap_Bsphimu"] = [=](const StandardModel& SM) { return new GammaPrime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_FB_Bsphimu"] = [=](const StandardModel& SM) { return new A_FB(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["BR_Bsphimu"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["Rphi_Bsphill"] = [=](const StandardModel& SM) { return new R_MVll(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RphiL_Bsphill"] = [=](const StandardModel& SM) { return new RL_MVll(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RphiT_Bsphill"] = [=](const StandardModel& SM) { return new RT_MVll(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["R6_Bsphill"] = [=](const StandardModel& SM) { return new R_6(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["ACP_Bsphimu"] = [=](const StandardModel& SM) { return new ACP_MVll(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P3CP_Bsphimu"] = [=](const StandardModel& SM) { return new P3CP(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["F_L_Bsphimu"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["M_1p_Bsphimu"] = [=](const StandardModel& SM) { return new M_1Prime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["M_2p_Bsphimu"] = [=](const StandardModel& SM) { return new M_2Prime(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_3_Bsphimu"] = [=](const StandardModel& SM) { return new S_3(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_4_Bsphimu"] = [=](const StandardModel& SM) { return new S_4(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_5_Bsphimu"] = [=](const StandardModel& SM) { return new S_5(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_7_Bsphimu"] = [=](const StandardModel& SM) { return new S_7(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_8_Bsphimu"] = [=](const StandardModel& SM) { return new S_8(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["S_9_Bsphimu"] = [=](const StandardModel& SM) { return new S_9(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_5_Bsphimu"] = [=](const StandardModel& SM) { return new A_5(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_6_Bsphimu"] = [=](const StandardModel& SM) { return new A_6(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_6c_Bsphimu"] = [=](const StandardModel& SM) { return new A_6c(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_8_Bsphimu"] = [=](const StandardModel& SM) { return new A_8(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["A_9_Bsphimu"] = [=](const StandardModel& SM) { return new A_9(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["P_1_Bsphie"] = [=](const StandardModel& SM) { return new P_1(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON); };
    obsThFactory["P_2_Bsphie"] = [=](const StandardModel& SM) { return new P_2(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON); };
    obsThFactory["P_3_Bsphie"] = [=](const StandardModel& SM) { return new P_3(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON); };
    obsThFactory["F_L_Bsphie"] = [=](const StandardModel& SM) { return new F_L(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON); };
    obsThFactory["S_3_Bsphie"] = [=](const StandardModel& SM) { return new S_3(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON); };

    obsThFactory["unitarity_bound_1_Bsphimu"] = [=](const StandardModel& SM) { return new unitarity_bound_f_F1(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["unitarity_bound_2_Bsphimu"] = [=](const StandardModel& SM) { return new unitarity_bound_g(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["unitarity_bound_3_Bsphimu"] = [=](const StandardModel& SM) { return new unitarity_bound_F2(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["unitarity_bound_4_Bsphimu"] = [=](const StandardModel& SM) { return new unitarity_bound_T1(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    obsThFactory["unitarity_bound_5_Bsphimu"] = [=](const StandardModel& SM) { return new unitarity_bound_T2_T0(SM, StandardModel::B_S, StandardModel::PHI, StandardModel::MU); };
    
    //----- B+ to K+ ll  -----
    obsThFactory["BR_BpKmu"] = [=](const StandardModel& SM) { return new BR_MPll(SM, StandardModel::B_P, StandardModel::K_P, StandardModel::MU); };
    obsThFactory["BR_BpKe"] = [=](const StandardModel& SM) { return new BR_MPll(SM, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON); };
    obsThFactory["dBR_BpKmu"] = [=](const StandardModel& SM) { return new dBR_MPll(SM, StandardModel::B_P, StandardModel::K_P, StandardModel::MU); };
    obsThFactory["dBR_BpKe"] = [=](const StandardModel& SM) { return new dBR_MPll(SM, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON); };
    obsThFactory["RK_BpKll"] = [=](const StandardModel& SM) { return new R_MPll(SM, StandardModel::B_P, StandardModel::K_P, StandardModel::MU, StandardModel::ELECTRON); };

    //----- B0 to K0 ll  -----
    obsThFactory["BR_B0Kmu"] = [=](const StandardModel& SM) { return new BR_MPll(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::MU); };
    obsThFactory["BR_B0Ke"] = [=](const StandardModel& SM) { return new BR_MPll(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON); };
    obsThFactory["dBR_B0Kmu"] = [=](const StandardModel& SM) { return new dBR_MPll(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::MU); };
    obsThFactory["dBR_B0Ke"] = [=](const StandardModel& SM) { return new dBR_MPll(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON); };
    obsThFactory["RK_B0Kll"] = [=](const StandardModel& SM) { return new R_MPll(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::MU, StandardModel::ELECTRON); };

    obsThFactory["DC9_hlambda"] = [=](const StandardModel& SM) { return new DC9_hlambda(SM, StandardModel::B_D, StandardModel::K_0, StandardModel::MU); };
    
    //----- b to Xc lnu -----
    obsThFactory["Q2moment1_BClnu"] = [=](const StandardModel& SM) { return new Q2moments_BClnu(SM, StandardModel::B_D, 1); };
    obsThFactory["Q2moment2_BClnu"] = [=](const StandardModel& SM) { return new Q2moments_BClnu(SM, StandardModel::B_D, 2); };
    obsThFactory["Q2moment3_BClnu"] = [=](const StandardModel& SM) { return new Q2moments_BClnu(SM, StandardModel::B_D, 3); };
    obsThFactory["Elmoment1_BClnu"] = [=](const StandardModel& SM) { return new Elmoments_BClnu(SM, StandardModel::B_D, 1); };
    obsThFactory["Elmoment2_BClnu"] = [=](const StandardModel& SM) { return new Elmoments_BClnu(SM, StandardModel::B_D, 2); };
    obsThFactory["Elmoment3_BClnu"] = [=](const StandardModel& SM) { return new Elmoments_BClnu(SM, StandardModel::B_D, 3); };
    obsThFactory["MXmoment1_BClnu"] = [=](const StandardModel& SM) { return new MXmoments_BClnu(SM, StandardModel::B_D, 1); };
    obsThFactory["MXmoment2_BClnu"] = [=](const StandardModel& SM) { return new MXmoments_BClnu(SM, StandardModel::B_D, 2); };
    obsThFactory["MXmoment3_BClnu"] = [=](const StandardModel& SM) { return new MXmoments_BClnu(SM, StandardModel::B_D, 3); };

    obsThFactory["PartialAverageBR_BClnu"] = [=](const StandardModel& SM) { return new PartialAverageBR_BClnu(SM, StandardModel::B_D); };
    obsThFactory["Vcb_BClnu"] = [=](const StandardModel& SM) { return new Vcb_BClnu(SM, StandardModel::B_D); };
    
    //----- B to K nunu-----
    obsThFactory["BR_BpKpnunu"] = [=](const StandardModel& SM) { return new BR_MPll(SM, StandardModel::B_P, StandardModel::K_P, QCD::NEUTRINO_1); }; // MAYBE WE SHOULD DEFINE ALL_NEUTRINO ?
    obsThFactory["BR_BpKstarpnunu"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_P, StandardModel::K_star_P, QCD::NEUTRINO_1); }; // MAYBE WE SHOULD DEFINE ALL_NEUTRINO ?
    obsThFactory["BR_BKstarnunu"] = [=](const StandardModel& SM) { return new BR_MVll(SM, StandardModel::B_D, StandardModel::K_star, QCD::NEUTRINO_1); }; // MAYBE WE SHOULD DEFINE ALL_NEUTRINO ?
    obsThFactory["AT1_BpKstarpnunu"] = [=](const StandardModel& SM) { return new A_T1(SM, StandardModel::B_P, StandardModel::K_star_P, QCD::NEUTRINO_1); }; // MAYBE WE SHOULD DEFINE ALL_NEUTRINO ?
    obsThFactory["AT1_BKstarnunu"] = [=](const StandardModel& SM) { return new A_T1(SM, StandardModel::B_D, StandardModel::K_star, QCD::NEUTRINO_1); }; // MAYBE WE SHOULD DEFINE ALL_NEUTRINO ?
 
    //----- B to D*lnu -----
    obsThFactory["Gammaw_MVlnu"] = [=](const StandardModel& SM) { return new Gammaw_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RDstar_MVlnu"] = [=](const StandardModel& SM) { return new RDstar_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["Gammacl_MVlnu"] = [=](const StandardModel& SM) { return new Gammacl_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["GammacV_MVlnu"] = [=](const StandardModel& SM) { return new GammacV_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["Gammachi_MVlnu"] = [=](const StandardModel& SM) { return new Gammachi_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["UnitarityV_MVlnu"] = [=](const StandardModel& SM) { return new UnitarityV_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON); };
    obsThFactory["UnitarityA_MVlnu"] = [=](const StandardModel& SM) { return new UnitarityA_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON); };
    obsThFactory["UnitarityV_D_Dst"] = [=](const StandardModel& SM) { return new UnitarityV_D_Dst(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["hA1_at_w1"] = [=](const StandardModel& SM) { return new FF_hA1atw1(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["hV_w"] = [=](const StandardModel& SM) { return new FF_hV(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["hA1_w"] = [=](const StandardModel& SM) { return new FF_hA1(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["hA2_w"] = [=](const StandardModel& SM) { return new FF_hA2(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["hA3_w"] = [=](const StandardModel& SM) { return new FF_hA3(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["R1_w"] = [=](const StandardModel& SM) { return new FF_R1(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["R2_w"] = [=](const StandardModel& SM) { return new FF_R2(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["R0_w"] = [=](const StandardModel& SM) { return new FF_R0(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU); };
    obsThFactory["FL_MVtaunu"] = [=](const StandardModel& SM) { return new FL_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU); };
    obsThFactory["Ptau_MVtaunu"] = [=](const StandardModel& SM) { return new Plep_MVlnu(SM, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU); };
    //----- B to Dlnu -----
    obsThFactory["Gammaw_MPlnu"] = [=](const StandardModel& SM) { return new Gammaw_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["RD_MPlnu"] = [=](const StandardModel& SM) { return new RD_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON); };
    obsThFactory["UnitarityV_MPlnu"] = [=](const StandardModel& SM) { return new UnitarityV_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["UnitarityA_MPlnu"] = [=](const StandardModel& SM) { return new UnitarityA_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["Unitarity_Strong_MPlnu"] = [=](const StandardModel& SM) { return new Unitarity_Strong_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["af0_0"] = [=](const StandardModel& SM) { return new af0_0(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["FF0_MPlnu"] = [=](const StandardModel& SM) { return new FF0_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };
    obsThFactory["FFplus_MPlnu"] = [=](const StandardModel& SM) { return new FFplus_MPlnu(SM, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON); };

    //----- B to tau nu  -----
    obsThFactory["btaunu"] = [=](const StandardModel& SM) { return new Btaunu(SM, StandardModel::B_P); };
    obsThFactory["bctaunu"] = [=](const StandardModel& SM) { return new Btaunu(SM, StandardModel::B_C); };


    //----- D to lepton nu  -----
    obsThFactory["Dmuonnu"] = [=](const StandardModel& SM) { return new Dleptonnu(SM, StandardModel::D_P, StandardModel::MU); };
    obsThFactory["Dtaunu"] = [=](const StandardModel& SM) { return new Dleptonnu(SM, StandardModel::D_P, StandardModel::TAU); };

    obsThFactory["Dsmuonnu"] = [=](const StandardModel& SM) { return new Dleptonnu(SM, StandardModel::D_S, StandardModel::MU); };
    obsThFactory["Dstaunu"] = [=](const StandardModel& SM) { return new Dleptonnu(SM, StandardModel::D_S, StandardModel::TAU); };


    //----- K to muon nu / Pi to muon nu  -----
    obsThFactory["Kmunu_o_Pmunu"] = [](const StandardModel& SM) { return new Kmunu_o_Pmunu(SM); };

    //----- P to muon nu -----
    obsThFactory["Pmunu"] = [](const StandardModel& SM) { return new Pmunu(SM); };


    //----- tau to K nu / tau to Pi nu  -----
    obsThFactory["tauKnu_o_tauPnu"] = [](const StandardModel& SM) { return new tauKnu_o_tauPnu(SM); };


    /* BEGIN: REMOVE FROM THE PACKAGE */
    obsThFactory["Deltaamu"] = [](const StandardModel& SM) { return new Deltaamu(SM); };
    /* END: REMOVE FROM THE PACKAGE */

    //-----  Lepton Flavour observables  -----
    obsThFactory["mu_e_gamma"] = [](const StandardModel& SM) { return new mu_e_gamma(SM); };
    obsThFactory["log_meg"] = [](const StandardModel& SM) { return new log_meg(SM); };
    obsThFactory["tau_mu_gamma"] = [](const StandardModel& SM) { return new tau_mu_gamma(SM); };
    obsThFactory["log_tmg"] = [](const StandardModel& SM) { return new log_tmg(SM); };
    obsThFactory["tau_e_gamma"] = [](const StandardModel& SM) { return new tau_e_gamma(SM); };
    obsThFactory["log_teg"] = [](const StandardModel& SM) { return new log_teg(SM); };
    obsThFactory["mu_3e"] = [](const StandardModel& SM) { return new mu_3e(SM); };
    obsThFactory["tau_3mu"] = [](const StandardModel& SM) { return new tau_3mu(SM); };
    obsThFactory["tau_3e"] = [](const StandardModel& SM) { return new tau_3e(SM); };
    obsThFactory["gminus2_mu"] = [](const StandardModel& SM) { return new gminus2_mu(SM); };
    obsThFactory["Robs_mu_e_gamma"] = [](const StandardModel& SM) { return new Robs_mu_e_gamma(SM); };
    obsThFactory["Robs_tau_mu_gamma"] = [](const StandardModel& SM) { return new Robs_tau_mu_gamma(SM); };
    obsThFactory["Robs_tau_mu_gamma_BelleII"] = [](const StandardModel& SM) { return new Robs_tau_mu_gamma_BelleII(SM); };
    obsThFactory["Robs_tau_e_gamma"] = [](const StandardModel& SM) { return new Robs_tau_e_gamma(SM); };
    obsThFactory["mueconversion_Ti"] = [](const StandardModel& SM) { return new mueconversion_Ti(SM); };

    obsThFactory["deltaRL_12_u"] = [](const StandardModel& SM) { return new deltaRL_12_u(SM); };
    obsThFactory["deltaRL_13_u"] = [](const StandardModel& SM) { return new deltaRL_13_u(SM); };
    obsThFactory["deltaRL_23_u"] = [](const StandardModel& SM) { return new deltaRL_23_u(SM); };
    obsThFactory["deltaRL_12_e"] = [](const StandardModel& SM) { return new deltaRL_12_e(SM); };
    obsThFactory["deltaRL_21_e"] = [](const StandardModel& SM) { return new deltaRL_21_e(SM); };
    obsThFactory["deltaRL_13_e"] = [](const StandardModel& SM) { return new deltaRL_13_e(SM); };
    obsThFactory["deltaRL_31_e"] = [](const StandardModel& SM) { return new deltaRL_31_e(SM); };
    obsThFactory["deltaRL_23_e"] = [](const StandardModel& SM) { return new deltaRL_23_e(SM); };
    obsThFactory["deltaRL_32_e"] = [](const StandardModel& SM) { return new deltaRL_32_e(SM); };

    obsThFactory["deltaLL1_q"] = [](const StandardModel& SM) { return new deltaLL1_q(SM); };
    obsThFactory["deltaLL2_q"] = [](const StandardModel& SM) { return new deltaLL2_q(SM); };
    obsThFactory["deltaLL3_q"] = [](const StandardModel& SM) { return new deltaLL3_q(SM); };
    obsThFactory["deltaRR1_u"] = [](const StandardModel& SM) { return new deltaRR1_u(SM); };
    obsThFactory["deltaRR2_u"] = [](const StandardModel& SM) { return new deltaRR2_u(SM); };
    obsThFactory["deltaRR3_u"] = [](const StandardModel& SM) { return new deltaRR3_u(SM); };
    obsThFactory["deltaRR1_d"] = [](const StandardModel& SM) { return new deltaRR1_d(SM); };
    obsThFactory["deltaRR2_d"] = [](const StandardModel& SM) { return new deltaRR2_d(SM); };
    obsThFactory["deltaRR3_d"] = [](const StandardModel& SM) { return new deltaRR3_d(SM); };
    obsThFactory["deltaLL1_l"] = [](const StandardModel& SM) { return new deltaLL1_l(SM); };
    obsThFactory["deltaLL2_l"] = [](const StandardModel& SM) { return new deltaLL2_l(SM); };
    obsThFactory["deltaLL3_l"] = [](const StandardModel& SM) { return new deltaLL3_l(SM); };
    obsThFactory["deltaRR1_e"] = [](const StandardModel& SM) { return new deltaRR1_e(SM); };
    obsThFactory["deltaRR2_e"] = [](const StandardModel& SM) { return new deltaRR2_e(SM); };
    obsThFactory["deltaRR3_e"] = [](const StandardModel& SM) { return new deltaRR3_e(SM); };

    obsThFactory["CCBu11"] = [](const StandardModel& SM) { return new CCBu11(SM); };
    obsThFactory["CCBu22"] = [](const StandardModel& SM) { return new CCBu22(SM); };
    obsThFactory["CCBu33"] = [](const StandardModel& SM) { return new CCBu33(SM); };
    obsThFactory["CCBu12"] = [](const StandardModel& SM) { return new CCBu12(SM); };
    obsThFactory["CCBu13"] = [](const StandardModel& SM) { return new CCBu13(SM); };
    obsThFactory["CCBu23"] = [](const StandardModel& SM) { return new CCBu23(SM); };
    obsThFactory["CCBd11"] = [](const StandardModel& SM) { return new CCBd11(SM); };
    obsThFactory["CCBd22"] = [](const StandardModel& SM) { return new CCBd22(SM); };
    obsThFactory["CCBd33"] = [](const StandardModel& SM) { return new CCBd33(SM); };
    obsThFactory["CCBd12"] = [](const StandardModel& SM) { return new CCBd12(SM); };
    obsThFactory["CCBd13"] = [](const StandardModel& SM) { return new CCBd13(SM); };
    obsThFactory["CCBd23"] = [](const StandardModel& SM) { return new CCBd23(SM); };
    obsThFactory["CCBe11"] = [](const StandardModel& SM) { return new CCBe11(SM); };
    obsThFactory["CCBe22"] = [](const StandardModel& SM) { return new CCBe22(SM); };
    obsThFactory["CCBe33"] = [](const StandardModel& SM) { return new CCBe33(SM); };
    obsThFactory["CCBe12"] = [](const StandardModel& SM) { return new CCBe12(SM); };
    obsThFactory["CCBe13"] = [](const StandardModel& SM) { return new CCBe13(SM); };
    obsThFactory["CCBe23"] = [](const StandardModel& SM) { return new CCBe23(SM); };

    obsThFactory["VacuumTunnelingRate"] = [](const StandardModel& SM) { return new FindAction(SM); };

    obsThFactory["logdeltaRL_13_e"] = [](const StandardModel& SM) { return new logdeltaRL_13_e(SM); };
    obsThFactory["logdeltaRL_23_e"] = [](const StandardModel& SM) { return new logdeltaRL_23_e(SM); };
    obsThFactory["logmslepton"] = [](const StandardModel& SM) { return new logmslepton(SM); };
    obsThFactory["mslepton"] = [](const StandardModel& SM) { return new mslepton(SM); };
    obsThFactory["deltaTEhat23"] = [](const StandardModel& SM) { return new deltaTEhat23(SM); };
    obsThFactory["deltaLLRR_l"] = [](const StandardModel& SM) { return new deltaLLRR_l(SM); };

    //-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["DmBd"] = [](const StandardModel& SM) { return new DmBd(SM); };
    obsThFactory["DmBs"] = [](const StandardModel& SM) { return new DmBs(SM); };
    obsThFactory["RmBd"] = [](const StandardModel& SM) { return new RmBd(SM); };
    obsThFactory["RmBs"] = [](const StandardModel& SM) { return new RmBs(SM); };
    obsThFactory["CBd"] = [](const StandardModel& SM) { return new CBd(SM); };
    obsThFactory["CBs"] = [](const StandardModel& SM) { return new CBs(SM); };
    obsThFactory["PhiBd"] = [](const StandardModel& SM) { return new PhiBd(SM); };
    obsThFactory["PhiBs"] = [](const StandardModel& SM) { return new PhiBs(SM); };
    obsThFactory["BBd"] = [](const StandardModel& SM) { return new BBd(SM); };
    obsThFactory["FBd"] = [](const StandardModel& SM) { return new FBd(SM); };
    obsThFactory["FBdSqrtBBd"] = [](const StandardModel& SM) { return new FBdSqrtBBd(SM); };
    obsThFactory["FBsSqrtBBs"] = [](const StandardModel& SM) { return new FBsSqrtBBs(SM); };
    obsThFactory["xi"] = [](const StandardModel& SM) { return new xi(SM); };
    obsThFactory["alpha"] = [](const StandardModel& SM) { return new Alpha(SM); };
    obsThFactory["alpha_2a"] = [](const StandardModel& SM) { return new Alpha_2a(SM); };
    obsThFactory["SJPsiK"] = [](const StandardModel& SM) { return new SJPsiK(SM); };
    obsThFactory["C2beta"] = [](const StandardModel& SM) { return new C2beta(SM); };
    obsThFactory["Phis_JPsiPhi"] = [](const StandardModel& SM) { return new Phis_JPsiPhi(SM); };
    obsThFactory["EpsilonK"] = [](const StandardModel& SM) { return new EpsilonK(SM); };
    obsThFactory["DmK"] = [](const StandardModel& SM) { return new DmK(SM); };
    obsThFactory["ImADC2"] = [](const StandardModel& SM) { return new ImADC2(SM); };
    /* BEGIN: REMOVE FROM THE PACKAGE */
    obsThFactory["M12D"] = [](const StandardModel& SM) { return new M12D(SM); };
    obsThFactory["ArgD"] = [](const StandardModel& SM) { return new ArgD(SM); };
    //semileptonic asymmetry
    obsThFactory["Asl_d_pole"] = [](const StandardModel& SM) { return new Asl_d_pole(SM); };
    obsThFactory["Asl_s_pole"] = [](const StandardModel& SM) { return new Asl_s_pole(SM); };
    obsThFactory["Asl_s_pole_NLO"] = [](const StandardModel& SM) { return new Asl_s_pole_NLO(SM); };
    obsThFactory["Asl_s_pole_LO"] = [](const StandardModel& SM) { return new Asl_s_pole_LO(SM); };
    obsThFactory["Asl_d_MSbar"] = [](const StandardModel& SM) { return new Asl_d_MSbar(SM); };
    obsThFactory["Asl_s_MSbar"] = [](const StandardModel& SM) { return new Asl_s_MSbar(SM); };
    obsThFactory["Asl_s_MSbar_NLO"] = [](const StandardModel& SM) { return new Asl_s_MSbar_NLO(SM); };
    obsThFactory["Asl_s_MSbar_LO"] = [](const StandardModel& SM) { return new Asl_s_MSbar_LO(SM); };
    obsThFactory["Asl_d_PS"] = [](const StandardModel& SM) { return new Asl_d_PS(SM); };
    obsThFactory["Asl_s_PS"] = [](const StandardModel& SM) { return new Asl_s_PS(SM); };
    obsThFactory["Asl_s_PS_NLO"] = [](const StandardModel& SM) { return new Asl_s_PS_NLO(SM); };
    obsThFactory["Asl_s_PS_LO"] = [](const StandardModel& SM) { return new Asl_s_PS_LO(SM); };
    obsThFactory["Asl_s_pole_fixmub"] = [](const StandardModel& SM) { return new Asl_s_pole_fixmub(SM); };
    obsThFactory["Asl_s_MSbar_fixmub"] = [](const StandardModel& SM) { return new Asl_s_MSbar_fixmub(SM); };
    obsThFactory["Asl_s_PS_fixmub"] = [](const StandardModel& SM) { return new Asl_s_PS_fixmub(SM); };
    obsThFactory["Asl_d_only1overmb"] = [](const StandardModel& SM) { return new Asl_d_only1overmb(SM); };
    obsThFactory["Asl_s_only1overmb"] = [](const StandardModel& SM) { return new Asl_s_only1overmb(SM); };
    obsThFactory["Asl_d_MSbar_NLO_tradBasis"] = [](const StandardModel& SM) { return new Asl_d_MSbar_NLO_tradBasis(SM); };
    obsThFactory["Asl_s_MSbar_NLO_tradBasis"] = [](const StandardModel& SM) { return new Asl_s_MSbar_NLO_tradBasis(SM); };
    obsThFactory["Asl_s_MSbar_NLO_RI"] = [](const StandardModel& SM) { return new Asl_s_MSbar_NLO_RI(SM); };
    obsThFactory["Asl_s_MSbar_NLO_RI_tradBasis"] = [](const StandardModel& SM) { return new Asl_s_MSbar_NLO_RI_tradBasis(SM); };
    obsThFactory["Asl_s_PS_NLO_RI"] = [](const StandardModel& SM) { return new Asl_s_PS_NLO_RI(SM); };
    obsThFactory["Asl_s_MSbar_partialNNLO"] = [](const StandardModel& SM) { return new Asl_s_MSbar_partialNNLO(SM); };
    obsThFactory["Asl_s_PS_partialNNLO"] = [](const StandardModel& SM) { return new Asl_s_PS_partialNNLO(SM); };
    obsThFactory["Asl_s_MSbar_partialN3LO"] = [](const StandardModel& SM) { return new Asl_s_MSbar_partialN3LO(SM); };
    obsThFactory["Asl_s_PS_partialN3LO"] = [](const StandardModel& SM) { return new Asl_s_PS_partialN3LO(SM); };
    //DeltaGamma
    obsThFactory["DGamma_d_pole"] = [](const StandardModel& SM) { return new DGamma_d_pole(SM); };
    obsThFactory["DGamma_s_pole"] = [](const StandardModel& SM) { return new DGamma_s_pole(SM); };
    obsThFactory["DGamma_s_pole_NLO"] = [](const StandardModel& SM) { return new DGamma_s_pole_NLO(SM); };
    obsThFactory["DGamma_s_pole_LO"] = [](const StandardModel& SM) { return new DGamma_s_pole_LO(SM); };
    obsThFactory["DGamma_d_MSbar"] = [](const StandardModel& SM) { return new DGamma_d_MSbar(SM); };
    obsThFactory["DGamma_s_MSbar"] = [](const StandardModel& SM) { return new DGamma_s_MSbar(SM); };
    obsThFactory["DGamma_s_MSbar_NLO"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_NLO(SM); };
    obsThFactory["DGamma_s_MSbar_LO"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_LO(SM); };
    obsThFactory["DGamma_d_PS"] = [](const StandardModel& SM) { return new DGamma_d_PS(SM); };
    obsThFactory["DGamma_s_PS"] = [](const StandardModel& SM) { return new DGamma_s_PS(SM); };
    obsThFactory["DGamma_s_PS_NLO"] = [](const StandardModel& SM) { return new DGamma_s_PS_NLO(SM); };
    obsThFactory["DGamma_s_PS_LO"] = [](const StandardModel& SM) { return new DGamma_s_PS_LO(SM); };
    obsThFactory["DGamma_s_pole_fixmub"] = [](const StandardModel& SM) { return new DGamma_s_pole_fixmub(SM); };
    obsThFactory["DGamma_s_MSbar_fixmub"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_fixmub(SM); };
    obsThFactory["DGamma_s_PS_fixmub"] = [](const StandardModel& SM) { return new DGamma_s_PS_fixmub(SM); };
    obsThFactory["DGamma_d_only1overmb"] = [](const StandardModel& SM) { return new DGamma_d_only1overmb(SM); };
    obsThFactory["DGamma_s_only1overmb"] = [](const StandardModel& SM) { return new DGamma_s_only1overmb(SM); };
    obsThFactory["DGamma_d_NLO_tradBasis"] = [](const StandardModel& SM) { return new DGamma_d_NLO_tradBasis(SM); };
    obsThFactory["DGamma_s_NLO_tradBasis"] = [](const StandardModel& SM) { return new DGamma_s_NLO_tradBasis(SM); };
    obsThFactory["DGamma_d_LO_tradBasis"] = [](const StandardModel& SM) { return new DGamma_d_LO_tradBasis(SM); };
    obsThFactory["DGamma_s_LO_tradBasis"] = [](const StandardModel& SM) { return new DGamma_s_LO_tradBasis(SM); };
    obsThFactory["DGamma_s_MSbar_NLO_RI"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_NLO_RI(SM); };
    obsThFactory["DGamma_s_MSbar_NLO_RI_tradBasis"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_NLO_RI_tradBasis(SM); };
    obsThFactory["DGamma_s_PS_NLO_RI"] = [](const StandardModel& SM) { return new DGamma_s_PS_NLO_RI(SM); };
    obsThFactory["DGamma_s_MSbar_partialNNLO"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_partialNNLO(SM); };
    obsThFactory["DGamma_s_PS_partialNNLO"] = [](const StandardModel& SM) { return new DGamma_s_PS_partialNNLO(SM); };
    obsThFactory["DGamma_s_MSbar_partialN3LO"] = [](const StandardModel& SM) { return new DGamma_s_MSbar_partialN3LO(SM); };
    obsThFactory["DGamma_s_PS_partialN3LO"] = [](const StandardModel& SM) { return new DGamma_s_PS_partialN3LO(SM); };
    //----- eps'/eps  -----
    obsThFactory["EpsilonP_O_Epsilon_ReA0EXP"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 0); };
    obsThFactory["EpsilonP_O_Epsilon_ReA2EXP"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 1); };
    obsThFactory["EpsilonP_O_Epsilon_ReA0LAT"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 2); };
    obsThFactory["EpsilonP_O_Epsilon_ReA2LAT"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 3); };
    obsThFactory["EpsilonP_O_Epsilon"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 4); };
    obsThFactory["EpsilonP_O_Epsilon_pureLAT"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 5); };
    obsThFactory["EpsilonP_O_Epsilon_ImA0"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 6); };
    obsThFactory["EpsilonP_O_Epsilon_ImA0LAT"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 7); };
    obsThFactory["EpsilonP_O_Epsilon_ImA2"] = [=](const StandardModel& SM) { return new EpsilonP_O_Epsilon(SM, 8); };
    obsThFactory["EpsilonP_O_Epsilon_TH"] = [](const StandardModel& SM) { return new EpsilonP_O_Epsilon_TH(SM); };
    // 1.3 GeV --> scale as in Table 5 of 1909.05610
    obsThFactory["EpsilonP_O_Epsilon_z1"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 0, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_z2"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 1, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y3"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 2, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y4"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 3, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y5"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 4, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y6"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 5, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y7"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 6, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y8"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 7, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y9"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 8, 1.3); };
    obsThFactory["EpsilonP_O_Epsilon_y10"] = [=](const StandardModel& SM) { return new WC_epspOeps(SM, 9, 1.3); };
    obsThFactory["BR_Kppnunu"] = [](const StandardModel& SM) { return new BR_Kppnunu(SM); };
    obsThFactory["BR_Kp0nunu"] = [](const StandardModel& SM) { return new BR_Kp0nunu(SM); };

    /* END: REMOVE FROM THE PACKAGE */
    //----- CKM  -----
    obsThFactory["Vud"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 1); };
    obsThFactory["Vus"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 2); };
    obsThFactory["Vub"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 3); };
    obsThFactory["Vcd"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 1); };
    obsThFactory["Vcs"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 2); };
    obsThFactory["Vcb"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 3); };
    obsThFactory["Vtd"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 1); };
    obsThFactory["Vts"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 2); };
    obsThFactory["Vtb"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 3); };
    obsThFactory["argVud"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 1, 1); };
    obsThFactory["argVus"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 2, 1); };
    obsThFactory["argVub"] = [=](const StandardModel& SM) { return new VCKM(SM, 1, 3, 1); };
    obsThFactory["argVcd"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 1, 1); };
    obsThFactory["argVcs"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 2, 1); };
    obsThFactory["argVcb"] = [=](const StandardModel& SM) { return new VCKM(SM, 2, 3, 1); };
    obsThFactory["argVtd"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 1, 1); };
    obsThFactory["argVts"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 2, 1); };
    obsThFactory["argVtb"] = [=](const StandardModel& SM) { return new VCKM(SM, 3, 3, 1); };
    obsThFactory["CKM_alpha"] = [](const StandardModel& SM) { return new CKM_Alpha(SM); };
    obsThFactory["CKM_gamma"] = [](const StandardModel& SM) { return new CKM_Gamma(SM); };
    obsThFactory["CKM_beta"] = [](const StandardModel& SM) { return new CKM_Beta(SM); };
    obsThFactory["CKM_betas"] = [](const StandardModel& SM) { return new CKM_Betas(SM); };
    obsThFactory["CKM_2betapgamma"] = [](const StandardModel& SM) { return new CKM_2BpG(SM); };
    obsThFactory["CKM_s2beta"] = [](const StandardModel& SM) { return new CKM_S2Beta(SM); };
    obsThFactory["CKM_c2beta"] = [](const StandardModel& SM) { return new CKM_C2Beta(SM); };
    obsThFactory["CKM_rho"] = [](const StandardModel& SM) { return new CKM_rho(SM); };
    obsThFactory["CKM_eta"] = [](const StandardModel& SM) { return new CKM_eta(SM); };
    obsThFactory["CKM_sintheta12"] = [](const StandardModel& SM) { return new CKM_SinTheta12(SM); };
    obsThFactory["CKM_sintheta13"] = [](const StandardModel& SM) { return new CKM_SinTheta13(SM); };
    obsThFactory["CKM_sintheta23"] = [](const StandardModel& SM) { return new CKM_SinTheta23(SM); };
    obsThFactory["CKM_delta"] = [](const StandardModel& SM) { return new CKM_Delta(SM); };
    obsThFactory["J_CP"] = [](const StandardModel& SM) { return new J_CP(SM); };
    obsThFactory["Rt"] = [](const StandardModel& SM) { return new CKM_Rt(SM); };
    obsThFactory["Rt_dms"] = [](const StandardModel& SM) { return new CKM_Rt_dms(SM); };
    obsThFactory["Rts"] = [](const StandardModel& SM) { return new CKM_Rts(SM); };
    obsThFactory["Rb"] = [](const StandardModel& SM) { return new CKM_Rb(SM); };
    obsThFactory["VtdoVts"] = [](const StandardModel& SM) { return new CKM_VtdoVts(SM); };
    obsThFactory["Abslam_t"] = [](const StandardModel& SM) { return new Abslam_t(SM); };
    obsThFactory["Abslam_c"] = [](const StandardModel& SM) { return new Abslam_c(SM); };
    obsThFactory["Abslam_u"] = [](const StandardModel& SM) { return new Abslam_u(SM); };
    obsThFactory["Abslam_td"] = [](const StandardModel& SM) { return new Abslam_td(SM); };
    obsThFactory["Abslam_cd"] = [](const StandardModel& SM) { return new Abslam_cd(SM); };
    obsThFactory["Abslam_ud"] = [](const StandardModel& SM) { return new Abslam_ud(SM); };
    obsThFactory["Abslam_ts"] = [](const StandardModel& SM) { return new Abslam_ts(SM); };
    obsThFactory["Abslam_cs"] = [](const StandardModel& SM) { return new Abslam_cs(SM); };
    obsThFactory["Abslam_us"] = [](const StandardModel& SM) { return new Abslam_us(SM); };
    obsThFactory["Relam_t"] = [](const StandardModel& SM) { return new Relam_t(SM); };
    obsThFactory["Relam_c"] = [](const StandardModel& SM) { return new Relam_c(SM); };
    obsThFactory["Relam_u"] = [](const StandardModel& SM) { return new Relam_u(SM); };
    obsThFactory["Relam_td"] = [](const StandardModel& SM) { return new Relam_td(SM); };
    obsThFactory["Relam_cd"] = [](const StandardModel& SM) { return new Relam_cd(SM); };
    obsThFactory["Relam_ud"] = [](const StandardModel& SM) { return new Relam_ud(SM); };
    obsThFactory["Relam_ts"] = [](const StandardModel& SM) { return new Relam_ts(SM); };
    obsThFactory["Relam_cs"] = [](const StandardModel& SM) { return new Relam_cs(SM); };
    obsThFactory["Relam_us"] = [](const StandardModel& SM) { return new Relam_us(SM); };
    obsThFactory["Imlam_t"] = [](const StandardModel& SM) { return new Imlam_t(SM); };
    obsThFactory["Imlam_c"] = [](const StandardModel& SM) { return new Imlam_c(SM); };
    obsThFactory["Imlam_u"] = [](const StandardModel& SM) { return new Imlam_u(SM); };
    obsThFactory["Imlam_td"] = [](const StandardModel& SM) { return new Imlam_td(SM); };
    obsThFactory["Imlam_cd"] = [](const StandardModel& SM) { return new Imlam_cd(SM); };
    obsThFactory["Imlam_ud"] = [](const StandardModel& SM) { return new Imlam_ud(SM); };
    obsThFactory["Imlam_ts"] = [](const StandardModel& SM) { return new Imlam_ts(SM); };
    obsThFactory["Imlam_cs"] = [](const StandardModel& SM) { return new Imlam_cs(SM); };
    obsThFactory["Imlam_us"] = [](const StandardModel& SM) { return new Imlam_us(SM); };

}
