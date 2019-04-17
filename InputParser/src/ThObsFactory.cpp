/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObservable.h"
#include "ParamObs.h"
#include "EWObservables.h"
#include "FlavourObservables.h"
#include "alpha_s.h"
#include "MtMSbar.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

ThObsFactory::ThObsFactory()
{
    //
    //-----  StandardModel observables  -----
    obsThFactory["MtMSbar"] = boost::factory<MtMSbar*>();
    obsThFactory["alpha_s_LO"] = boost::bind(boost::factory<alpha_s*>(), _1, LO);
    obsThFactory["alpha_s_NLO"] = boost::bind(boost::factory<alpha_s*>(), _1, NLO);
    obsThFactory["alpha_s_FULLNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, FULLNLO);
    obsThFactory["alpha_s_NNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, NNLO);
    obsThFactory["alpha_s_FULLNNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, FULLNNLO);
    //-----  Electroweak precision observables  -----
    obsThFactory["alphaMz"] = boost::factory<AlphaEmMz*>();
    obsThFactory["Mw"] = boost::factory<Mw*>();
    obsThFactory["GammaW"] = boost::factory<GammaW*>();
    obsThFactory["BrWlepton"] = boost::factory<BrWlepton*>();
    obsThFactory["BrWelectron"] = boost::factory<BrWelectron*>();
    obsThFactory["BrWmuon"] = boost::factory<BrWmuon*>();
    obsThFactory["BrWtau"] = boost::factory<BrWtau*>();
    obsThFactory["BrWhadrons"] = boost::factory<BrWhadrons*>();
    obsThFactory["RWc"] = boost::factory<RWc*>();
    obsThFactory["GammaZ"] = boost::factory<GammaZ*>();
    obsThFactory["sigmaHadron"] = boost::factory<sigmaHadron*>();
    obsThFactory["sin2thetaEff"] = boost::factory<sin2thetaEff*>();
    obsThFactory["sin2thetaEffelectron"] = boost::factory<sin2thetaEffel*>();
    obsThFactory["sin2thetaEffmuon"] = boost::factory<sin2thetaEffmu*>();
    obsThFactory["PtauPol"] = boost::factory<PtauPol*>();
    obsThFactory["Alepton"] = boost::factory<Alepton*>();
    obsThFactory["Aelectron"] = boost::factory<Aelectron*>();
    obsThFactory["Amuon"] = boost::factory<Amuon*>();
    obsThFactory["Atau"] = boost::factory<Atau*>();
    obsThFactory["Astrange"] = boost::factory<Astrange*>();
    obsThFactory["Acharm"] = boost::factory<Acharm*>();
    obsThFactory["Abottom"] = boost::factory<Abottom*>();
    obsThFactory["AFBlepton"] = boost::factory<AFBlepton*>();
    obsThFactory["AFBelectron"] = boost::factory<AFBelectron*>();
    obsThFactory["AFBmuon"] = boost::factory<AFBmuon*>();
    obsThFactory["AFBtau"] = boost::factory<AFBtau*>();
    obsThFactory["AFBstrange"] = boost::factory<AFBstrange*>();
    obsThFactory["AFBcharm"] = boost::factory<AFBcharm*>();
    obsThFactory["AFBbottom"] = boost::factory<AFBbottom*>();
    obsThFactory["Nneutrinos"] = boost::factory<Nneutrinos*>();
    obsThFactory["Rlepton"] = boost::factory<Rlepton*>();
    obsThFactory["Relectron"] = boost::factory<Relectron*>();
    obsThFactory["Rmuon"] = boost::factory<Rmuon*>();
    obsThFactory["Rtau"] = boost::factory<Rtau*>();
    obsThFactory["Rinv"] = boost::factory<Rinv*>();
    obsThFactory["Ruc"] = boost::factory<Ruc*>();
    obsThFactory["Rcharm"] = boost::factory<Rcharm*>();
    obsThFactory["Rbottom"] = boost::factory<Rbottom*>();

    //-----  Epsilon parameters  -----
    obsThFactory["epsilon1"] = boost::factory<Epsilon1*>();
    obsThFactory["epsilon2"] = boost::factory<Epsilon2*>();
    obsThFactory["epsilon3"] = boost::factory<Epsilon3*>();
    obsThFactory["epsilonb"] = boost::factory<Epsilonb*>();
          

    //-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["DmBd"] = boost::factory<DmBd*>();
    obsThFactory["DmBs"] = boost::factory<DmBs*>();
    obsThFactory["RmBs"] = boost::factory<RmBs*>();
    obsThFactory["SJPsiK"] = boost::factory<SJPsiK*>();
    obsThFactory["Betas_JPsiPhi"] = boost::factory<Betas_JPsiPhi*>();
    obsThFactory["EpsilonK"] = boost::factory<EpsilonK*>();
    obsThFactory["DmK"] = boost::factory<DmK*>();
    //----- CKM  -----
    obsThFactory["Vud"] = boost::bind(boost::factory<VCKM*>(), _1, 1, 1);
    obsThFactory["Vus"] = boost::bind(boost::factory<VCKM*>(), _1, 1, 2);
    obsThFactory["Vub"] = boost::bind(boost::factory<VCKM*>(), _1, 1, 3);
    obsThFactory["Vcd"] = boost::bind(boost::factory<VCKM*>(), _1, 2, 1);
    obsThFactory["Vcs"] = boost::bind(boost::factory<VCKM*>(), _1, 2, 2);
    obsThFactory["Vcb"] = boost::bind(boost::factory<VCKM*>(), _1, 2, 3);
    obsThFactory["Vtd"] = boost::bind(boost::factory<VCKM*>(), _1, 3, 1);
    obsThFactory["Vts"] = boost::bind(boost::factory<VCKM*>(), _1, 3, 2);
    obsThFactory["Vtb"] = boost::bind(boost::factory<VCKM*>(), _1, 3, 3);
    obsThFactory["alpha"] = boost::factory<CKM_Alpha*>();
    obsThFactory["alpha_2a"] = boost::factory<Alpha_2a*>();
    obsThFactory["gamma"] = boost::factory<CKM_Gamma*>();
    obsThFactory["beta"] = boost::factory<CKM_Beta*>();
    obsThFactory["betas"] = boost::factory<CKM_Betas*>();
    obsThFactory["2betapgamma"] = boost::factory<CKM_2BpG*>();
    obsThFactory["s2beta"] = boost::factory<CKM_S2Beta*>();
    obsThFactory["c2beta"] = boost::factory<CKM_C2Beta*>();
    obsThFactory["CKM_rho"] = boost::factory<CKM_rho*>();
    obsThFactory["CKM_eta"] = boost::factory<CKM_eta*>();
    obsThFactory["sintheta12"] = boost::factory<CKM_SinTheta12*>();
    obsThFactory["sintheta13"] = boost::factory<CKM_SinTheta13*>();
    obsThFactory["sintheta23"] = boost::factory<CKM_SinTheta23*>();
    obsThFactory["ckmdelta"] = boost::factory<CKM_Delta*>();
    obsThFactory["J_CP"] = boost::factory<J_CP*>();
    obsThFactory["Rt"] = boost::factory<CKM_Rt*>();
    obsThFactory["Rts"] = boost::factory<CKM_Rts*>();
    obsThFactory["Rb"] = boost::factory<CKM_Rb*>();
    obsThFactory["VtdoVts"] = boost::factory<CKM_VtdoVts*>();
    obsThFactory["Abslam_t"] = boost::factory<Abslam_t*>();
    obsThFactory["Abslam_c"] = boost::factory<Abslam_c*>();
    obsThFactory["Abslam_u"] = boost::factory<Abslam_u*>();
    obsThFactory["Abslam_td"] = boost::factory<Abslam_td*>();
    obsThFactory["Abslam_cd"] = boost::factory<Abslam_cd*>();
    obsThFactory["Abslam_ud"] = boost::factory<Abslam_ud*>();
    obsThFactory["Abslam_ts"] = boost::factory<Abslam_ts*>();
    obsThFactory["Abslam_cs"] = boost::factory<Abslam_cs*>();
    obsThFactory["Abslam_us"] = boost::factory<Abslam_us*>();
    obsThFactory["Relam_t"] = boost::factory<Relam_t*>();
    obsThFactory["Relam_c"] = boost::factory<Relam_c*>();
    obsThFactory["Relam_u"] = boost::factory<Relam_u*>();
    obsThFactory["Relam_td"] = boost::factory<Relam_td*>();
    obsThFactory["Relam_cd"] = boost::factory<Relam_cd*>();
    obsThFactory["Relam_ud"] = boost::factory<Relam_ud*>();
    obsThFactory["Relam_ts"] = boost::factory<Relam_ts*>();
    obsThFactory["Relam_cs"] = boost::factory<Relam_cs*>();
    obsThFactory["Relam_us"] = boost::factory<Relam_us*>();
    obsThFactory["Imlam_t"] = boost::factory<Imlam_t*>();
    obsThFactory["Imlam_c"] = boost::factory<Imlam_c*>();
    obsThFactory["Imlam_u"] = boost::factory<Imlam_u*>();
    obsThFactory["Imlam_td"] = boost::factory<Imlam_td*>();
    obsThFactory["Imlam_cd"] = boost::factory<Imlam_cd*>();
    obsThFactory["Imlam_ud"] = boost::factory<Imlam_ud*>();
    obsThFactory["Imlam_ts"] = boost::factory<Imlam_ts*>();
    obsThFactory["Imlam_cs"] = boost::factory<Imlam_cs*>();
    obsThFactory["Imlam_us"] = boost::factory<Imlam_us*>();
    //----- B(s) to mu mu  -----
    obsThFactory["BR_Bdmumu"] = boost::bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_D, StandardModel::MU);
    obsThFactory["BRbar_Bdmumu"] = boost::bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_D, StandardModel::MU);
    obsThFactory["Amumu_Bd"] = boost::bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_D, StandardModel::MU);
    obsThFactory["Smumu_Bd"] = boost::bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_D, StandardModel::MU);

    obsThFactory["BR_Bsmumu"] = boost::bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_S, StandardModel::MU);
    obsThFactory["BRbar_Bsmumu"] = boost::bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_S, StandardModel::MU);
    obsThFactory["Amumu_Bs"] = boost::bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_S, StandardModel::MU);
    obsThFactory["Smumu_Bs"] = boost::bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_S, StandardModel::MU);
    obsThFactory["BR_Bsee"] = boost::bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["BRbar_Bsee"] = boost::bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["Aee_Bs"] = boost::bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["See_Bs"] = boost::bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_S, StandardModel::ELECTRON);

    obsThFactory["BR_BdmumuOBR_Bsmumu"] = boost::factory<BdmumuOBsmumu*>();
   //----- b to q gamma  -----
    obsThFactory["BR_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 1);
    obsThFactory["ACP_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 2);
    obsThFactory["BR_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 1);
    obsThFactory["ACP_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 2);
    obsThFactory["BR_bqgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 1);
    obsThFactory["ACP_bqgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 2);

    //----- Wilson coefficients  -----
    obsThFactory["WC_real_C7g"] = boost::bind(boost::factory<WC_C7g*>(), _1, 0);
    obsThFactory["WC_imag_C7g"] = boost::bind(boost::factory<WC_C7g*>(), _1, 1);
    obsThFactory["WC_abs_C7g"] = boost::bind(boost::factory<WC_C7g*>(), _1, 2);
    obsThFactory["WC_arg_C7g"] = boost::bind(boost::factory<WC_C7g*>(), _1, 3);

    obsThFactory["WC_real_C9_mu"] = boost::bind(boost::factory<WC_C9*>(), _1, 0, StandardModel::MU);
    obsThFactory["WC_imag_C9_mu"] = boost::bind(boost::factory<WC_C9*>(), _1, 1, StandardModel::MU);
    obsThFactory["WC_abs_C9_mu"] = boost::bind(boost::factory<WC_C9*>(), _1, 2, StandardModel::MU);
    obsThFactory["WC_arg_C9_mu"] = boost::bind(boost::factory<WC_C9*>(), _1, 3, StandardModel::MU);

    obsThFactory["WC_real_C9_el"] = boost::bind(boost::factory<WC_C9*>(), _1, 0, StandardModel::ELECTRON);
    obsThFactory["WC_imag_C9_el"] = boost::bind(boost::factory<WC_C9*>(), _1, 1, StandardModel::ELECTRON);
    obsThFactory["WC_abs_C9_el"] = boost::bind(boost::factory<WC_C9*>(), _1, 2, StandardModel::ELECTRON);
    obsThFactory["WC_arg_C9_el"] = boost::bind(boost::factory<WC_C9*>(), _1, 3, StandardModel::ELECTRON);

    obsThFactory["WC_real_C10_mu"] = boost::bind(boost::factory<WC_C10*>(), _1, 0, StandardModel::MU);
    obsThFactory["WC_imag_C10_mu"] = boost::bind(boost::factory<WC_C10*>(), _1, 1, StandardModel::MU);
    obsThFactory["WC_abs_C10_mu"] = boost::bind(boost::factory<WC_C10*>(), _1, 2, StandardModel::MU);
    obsThFactory["WC_arg_C10_mu"] = boost::bind(boost::factory<WC_C10*>(), _1, 3, StandardModel::MU);

    obsThFactory["WC_real_C10_el"] = boost::bind(boost::factory<WC_C10*>(), _1, 0, StandardModel::ELECTRON);
    obsThFactory["WC_imag_C10_el"] = boost::bind(boost::factory<WC_C10*>(), _1, 1, StandardModel::ELECTRON);
    obsThFactory["WC_abs_C10_el"] = boost::bind(boost::factory<WC_C10*>(), _1, 2, StandardModel::ELECTRON);
    obsThFactory["WC_arg_C10_el"] = boost::bind(boost::factory<WC_C10*>(), _1, 3, StandardModel::ELECTRON);

    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_1_BdKste"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_2_BdKstmu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2_BdKste"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_3_BdKstmu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3_BdKste"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_4p_BdKstmu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4p_BdKste"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_5p_BdKstmu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5p_BdKste"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_6p_BdKstmu"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6p_BdKste"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_8p_BdKstmu"] = boost::bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_8p_BdKste"] = boost::bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["Gammap_BdKstmu"] = boost::bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_FB_BdKstmu"] = boost::bind(boost::factory<A_FB*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BR_BdKstmu"] = boost::bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BR_BdKste"] = boost::bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["RKst_BdKstll"] = boost::bind(boost::factory<R_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RKstL_BdKstll"] = boost::bind(boost::factory<RL_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RKstT_BdKstll"] = boost::bind(boost::factory<RT_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["R6_BdKstll"] = boost::bind(boost::factory<R_6*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["ACP_BdKstmu"] = boost::bind(boost::factory<ACP_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P3CP_BdKstmu"] = boost::bind(boost::factory<P3CP*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_L_BdKstmu"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_L_BdKste"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["M_1p_BdKstmu"] = boost::bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["M_2p_BdKstmu"] = boost::bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_3_BdKstmu"] = boost::bind(boost::factory<S_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_4_BdKstmu"] = boost::bind(boost::factory<S_4*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_5_BdKstmu"] = boost::bind(boost::factory<S_5*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_7_BdKstmu"] = boost::bind(boost::factory<S_7*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_8_BdKstmu"] = boost::bind(boost::factory<S_8*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_9_BdKstmu"] = boost::bind(boost::factory<S_9*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_5_BdKstmu"] = boost::bind(boost::factory<A_5*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_6_BdKstmu"] = boost::bind(boost::factory<A_6*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_6c_BdKstmu"] = boost::bind(boost::factory<A_6c*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_8_BdKstmu"] = boost::bind(boost::factory<A_8*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_9_BdKstmu"] = boost::bind(boost::factory<A_9*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["P_1f_BdKstmu"] = boost::bind(boost::factory<P_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2f_BdKstmu"] = boost::bind(boost::factory<P_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3f_BdKstmu"] = boost::bind(boost::factory<P_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4pf_BdKstmu"] = boost::bind(boost::factory<P_4Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5pf_BdKstmu"] = boost::bind(boost::factory<P_5Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6pf_BdKstmu"] = boost::bind(boost::factory<P_6Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_8pf_BdKstmu"] = boost::bind(boost::factory<P_8Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Gammapf_BdKstmu"] = boost::bind(boost::factory<GammaPrimef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BRf_BdKstmu"] = boost::bind(boost::factory<BRf_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_FBf_BdKstmu"] = boost::bind(boost::factory<A_FBf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_Lf_BdKstmu"] = boost::bind(boost::factory<F_Lf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_3f_BdKstmu"] = boost::bind(boost::factory<S_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_4f_BdKstmu"] = boost::bind(boost::factory<S_4f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_5f_BdKstmu"] = boost::bind(boost::factory<S_5f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_7f_BdKstmu"] = boost::bind(boost::factory<S_7f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_8f_BdKstmu"] = boost::bind(boost::factory<S_8f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_9f_BdKstmu"] = boost::bind(boost::factory<S_9f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_relationf"] = boost::bind(boost::factory<P_relationf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_relation_exactf"] = boost::bind(boost::factory<P_relation_exactf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["V0_BdKstmu"] = boost::bind(boost::factory<V0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Vp_BdKstmu"] = boost::bind(boost::factory<Vp*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Vm_BdKstmu"] = boost::bind(boost::factory<Vm*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["T0_BdKstmu"] = boost::bind(boost::factory<T0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Tp_BdKstmu"] = boost::bind(boost::factory<Tp*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Tm_BdKstmu"] = boost::bind(boost::factory<Tm*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_BdKstmu"] = boost::bind(boost::factory<S*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["QCDfC9_1f_BdKstmu"] = boost::bind(boost::factory<QCDfC9_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9_2f_BdKstmu"] = boost::bind(boost::factory<QCDfC9_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9_3f_BdKstmu"] = boost::bind(boost::factory<QCDfC9_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["QCDfC9p_1f_BdKstmu"] = boost::bind(boost::factory<QCDfC9p_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9p_2f_BdKstmu"] = boost::bind(boost::factory<QCDfC9p_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9p_3f_BdKstmu"] = boost::bind(boost::factory<QCDfC9p_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["Regtilde_1_BdKstmu"] = boost::bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Regtilde_2_BdKstmu"] = boost::bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Regtilde_3_BdKstmu"] = boost::bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);

    obsThFactory["Imgtilde_1_BdKstmu"] = boost::bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imgtilde_2_BdKstmu"] = boost::bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imgtilde_3_BdKstmu"] = boost::bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);

    obsThFactory["Absgtilde_1_BdKstmu"] = boost::bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absgtilde_2_BdKstmu"] = boost::bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absgtilde_3_BdKstmu"] = boost::bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);

    obsThFactory["Arggtilde_1_BdKstmu"] = boost::bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Arggtilde_2_BdKstmu"] = boost::bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Arggtilde_3_BdKstmu"] = boost::bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);

    obsThFactory["Reh_0_BdKstmu"] = boost::bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Reh_p_BdKstmu"] = boost::bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Reh_m_BdKstmu"] = boost::bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);

    obsThFactory["Imh_0_BdKstmu"] = boost::bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imh_p_BdKstmu"] = boost::bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imh_m_BdKstmu"] = boost::bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);

    obsThFactory["Absh_0_BdKstmu"] = boost::bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absh_p_BdKstmu"] = boost::bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absh_m_BdKstmu"] = boost::bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);

    obsThFactory["Argh_0_BdKstmu"] = boost::bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Argh_p_BdKstmu"] = boost::bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Argh_m_BdKstmu"] = boost::bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);

    //----- B+ to K*+ ll  -----
    obsThFactory["A_FB_BpKstmu"] = boost::bind(boost::factory<A_FB*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["F_L_BpKstmu"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["BR_BpKstmu"] = boost::bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);

    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["C_BKstgamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["S_BKstgamma"] = boost::bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ADG_BKstgamma"] = boost::bind(boost::factory<ADG_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_1"] = boost::bind(boost::factory<DC7_1*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_2"] = boost::bind(boost::factory<DC7_2*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_L"] = boost::bind(boost::factory<AbsDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_R"] = boost::bind(boost::factory<AbsDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_L_Bd"] = boost::bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_R_Bd"] = boost::bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_L_Bd"] = boost::bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_R_Bd"] = boost::bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["hp0_hm0"] = boost::bind(boost::factory<hp0_hm0*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_QCDF_Bd"] = boost::bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_QCDF_Bd_bar"] = boost::bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_QCDF_Bd"] = boost::bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_QCDF_Bd_bar"] = boost::bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_QCDF_Bd"] = boost::bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_QCDF_Bd_bar"] = boost::bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);


    //----- B+ to K*+ gamma  -----
    obsThFactory["BR_BpKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ACP_BpKstgamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_L_Bp"] = boost::bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_R_Bp"] = boost::bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_L_Bp"] = boost::bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_R_Bp"] = boost::bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["AbsDC7_QCDF_Bp"] = boost::bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["AbsDC7_QCDF_Bp_bar"] = boost::bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_QCDF_Bp"] = boost::bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_QCDF_Bp_bar"] = boost::bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_QCDF_Bp"] = boost::bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_QCDF_Bp_bar"] = boost::bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);

    //----- B to phi ll  -----
    obsThFactory["P_1_Bsphimu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_2_Bsphimu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_3_Bsphimu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_4p_Bsphimu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_5p_Bsphimu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_6p_Bsphimu"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_8p_Bsphimu"] = boost::bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::PHI, StandardModel::MU);
    obsThFactory["Gammap_Bsphimu"] = boost::bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_FB_Bsphimu"] = boost::bind(boost::factory<A_FB*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["BR_Bsphimu"] = boost::bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["Rphi_Bsphill"] = boost::bind(boost::factory<R_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RphiL_Bsphill"] = boost::bind(boost::factory<RL_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RphiT_Bsphill"] = boost::bind(boost::factory<RT_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["R6_Bsphill"] = boost::bind(boost::factory<R_6*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["ACP_Bsphimu"] = boost::bind(boost::factory<ACP_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P3CP_Bsphimu"] = boost::bind(boost::factory<P3CP*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["F_L_Bsphimu"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_1p_Bsphimu"] = boost::bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_2p_Bsphimu"] = boost::bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_3_Bsphimu"] = boost::bind(boost::factory<S_3*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_4_Bsphimu"] = boost::bind(boost::factory<S_4*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_5_Bsphimu"] = boost::bind(boost::factory<S_5*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_7_Bsphimu"] = boost::bind(boost::factory<S_7*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_8_Bsphimu"] = boost::bind(boost::factory<S_8*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_9_Bsphimu"] = boost::bind(boost::factory<S_9*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_5_Bsphimu"] = boost::bind(boost::factory<A_5*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_6_Bsphimu"] = boost::bind(boost::factory<A_6*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_6c_Bsphimu"] = boost::bind(boost::factory<A_6c*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_8_Bsphimu"] = boost::bind(boost::factory<A_8*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_9_Bsphimu"] = boost::bind(boost::factory<A_9*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);

    //----- B to PHI gamma  -----
    obsThFactory["BR_Bsphigamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["C_Bsphigamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["S_Bsphigamma"] = boost::bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ADG_Bsphigamma"] = boost::bind(boost::factory<ADG_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_L_Bs"] = boost::bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_R_Bs"] = boost::bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_L_Bs"] = boost::bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_R_Bs"] = boost::bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["AbsDC7_QCDF_Bs"] = boost::bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["AbsDC7_QCDF_Bs_bar"] = boost::bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_QCDF_Bs"] = boost::bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_QCDF_Bs_bar"] = boost::bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_QCDF_Bs"] = boost::bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_QCDF_Bs_bar"] = boost::bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);

    //----- B+ to K+ ll  -----
    obsThFactory["BR_BpKmu"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["BR_BpKe"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["dBR_BpKmu"] = boost::bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["dBR_BpKe"] = boost::bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["RK_BpKll"] = boost::bind(boost::factory<R_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU, StandardModel::ELECTRON);

    //----- B0 to K0 ll  -----
    obsThFactory["BR_B0Kmu"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);
    obsThFactory["BR_B0Ke"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON);
    obsThFactory["dBR_B0Kmu"] = boost::bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);
    obsThFactory["dBR_B0Ke"] = boost::bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON);
    obsThFactory["RK_B0Kll"] = boost::bind(boost::factory<R_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU, StandardModel::ELECTRON);

    //----- B to D*lnu -----
    obsThFactory["Gammaw_MVlnu"] = boost::bind(boost::factory<Gammaw_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RDstar_MVlnu"] = boost::bind(boost::factory<RDstar_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["Gammacl_MVlnu"] = boost::bind(boost::factory<Gammacl_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["GammacV_MVlnu"] = boost::bind(boost::factory<GammacV_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["Gammachi_MVlnu"] = boost::bind(boost::factory<Gammachi_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_MVlnu"] = boost::bind(boost::factory<UnitarityV_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityA_MVlnu"] = boost::bind(boost::factory<UnitarityA_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_D_Dst"] = boost::bind(boost::factory<UnitarityV_D_Dst*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["hA1_at_w1"] = boost::bind(boost::factory<FF_hA1atw1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hV_w"] = boost::bind(boost::factory<FF_hV*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA1_w"] = boost::bind(boost::factory<FF_hA1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA2_w"] = boost::bind(boost::factory<FF_hA2*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA3_w"] = boost::bind(boost::factory<FF_hA3*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R1_w"] = boost::bind(boost::factory<FF_R1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R2_w"] = boost::bind(boost::factory<FF_R2*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R0_w"] = boost::bind(boost::factory<FF_R0*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["FL_MVtaunu"] = boost::bind(boost::factory<FL_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU);
    obsThFactory["Ptau_MVtaunu"] = boost::bind(boost::factory<Plep_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU);
    //----- B to Dlnu -----
    obsThFactory["Gammaw_MPlnu"] = boost::bind(boost::factory<Gammaw_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RD_MPlnu"] = boost::bind(boost::factory<RD_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_MPlnu"] = boost::bind(boost::factory<UnitarityV_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityA_MPlnu"] = boost::bind(boost::factory<UnitarityA_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["Unitarity_Strong_MPlnu"] = boost::bind(boost::factory<Unitarity_Strong_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["af0_0"] = boost::bind(boost::factory<af0_0*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["FF0_MPlnu"] = boost::bind(boost::factory<FF0_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["FFplus_MPlnu"] = boost::bind(boost::factory<FFplus_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);

    //----- B to tau nu  -----
    obsThFactory["btaunu"] = boost::bind(boost::factory<Btaunu*>(), _1, StandardModel::B_P);
    obsThFactory["bctaunu"] = boost::bind(boost::factory<Btaunu*>(), _1, StandardModel::B_C);

}

void ThObsFactory::addObsToFactory(const std::string name, boost::function<ThObservable*(const StandardModel&) > funct)
{
    if (obsThFactory.find(name) == obsThFactory.end()) obsThFactory[name] = funct;
    else throw std::runtime_error("ERROR: Observable named: " + name + " already exists. Please give a different name to the user-defined observable " + name + ".");
}

ThObservable * ThObsFactory::CreateThMethod(const std::string& name, StandardModel& model) const
{
    if (model.isModelParam(name))
        return new ParamObs(model, name);
    if (obsThFactory.find(name) == obsThFactory.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThObsFactory.\nIf " + name + " is a parameter that is specific to an observable, please list it after the observable in the configuration file.\n");
    ThObservable * myThObs = obsThFactory.at(name)(model);
    if (!myThObs->getParametersForObservable().empty()) model.addParameters(myThObs->getParametersForObservable());
    return (myThObs);
}
