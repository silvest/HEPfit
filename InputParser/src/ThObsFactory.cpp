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
#include "HiggsThObservables.h"
#include "FlavourObservables.h"
#include "MtMSbar.h"
/** BEGIN: REMOVE FROM THE PACKAGE **/
#include "LeptonFlavourObservables.h"
#include "SUSYObservables.h"
#include "GeneralTHDMObservables.h"
#include "THDMObservables.h"
/** END: REMOVE FROM THE PACKAGE **/
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

ThObsFactory::ThObsFactory()
{
    //-----  StandardModel observables  -----
    obsThFactory["MtMSbar"] = boost::factory<MtMSbar*>();
    //-----  Electroweak precision observables  -----
    obsThFactory["Mw"] = boost::factory<Mw*>();
    obsThFactory["GammaW"] = boost::factory<GammaW*>();
    obsThFactory["GammaZ"] = boost::factory<GammaZ*>();
    obsThFactory["sigmaHadron"] = boost::factory<sigmaHadron*>();
    obsThFactory["sin2thetaEff"] = boost::factory<sin2thetaEff*>();
    obsThFactory["PtauPol"] = boost::factory<PtauPol*>();
    obsThFactory["Alepton"] = boost::factory<Alepton*>();
    obsThFactory["Acharm"] = boost::factory<Acharm*>();
    obsThFactory["Abottom"] = boost::factory<Abottom*>();
    obsThFactory["AFBlepton"] = boost::factory<AFBlepton*>();
    obsThFactory["AFBcharm"] = boost::factory<AFBcharm*>();
    obsThFactory["AFBbottom"] = boost::factory<AFBbottom*>();
    obsThFactory["Rlepton"] = boost::factory<Rlepton*>();
    obsThFactory["Rcharm"] = boost::factory<Rcharm*>();
    obsThFactory["Rbottom"] = boost::factory<Rbottom*>();

    //-----  Higgs Extension observables  ----------
    const double sqrt_s_LHC7 = 7.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC8 = 8.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC13 = 13.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC14 = 14.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCC100 = 100.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_TeV = 1.96;
    const double sqrt_s_TLEP = .24;
    const double sqrt_s_ILC250 = .25;
    const double sqrt_s_ILC500 = .5;
    const double sqrt_s_ILC1000 = 1.0;
    obsThFactory["ggH"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC8);
    obsThFactory["WH"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ZH"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH+ttH"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF+VH"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH7"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC7);
    obsThFactory["VBF7"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC7);
    obsThFactory["WH7"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ZH7"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC7);
    obsThFactory["VH7"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ttH7"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ggH8"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH+ttH8"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF8"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF+VH8"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH8"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["WH8"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ZH8"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH8"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH13"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH+ttH13"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF13"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF+VH13"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VH13"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["WH13"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ZH13"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ttH13"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH14"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH+ttH14"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF14"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF+VH14"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VH14"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["WH14"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ZH14"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ttH14"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH100"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH+ttH100"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF100"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF+VH100"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VH100"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["WH100"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ZH100"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ttH100"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH196"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_TeV);
    obsThFactory["VBF196"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_TeV);
    obsThFactory["VH196"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_TeV);
    obsThFactory["ttH196"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_TeV);
    obsThFactory["eeZH240"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_TLEP);
    obsThFactory["eeZH250"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC250);
    obsThFactory["eeZH500"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC500);
    obsThFactory["eeZH1000"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC1000);
    obsThFactory["eeWBF250"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC250);
    obsThFactory["eeWBF500"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC500);
    obsThFactory["eeWBF1000"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC1000);
    obsThFactory["eettH500"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_ILC500);
    obsThFactory["eettH1000"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_ILC1000);
    obsThFactory["BrHggRatio"] = boost::factory<BrHtoggRatio*>();
    obsThFactory["BrHWWRatio"] = boost::factory<BrHtoWWRatio*>();
    obsThFactory["BrHZZRatio"] = boost::factory<BrHtoZZRatio*>();
    obsThFactory["BrHZgaRatio"] = boost::factory<BrHtoZgaRatio*>();
    obsThFactory["BrHgagaRatio"] = boost::factory<BrHtogagaRatio*>();
    obsThFactory["BrHmumuRatio"] = boost::factory<BrHtomumuRatio*>();
    obsThFactory["BrHtautauRatio"] = boost::factory<BrHtotautauRatio*>();
    obsThFactory["BrHccRatio"] = boost::factory<BrHtoccRatio*>();
    obsThFactory["BrHbbRatio"] = boost::factory<BrHtobbRatio*>();

    //-----  Epsilon parameters  -----
    obsThFactory["epsilon1"] = boost::factory<Epsilon1*>();
    obsThFactory["epsilon2"] = boost::factory<Epsilon2*>();
    obsThFactory["epsilon3"] = boost::factory<Epsilon3*>();
    obsThFactory["epsilonb"] = boost::factory<Epsilonb*>();

    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  LEP-II two-fermion processes  -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189.,
        192., 196., 200., 202., 205., 207.};
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192.,
        196., 200., 202., 205., 207.};
    for (int i = 0; i < 12; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
        obsThFactory["sigmaqLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2sigmaHadron*>(), _1, sqrt_s[i]);
        obsThFactory["sigmamuLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2sigmaMu*>(), _1, sqrt_s[i]);
        obsThFactory["sigmatauLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2sigmaTau*>(), _1, sqrt_s[i]);
        obsThFactory["AFBmuLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2AFBmu*>(), _1, sqrt_s[i]);
        obsThFactory["AFBtauLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2AFBtau*>(), _1, sqrt_s[i]);
    }
    for (int i = 0; i < 10; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
        obsThFactory["AFBbottomLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2AFBbottom*>(), _1, sqrt_s_HF[i]);
        obsThFactory["AFBcharmLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2AFBcharm*>(), _1, sqrt_s_HF[i]);
        obsThFactory["RbottomLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2Rbottom*>(), _1, sqrt_s_HF[i]);
        obsThFactory["RcharmLEP2_" + sqrt_s_str] = boost::bind(boost::factory<LEP2Rcharm*>(), _1, sqrt_s_HF[i]);
    }
    /** END: REMOVE FROM THE PACKAGE **/

    //-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["DmBd"] = boost::factory<DmBd*>();
    obsThFactory["DmBs"] = boost::factory<DmBs*>();
    obsThFactory["SJPsiK"] = boost::factory<SJPsiK*>();
    obsThFactory["Betas_JPsiPhi"] = boost::factory<Betas_JPsiPhi*>();
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    obsThFactory["M12D"] = boost::factory<M12D*>();
    obsThFactory["ArgD"] = boost::factory<ArgD*>();
    /** END: REMOVE FROM THE PACKAGE **/
    obsThFactory["EpsilonK"] = boost::factory<EpsilonK*>();
    obsThFactory["DmK"] = boost::factory<DmK*>();
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //----- eps'/eps  -----
    obsThFactory["EpsiloP_o_Epsilon"] = boost::factory<EpsilonP_O_Epsilon*>();
    /** END: REMOVE FROM THE PACKAGE **/
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
    obsThFactory["BR_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 1);
    obsThFactory["BRbar_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 2);
    obsThFactory["Amumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 3);
    obsThFactory["Smumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 4);
    obsThFactory["BR_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 1);
    obsThFactory["BRbar_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 2);
    obsThFactory["Amumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 3);
    obsThFactory["Smumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 4);
    obsThFactory["BR_BdmumuOBR_Bsmumu"] = boost::factory<BdmumuOBsmumu*>();
   //----- b to q gamma  -----
    obsThFactory["BR_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 1);
    obsThFactory["BR_CPodd_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 2);
    obsThFactory["BR_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 1);
    obsThFactory["BR_CPodd_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 2);
    obsThFactory["BR_bqgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 1);
    parameterForObservable["BR_bsgamma"] = make_vector<std::string>() << "Gambino_mukin" << "Gambino_BRsem" << "Gambino_Mbkin" << "Gambino_Mcatmuc" 
                                                                      << "Gambino_mupi2" << "Gambino_rhoD3" << "Gambino_muG2" << "Gambino_rhoLS3";
    parameterForObservable["BR_CPodd_bsgamma"] = parameterForObservable["BR_bsgamma"];
    parameterForObservable["BR_bdgamma"] = parameterForObservable["BR_bsgamma"];
    parameterForObservable["BR_CPodd_bdgamma"] = parameterForObservable["BR_bsgamma"];
    parameterForObservable["BR_bqgamma"] = parameterForObservable["BR_bsgamma"];
    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_1_BdKste"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_2_BdKstmu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2_BdKste"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_3_BdKstmu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3_BdKste"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_4p_BdKstmu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5p_BdKstmu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6p_BdKstmu"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_8p_BdKstmu"] = boost::bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
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
    obsThFactory["A_6_BdKstmu"] = boost::bind(boost::factory<A_6*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
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
      
    obsThFactory["DC7_1"] = boost::bind(boost::factory<DC7_1*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_2"] = boost::bind(boost::factory<DC7_2*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["hp0_hm0"] = boost::bind(boost::factory<hp0_hm0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["hm0_h00"] = boost::bind(boost::factory<hm0_h00*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

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
    obsThFactory["A_6_Bsphimu"] = boost::bind(boost::factory<A_6*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_9_Bsphimu"] = boost::bind(boost::factory<A_9*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    //----- B to K ll  -----
    obsThFactory["BR_BKmu"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["BR_BKe"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["RK_BKll"] = boost::bind(boost::factory<R_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU, StandardModel::ELECTRON);
    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["C_BKstgamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["S_BKstgamma"] = boost::bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    //----- B to PHI gamma  -----
    obsThFactory["BR_Bsphigamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["C_Bsphigamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["S_Bsphigamma"] = boost::bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    
    //----- B to tau nu  -----
    obsThFactory["btaunu"] = boost::factory<Btaunu*>();
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  Lepton Flavour observables  -----
    obsThFactory["mu_e_gamma"] = boost::factory<mu_e_gamma*>();
    obsThFactory["log_meg"] = boost::factory<log_meg*>();
    obsThFactory["tau_mu_gamma"] = boost::factory<tau_mu_gamma*>();
    obsThFactory["log_tmg"] = boost::factory<log_tmg*>();
    obsThFactory["tau_e_gamma"] = boost::factory<tau_e_gamma*>();
    obsThFactory["log_teg"] = boost::factory<log_teg*>();
    obsThFactory["mu_3e"] = boost::factory<mu_3e*>();
    obsThFactory["tau_3mu"] = boost::factory<tau_3mu*>();
    obsThFactory["tau_3e"] = boost::factory<tau_3e*>();
    obsThFactory["gminus2_mu"] = boost::factory<gminus2_mu*>();
    obsThFactory["Robs_mu_e_gamma"] = boost::factory<Robs_mu_e_gamma*>();
    obsThFactory["Robs_tau_mu_gamma"] = boost::factory<Robs_tau_mu_gamma*>();
    obsThFactory["Robs_tau_e_gamma"] = boost::factory<Robs_tau_e_gamma*>();

    obsThFactory["deltaRL_12_u"] = boost::factory<deltaRL_12_u*>();
    obsThFactory["deltaRL_13_u"] = boost::factory<deltaRL_13_u*>();
    obsThFactory["deltaRL_23_u"] = boost::factory<deltaRL_23_u*>();
    obsThFactory["deltaRL_12_e"] = boost::factory<deltaRL_12_e*>();
    obsThFactory["deltaRL_21_e"] = boost::factory<deltaRL_21_e*>();
    obsThFactory["deltaRL_13_e"] = boost::factory<deltaRL_13_e*>();
    obsThFactory["deltaRL_31_e"] = boost::factory<deltaRL_31_e*>();
    obsThFactory["deltaRL_23_e"] = boost::factory<deltaRL_23_e*>();
    obsThFactory["deltaRL_32_e"] = boost::factory<deltaRL_32_e*>();

    obsThFactory["deltaLL1_q"] = boost::factory<deltaLL1_q*>();
    obsThFactory["deltaLL2_q"] = boost::factory<deltaLL2_q*>();
    obsThFactory["deltaLL3_q"] = boost::factory<deltaLL3_q*>();
    obsThFactory["deltaRR1_u"] = boost::factory<deltaRR1_u*>();
    obsThFactory["deltaRR2_u"] = boost::factory<deltaRR2_u*>();
    obsThFactory["deltaRR3_u"] = boost::factory<deltaRR3_u*>();
    obsThFactory["deltaRR1_d"] = boost::factory<deltaRR1_d*>();
    obsThFactory["deltaRR2_d"] = boost::factory<deltaRR2_d*>();
    obsThFactory["deltaRR3_d"] = boost::factory<deltaRR3_d*>();
    obsThFactory["deltaLL1_l"] = boost::factory<deltaLL1_l*>();
    obsThFactory["deltaLL2_l"] = boost::factory<deltaLL2_l*>();
    obsThFactory["deltaLL3_l"] = boost::factory<deltaLL3_l*>();
    obsThFactory["deltaRR1_e"] = boost::factory<deltaRR1_e*>();
    obsThFactory["deltaRR2_e"] = boost::factory<deltaRR2_e*>();
    obsThFactory["deltaRR3_e"] = boost::factory<deltaRR3_e*>();

    obsThFactory["CCBu11"] = boost::factory<CCBu11*>();
    obsThFactory["CCBu22"] = boost::factory<CCBu22*>();
    obsThFactory["CCBu33"] = boost::factory<CCBu33*>();
    obsThFactory["CCBu12"] = boost::factory<CCBu12*>();
    obsThFactory["CCBu13"] = boost::factory<CCBu13*>();
    obsThFactory["CCBu23"] = boost::factory<CCBu23*>();
    obsThFactory["CCBd11"] = boost::factory<CCBd11*>();
    obsThFactory["CCBd22"] = boost::factory<CCBd22*>();
    obsThFactory["CCBd33"] = boost::factory<CCBd33*>();
    obsThFactory["CCBd12"] = boost::factory<CCBd12*>();
    obsThFactory["CCBd13"] = boost::factory<CCBd13*>();
    obsThFactory["CCBd23"] = boost::factory<CCBd23*>();
    obsThFactory["CCBe11"] = boost::factory<CCBe11*>();
    obsThFactory["CCBe22"] = boost::factory<CCBe22*>();
    obsThFactory["CCBe33"] = boost::factory<CCBe33*>();
    obsThFactory["CCBe12"] = boost::factory<CCBe12*>();
    obsThFactory["CCBe13"] = boost::factory<CCBe13*>();
    obsThFactory["CCBe23"] = boost::factory<CCBe23*>();
    
    obsThFactory["logdeltaRL_13_e"] = boost::factory<logdeltaRL_13_e*>();
    obsThFactory["logdeltaRL_23_e"] = boost::factory<logdeltaRL_23_e*>();
    obsThFactory["logmslepton"] = boost::factory<logmslepton*>();
    obsThFactory["deltaTEhat23"] = boost::factory<deltaTEhat23*>();
    obsThFactory["deltaLLRR_l"] = boost::factory<deltaLLRR_l*>();

    /** END: REMOVE FROM THE PACKAGE **/
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  SUSY spectra and observables  -----
    obsThFactory["OutputSLHAfromFH"] = boost::factory<OutputSLHAfromFH*>(); // for debug
    obsThFactory["MHl"] = boost::bind(boost::factory<Mhiggs*>(), _1, 0);
    obsThFactory["MHh"] = boost::bind(boost::factory<Mhiggs*>(), _1, 1);
    obsThFactory["MHa"] = boost::bind(boost::factory<Mhiggs*>(), _1, 2);
    obsThFactory["MHp"] = boost::bind(boost::factory<Mhiggs*>(), _1, 3);
    obsThFactory["Msu1"] = boost::bind(boost::factory<Msup*>(), _1, 0);
    obsThFactory["Msu2"] = boost::bind(boost::factory<Msup*>(), _1, 1);
    obsThFactory["Msu3"] = boost::bind(boost::factory<Msup*>(), _1, 2);
    obsThFactory["Msu4"] = boost::bind(boost::factory<Msup*>(), _1, 3);
    obsThFactory["Msu5"] = boost::bind(boost::factory<Msup*>(), _1, 4);
    obsThFactory["Msu6"] = boost::bind(boost::factory<Msup*>(), _1, 5);
    obsThFactory["Msd1"] = boost::bind(boost::factory<Msdown*>(), _1, 0);
    obsThFactory["Msd2"] = boost::bind(boost::factory<Msdown*>(), _1, 1);
    obsThFactory["Msd3"] = boost::bind(boost::factory<Msdown*>(), _1, 2);
    obsThFactory["Msd4"] = boost::bind(boost::factory<Msdown*>(), _1, 3);
    obsThFactory["Msd5"] = boost::bind(boost::factory<Msdown*>(), _1, 4);
    obsThFactory["Msd6"] = boost::bind(boost::factory<Msdown*>(), _1, 5);
    obsThFactory["Mch1"] = boost::bind(boost::factory<Mchargino*>(), _1, 0);
    obsThFactory["Mch2"] = boost::bind(boost::factory<Mchargino*>(), _1, 1);
    obsThFactory["Mneu1"] = boost::bind(boost::factory<Mneutralino*>(), _1, 0);
    obsThFactory["Mneu2"] = boost::bind(boost::factory<Mneutralino*>(), _1, 1);
    obsThFactory["Mneu3"] = boost::bind(boost::factory<Mneutralino*>(), _1, 2);
    obsThFactory["Mneu4"] = boost::bind(boost::factory<Mneutralino*>(), _1, 3);
    obsThFactory["Mw_dRho"] = boost::factory<Mw_dRho*>();
    /** END: REMOVE FROM THE PACKAGE **/
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  THDM observables  -----
    obsThFactory["globalminimum"] = boost::factory<globalminimum*>();

    obsThFactory["ggF_tth_htobb"] = boost::factory<ggF_tth_htobb*>();
    obsThFactory["ggF_tth_htoWW"] = boost::factory<ggF_tth_htoWW*>();
    obsThFactory["ggF_tth_htotautau"] = boost::factory<ggF_tth_htotautau*>();
    obsThFactory["ggF_tth_htoZZ"] = boost::factory<ggF_tth_htoZZ*>();
    obsThFactory["ggF_tth_htogaga"] = boost::factory<ggF_tth_htogaga*>();
    obsThFactory["VBF_Vh_htobb"] = boost::factory<VBF_Vh_htobb*>();
    obsThFactory["VBF_Vh_htoWW"] = boost::factory<VBF_Vh_htoWW*>();
    obsThFactory["VBF_Vh_htotautau"] = boost::factory<VBF_Vh_htotautau*>();
    obsThFactory["VBF_Vh_htoZZ"] = boost::factory<VBF_Vh_htoZZ*>();
    obsThFactory["VBF_Vh_htogaga"] = boost::factory<VBF_Vh_htogaga*>();
    obsThFactory["Gamma_h_THDM"] = boost::factory<Gamma_h_THDM*>();
    obsThFactory["rh_gaga_THDM"] = boost::factory<rh_gaga_THDM*>();
    obsThFactory["rh_gg_THDM"] = boost::factory<rh_gg_THDM*>();

    obsThFactory["Hobs_ggF_H_tautau_ATLAS"] = boost::factory<Hobs_ggF_H_tautau_ATLAS*>();
    obsThFactory["Hobs_ggF_H_tautau_CMS"] = boost::factory<Hobs_ggF_H_tautau_CMS*>();
    obsThFactory["Hobs_bbF_H_tautau_ATLAS"] = boost::factory<Hobs_bbF_H_tautau_ATLAS*>();
    obsThFactory["Hobs_bbF_H_tautau_CMS"] = boost::factory<Hobs_bbF_H_tautau_CMS*>();
    obsThFactory["Hobs_pp_H_gaga_ATLAS"] = boost::factory<Hobs_pp_H_gaga_ATLAS*>();
    obsThFactory["Hobs_ggF_H_gaga_CMS"] = boost::factory<Hobs_ggF_H_gaga_CMS*>();
    obsThFactory["Hobs_mu_pp_H_VV_CMS"] = boost::factory<Hobs_mu_pp_H_VV_CMS*>();
    obsThFactory["Hobs_ggF_H_ZZ_ATLAS"] = boost::factory<Hobs_ggF_H_ZZ_ATLAS*>();
    obsThFactory["Hobs_VBF_H_ZZ_ATLAS"] = boost::factory<Hobs_VBF_H_ZZ_ATLAS*>();
    obsThFactory["Hobs_ggF_H_WW_ATLAS"] = boost::factory<Hobs_ggF_H_WW_ATLAS*>();
    obsThFactory["Hobs_VBF_H_WW_ATLAS"] = boost::factory<Hobs_VBF_H_WW_ATLAS*>();
    obsThFactory["Hobs_ggF_H_hh_ATLAS"] = boost::factory<Hobs_ggF_H_hh_ATLAS*>();
    obsThFactory["Hobs_ggF_H_hh_bbtautau_CMS"] = boost::factory<Hobs_ggF_H_hh_bbtautau_CMS*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_CMS"] = boost::factory<Hobs_pp_H_hh_bbbb_CMS*>();
    obsThFactory["Hobs_pp_H_hh_gagabb_CMS"] = boost::factory<Hobs_pp_H_hh_gagabb_CMS*>();
    obsThFactory["Hobs_ggF_H_tt_ATLAS"] = boost::factory<Hobs_ggF_H_tt_ATLAS*>();
    obsThFactory["Hobs_bbF_H_bb_CMS"] = boost::factory<Hobs_bbF_H_bb_CMS*>();
    obsThFactory["Robs_ggF_H_tautau_ATLAS"] = boost::factory<Robs_ggF_H_tautau_ATLAS*>();
    obsThFactory["Robs_ggF_H_tautau_CMS"] = boost::factory<Robs_ggF_H_tautau_CMS*>();
    obsThFactory["Robs_bbF_H_tautau_ATLAS"] = boost::factory<Robs_bbF_H_tautau_ATLAS*>();
    obsThFactory["Robs_bbF_H_tautau_CMS"] = boost::factory<Robs_bbF_H_tautau_CMS*>();
    obsThFactory["Robs_pp_H_gaga_ATLAS"] = boost::factory<Robs_pp_H_gaga_ATLAS*>();
    obsThFactory["Robs_ggF_H_gaga_CMS"] = boost::factory<Robs_ggF_H_gaga_CMS*>();
    obsThFactory["Robs_mu_pp_H_VV_CMS"] = boost::factory<Robs_mu_pp_H_VV_CMS*>();
    obsThFactory["Robs_ggF_H_ZZ_ATLAS"] = boost::factory<Robs_ggF_H_ZZ_ATLAS*>();
    obsThFactory["Robs_VBF_H_ZZ_ATLAS"] = boost::factory<Robs_VBF_H_ZZ_ATLAS*>();
    obsThFactory["Robs_ggF_H_WW_ATLAS"] = boost::factory<Robs_ggF_H_WW_ATLAS*>();
    obsThFactory["Robs_VBF_H_WW_ATLAS"] = boost::factory<Robs_VBF_H_WW_ATLAS*>();
    obsThFactory["Robs_ggF_H_hh_ATLAS"] = boost::factory<Robs_ggF_H_hh_ATLAS*>();
    obsThFactory["Robs_ggF_H_hh_bbtautau_CMS"] = boost::factory<Robs_ggF_H_hh_bbtautau_CMS*>();
    obsThFactory["Robs_pp_H_hh_bbbb_CMS"] = boost::factory<Robs_pp_H_hh_bbbb_CMS*>();
    obsThFactory["Robs_pp_H_hh_gagabb_CMS"] = boost::factory<Robs_pp_H_hh_gagabb_CMS*>();
    obsThFactory["Robs_ggF_H_tt_ATLAS"] = boost::factory<Robs_ggF_H_tt_ATLAS*>();
    obsThFactory["Robs_bbF_H_bb_CMS"] = boost::factory<Robs_bbF_H_bb_CMS*>();
    obsThFactory["log10_ggF_H_tautau_TH"] = boost::factory<log10_ggF_H_tautau_TH*>();
    obsThFactory["log10_bbF_H_tautau_TH"] = boost::factory<log10_bbF_H_tautau_TH*>();
    obsThFactory["log10_pp_H_gaga_TH"] = boost::factory<log10_pp_H_gaga_TH*>();
    obsThFactory["log10_ggF_H_gaga_TH"] = boost::factory<log10_ggF_H_gaga_TH*>();
    obsThFactory["log10_mu_pp_H_VV_TH"] = boost::factory<log10_mu_pp_H_VV_TH*>();
    obsThFactory["log10_ggF_H_ZZ_TH"] = boost::factory<log10_ggF_H_ZZ_TH*>();
    obsThFactory["log10_VBF_H_ZZ_TH"] = boost::factory<log10_VBF_H_ZZ_TH*>();
    obsThFactory["log10_ggF_H_WW_TH"] = boost::factory<log10_ggF_H_WW_TH*>();
    obsThFactory["log10_VBF_H_WW_TH"] = boost::factory<log10_VBF_H_WW_TH*>();
    obsThFactory["log10_ggF_H_hh_TH"] = boost::factory<log10_ggF_H_hh_TH*>();
    obsThFactory["log10_ggF_H_hh_bbtautau_TH"] = boost::factory<log10_ggF_H_hh_bbtautau_TH*>();
    obsThFactory["log10_pp_H_hh_bbbb_TH"] = boost::factory<log10_pp_H_hh_bbbb_TH*>();
    obsThFactory["log10_pp_H_hh_gagabb_TH"] = boost::factory<log10_pp_H_hh_gagabb_TH*>();
    obsThFactory["log10_ggF_H_tt_TH"] = boost::factory<log10_ggF_H_tt_TH*>();
    obsThFactory["log10_bbF_H_bb_TH"] = boost::factory<log10_bbF_H_bb_TH*>();
    obsThFactory["Gamma_HH_THDM"] = boost::factory<Gamma_HH_THDM*>();
    obsThFactory["rHH_gg_THDM"] = boost::factory<rHH_gg_THDM*>();
    obsThFactory["BR_HH_hh_THDM"] = boost::factory<BR_HH_hh_THDM*>();
    obsThFactory["BR_HH_AA_THDM"] = boost::factory<BR_HH_AA_THDM*>();
    obsThFactory["BR_HH_HpHm_THDM"] = boost::factory<BR_HH_HpHm_THDM*>();
    obsThFactory["BR_HH_AZ_THDM"] = boost::factory<BR_HH_AZ_THDM*>();
    obsThFactory["BR_HH_HpW_THDM"] = boost::factory<BR_HH_HpW_THDM*>();

    obsThFactory["Hobs_ggF_A_tautau_ATLAS"] = boost::factory<Hobs_ggF_A_tautau_ATLAS*>();
    obsThFactory["Hobs_ggF_A_tautau_CMS"] = boost::factory<Hobs_ggF_A_tautau_CMS*>();
    obsThFactory["Hobs_bbF_A_tautau_ATLAS"] = boost::factory<Hobs_bbF_A_tautau_ATLAS*>();
    obsThFactory["Hobs_bbF_A_tautau_CMS"] = boost::factory<Hobs_bbF_A_tautau_CMS*>();
    obsThFactory["Hobs_pp_A_gaga_ATLAS"] = boost::factory<Hobs_pp_A_gaga_ATLAS*>();
    obsThFactory["Hobs_ggF_A_gaga_CMS"] = boost::factory<Hobs_ggF_A_gaga_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_bbll_CMS"] = boost::factory<Hobs_ggF_A_hZ_bbll_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_bbZ_ATLAS"] = boost::factory<Hobs_ggF_A_hZ_bbZ_ATLAS*>();
    obsThFactory["Hobs_ggF_A_hZ_tautaull_CMS"] = boost::factory<Hobs_ggF_A_hZ_tautaull_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_tautauZ_ATLAS"] = boost::factory<Hobs_ggF_A_hZ_tautauZ_ATLAS*>();
    obsThFactory["Hobs_ggF_A_tt_ATLAS"] = boost::factory<Hobs_ggF_A_tt_ATLAS*>();
    obsThFactory["Hobs_bbF_A_bb_CMS"] = boost::factory<Hobs_bbF_A_bb_CMS*>();
    obsThFactory["Robs_ggF_A_tautau_ATLAS"] = boost::factory<Robs_ggF_A_tautau_ATLAS*>();
    obsThFactory["Robs_ggF_A_tautau_CMS"] = boost::factory<Robs_ggF_A_tautau_CMS*>();
    obsThFactory["Robs_bbF_A_tautau_ATLAS"] = boost::factory<Robs_bbF_A_tautau_ATLAS*>();
    obsThFactory["Robs_bbF_A_tautau_CMS"] = boost::factory<Robs_bbF_A_tautau_CMS*>();
    obsThFactory["Robs_pp_A_gaga_ATLAS"] = boost::factory<Robs_pp_A_gaga_ATLAS*>();
    obsThFactory["Robs_ggF_A_gaga_CMS"] = boost::factory<Robs_ggF_A_gaga_CMS*>();
    obsThFactory["Robs_ggF_A_hZ_bbll_CMS"] = boost::factory<Robs_ggF_A_hZ_bbll_CMS*>();
    obsThFactory["Robs_ggF_A_hZ_bbZ_ATLAS"] = boost::factory<Robs_ggF_A_hZ_bbZ_ATLAS*>();
    obsThFactory["Robs_ggF_A_hZ_tautaull_CMS"] = boost::factory<Robs_ggF_A_hZ_tautaull_CMS*>();
    obsThFactory["Robs_ggF_A_hZ_tautauZ_ATLAS"] = boost::factory<Robs_ggF_A_hZ_tautauZ_ATLAS*>();
    obsThFactory["Robs_ggF_A_tt_ATLAS"] = boost::factory<Robs_ggF_A_tt_ATLAS*>();
    obsThFactory["Robs_bbF_A_bb_CMS"] = boost::factory<Robs_bbF_A_bb_CMS*>();
    obsThFactory["log10_ggF_A_tautau_TH"] = boost::factory<log10_ggF_A_tautau_TH*>();
    obsThFactory["log10_bbF_A_tautau_TH"] = boost::factory<log10_bbF_A_tautau_TH*>();
    obsThFactory["log10_pp_A_gaga_TH"] = boost::factory<log10_pp_A_gaga_TH*>();
    obsThFactory["log10_ggF_A_gaga_TH"] = boost::factory<log10_ggF_A_gaga_TH*>();
    obsThFactory["log10_ggF_A_hZ_bbll_TH"] = boost::factory<log10_ggF_A_hZ_bbll_TH*>();
    obsThFactory["log10_ggF_A_hZ_bbZ_TH"] = boost::factory<log10_ggF_A_hZ_bbZ_TH*>();
    obsThFactory["log10_ggF_A_hZ_tautaull_TH"] = boost::factory<log10_ggF_A_hZ_tautaull_TH*>();
    obsThFactory["log10_ggF_A_hZ_tautauZ_TH"] = boost::factory<log10_ggF_A_hZ_tautauZ_TH*>();
    obsThFactory["log10_ggF_A_tt_TH"] = boost::factory<log10_ggF_A_tt_TH*>();
    obsThFactory["log10_bbF_A_bb_TH"] = boost::factory<log10_bbF_A_bb_TH*>();
    obsThFactory["Gamma_A_THDM"] = boost::factory<Gamma_A_THDM*>();
    obsThFactory["rA_gg_THDM"] = boost::factory<rA_gg_THDM*>();
    obsThFactory["BR_A_HZ_THDM"] = boost::factory<BR_A_HZ_THDM*>();
    obsThFactory["BR_A_hZ_THDM"] = boost::factory<BR_A_hZ_THDM*>();
    obsThFactory["BR_A_HpW_THDM"] = boost::factory<BR_A_HpW_THDM*>();

    obsThFactory["mHh"] = boost::factory<mass_mHh*>();
    obsThFactory["mA"] = boost::factory<mass_mA*>();
    obsThFactory["mHp"] = boost::factory<mass_mHp*>();
    obsThFactory["mHhmmA"] = boost::factory<massdifference_mHhmmA*>();
    obsThFactory["mAmmHh"] = boost::factory<massdifference_mAmmHh*>();
    obsThFactory["mHhmmHp"] = boost::factory<massdifference_mHhmmHp*>();
    obsThFactory["mHpmmHh"] = boost::factory<massdifference_mHpmmHh*>();
    obsThFactory["mAmmHp"] = boost::factory<massdifference_mAmmHp*>();
    obsThFactory["mHpmmA"] = boost::factory<massdifference_mHpmmA*>();
    obsThFactory["m11_2"] = boost::factory<m11_2*>();
    obsThFactory["m22_2"] = boost::factory<m22_2*>();
    obsThFactory["lambda1"] = boost::factory<lambda1*>();
    obsThFactory["lambda2"] = boost::factory<lambda2*>();
    obsThFactory["lambda3"] = boost::factory<lambda3*>();
    obsThFactory["lambda4"] = boost::factory<lambda4*>();
    obsThFactory["lambda5"] = boost::factory<lambda5*>();
    obsThFactory["g_hhh"] = boost::factory<g_hhh*>();
    obsThFactory["g_hhHh"] = boost::factory<g_hhHh*>();
    obsThFactory["g_hHhHh"] = boost::factory<g_hHhHh*>();
    obsThFactory["g_HhHhHh"] = boost::factory<g_HhHhHh*>();
    obsThFactory["g_hAA"] = boost::factory<g_hAA*>();
    obsThFactory["g_HhAA"] = boost::factory<g_HhAA*>();
    obsThFactory["g_hHpHm"] = boost::factory<g_hHpHm*>();
    obsThFactory["g_HhHpHm"] = boost::factory<g_HhHpHm*>();

    obsThFactory["positivity1"] = boost::factory<positivity1*>();
    obsThFactory["positivity2"] = boost::factory<positivity2*>();

    obsThFactory["unitarity1"] = boost::factory<unitarity1*>();
    obsThFactory["unitarity2"] = boost::factory<unitarity2*>();
    obsThFactory["unitarity3"] = boost::factory<unitarity3*>();
    obsThFactory["unitarity4"] = boost::factory<unitarity4*>();
    obsThFactory["unitarity5"] = boost::factory<unitarity5*>();
    obsThFactory["unitarity6"] = boost::factory<unitarity6*>();
    obsThFactory["unitarity7"] = boost::factory<unitarity7*>();
    obsThFactory["unitarity8"] = boost::factory<unitarity8*>();
    obsThFactory["unitarity9"] = boost::factory<unitarity9*>();
    obsThFactory["unitarity10"] = boost::factory<unitarity10*>();
    obsThFactory["unitarity11"] = boost::factory<unitarity11*>();
    obsThFactory["unitarity12"] = boost::factory<unitarity12*>();

    obsThFactory["DeltaS"] = boost::factory<DeltaS*>();
    obsThFactory["DeltaT"] = boost::factory<DeltaT*>();
    obsThFactory["DeltaU"] = boost::factory<DeltaU*>();

    obsThFactory["B_BtoXsgammaTHDM"] = boost::factory<bsgammaTHDM*>();
//    obsThFactory["BR_BsmumuTHDM"] = boost::factory<BR_BsmumuTHDM*>();
//    obsThFactory["BR_BdmumuTHDM"] = boost::factory<BR_BdmumuTHDM*>();

    obsThFactory["Q_st"] = boost::factory<Q_st*>();
    obsThFactory["DeltaQ_THDM"] = boost::factory<DeltaQ_THDM*>();
    obsThFactory["g1atQ"] = boost::factory<g1atQ*>();
    obsThFactory["g2atQ"] = boost::factory<g2atQ*>();
    obsThFactory["g3atQ"] = boost::factory<g3atQ*>();
    obsThFactory["YtopatQ"] = boost::factory<YtopatQ*>();
    obsThFactory["YbottomatQ"] = boost::factory<YbottomatQ*>();
    obsThFactory["YtauatQ"] = boost::factory<YtauatQ*>();
    obsThFactory["m11_2atQ"] = boost::factory<m11_2atQ*>();
    obsThFactory["m22_2atQ"] = boost::factory<m22_2atQ*>();
    obsThFactory["m12_2atQ"] = boost::factory<m12_2atQ*>();
    obsThFactory["lambda1atQ"] = boost::factory<lambda1atQ*>();
    obsThFactory["lambda2atQ"] = boost::factory<lambda2atQ*>();
    obsThFactory["lambda3atQ"] = boost::factory<lambda3atQ*>();
    obsThFactory["lambda4atQ"] = boost::factory<lambda4atQ*>();
    obsThFactory["lambda5atQ"] = boost::factory<lambda5atQ*>();

    obsThFactory["unitarityNLO1"] = boost::factory<unitarityNLO1*>();
    obsThFactory["unitarityNLO2"] = boost::factory<unitarityNLO2*>();
    obsThFactory["unitarityNLO3"] = boost::factory<unitarityNLO3*>();
    obsThFactory["unitarityNLO4"] = boost::factory<unitarityNLO4*>();
    obsThFactory["unitarityNLO5"] = boost::factory<unitarityNLO5*>();
    obsThFactory["unitarityNLO6"] = boost::factory<unitarityNLO6*>();
    obsThFactory["unitarityNLO7"] = boost::factory<unitarityNLO7*>();
    obsThFactory["unitarityNLO8"] = boost::factory<unitarityNLO8*>();
    obsThFactory["unitarityNLO9"] = boost::factory<unitarityNLO9*>();
    obsThFactory["unitarityNLO10"] = boost::factory<unitarityNLO10*>();
    obsThFactory["unitarityNLO11"] = boost::factory<unitarityNLO11*>();
    obsThFactory["unitarityNLO12"] = boost::factory<unitarityNLO12*>();
    obsThFactory["unitarityNLO13"] = boost::factory<unitarityNLO13*>();
    obsThFactory["unitarityNLO14"] = boost::factory<unitarityNLO14*>();
    obsThFactory["unitarityNLO15"] = boost::factory<unitarityNLO15*>();
    obsThFactory["unitarityNLO16"] = boost::factory<unitarityNLO16*>();
    obsThFactory["unitarityNLO17"] = boost::factory<unitarityNLO17*>();
    obsThFactory["unitarityNLO18"] = boost::factory<unitarityNLO18*>();
    obsThFactory["unitarityNLO19"] = boost::factory<unitarityNLO19*>();
    obsThFactory["unitarityNLO20"] = boost::factory<unitarityNLO20*>();
    obsThFactory["unitarityNLO21"] = boost::factory<unitarityNLO21*>();
    obsThFactory["unitarityNLO22"] = boost::factory<unitarityNLO22*>();
    obsThFactory["unitarityNLO23"] = boost::factory<unitarityNLO23*>();
    obsThFactory["unitarityNLO24"] = boost::factory<unitarityNLO24*>();
    obsThFactory["unitarityNLO25"] = boost::factory<unitarityNLO25*>();
    obsThFactory["unitarityNLO26"] = boost::factory<unitarityNLO26*>();
    obsThFactory["unitarityNLOev1"] = boost::factory<unitarityNLOev1*>();
    obsThFactory["unitarityNLOev2"] = boost::factory<unitarityNLOev2*>();
    obsThFactory["unitarityNLOev3"] = boost::factory<unitarityNLOev3*>();
    obsThFactory["unitarityNLOev4"] = boost::factory<unitarityNLOev4*>();
    obsThFactory["unitarityNLOev5"] = boost::factory<unitarityNLOev5*>();
    obsThFactory["unitarityNLOev6"] = boost::factory<unitarityNLOev6*>();
    obsThFactory["unitarityNLOev7"] = boost::factory<unitarityNLOev7*>();
    obsThFactory["unitarityNLOev8"] = boost::factory<unitarityNLOev8*>();
    obsThFactory["unitarityNLOev9"] = boost::factory<unitarityNLOev9*>();
    obsThFactory["unitarityNLOev10"] = boost::factory<unitarityNLOev10*>();
    obsThFactory["unitarityNLOev11"] = boost::factory<unitarityNLOev11*>();
    obsThFactory["unitarityNLOev12"] = boost::factory<unitarityNLOev12*>();
    obsThFactory["unitarityNLOev13"] = boost::factory<unitarityNLOev13*>();
    obsThFactory["unitarityNLOev14"] = boost::factory<unitarityNLOev14*>();
    obsThFactory["unitarityNLOev15"] = boost::factory<unitarityNLOev15*>();
    obsThFactory["unitarityNLOev16"] = boost::factory<unitarityNLOev16*>();
    obsThFactory["unitarityNLOev17"] = boost::factory<unitarityNLOev17*>();
    obsThFactory["unitarityNLOev18"] = boost::factory<unitarityNLOev18*>();
    obsThFactory["unitarityRp1"] = boost::factory<unitarityRp1*>();
    obsThFactory["unitarityRp2"] = boost::factory<unitarityRp2*>();
    obsThFactory["unitarityRp3"] = boost::factory<unitarityRp3*>();
    obsThFactory["unitarityRp4"] = boost::factory<unitarityRp4*>();
    obsThFactory["unitarityRp5"] = boost::factory<unitarityRp5*>();
    obsThFactory["unitarityRp6"] = boost::factory<unitarityRp6*>();
    obsThFactory["unitarityRp7"] = boost::factory<unitarityRp7*>();
    obsThFactory["unitarityRp8"] = boost::factory<unitarityRp8*>();
    obsThFactory["unitarityRp9"] = boost::factory<unitarityRp9*>();
    obsThFactory["unitarityRp10"] = boost::factory<unitarityRp10*>();
    obsThFactory["unitarityRp11"] = boost::factory<unitarityRp11*>();
    obsThFactory["unitarityRp12"] = boost::factory<unitarityRp12*>();
    obsThFactory["unitarityRp13"] = boost::factory<unitarityRp13*>();
    obsThFactory["unitarityRp14"] = boost::factory<unitarityRp14*>();
    obsThFactory["unitarityRp15"] = boost::factory<unitarityRp15*>();
    obsThFactory["unitarityRp16"] = boost::factory<unitarityRp16*>();
    obsThFactory["unitarityRp17"] = boost::factory<unitarityRp17*>();
    obsThFactory["unitarityRp18"] = boost::factory<unitarityRp18*>();
    obsThFactory["unitarityRp19"] = boost::factory<unitarityRp19*>();
    obsThFactory["unitarityRp20"] = boost::factory<unitarityRp20*>();
    obsThFactory["unitarityRp21"] = boost::factory<unitarityRp21*>();
    obsThFactory["unitarityRp22"] = boost::factory<unitarityRp22*>();
    /** END: REMOVE FROM THE PACKAGE **/

    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  GeneralTHDM observables  -----
    //obsThFactory["Re_sigma_u"] = boost::factory<Re_sigma_u*>();
    /** END: REMOVE FROM THE PACKAGE **/
}

void ThObsFactory::addObsToFactory(const std::string name, boost::function<ThObservable*(const StandardModel&) > funct)
{
    obsThFactory[name] = funct;
}

ThObservable * ThObsFactory::CreateThMethod(const std::string& name, StandardModel& model) const
{
    if (model.isModelParam(name))
        return new ParamObs(model, name);
    if (obsThFactory.find(name) == obsThFactory.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThObsFactory");
    if (parameterForObservable.find(name) != parameterForObservable.end()) model.addParameters(parameterForObservable.at(name));
    return (obsThFactory.at(name)(model));
}
