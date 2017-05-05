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
#include "alpha_s.h"
#include "LeptonFlavourObservables.h"
#include "SUSYObservables.h"
/** BEGIN: REMOVE FROM THE PACKAGE **/
#include "GeneralTHDMObservables.h"
#include "LRSMObservables.h"
#include "THDMWObservables.h"
/** END: REMOVE FROM THE PACKAGE **/
#include "THDMObservables.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

ThObsFactory::ThObsFactory()
{
    //-----  StandardModel observables  -----
    obsThFactory["MtMSbar"] = boost::factory<MtMSbar*>();
    obsThFactory["alpha_s_LO"] = boost::bind(boost::factory<alpha_s*>(), _1, LO);
    obsThFactory["alpha_s_NLO"] = boost::bind(boost::factory<alpha_s*>(), _1, NLO);
    obsThFactory["alpha_s_FULLNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, FULLNLO);
    obsThFactory["alpha_s_NNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, NNLO);
    obsThFactory["alpha_s_FULLNNLO"] = boost::bind(boost::factory<alpha_s*>(), _1, FULLNNLO);
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
    const double sqrt_s_FCCee240 = .24;
    const double sqrt_s_FCCee350 = .35;
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
    obsThFactory["eeZH240"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_FCCee240);
    obsThFactory["eeZH250"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC250);
    obsThFactory["eeZH500"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC500);
    obsThFactory["eeZH1000"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_ILC1000);
    obsThFactory["eeWBF250"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC250);
    obsThFactory["eeWBF350"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_FCCee350);
    obsThFactory["eeWBF500"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC500);
    obsThFactory["eeWBF1000"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_ILC1000);
    obsThFactory["eettH500"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_ILC500);
    obsThFactory["eettH1000"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_ILC1000);
    obsThFactory["GammaHRatio"] = boost::factory<GammaHRatio*>();
    obsThFactory["BrHggRatio"] = boost::factory<BrHtoggRatio*>();
    obsThFactory["BrHWWRatio"] = boost::factory<BrHtoWWRatio*>();
    obsThFactory["BrHZZRatio"] = boost::factory<BrHtoZZRatio*>();
    obsThFactory["BrHZgaRatio"] = boost::factory<BrHtoZgaRatio*>();
    obsThFactory["BrHgagaRatio"] = boost::factory<BrHtogagaRatio*>();
    obsThFactory["BrHmumuRatio"] = boost::factory<BrHtomumuRatio*>();
    obsThFactory["BrHtautauRatio"] = boost::factory<BrHtotautauRatio*>();
    obsThFactory["BrHccRatio"] = boost::factory<BrHtoccRatio*>();
    obsThFactory["BrHbbRatio"] = boost::factory<BrHtobbRatio*>();
    obsThFactory["UpperLimit_ppHZgammaA"] = boost::bind(boost::factory<UpperLimit_ppHZgammaA*>(), _1, sqrt_s_LHC8);
    obsThFactory["UpperLimit_ppHZgammaC"] = boost::bind(boost::factory<UpperLimit_ppHZgammaC*>(), _1, sqrt_s_LHC8);

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
    obsThFactory["ACP_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 2);
    obsThFactory["BR_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 1);
    obsThFactory["ACP_bdgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 2);
    obsThFactory["BR_bqgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 1);
    obsThFactory["ACP_bqgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 2);

    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_1_BdKste"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_2_BdKstmu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2_BdKste"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_3_BdKstmu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3_BdKste"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_4p_BdKstmu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5p_BdKstmu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4p_BdKste"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_5p_BdKste"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
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
    
/** BEGIN: REMOVE FROM THE PACKAGE **/
    //----- B to X_q ll -----
    obsThFactory["R_BXsee"] = boost::bind(boost::factory<R_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["Rlow_BXsee"] = boost::bind(boost::factory<Rlow_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["Rhigh_BXsee"] = boost::bind(boost::factory<Rhigh_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
/** END: REMOVE FROM THE PACKAGE **/
    
    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["C_BKstgamma"] = boost::bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["S_BKstgamma"] = boost::bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ADG_BKstgamma"] = boost::bind(boost::factory<ADG_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_1"] = boost::bind(boost::factory<DC7_1*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_2"] = boost::bind(boost::factory<DC7_2*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_L"] = boost::bind(boost::factory<AbsDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_R"] = boost::bind(boost::factory<AbsDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_L"] = boost::bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_R"] = boost::bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_L"] = boost::bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_R"] = boost::bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["hp0_hm0"] = boost::bind(boost::factory<hp0_hm0*>(), _1, StandardModel::B_D, StandardModel::K_star);
    
    //----- B+ to K*+ gamma  -----
    obsThFactory["BR_BpKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star);/*** NOTE: K_star+ is NOT distinguished from K_star ***/
    
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

    //----- B to K ll  -----
    obsThFactory["BR_BKmu"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["BR_BKe"] = boost::bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["RK_BKll"] = boost::bind(boost::factory<R_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU, StandardModel::ELECTRON);
    
    //----- B to tau nu  -----
    obsThFactory["btaunu"] = boost::factory<Btaunu*>();
    
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
    obsThFactory["Robs_tau_mu_gamma_BelleII"] = boost::factory<Robs_tau_mu_gamma_BelleII*>();
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
    
    obsThFactory["VacuumTunnelingRate"] = boost::factory<FindAction*>();
    
    obsThFactory["logdeltaRL_13_e"] = boost::factory<logdeltaRL_13_e*>();
    obsThFactory["logdeltaRL_23_e"] = boost::factory<logdeltaRL_23_e*>();
    obsThFactory["logmslepton"] = boost::factory<logmslepton*>();
    obsThFactory["mslepton"] = boost::factory<mslepton*>();
    obsThFactory["deltaTEhat23"] = boost::factory<deltaTEhat23*>();
    obsThFactory["deltaLLRR_l"] = boost::factory<deltaLLRR_l*>();
    
    //-----  SUSY spectra and observables  -----
/** BEGIN: REMOVE FROM THE PACKAGE **/
    obsThFactory["OutputSLHAfromFH"] = boost::factory<OutputSLHAfromFH*>(); // for debug
/** END: REMOVE FROM THE PACKAGE **/
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
    obsThFactory["Msl1"] = boost::bind(boost::factory<Mslepton*>(), _1, 0);
    obsThFactory["Msl2"] = boost::bind(boost::factory<Mslepton*>(), _1, 1);
    obsThFactory["Msl3"] = boost::bind(boost::factory<Mslepton*>(), _1, 2);
    obsThFactory["Msl4"] = boost::bind(boost::factory<Mslepton*>(), _1, 3);
    obsThFactory["Msl5"] = boost::bind(boost::factory<Mslepton*>(), _1, 4);
    obsThFactory["Msl6"] = boost::bind(boost::factory<Mslepton*>(), _1, 5);
    obsThFactory["Msnu1"] = boost::bind(boost::factory<Msneutrino*>(), _1, 0);
    obsThFactory["Msnu2"] = boost::bind(boost::factory<Msneutrino*>(), _1, 1);
    obsThFactory["Msnu3"] = boost::bind(boost::factory<Msneutrino*>(), _1, 2);
    obsThFactory["Mch1"] = boost::bind(boost::factory<Mchargino*>(), _1, 0);
    obsThFactory["Mch2"] = boost::bind(boost::factory<Mchargino*>(), _1, 1);
    obsThFactory["Mneu1"] = boost::bind(boost::factory<Mneutralino*>(), _1, 0);
    obsThFactory["Mneu2"] = boost::bind(boost::factory<Mneutralino*>(), _1, 1);
    obsThFactory["Mneu3"] = boost::bind(boost::factory<Mneutralino*>(), _1, 2);
    obsThFactory["Mneu4"] = boost::bind(boost::factory<Mneutralino*>(), _1, 3);
    obsThFactory["Mw_dRho"] = boost::factory<Mw_dRho*>();
    
    //-----  THDM observables  -----
    obsThFactory["globalminimum"] = boost::factory<globalminimum*>();

    obsThFactory["mu_ggF_tth_htobb8"] = boost::factory<ggF_tth_htobb8*>();
    obsThFactory["mu_ggF_tth_htoWW8"] = boost::factory<ggF_tth_htoWW8*>();
    obsThFactory["mu_ggF_tth_htotautau8"] = boost::factory<ggF_tth_htotautau8*>();
    obsThFactory["mu_ggF_tth_htoZZ8"] = boost::factory<ggF_tth_htoZZ8*>();
    obsThFactory["mu_ggF_tth_htogaga8"] = boost::factory<ggF_tth_htogaga8*>();
    obsThFactory["mu_ggF_tth_htobb13"] = boost::factory<ggF_tth_htobb13*>();
    obsThFactory["mu_ggF_tth_htoWW13"] = boost::factory<ggF_tth_htoWW13*>();
    obsThFactory["mu_ggF_tth_htotautau13"] = boost::factory<ggF_tth_htotautau13*>();
    obsThFactory["mu_ggF_tth_htoZZ13"] = boost::factory<ggF_tth_htoZZ13*>();
    obsThFactory["mu_ggF_tth_htogaga13"] = boost::factory<ggF_tth_htogaga13*>();
    obsThFactory["mu_VBF_Vh_htobb"] = boost::factory<VBF_Vh_htobb*>();
    obsThFactory["mu_VBF_Vh_htoWW"] = boost::factory<VBF_Vh_htoWW*>();
    obsThFactory["mu_VBF_Vh_htotautau"] = boost::factory<VBF_Vh_htotautau*>();
    obsThFactory["mu_VBF_Vh_htoZZ"] = boost::factory<VBF_Vh_htoZZ*>();
    obsThFactory["mu_VBF_Vh_htogaga"] = boost::factory<VBF_Vh_htogaga*>();
    obsThFactory["mu_ggF_htobb"] = boost::factory<ggF_htobb*>();
    obsThFactory["mu_ggF_htoWW"] = boost::factory<ggF_htoWW*>();
    obsThFactory["mu_ggF_htotautau"] = boost::factory<ggF_htotautau*>();
    obsThFactory["mu_ggF_htoZZ"] = boost::factory<ggF_htoZZ*>();
    obsThFactory["mu_ggF_htogaga"] = boost::factory<ggF_htogaga*>();
    obsThFactory["mu_tth_htobb"] = boost::factory<tth_htobb*>();
    obsThFactory["mu_tth_htoWW"] = boost::factory<tth_htoWW*>();
    obsThFactory["mu_tth_htotautau"] = boost::factory<tth_htotautau*>();
    obsThFactory["mu_tth_htoZZ"] = boost::factory<tth_htoZZ*>();
    obsThFactory["mu_tth_htogaga"] = boost::factory<tth_htogaga*>();
    obsThFactory["mu_htoWW"] = boost::factory<mu_htoWW*>();
    obsThFactory["mu_htotautau"] = boost::factory<mu_htotautau*>();
    obsThFactory["Gamma_h_THDM"] = boost::factory<Gamma_h_THDM*>();
    obsThFactory["rh_gaga_THDM"] = boost::factory<rh_gaga_THDM*>();
    obsThFactory["rh_gg_THDM"] = boost::factory<rh_gg_THDM*>();

    obsThFactory["Hobs_ggF_H_tautau_ATLAS8"] = boost::factory<Hobs_ggF_H_tautau_ATLAS8*>();
    obsThFactory["Hobs_ggF_H_tautau_CMS8"] = boost::factory<Hobs_ggF_H_tautau_CMS8*>();
    obsThFactory["Hobs_bbF_H_tautau_ATLAS8"] = boost::factory<Hobs_bbF_H_tautau_ATLAS8*>();
    obsThFactory["Hobs_bbF_H_tautau_CMS8"] = boost::factory<Hobs_bbF_H_tautau_CMS8*>();
    obsThFactory["Hobs_pp_H_gaga_ATLAS8"] = boost::factory<Hobs_pp_H_gaga_ATLAS8*>();
    obsThFactory["Hobs_ggF_H_gaga_CMS8"] = boost::factory<Hobs_ggF_H_gaga_CMS8*>();
    obsThFactory["Hobs_pp_H_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_H_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_H_Zga_llga_CMS8"] = boost::factory<Hobs_pp_H_Zga_llga_CMS8*>();
    obsThFactory["Hobs_mu_pp_H_VV_CMS8"] = boost::factory<Hobs_mu_pp_H_VV_CMS8*>();
    obsThFactory["Hobs_ggF_H_ZZ_ATLAS8"] = boost::factory<Hobs_ggF_H_ZZ_ATLAS8*>();
    obsThFactory["Hobs_VBF_H_ZZ_ATLAS8"] = boost::factory<Hobs_VBF_H_ZZ_ATLAS8*>();
    obsThFactory["Hobs_ggF_H_WW_ATLAS8"] = boost::factory<Hobs_ggF_H_WW_ATLAS8*>();
    obsThFactory["Hobs_VBF_H_WW_ATLAS8"] = boost::factory<Hobs_VBF_H_WW_ATLAS8*>();
    obsThFactory["Hobs_ggF_H_hh_ATLAS8"] = boost::factory<Hobs_ggF_H_hh_ATLAS8*>();
    obsThFactory["Hobs_pp_H_hh_CMS8"] = boost::factory<Hobs_pp_H_hh_CMS8*>();
    obsThFactory["Hobs_ggF_H_hh_bbtautau_CMS8"] = boost::factory<Hobs_ggF_H_hh_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_CMS8"] = boost::factory<Hobs_pp_H_hh_bbbb_CMS8*>();
    obsThFactory["Hobs_pp_H_hh_gagabb_CMS8"] = boost::factory<Hobs_pp_H_hh_gagabb_CMS8*>();
    obsThFactory["Hobs_ggF_H_tt_ATLAS8"] = boost::factory<Hobs_ggF_H_tt_ATLAS8*>();
    obsThFactory["Hobs_bbF_H_bb_CMS8"] = boost::factory<Hobs_bbF_H_bb_CMS8*>();
    obsThFactory["Robs_ggF_H_tautau_ATLAS8"] = boost::factory<Robs_ggF_H_tautau_ATLAS8*>();
    obsThFactory["Robs_ggF_H_tautau_CMS8"] = boost::factory<Robs_ggF_H_tautau_CMS8*>();
    obsThFactory["Robs_bbF_H_tautau_ATLAS8"] = boost::factory<Robs_bbF_H_tautau_ATLAS8*>();
    obsThFactory["Robs_bbF_H_tautau_CMS8"] = boost::factory<Robs_bbF_H_tautau_CMS8*>();
    obsThFactory["Robs_pp_H_gaga_ATLAS8"] = boost::factory<Robs_pp_H_gaga_ATLAS8*>();
    obsThFactory["Robs_ggF_H_gaga_CMS8"] = boost::factory<Robs_ggF_H_gaga_CMS8*>();
    obsThFactory["Robs_pp_H_Zga_llga_ATLAS8"] = boost::factory<Robs_pp_H_Zga_llga_ATLAS8*>();
    obsThFactory["Robs_pp_H_Zga_llga_CMS8"] = boost::factory<Robs_pp_H_Zga_llga_CMS8*>();
    obsThFactory["Robs_mu_pp_H_VV_CMS8"] = boost::factory<Robs_mu_pp_H_VV_CMS8*>();
    obsThFactory["Robs_ggF_H_ZZ_ATLAS8"] = boost::factory<Robs_ggF_H_ZZ_ATLAS8*>();
    obsThFactory["Robs_VBF_H_ZZ_ATLAS8"] = boost::factory<Robs_VBF_H_ZZ_ATLAS8*>();
    obsThFactory["Robs_ggF_H_WW_ATLAS8"] = boost::factory<Robs_ggF_H_WW_ATLAS8*>();
    obsThFactory["Robs_VBF_H_WW_ATLAS8"] = boost::factory<Robs_VBF_H_WW_ATLAS8*>();
    obsThFactory["Robs_ggF_H_hh_ATLAS8"] = boost::factory<Robs_ggF_H_hh_ATLAS8*>();
    obsThFactory["Robs_pp_H_hh_CMS8"] = boost::factory<Robs_pp_H_hh_CMS8*>();
    obsThFactory["Robs_ggF_H_hh_bbtautau_CMS8"] = boost::factory<Robs_ggF_H_hh_bbtautau_CMS8*>();
    obsThFactory["Robs_pp_H_hh_bbbb_CMS8"] = boost::factory<Robs_pp_H_hh_bbbb_CMS8*>();
    obsThFactory["Robs_pp_H_hh_gagabb_CMS8"] = boost::factory<Robs_pp_H_hh_gagabb_CMS8*>();
    obsThFactory["Robs_ggF_H_tt_ATLAS8"] = boost::factory<Robs_ggF_H_tt_ATLAS8*>();
    obsThFactory["Robs_bbF_H_bb_CMS8"] = boost::factory<Robs_bbF_H_bb_CMS8*>();
    obsThFactory["log10_ggF_H_tautau_TH8"] = boost::factory<log10_ggF_H_tautau_TH8*>();
    obsThFactory["log10_bbF_H_tautau_TH8"] = boost::factory<log10_bbF_H_tautau_TH8*>();
    obsThFactory["log10_pp_H_gaga_TH8"] = boost::factory<log10_pp_H_gaga_TH8*>();
    obsThFactory["log10_ggF_H_gaga_TH8"] = boost::factory<log10_ggF_H_gaga_TH8*>();
    obsThFactory["log10_pp_H_Zga_llga_TH8"] = boost::factory<log10_pp_H_Zga_llga_TH8*>();
    obsThFactory["log10_mu_pp_H_VV_TH8"] = boost::factory<log10_mu_pp_H_VV_TH8*>();
    obsThFactory["log10_ggF_H_ZZ_TH8"] = boost::factory<log10_ggF_H_ZZ_TH8*>();
    obsThFactory["log10_VBF_H_ZZ_TH8"] = boost::factory<log10_VBF_H_ZZ_TH8*>();
    obsThFactory["log10_ggF_H_WW_TH8"] = boost::factory<log10_ggF_H_WW_TH8*>();
    obsThFactory["log10_VBF_H_WW_TH8"] = boost::factory<log10_VBF_H_WW_TH8*>();
    obsThFactory["log10_ggF_H_hh_TH8"] = boost::factory<log10_ggF_H_hh_TH8*>();
    obsThFactory["log10_pp_H_hh_TH8"] = boost::factory<log10_pp_H_hh_TH8*>();
    obsThFactory["log10_ggF_H_hh_bbtautau_TH8"] = boost::factory<log10_ggF_H_hh_bbtautau_TH8*>();
    obsThFactory["log10_pp_H_hh_bbbb_TH8"] = boost::factory<log10_pp_H_hh_bbbb_TH8*>();
    obsThFactory["log10_pp_H_hh_gagabb_TH8"] = boost::factory<log10_pp_H_hh_gagabb_TH8*>();
    obsThFactory["log10_ggF_H_tt_TH8"] = boost::factory<log10_ggF_H_tt_TH8*>();
    obsThFactory["log10_bbF_H_bb_TH8"] = boost::factory<log10_bbF_H_bb_TH8*>();
    obsThFactory["Hobs_ttF_H_tt_ATLAS13"] = boost::factory<Hobs_ttF_H_tt_ATLAS13*>();
    obsThFactory["Hobs_bbF_H_tt_ATLAS13"] = boost::factory<Hobs_bbF_H_tt_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_tautau_ATLAS13"] = boost::factory<Hobs_ggF_H_tautau_ATLAS13*>();
    obsThFactory["Hobs_bbF_H_tautau_ATLAS13"] = boost::factory<Hobs_bbF_H_tautau_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_tautau_CMS13"] = boost::factory<Hobs_ggF_H_tautau_CMS13*>();
    obsThFactory["Hobs_bbF_H_tautau_CMS13"] = boost::factory<Hobs_bbF_H_tautau_CMS13*>();
    obsThFactory["Hobs_pp_H_gaga_ATLAS13"] = boost::factory<Hobs_pp_H_gaga_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_gaga_CMS13"] = boost::factory<Hobs_ggF_H_gaga_CMS13*>();
    obsThFactory["Hobs_pp_H_Zga_llga_ATLAS13"] = boost::factory<Hobs_pp_H_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_pp_H_Zga_llga_CMS13"] = boost::factory<Hobs_pp_H_Zga_llga_CMS13*>();
    obsThFactory["Hobs_pp_H_Zga_qqga_CMS13"] = boost::factory<Hobs_pp_H_Zga_qqga_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llnunu_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llnunu_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llll_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llll_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Hobs_pp_H_ZZ_llll_CMS13"] = boost::factory<Hobs_pp_H_ZZ_llll_CMS13*>();
    obsThFactory["Hobs_VBF_VH_H_ZZ_llll_CMS13"] = boost::factory<Hobs_VBF_VH_H_ZZ_llll_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llqq_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llqq_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_nunuqq_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_nunuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_H_ZZ_llqq_CMS13"] = boost::factory<Hobs_pp_H_ZZ_llqq_CMS13*>();
    obsThFactory["Hobs_ggF_H_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_ggF_H_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_WW_enumunu_ATLAS13"] = boost::factory<Hobs_ggF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VBF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_ggF_VBF_H_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggF_VBF_H_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bbgaga_ATLAS13"] = boost::factory<Hobs_pp_H_hh_bbgaga_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbgaga_CMS13"] = boost::factory<Hobs_pp_H_hh_bbgaga_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_ATLAS13"] = boost::factory<Hobs_pp_H_hh_bbbb_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_CMS13"] = boost::factory<Hobs_pp_H_hh_bbbb_CMS13*>();
    obsThFactory["Hobs_ggF_H_hh_gagaWW_ATLAS13"] = boost::factory<Hobs_ggF_H_hh_gagaWW_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbtautau_CMS13"] = boost::factory<Hobs_pp_H_hh_bbtautau_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bbtautau1_CMS13"] = boost::factory<Hobs_pp_H_hh_bbtautau1_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bblnulnu_CMS13"] = boost::factory<Hobs_pp_H_hh_bblnulnu_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bbVV_CMS13"] = boost::factory<Hobs_pp_H_hh_bbVV_CMS13*>();
    obsThFactory["Hobs_pp_H_bb_CMS13"] = boost::factory<Hobs_pp_H_bb_CMS13*>();
    obsThFactory["Robs_ttF_H_tt_ATLAS13"] = boost::factory<Robs_ttF_H_tt_ATLAS13*>();
    obsThFactory["Robs_bbF_H_tt_ATLAS13"] = boost::factory<Robs_bbF_H_tt_ATLAS13*>();
    obsThFactory["Robs_ggF_H_tautau_ATLAS13"] = boost::factory<Robs_ggF_H_tautau_ATLAS13*>();
    obsThFactory["Robs_bbF_H_tautau_ATLAS13"] = boost::factory<Robs_bbF_H_tautau_ATLAS13*>();
    obsThFactory["Robs_ggF_H_tautau_CMS13"] = boost::factory<Robs_ggF_H_tautau_CMS13*>();
    obsThFactory["Robs_bbF_H_tautau_CMS13"] = boost::factory<Robs_bbF_H_tautau_CMS13*>();
    obsThFactory["Robs_pp_H_gaga_ATLAS13"] = boost::factory<Robs_pp_H_gaga_ATLAS13*>();
    obsThFactory["Robs_ggF_H_gaga_CMS13"] = boost::factory<Robs_ggF_H_gaga_CMS13*>();
    obsThFactory["Robs_pp_H_Zga_llga_ATLAS13"] = boost::factory<Robs_pp_H_Zga_llga_ATLAS13*>();
    obsThFactory["Robs_pp_H_Zga_llga_CMS13"] = boost::factory<Robs_pp_H_Zga_llga_CMS13*>();
    obsThFactory["Robs_pp_H_Zga_qqga_CMS13"] = boost::factory<Robs_pp_H_Zga_qqga_CMS13*>();
    obsThFactory["Robs_ggF_H_ZZ_llnunu_ATLAS13"] = boost::factory<Robs_ggF_H_ZZ_llnunu_ATLAS13*>();
    obsThFactory["Robs_ggF_H_ZZ_llll_ATLAS13"] = boost::factory<Robs_ggF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Robs_VBF_H_ZZ_llll_ATLAS13"] = boost::factory<Robs_VBF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Robs_pp_H_ZZ_llll_CMS13"] = boost::factory<Robs_pp_H_ZZ_llll_CMS13*>();
    obsThFactory["Robs_VBF_VH_H_ZZ_llll_CMS13"] = boost::factory<Robs_VBF_VH_H_ZZ_llll_CMS13*>();
    obsThFactory["Robs_ggF_H_ZZ_llqq_ATLAS13"] = boost::factory<Robs_ggF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Robs_VBF_H_ZZ_llqq_ATLAS13"] = boost::factory<Robs_VBF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Robs_ggF_H_ZZ_nunuqq_ATLAS13"] = boost::factory<Robs_ggF_H_ZZ_nunuqq_ATLAS13*>();
    obsThFactory["Robs_pp_H_ZZ_llqq_CMS13"] = boost::factory<Robs_pp_H_ZZ_llqq_CMS13*>();
    obsThFactory["Robs_ggF_H_WW_lnuqq_ATLAS13"] = boost::factory<Robs_ggF_H_WW_lnuqq_ATLAS13*>();
    obsThFactory["Robs_ggF_H_WW_enumunu_ATLAS13"] = boost::factory<Robs_ggF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Robs_VBF_H_WW_enumunu_ATLAS13"] = boost::factory<Robs_VBF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Robs_ggF_VBF_H_WW_lnulnu_CMS13"] = boost::factory<Robs_ggF_VBF_H_WW_lnulnu_CMS13*>();
    obsThFactory["Robs_pp_H_hh_bbgaga_ATLAS13"] = boost::factory<Robs_pp_H_hh_bbgaga_ATLAS13*>();
    obsThFactory["Robs_pp_H_hh_bbgaga_CMS13"] = boost::factory<Robs_pp_H_hh_bbgaga_CMS13*>();
    obsThFactory["Robs_pp_H_hh_bbbb_ATLAS13"] = boost::factory<Robs_pp_H_hh_bbbb_ATLAS13*>();
    obsThFactory["Robs_pp_H_hh_bbbb_CMS13"] = boost::factory<Robs_pp_H_hh_bbbb_CMS13*>();
    obsThFactory["Robs_ggF_H_hh_gagaWW_ATLAS13"] = boost::factory<Robs_ggF_H_hh_gagaWW_ATLAS13*>();
    obsThFactory["Robs_pp_H_hh_bbtautau_CMS13"] = boost::factory<Robs_pp_H_hh_bbtautau_CMS13*>();
    obsThFactory["Robs_pp_H_hh_bbtautau1_CMS13"] = boost::factory<Robs_pp_H_hh_bbtautau1_CMS13*>();
    obsThFactory["Robs_pp_H_hh_bblnulnu_CMS13"] = boost::factory<Robs_pp_H_hh_bblnulnu_CMS13*>();
    obsThFactory["Robs_pp_H_hh_bbVV_CMS13"] = boost::factory<Robs_pp_H_hh_bbVV_CMS13*>();
    obsThFactory["Robs_pp_H_bb_CMS13"] = boost::factory<Robs_pp_H_bb_CMS13*>();
    obsThFactory["log10_ggF_H_tautau_TH13"] = boost::factory<log10_ggF_H_tautau_TH13*>();
    obsThFactory["log10_bbF_H_tautau_TH13"] = boost::factory<log10_bbF_H_tautau_TH13*>();
    obsThFactory["log10_pp_H_gaga_TH13"] = boost::factory<log10_pp_H_gaga_TH13*>();
    obsThFactory["log10_ggF_H_gaga_TH13"] = boost::factory<log10_ggF_H_gaga_TH13*>();
    obsThFactory["log10_pp_H_Zga_TH13"] = boost::factory<log10_pp_H_Zga_TH13*>();
    obsThFactory["log10_pp_H_ZZ_TH13"] = boost::factory<log10_pp_H_ZZ_TH13*>();
    obsThFactory["log10_ggF_H_ZZ_TH13"] = boost::factory<log10_ggF_H_ZZ_TH13*>();
    obsThFactory["log10_VBF_H_ZZ_TH13"] = boost::factory<log10_VBF_H_ZZ_TH13*>();
    obsThFactory["log10_ggF_H_ZZ_llll_TH13"] = boost::factory<log10_ggF_H_ZZ_llll_TH13*>();
    obsThFactory["log10_VBF_H_ZZ_llll_TH13"] = boost::factory<log10_VBF_H_ZZ_llll_TH13*>();
    obsThFactory["log10_pp_H_ZZ_llll_TH13"] = boost::factory<log10_pp_H_ZZ_llll_TH13*>();
    obsThFactory["log10_VBF_VH_H_ZZ_llll_TH13"] = boost::factory<log10_VBF_VH_H_ZZ_llll_TH13*>();
    obsThFactory["log10_ggF_H_WW_TH13"] = boost::factory<log10_ggF_H_WW_TH13*>();
    obsThFactory["log10_VBF_H_WW_TH13"] = boost::factory<log10_VBF_H_WW_TH13*>();
    obsThFactory["log10_ggF_VBF_H_WW_lnulnu_TH13"] = boost::factory<log10_ggF_VBF_H_WW_lnulnu_TH13*>();
    obsThFactory["log10_ggF_H_hh_TH13"] = boost::factory<log10_ggF_H_hh_TH13*>();
    obsThFactory["log10_pp_H_hh_TH13"] = boost::factory<log10_pp_H_hh_TH13*>();
    obsThFactory["log10_pp_H_hh_bbbb_TH13"] = boost::factory<log10_pp_H_hh_bbbb_TH13*>();
    obsThFactory["log10_pp_H_hh_gagabb_TH13"] = boost::factory<log10_pp_H_hh_gagabb_TH13*>();
    obsThFactory["log10_pp_H_hh_bbtautau_TH13"] = boost::factory<log10_pp_H_hh_bbtautau_TH13*>();
    obsThFactory["log10_pp_H_hh_bblnulnu_TH13"] = boost::factory<log10_pp_H_hh_bblnulnu_TH13*>();
    obsThFactory["log10_pp_H_hh_bbVV_TH13"] = boost::factory<log10_pp_H_hh_bbVV_TH13*>();
    obsThFactory["log10_tt_H_tt_TH13"] = boost::factory<log10_tt_H_tt_TH13*>();
    obsThFactory["log10_bb_H_tt_TH13"] = boost::factory<log10_bb_H_tt_TH13*>();
    obsThFactory["log10_pp_H_bb_TH13"] = boost::factory<log10_pp_H_bb_TH13*>();
    obsThFactory["Gamma_HH_THDM"] = boost::factory<Gamma_HH_THDM*>();
    obsThFactory["rHH_gg_THDM"] = boost::factory<rHH_gg_THDM*>();
    obsThFactory["BR_HH_hh_THDM"] = boost::factory<BR_HH_hh_THDM*>();
    obsThFactory["BR_HH_AA_THDM"] = boost::factory<BR_HH_AA_THDM*>();
    obsThFactory["BR_HH_HpHm_THDM"] = boost::factory<BR_HH_HpHm_THDM*>();
    obsThFactory["BR_HH_AZ_THDM"] = boost::factory<BR_HH_AZ_THDM*>();
    obsThFactory["BR_HH_HpW_THDM"] = boost::factory<BR_HH_HpW_THDM*>();

    obsThFactory["Hobs_ggF_A_tautau_ATLAS8"] = boost::factory<Hobs_ggF_A_tautau_ATLAS8*>();
    obsThFactory["Hobs_ggF_A_tautau_CMS8"] = boost::factory<Hobs_ggF_A_tautau_CMS8*>();
    obsThFactory["Hobs_bbF_A_tautau_ATLAS8"] = boost::factory<Hobs_bbF_A_tautau_ATLAS8*>();
    obsThFactory["Hobs_bbF_A_tautau_CMS8"] = boost::factory<Hobs_bbF_A_tautau_CMS8*>();
    obsThFactory["Hobs_pp_A_gaga_ATLAS8"] = boost::factory<Hobs_pp_A_gaga_ATLAS8*>();
    obsThFactory["Hobs_ggF_A_gaga_CMS8"] = boost::factory<Hobs_ggF_A_gaga_CMS8*>();
    obsThFactory["Hobs_pp_A_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_A_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_A_Zga_llga_CMS8"] = boost::factory<Hobs_pp_A_Zga_llga_CMS8*>();
    obsThFactory["Hobs_ggF_A_hZ_bbll_CMS8"] = boost::factory<Hobs_ggF_A_hZ_bbll_CMS8*>();
    obsThFactory["Hobs_ggF_A_hZ_bbZ_ATLAS8"] = boost::factory<Hobs_ggF_A_hZ_bbZ_ATLAS8*>();
    obsThFactory["Hobs_ggF_A_hZ_tautaull_CMS8"] = boost::factory<Hobs_ggF_A_hZ_tautaull_CMS8*>();
    obsThFactory["Hobs_ggF_A_hZ_tautauZ_ATLAS8"] = boost::factory<Hobs_ggF_A_hZ_tautauZ_ATLAS8*>();
    obsThFactory["Hobs_ggF_A_tt_ATLAS8"] = boost::factory<Hobs_ggF_A_tt_ATLAS8*>();
    obsThFactory["Hobs_bbF_A_bb_CMS8"] = boost::factory<Hobs_bbF_A_bb_CMS8*>();
    obsThFactory["Robs_ggF_A_tautau_ATLAS8"] = boost::factory<Robs_ggF_A_tautau_ATLAS8*>();
    obsThFactory["Robs_ggF_A_tautau_CMS8"] = boost::factory<Robs_ggF_A_tautau_CMS8*>();
    obsThFactory["Robs_bbF_A_tautau_ATLAS8"] = boost::factory<Robs_bbF_A_tautau_ATLAS8*>();
    obsThFactory["Robs_bbF_A_tautau_CMS8"] = boost::factory<Robs_bbF_A_tautau_CMS8*>();
    obsThFactory["Robs_pp_A_gaga_ATLAS8"] = boost::factory<Robs_pp_A_gaga_ATLAS8*>();
    obsThFactory["Robs_ggF_A_gaga_CMS8"] = boost::factory<Robs_ggF_A_gaga_CMS8*>();
    obsThFactory["Robs_pp_A_Zga_llga_ATLAS8"] = boost::factory<Robs_pp_A_Zga_llga_ATLAS8*>();
    obsThFactory["Robs_pp_A_Zga_llga_CMS8"] = boost::factory<Robs_pp_A_Zga_llga_CMS8*>();
    obsThFactory["Robs_ggF_A_hZ_bbll_CMS8"] = boost::factory<Robs_ggF_A_hZ_bbll_CMS8*>();
    obsThFactory["Robs_ggF_A_hZ_bbZ_ATLAS8"] = boost::factory<Robs_ggF_A_hZ_bbZ_ATLAS8*>();
    obsThFactory["Robs_ggF_A_hZ_tautaull_CMS8"] = boost::factory<Robs_ggF_A_hZ_tautaull_CMS8*>();
    obsThFactory["Robs_ggF_A_hZ_tautauZ_ATLAS8"] = boost::factory<Robs_ggF_A_hZ_tautauZ_ATLAS8*>();
    obsThFactory["Robs_ggF_A_tt_ATLAS8"] = boost::factory<Robs_ggF_A_tt_ATLAS8*>();
    obsThFactory["Robs_bbF_A_bb_CMS8"] = boost::factory<Robs_bbF_A_bb_CMS8*>();
    obsThFactory["log10_ggF_A_tautau_TH8"] = boost::factory<log10_ggF_A_tautau_TH8*>();
    obsThFactory["log10_bbF_A_tautau_TH8"] = boost::factory<log10_bbF_A_tautau_TH8*>();
    obsThFactory["log10_pp_A_gaga_TH8"] = boost::factory<log10_pp_A_gaga_TH8*>();
    obsThFactory["log10_ggF_A_gaga_TH8"] = boost::factory<log10_ggF_A_gaga_TH8*>();
    obsThFactory["log10_pp_A_Zga_llga_TH8"] = boost::factory<log10_pp_A_Zga_llga_TH8*>();
    obsThFactory["log10_ggF_A_hZ_bbll_TH8"] = boost::factory<log10_ggF_A_hZ_bbll_TH8*>();
    obsThFactory["log10_ggF_A_hZ_bbZ_TH8"] = boost::factory<log10_ggF_A_hZ_bbZ_TH8*>();
    obsThFactory["log10_ggF_A_hZ_tautaull_TH8"] = boost::factory<log10_ggF_A_hZ_tautaull_TH8*>();
    obsThFactory["log10_ggF_A_hZ_tautauZ_TH8"] = boost::factory<log10_ggF_A_hZ_tautauZ_TH8*>();
    obsThFactory["log10_ggF_A_tt_TH8"] = boost::factory<log10_ggF_A_tt_TH8*>();
    obsThFactory["log10_bbF_A_bb_TH8"] = boost::factory<log10_bbF_A_bb_TH8*>();
    obsThFactory["Hobs_ttF_A_tt_ATLAS13"] = boost::factory<Hobs_ttF_A_tt_ATLAS13*>();
    obsThFactory["Hobs_bbF_A_tt_ATLAS13"] = boost::factory<Hobs_bbF_A_tt_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_tautau_ATLAS13"] = boost::factory<Hobs_ggF_A_tautau_ATLAS13*>();
    obsThFactory["Hobs_bbF_A_tautau_ATLAS13"] = boost::factory<Hobs_bbF_A_tautau_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_tautau_CMS13"] = boost::factory<Hobs_ggF_A_tautau_CMS13*>();
    obsThFactory["Hobs_bbF_A_tautau_CMS13"] = boost::factory<Hobs_bbF_A_tautau_CMS13*>();
    obsThFactory["Hobs_pp_A_gaga_ATLAS13"] = boost::factory<Hobs_pp_A_gaga_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_gaga_CMS13"] = boost::factory<Hobs_ggF_A_gaga_CMS13*>();
    obsThFactory["Hobs_pp_A_Zga_llga_ATLAS13"] = boost::factory<Hobs_pp_A_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_pp_A_Zga_llga_CMS13"] = boost::factory<Hobs_pp_A_Zga_llga_CMS13*>();
    obsThFactory["Hobs_pp_A_Zga_qqga_CMS13"] = boost::factory<Hobs_pp_A_Zga_qqga_CMS13*>();
    obsThFactory["Hobs_ggF_A_hZ_bbZ_ATLAS13"] = boost::factory<Hobs_ggF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bbF_A_hZ_bbZ_ATLAS13"] = boost::factory<Hobs_bbF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Hobs_pp_A_bb_CMS13"] = boost::factory<Hobs_pp_A_bb_CMS13*>();
    obsThFactory["Robs_ttF_A_tt_ATLAS13"] = boost::factory<Robs_ttF_A_tt_ATLAS13*>();
    obsThFactory["Robs_bbF_A_tt_ATLAS13"] = boost::factory<Robs_bbF_A_tt_ATLAS13*>();
    obsThFactory["Robs_ggF_A_tautau_ATLAS13"] = boost::factory<Robs_ggF_A_tautau_ATLAS13*>();
    obsThFactory["Robs_bbF_A_tautau_ATLAS13"] = boost::factory<Robs_bbF_A_tautau_ATLAS13*>();
    obsThFactory["Robs_ggF_A_tautau_CMS13"] = boost::factory<Robs_ggF_A_tautau_CMS13*>();
    obsThFactory["Robs_bbF_A_tautau_CMS13"] = boost::factory<Robs_bbF_A_tautau_CMS13*>();
    obsThFactory["Robs_pp_A_gaga_ATLAS13"] = boost::factory<Robs_pp_A_gaga_ATLAS13*>();
    obsThFactory["Robs_ggF_A_gaga_CMS13"] = boost::factory<Robs_ggF_A_gaga_CMS13*>();
    obsThFactory["Robs_pp_A_Zga_llga_ATLAS13"] = boost::factory<Robs_pp_A_Zga_llga_ATLAS13*>();
    obsThFactory["Robs_pp_A_Zga_llga_CMS13"] = boost::factory<Robs_pp_A_Zga_llga_CMS13*>();
    obsThFactory["Robs_pp_A_Zga_qqga_CMS13"] = boost::factory<Robs_pp_A_Zga_qqga_CMS13*>();
    obsThFactory["Robs_ggF_A_hZ_bbZ_ATLAS13"] = boost::factory<Robs_ggF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Robs_bbF_A_hZ_bbZ_ATLAS13"] = boost::factory<Robs_bbF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Robs_pp_A_bb_CMS13"] = boost::factory<Robs_pp_A_bb_CMS13*>();
    obsThFactory["log10_ggF_A_tautau_TH13"] = boost::factory<log10_ggF_A_tautau_TH13*>();
    obsThFactory["log10_bbF_A_tautau_TH13"] = boost::factory<log10_bbF_A_tautau_TH13*>();
    obsThFactory["log10_pp_A_gaga_TH13"] = boost::factory<log10_pp_A_gaga_TH13*>();
    obsThFactory["log10_ggF_A_gaga_TH13"] = boost::factory<log10_ggF_A_gaga_TH13*>();
    obsThFactory["log10_pp_A_Zga_TH13"] = boost::factory<log10_pp_A_Zga_TH13*>();
    obsThFactory["log10_ggF_A_hZ_bbZ_TH13"] = boost::factory<log10_ggF_A_hZ_bbZ_TH13*>();
    obsThFactory["log10_bbF_A_hZ_bbZ_TH13"] = boost::factory<log10_bbF_A_hZ_bbZ_TH13*>();
    obsThFactory["log10_ttF_A_tt_TH13"] = boost::factory<log10_ttF_A_tt_TH13*>();
    obsThFactory["log10_bbF_A_tt_TH13"] = boost::factory<log10_bbF_A_tt_TH13*>();
    obsThFactory["log10_pp_A_bb_TH13"] = boost::factory<log10_pp_A_bb_TH13*>();
    obsThFactory["Gamma_A_THDM"] = boost::factory<Gamma_A_THDM*>();
    obsThFactory["rA_gg_THDM"] = boost::factory<rA_gg_THDM*>();
    obsThFactory["BR_A_HZ_THDM"] = boost::factory<BR_A_HZ_THDM*>();
    obsThFactory["BR_A_hZ_THDM"] = boost::factory<BR_A_hZ_THDM*>();
    obsThFactory["BR_A_HpW_THDM"] = boost::factory<BR_A_HpW_THDM*>();

    obsThFactory["Hobs_pp_Hpm_taunu_ATLAS8"] = boost::factory<Hobs_pp_Hpm_taunu_ATLAS8*>();
    obsThFactory["Hobs_pp_Hp_taunu_CMS8"] = boost::factory<Hobs_pp_Hp_taunu_CMS8*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS8"] = boost::factory<Hobs_pp_Hp_tb_ATLAS8*>();
    obsThFactory["Hobs_pp_Hp_tb_CMS8"] = boost::factory<Hobs_pp_Hp_tb_CMS8*>();
    obsThFactory["Hobs_pp_Hpm_taunu_ATLAS13"] = boost::factory<Hobs_pp_Hpm_taunu_ATLAS13*>();
    obsThFactory["Hobs_pp_Hpm_taunu_CMS13"] = boost::factory<Hobs_pp_Hpm_taunu_CMS13*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13_1"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13_1*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13_2"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13_2*>();
    obsThFactory["Robs_pp_Hpm_taunu_ATLAS8"] = boost::factory<Robs_pp_Hpm_taunu_ATLAS8*>();
    obsThFactory["Robs_pp_Hp_taunu_CMS8"] = boost::factory<Robs_pp_Hp_taunu_CMS8*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS8"] = boost::factory<Robs_pp_Hp_tb_ATLAS8*>();
    obsThFactory["Robs_pp_Hp_tb_CMS8"] = boost::factory<Robs_pp_Hp_tb_CMS8*>();
    obsThFactory["Robs_pp_Hpm_taunu_ATLAS13"] = boost::factory<Robs_pp_Hpm_taunu_ATLAS13*>();
    obsThFactory["Robs_pp_Hpm_taunu_CMS13"] = boost::factory<Robs_pp_Hpm_taunu_CMS13*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13"] = boost::factory<Robs_pp_Hp_tb_ATLAS13*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13_1"] = boost::factory<Robs_pp_Hp_tb_ATLAS13_1*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13_2"] = boost::factory<Robs_pp_Hp_tb_ATLAS13_2*>();
    obsThFactory["log10_pp_Hpm_taunu_TH8"] = boost::factory<log10_pp_Hpm_taunu_TH8*>();
    obsThFactory["log10_pp_Hp_tb_TH8"] = boost::factory<log10_pp_Hp_tb_TH8*>();
    obsThFactory["log10_pp_Hp_TH8"] = boost::factory<log10_pp_Hp_TH8*>();
    obsThFactory["log10_pp_Hpm_taunu_TH13"] = boost::factory<log10_pp_Hpm_taunu_TH13*>();
    obsThFactory["log10_pp_Hp_tb_TH13"] = boost::factory<log10_pp_Hp_tb_TH13*>();
    obsThFactory["Gamma_Hp_THDM"] = boost::factory<Gamma_Hp_THDM*>();

    obsThFactory["tanbeta"] = boost::factory<tanbeta*>();
    obsThFactory["mHl_THDM"] = boost::factory<mass_mHl*>();
    obsThFactory["mHh"] = boost::factory<mass_mHh*>();
    obsThFactory["mA"] = boost::factory<mass_mA*>();
    obsThFactory["mHp"] = boost::factory<mass_mHp*>();
    obsThFactory["mHlmmA"] = boost::factory<massdifference_mHlmmA*>();
    obsThFactory["mAmmHl"] = boost::factory<massdifference_mAmmHl*>();
    obsThFactory["mHlmmHp"] = boost::factory<massdifference_mHlmmHp*>();
    obsThFactory["mHpmmHl"] = boost::factory<massdifference_mHpmmHl*>();
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
    obsThFactory["lambda345"] = boost::factory<lambda345*>();
    obsThFactory["g_hhh"] = boost::factory<g_hhh*>();
    obsThFactory["g_hhHh"] = boost::factory<g_hhHh*>();
    obsThFactory["g_hHhHh"] = boost::factory<g_hHhHh*>();
    obsThFactory["g_HhHhHh"] = boost::factory<g_HhHhHh*>();
    obsThFactory["g_hAA"] = boost::factory<g_hAA*>();
    obsThFactory["g_HhAA"] = boost::factory<g_HhAA*>();
    obsThFactory["g_hHpHm"] = boost::factory<g_hHpHm*>();
    obsThFactory["g_HhHpHm"] = boost::factory<g_HhHpHm*>();
    obsThFactory["Y1_THDM"] = boost::factory<Y1_THDM*>();
    obsThFactory["Y2_THDM"] = boost::factory<Y2_THDM*>();
    obsThFactory["Y3_THDM"] = boost::factory<Y3_THDM*>();
    obsThFactory["Z1_THDM"] = boost::factory<Z1_THDM*>();
    obsThFactory["Z2_THDM"] = boost::factory<Z2_THDM*>();
    obsThFactory["Z3_THDM"] = boost::factory<Z3_THDM*>();
    obsThFactory["Z4_THDM"] = boost::factory<Z4_THDM*>();
    obsThFactory["Z5_THDM"] = boost::factory<Z5_THDM*>();
    obsThFactory["Z6_THDM"] = boost::factory<Z6_THDM*>();
    obsThFactory["Z7_THDM"] = boost::factory<Z7_THDM*>();
    obsThFactory["xi0_THDM"] = boost::factory<xi0_THDM*>();
    obsThFactory["xi1_THDM"] = boost::factory<xi1_THDM*>();
    obsThFactory["xi3_THDM"] = boost::factory<xi3_THDM*>();
    obsThFactory["eta00_THDM"] = boost::factory<eta00_THDM*>();
    obsThFactory["eta3_THDM"] = boost::factory<eta3_THDM*>();
    obsThFactory["E11_THDM"] = boost::factory<E11_THDM*>();
    obsThFactory["E22_THDM"] = boost::factory<E22_THDM*>();
    obsThFactory["E33_THDM"] = boost::factory<E33_THDM*>();
    obsThFactory["HHlambda1"] = boost::factory<HHlambda1*>();
    obsThFactory["HHlambda2"] = boost::factory<HHlambda2*>();
    obsThFactory["HHlambda3"] = boost::factory<HHlambda3*>();
    obsThFactory["HHlambda4"] = boost::factory<HHlambda4*>();
    obsThFactory["HHlambda5"] = boost::factory<HHlambda5*>();
    obsThFactory["HHlambda6"] = boost::factory<HHlambda6*>();

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
    obsThFactory["RBDstartaunu"] = boost::factory<RBDstartaunu*>();
    obsThFactory["RBDtaunu"] = boost::factory<RBDtaunu*>();
    obsThFactory["obsBDtaunu_SM"] = boost::factory<obsBDtaunu_SM*>();
    obsThFactory["obsBDtaunu_A"] = boost::factory<obsBDtaunu_A*>();
    obsThFactory["obsBDtaunu_B"] = boost::factory<obsBDtaunu_B*>();
    obsThFactory["obsBDstartaunu_SM"] = boost::factory<obsBDstartaunu_SM*>();
    obsThFactory["obsBDstartaunu_A"] = boost::factory<obsBDstartaunu_A*>();
    obsThFactory["obsBDstartaunu_B"] = boost::factory<obsBDstartaunu_B*>();
    obsThFactory["THDMgminus2_mu"] = boost::factory<THDMgminus2_mu*>();

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

    obsThFactory["unitaritya10odd"] = boost::factory<unitarityNLO14*>();
    obsThFactory["unitaritya11odd"] = boost::factory<unitarityNLO24*>();
    obsThFactory["unitaritya00evenp"] = boost::factory<unitarityNLOev1*>();
    obsThFactory["unitaritya00evenm"] = boost::factory<unitarityNLOev2*>();
    obsThFactory["unitaritya00oddp"] = boost::factory<unitarityNLOev3*>();
    obsThFactory["unitaritya00oddm"] = boost::factory<unitarityNLOev4*>();
    obsThFactory["unitaritya01evenp"] = boost::factory<unitarityNLOev5*>();
    obsThFactory["unitaritya01evenm"] = boost::factory<unitarityNLOev6*>();
    obsThFactory["unitaritya01oddp"] = boost::factory<unitarityNLOev9*>();
    obsThFactory["unitaritya01oddm"] = boost::factory<unitarityNLOev10*>();
    obsThFactory["unitaritya11evenp"] = boost::factory<unitarityNLOev13*>();
    obsThFactory["unitaritya11evenm"] = boost::factory<unitarityNLOev14*>();
    obsThFactory["unitarityRp1"] = boost::factory<unitarityRp1*>();
    obsThFactory["unitarityRp2"] = boost::factory<unitarityRp2*>();
    obsThFactory["unitarityRp3"] = boost::factory<unitarityRp3*>();
    obsThFactory["unitarityRp4"] = boost::factory<unitarityRp4*>();
    obsThFactory["unitarityRp5"] = boost::factory<unitarityRp5*>();
    obsThFactory["unitarityRp6"] = boost::factory<unitarityRp6*>();
    obsThFactory["unitarityRp9"] = boost::factory<unitarityRp9*>();
    obsThFactory["unitarityRp10"] = boost::factory<unitarityRp10*>();
    obsThFactory["unitarityRp13"] = boost::factory<unitarityRp13*>();
    obsThFactory["unitarityRp14"] = boost::factory<unitarityRp14*>();
    obsThFactory["unitarityRp19"] = boost::factory<unitarityRp19*>();
    obsThFactory["unitarityRp20"] = boost::factory<unitarityRp20*>();
    obsThFactory["unitarityR1"] = boost::factory<unitarityR1*>();
    obsThFactory["unitarityR2"] = boost::factory<unitarityR2*>();
    obsThFactory["unitarityR3"] = boost::factory<unitarityR3*>();
    obsThFactory["unitarityR4"] = boost::factory<unitarityR4*>();
    obsThFactory["unitarityR5"] = boost::factory<unitarityR5*>();
    obsThFactory["unitarityR6"] = boost::factory<unitarityR6*>();
    obsThFactory["unitarityR9"] = boost::factory<unitarityR9*>();
    obsThFactory["unitarityR10"] = boost::factory<unitarityR10*>();
    obsThFactory["unitarityR13"] = boost::factory<unitarityR13*>();
    obsThFactory["unitarityR14"] = boost::factory<unitarityR14*>();
    obsThFactory["unitarityR19"] = boost::factory<unitarityR19*>();
    obsThFactory["unitarityR20"] = boost::factory<unitarityR20*>();
    obsThFactory["unitaritya00evenpRe"] = boost::factory<unitaritya00evenpRe*>();
    obsThFactory["unitaritya00evenpIm"] = boost::factory<unitaritya00evenpIm*>();
    obsThFactory["unitaritya00evenmRe"] = boost::factory<unitaritya00evenmRe*>();
    obsThFactory["unitaritya00evenmIm"] = boost::factory<unitaritya00evenmIm*>();
    obsThFactory["unitaritya00oddpRe"] = boost::factory<unitaritya00oddpRe*>();
    obsThFactory["unitaritya00oddpIm"] = boost::factory<unitaritya00oddpIm*>();
    obsThFactory["unitaritya00oddmRe"] = boost::factory<unitaritya00oddmRe*>();
    obsThFactory["unitaritya00oddmIm"] = boost::factory<unitaritya00oddmIm*>();
    obsThFactory["unitaritya01evenpRe"] = boost::factory<unitaritya01evenpRe*>();
    obsThFactory["unitaritya01evenpIm"] = boost::factory<unitaritya01evenpIm*>();
    obsThFactory["unitaritya01evenmRe"] = boost::factory<unitaritya01evenmRe*>();
    obsThFactory["unitaritya01evenmIm"] = boost::factory<unitaritya01evenmIm*>();
    obsThFactory["unitaritya01oddpRe"] = boost::factory<unitaritya01oddpRe*>();
    obsThFactory["unitaritya01oddpIm"] = boost::factory<unitaritya01oddpIm*>();
    obsThFactory["unitaritya01oddmRe"] = boost::factory<unitaritya01oddmRe*>();
    obsThFactory["unitaritya01oddmIm"] = boost::factory<unitaritya01oddmIm*>();
    obsThFactory["unitaritya10oddRe"] = boost::factory<unitaritya10oddRe*>();
    obsThFactory["unitaritya10oddIm"] = boost::factory<unitaritya10oddIm*>();
    obsThFactory["unitaritya11evenpRe"] = boost::factory<unitaritya11evenpRe*>();
    obsThFactory["unitaritya11evenpIm"] = boost::factory<unitaritya11evenpIm*>();
    obsThFactory["unitaritya11evenmRe"] = boost::factory<unitaritya11evenmRe*>();
    obsThFactory["unitaritya11evenmIm"] = boost::factory<unitaritya11evenmIm*>();
    obsThFactory["unitaritya11oddRe"] = boost::factory<unitaritya11oddRe*>();
    obsThFactory["unitaritya11oddIm"] = boost::factory<unitaritya11oddIm*>();
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  GeneralTHDM observables  -----
    obsThFactory["mH1"] = boost::factory<mH1_GTHDM*>();
    obsThFactory["mH2"] = boost::factory<mH2_GTHDM*>();
    obsThFactory["mH3"] = boost::factory<mH3_GTHDM*>();
    obsThFactory["mHlight"] = boost::factory<mHlight_GTHDM*>();
    obsThFactory["mHmedium"] = boost::factory<mHmedium_GTHDM*>();
    obsThFactory["mHheavy"] = boost::factory<mHheavy_GTHDM*>();
    obsThFactory["mHp"] = boost::factory<mHp_GTHDM*>();
    obsThFactory["mH3mmH2"] = boost::factory<mH3mmH2_GTHDM*>();
    obsThFactory["mH3mmHp"] = boost::factory<mH3mmHp_GTHDM*>();
    obsThFactory["mH3mmH1"] = boost::factory<mH3mmH1_GTHDM*>();
    obsThFactory["mH2mmHp"] = boost::factory<mH2mmHp_GTHDM*>();
    obsThFactory["mH2mmH1"] = boost::factory<mH2mmH1_GTHDM*>();
    obsThFactory["mHpmmH1"] = boost::factory<mHpmmH1_GTHDM*>();
    obsThFactory["mH1sq"] = boost::factory<mH1sq_GTHDM*>();
    obsThFactory["mH2sq"] = boost::factory<mH2sq_GTHDM*>();
    obsThFactory["mH3sq"] = boost::factory<mH3sq_GTHDM*>();
    obsThFactory["Msq11"] = boost::factory<Msq11_GTHDM*>();
    obsThFactory["Msq12"] = boost::factory<Msq12_GTHDM*>();
    obsThFactory["Msq13"] = boost::factory<Msq13_GTHDM*>();
    obsThFactory["Msq22"] = boost::factory<Msq22_GTHDM*>();
    obsThFactory["Msq23"] = boost::factory<Msq23_GTHDM*>();
    obsThFactory["Msq33"] = boost::factory<Msq33_GTHDM*>();
    obsThFactory["M2_GTHDM"] = boost::factory<M2_GTHDM*>();
    obsThFactory["m11_2_GTHDM"] = boost::factory<m11_2_GTHDM*>();
    obsThFactory["m22_2_GTHDM"] = boost::factory<m22_2_GTHDM*>();
    obsThFactory["Imm12_2_GTHDM"] = boost::factory<Imm12_2_GTHDM*>();
    obsThFactory["lambda1_GTHDM"] = boost::factory<lambda1_GTHDM*>();
    obsThFactory["lambda2_GTHDM"] = boost::factory<lambda2_GTHDM*>();
    obsThFactory["lambda3_GTHDM"] = boost::factory<lambda3_GTHDM*>();
    obsThFactory["lambda4_GTHDM"] = boost::factory<lambda4_GTHDM*>();
    obsThFactory["Relambda5_GTHDM"] = boost::factory<Relambda5_GTHDM*>();
    obsThFactory["v1_GTHDM"] = boost::factory<v1_GTHDM*>();
    obsThFactory["v2_GTHDM"] = boost::factory<v2_GTHDM*>();

    obsThFactory["unitarity1_GTHDM"] = boost::factory<unitarity1_GTHDM*>();
    obsThFactory["unitarity2_GTHDM"] = boost::factory<unitarity2_GTHDM*>();
    obsThFactory["unitarity3_GTHDM"] = boost::factory<unitarity3_GTHDM*>();
    obsThFactory["unitarity4_GTHDM"] = boost::factory<unitarity4_GTHDM*>();
    obsThFactory["unitarity5_GTHDM"] = boost::factory<unitarity5_GTHDM*>();
    obsThFactory["unitarity6_GTHDM"] = boost::factory<unitarity6_GTHDM*>();
    obsThFactory["unitarity7_GTHDM"] = boost::factory<unitarity7_GTHDM*>();
    obsThFactory["unitarity8_GTHDM"] = boost::factory<unitarity8_GTHDM*>();
    obsThFactory["unitarity9_GTHDM"] = boost::factory<unitarity9_GTHDM*>();
    obsThFactory["unitarity10_GTHDM"] = boost::factory<unitarity10_GTHDM*>();
    obsThFactory["unitarity11_GTHDM"] = boost::factory<unitarity11_GTHDM*>();
    obsThFactory["unitarity12_GTHDM"] = boost::factory<unitarity12_GTHDM*>();

    obsThFactory["stability1_GTHDM"] = boost::factory<stability1_GTHDM*>();
    obsThFactory["stability2_GTHDM"] = boost::factory<stability2_GTHDM*>();
    obsThFactory["stability3_GTHDM"] = boost::factory<stability3_GTHDM*>();
    obsThFactory["stability4_GTHDM"] = boost::factory<stability4_GTHDM*>();

//    obsThFactory["GTHDMgminus2_mu"] = boost::factory<GeneralTHDMgminus2_mu*>();

    /** END: REMOVE FROM THE PACKAGE **/
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  LeftRightSymmetric model observables  -----
    obsThFactory["mu1_2_LRSM"] = boost::factory<mu1_2_LRSM*>();
    obsThFactory["mu2_2_LRSM"] = boost::factory<mu2_2_LRSM*>();
    obsThFactory["mu3_2_LRSM"] = boost::factory<mu3_2_LRSM*>();
    obsThFactory["rho2_LRSM"] = boost::factory<rho2_LRSM*>();
    obsThFactory["rho3_LRSM"] = boost::factory<rho3_LRSM*>();
    obsThFactory["alpha3_LRSM"] = boost::factory<alpha3_LRSM*>();
    obsThFactory["mH00_LRSM"] = boost::bind(boost::factory<MH0_LRSM*>(), _1, 0);
    obsThFactory["mH01_LRSM"] = boost::bind(boost::factory<MH0_LRSM*>(), _1, 1);
    obsThFactory["mH02_LRSM"] = boost::bind(boost::factory<MH0_LRSM*>(), _1, 2);
    obsThFactory["mH03_LRSM"] = boost::bind(boost::factory<MH0_LRSM*>(), _1, 3);
    obsThFactory["mH04_LRSM"] = boost::bind(boost::factory<MH0_LRSM*>(), _1, 4);
    obsThFactory["mH05_LRSM"] = boost::factory<MH05_LRSM*>();
    obsThFactory["mH06_LRSM"] = boost::factory<MH06_LRSM*>();
    obsThFactory["MH01_app1"] = boost::factory<MH01_app1*>();
    obsThFactory["MH01_app"] = boost::bind(boost::factory<MH0_app*>(), _1, 0);
    obsThFactory["MH02_app"] = boost::bind(boost::factory<MH0_app*>(), _1, 1);
    obsThFactory["MH03_app"] = boost::bind(boost::factory<MH0_app*>(), _1, 2);
    obsThFactory["MH04_app"] = boost::bind(boost::factory<MH0_app*>(), _1, 3);

    /** END: REMOVE FROM THE PACKAGE **/

    /** BEGIN: REMOVE FROM THE PACKAGE **/
    obsThFactory["unitarity1_GTHDM"] = boost::factory<unitarity1_GTHDM*>();
    /** END: REMOVE FROM THE PACKAGE **/
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
