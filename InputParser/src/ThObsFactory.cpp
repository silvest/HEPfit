/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObservable.h"
#include "ParamObs.h"
#include "EWObservables.h"
#include "HiggsThObservables.h"
#include "FlavourObservables.h"
/** BEGIN: REMOVE FROM THE PACKAGE **/
#include "LeptonFlavourObservables.h"
#include "SUSYObservables.h"
#include "THDMObservables.h"
/** END: REMOVE FROM THE PACKAGE **/
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

ThObsFactory::ThObsFactory()
{
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
    const double sqrt_s_TeV = 1.96;
    const double sqrt_s_TLEP = .24;
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
    obsThFactory["ggH196"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_TeV);
    obsThFactory["VBF196"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_TeV);
    obsThFactory["VH196"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_TeV);
    obsThFactory["ttH196"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_TeV);
    obsThFactory["eeZH240"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_TLEP);
    obsThFactory["BrHggRatio"] = boost::factory<BrHtoggRatio*>();
    obsThFactory["BrHWWRatio"] = boost::factory<BrHtoWWRatio*>();
    obsThFactory["BrHZZRatio"] = boost::factory<BrHtoZZRatio*>();
    obsThFactory["BrHZgaRatio"] = boost::factory<BrHtoZgaRatio*>();
    obsThFactory["BrHgagaRatio"] = boost::factory<BrHtogagaRatio*>();
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
    obsThFactory["alpha"] = boost::factory<Alpha*>();
    obsThFactory["alpha_2a"] = boost::factory<Alpha_2a*>();
    obsThFactory["gamma"] = boost::factory<CKMGamma*>();
    obsThFactory["Abslam_t"] = boost::factory<Abslam_t*>();
    obsThFactory["Abslam_c"] = boost::factory<Abslam_c*>();
    obsThFactory["Abslam_u"] = boost::factory<Abslam_u*>();
    obsThFactory["Abslam_td"] = boost::factory<Abslam_td*>();
    obsThFactory["Abslam_cd"] = boost::factory<Abslam_cd*>();
    obsThFactory["Abslam_ud"] = boost::factory<Abslam_ud*>();
    obsThFactory["Abslam_ts"] = boost::factory<Abslam_ts*>();
    obsThFactory["Abslam_cs"] = boost::factory<Abslam_cs*>();
    obsThFactory["Abslam_us"] = boost::factory<Abslam_us*>();
    //----- B(s) to mu mu  -----
    obsThFactory["BR_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 1);
    obsThFactory["BRbar_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 2);
    obsThFactory["Amumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 3);
    obsThFactory["Smumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 4);
    obsThFactory["BR_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 1);
    obsThFactory["BRbar_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 2);
    obsThFactory["Amumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 3);
    obsThFactory["Smumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 4);
    //----- b to s gamma  -----
    obsThFactory["BR_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 1);
    obsThFactory["BR_CPodd_bsgamma"] = boost::bind(boost::factory<Bsgamma*>(), _1, 2);
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
    obsThFactory["ACP_BKstgamma"] = boost::bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    //----- B to PHI gamma  -----
    obsThFactory["BR_Bsphigamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ACP_Bsphigamma"] = boost::bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    
    //----- B to tau nu  -----
    obsThFactory["btaunu"] = boost::factory<Btaunu*>();
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  Lepton Flavour observables  -----
    obsThFactory["mu_e_gamma"] = boost::factory<mu_e_gamma*>();
    obsThFactory["tau_mu_gamma"] = boost::factory<tau_mu_gamma*>();
    obsThFactory["tau_e_gamma"] = boost::factory<tau_e_gamma*>();
    obsThFactory["mu_3e"] = boost::factory<mu_3e*>();
    obsThFactory["tau_3mu"] = boost::factory<tau_3mu*>();
    obsThFactory["tau_3e"] = boost::factory<tau_3e*>();
    /** END: REMOVE FROM THE PACKAGE **/
    
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    //-----  SUSY spectra and observables  -----
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

    obsThFactory["Hobs_ggF_H_tautau_ATLAS"] = boost::factory<Hobs_ggF_H_tautau_ATLAS*>();
    obsThFactory["Hobs_ggF_H_tautau_CMS"] = boost::factory<Hobs_ggF_H_tautau_CMS*>();
    obsThFactory["Hobs_bbF_H_tautau_ATLAS"] = boost::factory<Hobs_bbF_H_tautau_ATLAS*>();
    obsThFactory["Hobs_bbF_H_tautau_CMS"] = boost::factory<Hobs_bbF_H_tautau_CMS*>();
    obsThFactory["Hobs_ggF_H_gaga_ATLAS"] = boost::factory<Hobs_ggF_H_gaga_ATLAS*>();
    obsThFactory["Hobs_ggF_H_gaga_CMS"] = boost::factory<Hobs_ggF_H_gaga_CMS*>();
    obsThFactory["Hobs_pp_H_ZZ_CMS"] = boost::factory<Hobs_pp_H_ZZ_CMS*>();
    obsThFactory["Hobs_ggF_H_WW_ATLAS"] = boost::factory<Hobs_ggF_H_WW_ATLAS*>();
    obsThFactory["Hobs_VBF_H_WW_ATLAS"] = boost::factory<Hobs_VBF_H_WW_ATLAS*>();
    obsThFactory["Hobs_ggF_H_hh_ATLAS"] = boost::factory<Hobs_ggF_H_hh_ATLAS*>();
    obsThFactory["Hobs_ggF_H_hh_bbtautau_CMS"] = boost::factory<Hobs_ggF_H_hh_bbtautau_CMS*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_CMS"] = boost::factory<Hobs_pp_H_hh_bbbb_CMS*>();
    obsThFactory["Hobs_pp_H_hh_gagabb_CMS"] = boost::factory<Hobs_pp_H_hh_gagabb_CMS*>();
    obsThFactory["Hobs_pp_H_tt_ATLAS"] = boost::factory<Hobs_pp_H_tt_ATLAS*>();
    obsThFactory["Hobs_bbF_H_bb_CMS"] = boost::factory<Hobs_bbF_H_bb_CMS*>();

    obsThFactory["log10_ggF_H_tautau_TH"] = boost::factory<log10_ggF_H_tautau_TH*>();
    obsThFactory["log10_bbF_H_tautau_TH"] = boost::factory<log10_bbF_H_tautau_TH*>();
    obsThFactory["log10_ggF_H_gaga_TH"] = boost::factory<log10_ggF_H_gaga_TH*>();
    obsThFactory["log10_pp_H_ZZ_TH"] = boost::factory<log10_pp_H_ZZ_TH*>();
    obsThFactory["log10_ggF_H_WW_TH"] = boost::factory<log10_ggF_H_WW_TH*>();
    obsThFactory["log10_VBF_H_WW_TH"] = boost::factory<log10_VBF_H_WW_TH*>();
    obsThFactory["log10_ggF_H_hh_TH"] = boost::factory<log10_ggF_H_hh_TH*>();
    obsThFactory["log10_ggF_H_hh_bbtautau_TH"] = boost::factory<log10_ggF_H_hh_bbtautau_TH*>();
    obsThFactory["log10_pp_H_hh_bbbb_TH"] = boost::factory<log10_pp_H_hh_bbbb_TH*>();
    obsThFactory["log10_pp_H_hh_gagabb_TH"] = boost::factory<log10_pp_H_hh_gagabb_TH*>();
    obsThFactory["log10_pp_H_tt_TH"] = boost::factory<log10_pp_H_tt_TH*>();
    obsThFactory["log10_bbF_H_bb_TH"] = boost::factory<log10_bbF_H_bb_TH*>();

    obsThFactory["Hobs_ggF_A_tautau_ATLAS"] = boost::factory<Hobs_ggF_A_tautau_ATLAS*>();
    obsThFactory["Hobs_ggF_A_tautau_CMS"] = boost::factory<Hobs_ggF_A_tautau_CMS*>();
    obsThFactory["Hobs_bbF_A_tautau_ATLAS"] = boost::factory<Hobs_bbF_A_tautau_ATLAS*>();
    obsThFactory["Hobs_bbF_A_tautau_CMS"] = boost::factory<Hobs_bbF_A_tautau_CMS*>();
    obsThFactory["Hobs_ggF_A_gaga_ATLAS"] = boost::factory<Hobs_ggF_A_gaga_ATLAS*>();
    obsThFactory["Hobs_ggF_A_gaga_CMS"] = boost::factory<Hobs_ggF_A_gaga_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_bbll_CMS"] = boost::factory<Hobs_ggF_A_hZ_bbll_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_bbZ_ATLAS"] = boost::factory<Hobs_ggF_A_hZ_bbZ_ATLAS*>();
    obsThFactory["Hobs_ggF_A_hZ_tautaull_CMS"] = boost::factory<Hobs_ggF_A_hZ_tautaull_CMS*>();
    obsThFactory["Hobs_ggF_A_hZ_tautauZ_ATLAS"] = boost::factory<Hobs_ggF_A_hZ_tautauZ_ATLAS*>();
    obsThFactory["Hobs_pp_A_tt_ATLAS"] = boost::factory<Hobs_pp_A_tt_ATLAS*>();
    obsThFactory["Hobs_bbF_A_bb_CMS"] = boost::factory<Hobs_bbF_A_bb_CMS*>();

    obsThFactory["log10_ggF_A_tautau_TH"] = boost::factory<log10_ggF_A_tautau_TH*>();
    obsThFactory["log10_bbF_A_tautau_TH"] = boost::factory<log10_bbF_A_tautau_TH*>();
    obsThFactory["log10_ggF_A_gaga_TH"] = boost::factory<log10_ggF_A_gaga_TH*>();
    obsThFactory["log10_ggF_A_hZ_bbll_TH"] = boost::factory<log10_ggF_A_hZ_bbll_TH*>();
    obsThFactory["log10_ggF_A_hZ_bbZ_TH"] = boost::factory<log10_ggF_A_hZ_bbZ_TH*>();
    obsThFactory["log10_ggF_A_hZ_tautaull_TH"] = boost::factory<log10_ggF_A_hZ_tautaull_TH*>();
    obsThFactory["log10_ggF_A_hZ_tautauZ_TH"] = boost::factory<log10_ggF_A_hZ_tautauZ_TH*>();
    obsThFactory["log10_pp_A_tt_TH"] = boost::factory<log10_pp_A_tt_TH*>();
    obsThFactory["log10_bbF_A_bb_TH"] = boost::factory<log10_bbF_A_bb_TH*>();

    obsThFactory["mHh"] = boost::factory<mass_mHh*>();
    obsThFactory["mA"] = boost::factory<mass_mA*>();
    obsThFactory["mHp"] = boost::factory<mass_mHp*>();
    obsThFactory["mHhmmA"] = boost::factory<massdifference_mHhmmA*>();
    obsThFactory["mAmmHh"] = boost::factory<massdifference_mAmmHh*>();
    obsThFactory["mHhmmHp"] = boost::factory<massdifference_mHhmmHp*>();
    obsThFactory["mHpmmHh"] = boost::factory<massdifference_mHpmmHh*>();
    obsThFactory["mAmmHp"] = boost::factory<massdifference_mAmmHp*>();
    obsThFactory["mHpmmA"] = boost::factory<massdifference_mHpmmA*>();
    obsThFactory["lambda1"] = boost::factory<lambda1*>();
    obsThFactory["lambda2"] = boost::factory<lambda2*>();
    obsThFactory["lambda3"] = boost::factory<lambda3*>();
    obsThFactory["lambda4"] = boost::factory<lambda4*>();
    obsThFactory["lambda5"] = boost::factory<lambda5*>();

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
    /** END: REMOVE FROM THE PACKAGE **/
}

void ThObsFactory::addObsToFactory(const std::string name, boost::function<ThObservable*(const StandardModel&) > funct)
{
    obsThFactory[name] = funct;
}

ThObservable * ThObsFactory::CreateThMethod(const std::string& name, const StandardModel& model) const
{
    if (model.isModelParam(name))
        return new ParamObs(model, name);
    if (obsThFactory.find(name) == obsThFactory.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThObsFactory");
    return (obsThFactory.at(name)(model));
}
