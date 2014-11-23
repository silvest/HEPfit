#include "ThObsFactory.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <EWObservables.h>
#include <HiggsThObservables.h>
#include <FlavourObservables.h>
#include <LeptonFlavourObservables.h>
#include <SUSYObservables.h>
#include <ParamObs.h>

ThObsFactory::ThObsFactory(){
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
    obsThFactory["ggH7"] = boost::bind(boost::factory<muggH*>(),  _1, sqrt_s_LHC7);
    obsThFactory["VBF7"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC7);
    obsThFactory["VH7"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ttH7"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ggH8"] = boost::bind(boost::factory<muggH*>(),  _1, sqrt_s_LHC8);
    obsThFactory["VBF8"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH8"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH8"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH196"] = boost::bind(boost::factory<muggH*>(),  _1, sqrt_s_TeV);
    obsThFactory["VBF196"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_TeV);
    obsThFactory["VH196"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_TeV);
    obsThFactory["ttH196"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_TeV);
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

//-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["Dmd1"] = boost::factory<DmBd*>();
    obsThFactory["Dms1"] = boost::factory<DmBs*>();
    obsThFactory["SJPsiK"] = boost::factory<SJPsiK*>();
    obsThFactory["SJPsiPhi"] = boost::factory<SJPsiPhi*>();
    obsThFactory["M12D"] = boost::factory<M12D*>();
    obsThFactory["ArgD"] = boost::factory<ArgD*>();
    obsThFactory["EpsilonK"] = boost::factory<EpsilonK*>();
    obsThFactory["DmK"] = boost::factory<DmK*>();
    //----- eps'/eps  -----
    obsThFactory["EpsiloP_o_Epsilon"] = boost::factory<EpsilonP_O_Epsilon*>();
    //----- CKM  -----
    obsThFactory["Vud"] = boost::factory<Vud*>();
    obsThFactory["Vus"] = boost::factory<Vus*>();
    obsThFactory["Vub"] = boost::factory<Vub*>();
    obsThFactory["Vcb"] = boost::factory<Vcb*>();
    obsThFactory["alpha"] = boost::factory<Alpha*>();
    obsThFactory["alpha_2a"] = boost::factory<Alpha_2a*>();
    obsThFactory["gamma"] = boost::factory<CKMGamma*>();
    //----- B(s) to mu mu  -----
    obsThFactory["BR_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 1);
    obsThFactory["BRbar_Bdmumu"] = boost::bind(boost::factory<Bdmumu*>(), _1, 2);
    obsThFactory["Amumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 3);
    obsThFactory["Smumu_Bd"] = boost::bind(boost::factory<Bdmumu*>(), _1, 4);
    obsThFactory["BR_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 1);
    obsThFactory["BRbar_Bsmumu"] = boost::bind(boost::factory<Bsmumu*>(), _1, 2);
    obsThFactory["Amumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 3);
    obsThFactory["Smumu_Bs"] = boost::bind(boost::factory<Bsmumu*>(), _1, 4);
    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2_BdKstmu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3_BdKstmu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4p_BdKstmu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5p_BdKstmu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6p_BdKstmu"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Gammap_BdKstmu"] = boost::bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_FB_BdKstmu"] = boost::bind(boost::factory<A_FB*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BR_BdKstmu"] = boost::bind(boost::factory<BR_BKstarll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["ACP_BdKstmu"] = boost::bind(boost::factory<ACP*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P3CP_BdKstmu"] = boost::bind(boost::factory<P3CP*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_L_BdKstmu"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["M_1p_BdKstmu"] = boost::bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["M_2p_BdKstmu"] = boost::bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    //----- B to phi ll  -----
    obsThFactory["P_1_Bsphimu"] = boost::bind(boost::factory<P_1*>(), _1, StandardModel::B_S , StandardModel::PHI , StandardModel::MU);
    obsThFactory["P_2_Bsphimu"] = boost::bind(boost::factory<P_2*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_3_Bsphimu"] = boost::bind(boost::factory<P_3*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_4p_Bsphimu"] = boost::bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_5p_Bsphimu"] = boost::bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_6p_Bsphimu"] = boost::bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["Gammap_Bsphimu"] = boost::bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_FB_Bsphimu"] = boost::bind(boost::factory<A_FB*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["BR_Bsphimu"] = boost::bind(boost::factory<BR_BKstarll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["ACP_Bsphimu"] = boost::bind(boost::factory<ACP*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P3CP_Bsphimu"] = boost::bind(boost::factory<P3CP*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["F_L_Bsphimu"] = boost::bind(boost::factory<F_L*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_1p_Bsphimu"] = boost::bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_2p_Bsphimu"] = boost::bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    //----- B to K ll  -----
    obsThFactory["BR_BKmu"] = boost::bind(boost::factory<BR_MPll_mu*>(), _1);
    obsThFactory["BR_BKe"] = boost::bind(boost::factory<BR_MPll_e*>(), _1);
    obsThFactory["R_BKll"] = boost::bind(boost::factory<R_MPll*>(), _1);
    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = boost::bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);

//
//    //-----  Lepton Flavour observables  -----
    obsThFactory["li_lj_gamma"] = boost::factory<li_lj_gamma*>();

    //-----  SM input parameters  -----
//    obsThFactory["Mz"] = boost::factory<StandardModelParams*>(myMO, "Mz");
//    obsThFactory["mHl"] = boost::factory<StandardModelParams*>(myMO, "mHl");

    //-----  NP input parameters, etc.   -----
//    if (myModel.ModelName().compare("NPEffective1") == 0
//            || myModel.ModelName().compare("NPEffective2") == 0) {
//        obsThFactory["cHQ1pPLUScHQ2p_NP"] = boost::factory<NewPhysicsParams*>(myMO, "cHQ1pPLUScHQ2p_NP");
//        obsThFactory["cHQ2pMINUScHQ2_NP"] = boost::factory<NewPhysicsParams*>(myMO, "cHQ2pMINUScHQ2_NP");
//        obsThFactory["cHQ3pPLUScHQ3_NP"] = boost::factory<NewPhysicsParams*>(myMO, "cHQ3pPLUScHQ3_NP");
//        obsThFactory["c_Ae_NP"] = boost::factory<NewPhysicsParams*>(myMO, "c_Ae_NP");
//        obsThFactory["c_GammaZ_uds_NP"] = boost::factory<NewPhysicsParams*>(myMO, "c_GammaZ_uds_NP");
//    }
//    if ((myModel.ModelName().compare("NPZbbbar") == 0)
//            || (myModel.ModelName().compare("NPZbbbarLR") == 0)) {
//        obsThFactory["deltaGVb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaGVb");
//        obsThFactory["deltaGAb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaGAb");
//        obsThFactory["deltaGLb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaGLb");
//        obsThFactory["deltaGRb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaGRb");
//        obsThFactory["deltaRhoZb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaRhoZb");
//        obsThFactory["deltaKappaZb"] = boost::factory<NewPhysicsParams*>(myMO, "deltaKappaZb");
//    }
//
//    //-----  SUSY spectra and observables  -----
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

}

void ThObsFactory::addObsToFactory (const std::string name, boost::function<ThObservable*(const StandardModel&) > funct)
{
    obsThFactory[name] = funct;
}

ThObservable * ThObsFactory::CreateThMethod(const std::string& name, const StandardModel& model) const
{
    if(model.isModelParam(name))
        return new ParamObs(model, name);
    if (obsThFactory.find(name) == obsThFactory.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThObsFactory");
    return (obsThFactory.at(name)(model));
}
