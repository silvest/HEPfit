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
#include "NP_couplings.h"
#include "OptimizedObservablesSMEFTd6.h"
#include "HiggsThObservables.h"
#include "DiBosonThObservables.h"
#include "FlavourObservables.h"
#include "NPSMEFT6dtopquark.h"
#include "MtMSbar.h"
#include "alpha_s.h"
#include "LeptonFlavourObservables.h"
#include "SUSYObservables.h"
#include "GeorgiMachacekObservables.h"
#include "LoopMediators.h"
#include "THDMObservables.h"
#include "LRSMObservables.h"
/* BEGIN: REMOVE FROM THE PACKAGE */
#include "GeneralTHDMObservables.h"
#include "THDMWObservables.h"
/* END: REMOVE FROM THE PACKAGE */
#include <boost/lexical_cast.hpp>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

ThObsFactory::ThObsFactory()
{
    //-----  Energies of different colliders  -----
    const double sqrt_s_LEP2_WWcos1 = 0.18266; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos2 = 0.18909; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos3 = 0.19838; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos4 = 0.20592; ///< the center-of-mass energy in TeV
    //
    const double cos1_LEP2_WW = -0.95;
    const double cos2_LEP2_WW = -0.76;
    const double cos3_LEP2_WW = -0.57;
    const double cos4_LEP2_WW = -0.38;
    const double cos5_LEP2_WW = -0.19;
    const double cos6_LEP2_WW = 0.;
    const double cos7_LEP2_WW = 0.19;
    const double cos8_LEP2_WW = 0.38;
    const double cos9_LEP2_WW = 0.57;
    const double cos10_LEP2_WW = 0.76;
    const double cos11_LEP2_WW = 0.95;
    //
    const double cos1_ee_WW = -1.0;
    const double cos2_ee_WW = -0.8;
    const double cos3_ee_WW = -0.6;
    const double cos4_ee_WW = -0.4;
    const double cos5_ee_WW = -0.2;
    const double cos6_ee_WW = 0.;
    const double cos7_ee_WW = 0.2;
    const double cos8_ee_WW = 0.4;
    const double cos9_ee_WW = 0.6;
    const double cos10_ee_WW = 0.8;
    const double cos11_ee_WW = 1.0;
    //
    const double sqrt_s_LEP2_161 = 0.1613; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_172 = 0.1721; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_183 = 0.1827; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_189 = 0.1886; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_192 = 0.1916; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_196 = 0.1955; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_200 = 0.1995; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_202 = 0.2016; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_205 = 0.2049; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_206 = 0.2059; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_207 = 0.2066; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_208 = 0.208; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_LHC7 = 7.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC8 = 8.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC13 = 13.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC14 = 14.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC27 = 27.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCC100 = 100.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_TeV = 1.96; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_leptcoll_125 = .125; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_161 = .161; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_240 = .240; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_250 = .250; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_350 = .350; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_365 = .365; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_380 = .380; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_500 = .500; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1000 = 1.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1400 = 1.4; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1500 = 1.5; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_3000 = 3.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_10000 = 10.0; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_LHeC_1_3 = 1.3; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHeC_1_8 = 1.8; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_FCCep_3_5 = 3.5; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCCep_5 = 5.0; ///< the center-of-mass energy in TeV
    // Polarizations at lepton colliders
    const double pol_0 = 0.0;
    const double pol_20 = 20.0;
    const double pol_30 = 30.0;
    const double pol_80 = 80.0;
    //
    //-----  StandardModel observables  -----
    obsThFactory["MtMSbar"] = boost::factory<MtMSbar*>();
    obsThFactory["alpha_s_LO"] = bind(boost::factory<alpha_s*>(), _1, LO);
    obsThFactory["alpha_s_NLO"] = bind(boost::factory<alpha_s*>(), _1, NLO);
    obsThFactory["alpha_s_FULLNLO"] = bind(boost::factory<alpha_s*>(), _1, FULLNLO);
    obsThFactory["alpha_s_NNLO"] = bind(boost::factory<alpha_s*>(), _1, NNLO);
    obsThFactory["alpha_s_FULLNNLO"] = bind(boost::factory<alpha_s*>(), _1, FULLNNLO);
    //-----  Electroweak precision observables  -----
    obsThFactory["alphaMz"] = boost::factory<AlphaEmMz*>();
    obsThFactory["Dalpha_5h_Mz"] = boost::factory<DAlpha5hadMz*>();
    obsThFactory["Mw"] = boost::factory<Mw*>();
    obsThFactory["GammaW"] = boost::factory<GammaW*>();
    obsThFactory["BrWlepton"] = boost::factory<BrWlepton*>();
    obsThFactory["BrWelectron"] = boost::factory<BrWelectron*>();
    obsThFactory["BrWmuon"] = boost::factory<BrWmuon*>();
    obsThFactory["BrWtau"] = boost::factory<BrWtau*>();
    obsThFactory["BrWhadrons"] = boost::factory<BrWhadrons*>();
    obsThFactory["RWc"] = boost::factory<RWc*>();
    obsThFactory["RW_mu_e"] = boost::factory<RWmue*>();
    obsThFactory["RW_tau_e"] = boost::factory<RWtaue*>();
    obsThFactory["RW_tau_mu"] = boost::factory<RWtaumu*>();
    obsThFactory["RZ_mu_e"] = boost::factory<RZmue*>();
    obsThFactory["RZ_tau_e"] = boost::factory<RZtaue*>();
    obsThFactory["RZ_tau_mu"] = boost::factory<RZtaumu*>();
    obsThFactory["GammaZ"] = boost::factory<GammaZ*>();
    obsThFactory["GammaZhad"] = boost::factory<GammaZhad*>();
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
    obsThFactory["Rneutrinos"] = boost::factory<Rneutrinos*>();
    obsThFactory["Rinv"] = boost::factory<Rinv*>();
    obsThFactory["Ruc"] = boost::factory<Ruc*>();
    obsThFactory["Rcharm"] = boost::factory<Rcharm*>();
    obsThFactory["Rbottom"] = boost::factory<Rbottom*>();
    //-----  Triple gauge coupling observables  -----
    obsThFactory["deltag1Z"] = boost::factory<deltag1Z*>();
    obsThFactory["deltaKgamma"] = boost::factory<deltaKgamma*>();
    obsThFactory["lambdaZ"] = boost::factory<lambdaZ*>();
    //-----  ee -> WW observables: LEP2 total cross section  -----
    obsThFactory["eeWW_LEP2_161"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_161);
    obsThFactory["eeWW_LEP2_172"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_172);
    obsThFactory["eeWW_LEP2_183"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_183);
    obsThFactory["eeWW_LEP2_189"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWW_LEP2_192"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWW_LEP2_196"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWW_LEP2_200"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWW_LEP2_202"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWW_LEP2_205"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWW_LEP2_206"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["eeWW_LEP2_207"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["eeWW_LEP2_208"] = bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_208);
    // Similar observables, defined only for the d 6 SMEFT, from arXiv: 1606.06693 [hep-ph].
    obsThFactory["eeWWlept_LEP2_189"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWWlept_LEP2_192"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWWlept_LEP2_196"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWWlept_LEP2_200"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWWlept_LEP2_202"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWWlept_LEP2_205"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWWlept_LEP2_206"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["eeWWlept_LEP2_207"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["eeWWlept_LEP2_208"] = bind(boost::factory<xseeWWlept*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["eeWWsemil_LEP2_189"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWWsemil_LEP2_192"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWWsemil_LEP2_196"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWWsemil_LEP2_200"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWWsemil_LEP2_202"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWWsemil_LEP2_205"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWWsemil_LEP2_206"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["eeWWsemil_LEP2_207"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["eeWWsemil_LEP2_208"] = bind(boost::factory<xseeWWsemil*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["eeWWhad_LEP2_189"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWWhad_LEP2_192"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWWhad_LEP2_196"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWWhad_LEP2_200"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWWhad_LEP2_202"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWWhad_LEP2_205"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWWhad_LEP2_206"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["eeWWhad_LEP2_207"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["eeWWhad_LEP2_208"] = bind(boost::factory<xseeWWhad*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["eeWWtot_LEP2_189"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWWtot_LEP2_192"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWWtot_LEP2_196"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWWtot_LEP2_200"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWWtot_LEP2_202"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWWtot_LEP2_205"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWWtot_LEP2_206"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["eeWWtot_LEP2_207"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["eeWWtot_LEP2_208"] = bind(boost::factory<xseeWWtot*>(), _1, sqrt_s_LEP2_208);    
    // The same, but only the new physics contribution
    obsThFactory["deltaeeWWlept_LEP2_189"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["deltaeeWWlept_LEP2_192"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["deltaeeWWlept_LEP2_196"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["deltaeeWWlept_LEP2_200"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["deltaeeWWlept_LEP2_202"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["deltaeeWWlept_LEP2_205"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["deltaeeWWlept_LEP2_206"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["deltaeeWWlept_LEP2_207"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["deltaeeWWlept_LEP2_208"] = bind(boost::factory<deltaxseeWWlept*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["deltaeeWWsemil_LEP2_189"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["deltaeeWWsemil_LEP2_192"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["deltaeeWWsemil_LEP2_196"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["deltaeeWWsemil_LEP2_200"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["deltaeeWWsemil_LEP2_202"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["deltaeeWWsemil_LEP2_205"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["deltaeeWWsemil_LEP2_206"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["deltaeeWWsemil_LEP2_207"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["deltaeeWWsemil_LEP2_208"] = bind(boost::factory<deltaxseeWWsemil*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["deltaeeWWhad_LEP2_189"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["deltaeeWWhad_LEP2_192"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["deltaeeWWhad_LEP2_196"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["deltaeeWWhad_LEP2_200"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["deltaeeWWhad_LEP2_202"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["deltaeeWWhad_LEP2_205"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["deltaeeWWhad_LEP2_206"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["deltaeeWWhad_LEP2_207"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["deltaeeWWhad_LEP2_208"] = bind(boost::factory<deltaxseeWWhad*>(), _1, sqrt_s_LEP2_208);
    //
    obsThFactory["deltaeeWWtot_LEP2_189"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["deltaeeWWtot_LEP2_192"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["deltaeeWWtot_LEP2_196"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["deltaeeWWtot_LEP2_200"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["deltaeeWWtot_LEP2_202"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["deltaeeWWtot_LEP2_205"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["deltaeeWWtot_LEP2_206"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_206);
    obsThFactory["deltaeeWWtot_LEP2_207"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_207);
    obsThFactory["deltaeeWWtot_LEP2_208"] = bind(boost::factory<deltaxseeWWtot*>(), _1, sqrt_s_LEP2_208);
    //-----  ee -> WW observables: LEP2 differential cross section  -----
    obsThFactory["deeWWdcos_LEP2_183_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_189_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_198_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_206_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos10_LEP2_WW, cos11_LEP2_WW);
    // Similar observables, defined only for the d 6 SMEFT, from arXiv: 1606.06693 [hep-ph].
    obsThFactory["deeWW_LEP2_183_Bin1"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 1);
    obsThFactory["deeWW_LEP2_183_Bin2"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 2);
    obsThFactory["deeWW_LEP2_183_Bin3"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 3);
    obsThFactory["deeWW_LEP2_183_Bin4"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 4);    
    //
    obsThFactory["deeWW_LEP2_206_Bin1"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 1);
    obsThFactory["deeWW_LEP2_206_Bin2"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 2);
    obsThFactory["deeWW_LEP2_206_Bin3"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 3);
    obsThFactory["deeWW_LEP2_206_Bin4"] = bind(boost::factory<dxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 4);
    // The same but only the NP contribution
    obsThFactory["deltadeeWW_LEP2_183_Bin1"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 1);
    obsThFactory["deltadeeWW_LEP2_183_Bin2"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 2);
    obsThFactory["deltadeeWW_LEP2_183_Bin3"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 3);
    obsThFactory["deltadeeWW_LEP2_183_Bin4"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_183, 4);    
    //
    obsThFactory["deltadeeWW_LEP2_206_Bin1"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 1);
    obsThFactory["deltadeeWW_LEP2_206_Bin2"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 2);
    obsThFactory["deltadeeWW_LEP2_206_Bin3"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 3);
    obsThFactory["deltadeeWW_LEP2_206_Bin4"] = bind(boost::factory<deltadxseeWWLEP2Bin*>(), _1, sqrt_s_LEP2_206, 4);
    //-----  ee -> WW observables: Future colliders differential cross section  -----
    obsThFactory["deeWWdcos_161_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_161_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_161_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_161_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_161_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_161_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_161_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_161_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_161_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_161_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_240_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_240_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_240_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_240_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_240_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_240_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_240_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_240_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_240_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_240_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_250_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_250_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_250_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_250_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_250_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_250_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_250_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_250_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_250_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_250_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_350_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_350_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_350_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_350_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_350_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_350_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_350_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_350_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_350_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_350_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_365_Bin1"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_365_Bin2"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_365_Bin3"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_365_Bin4"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_365_Bin5"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_365_Bin6"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_365_Bin7"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_365_Bin8"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_365_Bin9"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_365_Bin10"] = bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos10_ee_WW, cos11_ee_WW);
    //-----  ee -> WW observables: total rates (ratio with the SM)  -----
    obsThFactory["eeWW161"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_161);
    //
    obsThFactory["eeWW240"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeWW240_p80_m30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_240, pol_80, -pol_30);
    obsThFactory["eeWW240_m80_p30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_240, -pol_80, pol_30);
    //
    obsThFactory["eeWW250"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeWW250_p80_m30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeWW250_m80_p30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeWW250_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeWW250_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeWW350"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeWW350_p80_m30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeWW350_m80_p30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeWW350_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeWW350_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeWW365"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeWW365_p80_m30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeWW365_m80_p30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    //
    obsThFactory["eeWW380_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeWW380_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeWW500"] = bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeWW500_p80_m30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeWW500_m80_p30"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeWW500_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeWW500_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);    
    //
    obsThFactory["eeWW1000_p80_m20"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeWW1000_m80_p20"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    //
    obsThFactory["eeWW1500_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeWW1500_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeWW3000_p80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeWW3000_m80_0"] = bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //----- High Energy diboson observables at hadron colliders
    obsThFactory["ppZHprobe14"] = bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_LHC14);
    obsThFactory["ppZHprobe27"] = bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_LHC27);
    obsThFactory["ppZHprobe100"] = bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["mupTVppWZ_14_Bin1"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 100., 150.);
    obsThFactory["mupTVppWZ_14_Bin2"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 150., 220.);
    obsThFactory["mupTVppWZ_14_Bin3"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 220., 300.);
    obsThFactory["mupTVppWZ_14_Bin4"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 300., 500.);
    obsThFactory["mupTVppWZ_14_Bin5"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 500., 750.);
    obsThFactory["mupTVppWZ_14_Bin6"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 750., 1200.);
    //
    obsThFactory["mupTVppWZ_27_Bin1"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 150., 220.);
    obsThFactory["mupTVppWZ_27_Bin2"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 220., 300.);
    obsThFactory["mupTVppWZ_27_Bin3"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 300., 500.);
    obsThFactory["mupTVppWZ_27_Bin4"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 500., 750.);
    obsThFactory["mupTVppWZ_27_Bin5"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 750., 1200.);
    obsThFactory["mupTVppWZ_27_Bin6"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 1200., 1800.);
    //
    obsThFactory["mupTVppWZ_100_Bin1"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 220., 300.);
    obsThFactory["mupTVppWZ_100_Bin2"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 300., 500.);
    obsThFactory["mupTVppWZ_100_Bin3"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 500., 750.);
    obsThFactory["mupTVppWZ_100_Bin4"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 750., 1200.);
    obsThFactory["mupTVppWZ_100_Bin5"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 1200., 1800.);
    obsThFactory["mupTVppWZ_100_Bin6"] = bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 1800., 2400.);    
    //-----  Observables for particle couplings -----
    //-----  Zff couplings observables  ----------
    obsThFactory["deltagZveveL"] = boost::factory<deltagZveveL*>();
    obsThFactory["deltagZvmuvmuL"] = boost::factory<deltagZvmuvmuL*>();
    obsThFactory["deltagZvtavtaL"] = boost::factory<deltagZvtavtaL*>();
    obsThFactory["deltagZeeL"] = boost::factory<deltagZeeL*>();
    obsThFactory["deltagZeeR"] = boost::factory<deltagZeeR*>();
    obsThFactory["deltagZmumuL"] = boost::factory<deltagZmumuL*>();
    obsThFactory["deltagZmumuR"] = boost::factory<deltagZmumuR*>();
    obsThFactory["deltagZtataL"] = boost::factory<deltagZtataL*>();
    obsThFactory["deltagZtataR"] = boost::factory<deltagZtataR*>();
    obsThFactory["deltagZuuL"] = boost::factory<deltagZuuL*>();
    obsThFactory["deltagZuuR"] = boost::factory<deltagZuuR*>();
    obsThFactory["deltagZuuV"] = boost::factory<deltagZuuV*>();
    obsThFactory["deltagZuuA"] = boost::factory<deltagZuuA*>();
    obsThFactory["deltagZccL"] = boost::factory<deltagZccL*>();
    obsThFactory["deltagZccR"] = boost::factory<deltagZccR*>();
    obsThFactory["deltagZttL"] = boost::factory<deltagZttL*>();
    obsThFactory["deltagZttR"] = boost::factory<deltagZttR*>();
    obsThFactory["deltagZttV"] = boost::factory<deltagZttV*>();
    obsThFactory["deltagZttA"] = boost::factory<deltagZttA*>();
    obsThFactory["deltagZddL"] = boost::factory<deltagZddL*>();
    obsThFactory["deltagZddR"] = boost::factory<deltagZddR*>();
    obsThFactory["deltagZddV"] = boost::factory<deltagZddV*>();
    obsThFactory["deltagZddA"] = boost::factory<deltagZddA*>();
    obsThFactory["deltagZssL"] = boost::factory<deltagZssL*>();
    obsThFactory["deltagZssR"] = boost::factory<deltagZssR*>();
    obsThFactory["deltagZbbL"] = boost::factory<deltagZbbL*>();
    obsThFactory["deltagZbbR"] = boost::factory<deltagZbbR*>();
    //-----  Wff couplings observables  ----------
    obsThFactory["deltaUWeve"] = boost::factory<deltaUWeve*>();
    obsThFactory["deltaUWmuvmu"] = boost::factory<deltaUWmuvmu*>();
    obsThFactory["deltaUWtavta"] = boost::factory<deltaUWtavta*>();
    obsThFactory["deltaVudL"] = boost::factory<deltaVudL*>();
    obsThFactory["deltaVudR"] = boost::factory<deltaVudR*>();
    obsThFactory["deltaVcsL"] = boost::factory<deltaVcsL*>();
    obsThFactory["deltaVcsR"] = boost::factory<deltaVcsR*>();
    obsThFactory["deltaVtbL"] = boost::factory<deltaVtbL*>();
    obsThFactory["deltaVtbR"] = boost::factory<deltaVtbR*>();
    //-----  Hff couplings observables  ----------
    obsThFactory["gHmumueff"] = boost::factory<gHmumueff*>();
    obsThFactory["gHtataeff"] = boost::factory<gHtataeff*>();
    obsThFactory["gHcceff"] = boost::factory<gHcceff*>();
    obsThFactory["gHbbeff"] = boost::factory<gHbbeff*>();
    obsThFactory["deltagHee"] = boost::factory<deltagHee*>();
    obsThFactory["deltagHmumu"] = boost::factory<deltagHmumu*>();
    obsThFactory["deltagHtata"] = boost::factory<deltagHtata*>();
    obsThFactory["deltagHuu"] = boost::factory<deltagHuu*>();
    obsThFactory["deltagHcc"] = boost::factory<deltagHcc*>();
    obsThFactory["deltagHtt"] = boost::factory<deltagHtt*>();
    obsThFactory["deltagHdd"] = boost::factory<deltagHdd*>();
    obsThFactory["deltagHss"] = boost::factory<deltagHss*>();
    obsThFactory["deltagHbb"] = boost::factory<deltagHbb*>();
    //-----  HGG couplings observables  ----------
    obsThFactory["gHGGeff"] = boost::factory<gHGGeff*>();
    obsThFactory["deltagHGG"] = boost::factory<deltagHGG*>();
    //-----  HZZ couplings observables  ----------
    obsThFactory["gHZZeff"] = boost::factory<gHZZeff*>();
    obsThFactory["deltagHZZ"] = boost::factory<deltagHZZ*>();
    obsThFactory["gHZZ1"] = boost::factory<gHZZ1*>();
    obsThFactory["gHZZ2"] = boost::factory<gHZZ2*>();
    //-----  HAA couplings observables  ----------
    obsThFactory["gHAAeff"] = boost::factory<gHAAeff*>();
    obsThFactory["deltagHAA"] = boost::factory<deltagHAA*>();
    //-----  HZA couplings observables  ----------
    obsThFactory["gHZAeff"] = boost::factory<gHZAeff*>();
    obsThFactory["deltagHZA"] = boost::factory<deltagHZA*>();
    obsThFactory["gHZA2"] = boost::factory<gHZA2*>();
    //-----  HWW couplings observables  ----------
    obsThFactory["gHWWeff"] = boost::factory<gHWWeff*>();
    obsThFactory["deltagHWW"] = boost::factory<deltagHWW*>();
    obsThFactory["gHWW1"] = boost::factory<gHWW1*>();
    obsThFactory["gHWW2"] = boost::factory<gHWW2*>();
    //-----  HHH couplings observables  ----------
    obsThFactory["deltalHHH"] = boost::factory<deltalHHH*>();
    //-----  Other Higgs couplings observables  ----------
    obsThFactory["gHWZeff_Ratio"] = boost::factory<gHWZeff*>();
    obsThFactory["gHbWeff_Ratio"] = boost::factory<gHbWeff*>();
    obsThFactory["gHtaWeff_Ratio"] = boost::factory<gHtaWeff*>();
    //-----  VVV couplings observables  ----------
    obsThFactory["deltag1ZEff"] = boost::factory<deltag1ZEff*>();
    obsThFactory["deltaKgammaEff"] = boost::factory<deltaKgammaEff*>();
    //-----  Basic interactions of the so-called Higgs basis  ----------
    obsThFactory["deltayt_HB"] = boost::factory<deltaytHB*>();
    obsThFactory["deltayb_HB"] = boost::factory<deltaybHB*>();
    obsThFactory["deltaytau_HB"] = boost::factory<deltaytauHB*>();
    obsThFactory["deltayc_HB"] = boost::factory<deltaycHB*>();
    obsThFactory["deltaymu_HB"] = boost::factory<deltaymuHB*>();
    obsThFactory["deltacZ_HB"] = boost::factory<deltacZHB*>();
    obsThFactory["cZBox_HB"] = boost::factory<cZBoxHB*>();
    obsThFactory["cZZ_HB"] = boost::factory<cZZHB*>();
    obsThFactory["cZga_HB"] = boost::factory<cZgaHB*>();
    obsThFactory["cgaga_HB"] = boost::factory<cgagaHB*>();
    obsThFactory["cgg_HB"] = boost::factory<cggHB*>();
    obsThFactory["cggEff_HB"] = boost::factory<cggEffHB*>();
    obsThFactory["lambz_HB"] = boost::factory<lambzHB*>();
    //-----  Other useful observables to work with new physics  ----------
    //-----  Z couplings with leptons ---------
    obsThFactory["delgZeL"] = bind(boost::factory<delgZlL*>(), _1, StandardModel::ELECTRON);
    obsThFactory["delgZeR"] = bind(boost::factory<delgZlR*>(), _1, StandardModel::ELECTRON);
    obsThFactory["delgZmuL"] = bind(boost::factory<delgZlL*>(), _1, StandardModel::MU);
    obsThFactory["delgZmuR"] = bind(boost::factory<delgZlR*>(), _1, StandardModel::MU);
    obsThFactory["delgZtaL"] = bind(boost::factory<delgZlL*>(), _1, StandardModel::TAU);
    obsThFactory["delgZtaR"] = bind(boost::factory<delgZlR*>(), _1, StandardModel::TAU);
    //-----  Z couplings with up sector quarks ---------
    obsThFactory["delgZuL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::UP);
    obsThFactory["delgZuR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::UP);
    obsThFactory["delgZcL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::CHARM);
    obsThFactory["delgZcR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::CHARM);
    obsThFactory["delgZtL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::TOP);
    obsThFactory["delgZtR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::TOP);
    //-----  Z couplings with down sector quarks ---------
    obsThFactory["delgZdL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::DOWN);
    obsThFactory["delgZdR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::DOWN);
    obsThFactory["delgZsL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::STRANGE);
    obsThFactory["delgZsR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::STRANGE);
    obsThFactory["delgZbL"] = bind(boost::factory<delgZqL*>(), _1, StandardModel::BOTTOM);
    obsThFactory["delgZbR"] = bind(boost::factory<delgZqR*>(), _1, StandardModel::BOTTOM);
    obsThFactory["deltaMW"] = boost::factory<deltaMW*>();
    //-----  Oblique Parameters ---------    
    obsThFactory["oblSpar"] = boost::factory<oblS*>();
    obsThFactory["oblTpar"] = boost::factory<oblT*>();
    obsThFactory["oblWpar"] = boost::factory<oblW*>();
    obsThFactory["oblYpar"] = boost::factory<oblY*>();


    //-----  Combinations of Warsaw basis coefficients constrained by EWPO  ----------
    obsThFactory["CEWHL1_11"] = boost::factory<CEWHL111*>();
    obsThFactory["CEWHL1_22"] = boost::factory<CEWHL122*>();
    obsThFactory["CEWHL1_33"] = boost::factory<CEWHL133*>();
    obsThFactory["CEWHL3_11"] = boost::factory<CEWHL311*>();
    obsThFactory["CEWHL3_22"] = boost::factory<CEWHL322*>();
    obsThFactory["CEWHL3_33"] = boost::factory<CEWHL333*>();
    obsThFactory["CEWHQ1_11"] = boost::factory<CEWHQ111*>();
    obsThFactory["CEWHQ1_22"] = boost::factory<CEWHQ122*>();
    obsThFactory["CEWHQ1_33"] = boost::factory<CEWHQ133*>();
    obsThFactory["CEWHQ3_11"] = boost::factory<CEWHQ311*>();
    obsThFactory["CEWHQ3_22"] = boost::factory<CEWHQ322*>();
    obsThFactory["CEWHQ3_33"] = boost::factory<CEWHQ333*>();
    obsThFactory["CEWHQd_33"] = boost::factory<CEWHQd33*>();
    obsThFactory["CEWHe_11"] = boost::factory<CEWHe11*>();
    obsThFactory["CEWHe_22"] = boost::factory<CEWHe22*>();
    obsThFactory["CEWHe_33"] = boost::factory<CEWHe33*>();
    obsThFactory["CEWHu_11"] = boost::factory<CEWHu11*>();
    obsThFactory["CEWHu_22"] = boost::factory<CEWHu22*>();
    obsThFactory["CEWHu_33"] = boost::factory<CEWHu33*>();
    obsThFactory["CEWHd_11"] = boost::factory<CEWHd11*>();
    obsThFactory["CEWHd_22"] = boost::factory<CEWHd22*>();
    obsThFactory["CEWHd_33"] = boost::factory<CEWHd33*>();




    //-----  Auxiliary observables to work with new physics  ----------
    obsThFactory["AuxObsNP1"] = boost::factory<AuxObsNP1*>();
    obsThFactory["AuxObsNP2"] = boost::factory<AuxObsNP2*>();
    obsThFactory["AuxObsNP3"] = boost::factory<AuxObsNP3*>();
    obsThFactory["AuxObsNP4"] = boost::factory<AuxObsNP4*>();
    obsThFactory["AuxObsNP5"] = boost::factory<AuxObsNP5*>();
    obsThFactory["AuxObsNP6"] = boost::factory<AuxObsNP6*>();
    obsThFactory["AuxObsNP7"] = boost::factory<AuxObsNP7*>();
    obsThFactory["AuxObsNP8"] = boost::factory<AuxObsNP8*>();
    obsThFactory["AuxObsNP9"] = boost::factory<AuxObsNP9*>();
    obsThFactory["AuxObsNP10"] = boost::factory<AuxObsNP10*>();
    obsThFactory["AuxObsNP11"] = boost::factory<AuxObsNP11*>();
    obsThFactory["AuxObsNP12"] = boost::factory<AuxObsNP12*>();
    obsThFactory["AuxObsNP13"] = boost::factory<AuxObsNP13*>();
    obsThFactory["AuxObsNP14"] = boost::factory<AuxObsNP14*>();
    obsThFactory["AuxObsNP15"] = boost::factory<AuxObsNP15*>();
    obsThFactory["AuxObsNP16"] = boost::factory<AuxObsNP16*>();
    obsThFactory["AuxObsNP17"] = boost::factory<AuxObsNP17*>();
    obsThFactory["AuxObsNP18"] = boost::factory<AuxObsNP18*>();
    obsThFactory["AuxObsNP19"] = boost::factory<AuxObsNP19*>();
    obsThFactory["AuxObsNP20"] = boost::factory<AuxObsNP20*>();
    obsThFactory["AuxObsNP21"] = boost::factory<AuxObsNP21*>();
    obsThFactory["AuxObsNP22"] = boost::factory<AuxObsNP22*>();
    obsThFactory["AuxObsNP23"] = boost::factory<AuxObsNP23*>();
    obsThFactory["AuxObsNP24"] = boost::factory<AuxObsNP24*>();
    obsThFactory["AuxObsNP25"] = boost::factory<AuxObsNP25*>();
    obsThFactory["AuxObsNP26"] = boost::factory<AuxObsNP26*>();
    obsThFactory["AuxObsNP27"] = boost::factory<AuxObsNP27*>();
    obsThFactory["AuxObsNP28"] = boost::factory<AuxObsNP28*>();
    obsThFactory["AuxObsNP29"] = boost::factory<AuxObsNP29*>();
    obsThFactory["AuxObsNP30"] = boost::factory<AuxObsNP30*>();

    //-----  Higgs observables  ----------

    //-----  Production cross sections (ratios with SM)  ----------
    obsThFactory["ggH"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC8);
    obsThFactory["WH"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ZH"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH+ttH"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF+VH"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["tHq"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH7"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC7);
    obsThFactory["VBF7"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC7);
    obsThFactory["WH7"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ZH7"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC7);
    obsThFactory["VH7"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC7);
    obsThFactory["ttH7"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC7);
    obsThFactory["tHq7"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC7);
    obsThFactory["ggH8"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH+ttH8"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF8"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBF+VH8"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["VBFgamma8"] = bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH8"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["WH8"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ZH8"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH8"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["tHq8"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH13"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH+ttH13"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF13"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF+VH13"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBFgamma13"] = bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC13);
    obsThFactory["VH13"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["WH13"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ZH13"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC13);   
    obsThFactory["VHpT25013"] = bind(boost::factory<muVHpT250*>(), _1, sqrt_s_LHC13);
    obsThFactory["WHpT25013"] = bind(boost::factory<muWHpT250*>(), _1, sqrt_s_LHC13);
    obsThFactory["ZHpT25013"] = bind(boost::factory<muZHpT250*>(), _1, sqrt_s_LHC13);
    obsThFactory["ttH13"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["tHq13"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH14"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH+ttH14"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF14"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF+VH14"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBFgamma14"] = bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC14);
    obsThFactory["VH14"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["WH14"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ZH14"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ttH14"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["tHq14"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH27"] = bind(boost::factory<muggH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ggH+ttH27"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBF27"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBF+VH27"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBFgamma27"] = bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC27);
    obsThFactory["VH27"] = bind(boost::factory<muVH*>(), _1, sqrt_s_LHC27);
    obsThFactory["WH27"] = bind(boost::factory<muWH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ZH27"] = bind(boost::factory<muZH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ttH27"] = bind(boost::factory<muttH*>(), _1, sqrt_s_LHC27);
    obsThFactory["tHq27"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_LHC27);
    obsThFactory["ggH100"] = bind(boost::factory<muggH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH+ttH100"] = bind(boost::factory<muggHpttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF100"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF+VH100"] = bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBFgamma100"] = bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_FCC100);
    obsThFactory["VH100"] = bind(boost::factory<muVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["WH100"] = bind(boost::factory<muWH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ZH100"] = bind(boost::factory<muZH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ttH100"] = bind(boost::factory<muttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["tHq100"] = bind(boost::factory<mutHq*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH196"] = bind(boost::factory<muggH*>(), _1, sqrt_s_TeV);
    obsThFactory["VBF196"] = bind(boost::factory<muVBF*>(), _1, sqrt_s_TeV);
    obsThFactory["VH196"] = bind(boost::factory<muVH*>(), _1, sqrt_s_TeV);
    obsThFactory["ttH196"] = bind(boost::factory<muttH*>(), _1, sqrt_s_TeV);
    //
    obsThFactory["eeZH240"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZH250"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZH350"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZH365"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZH380"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZH500"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZH1000"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZH1400"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZH1500"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZH3000"] = bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mumuZH3000"] = bind(boost::factory<mummZH*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumuZH10000"] = bind(boost::factory<mummZH*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["eeZllH240"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZllH250"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZllH350"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZllH365"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZllH380"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZllH500"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZllH1000"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZllH1400"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZllH1500"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZllH3000"] = bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZqqH240"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZqqH250"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZqqH350"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZqqH365"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZqqH380"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZqqH500"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZqqH1000"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZqqH1400"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZqqH1500"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZqqH3000"] = bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZH250_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZH250_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZH250_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZH250_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZH350_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZH350_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZH350_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZH350_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZH365_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZH365_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZH365_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZH365_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZH380_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZH380_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZH380_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZH380_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZH500_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZH500_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZH500_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZH500_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZH1000_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZH1000_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZH1000_p80_m20"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeZH1000_m80_p20"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeZH1000_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZH1000_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZH1400_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZH1400_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZH1400_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZH1400_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZH1500_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZH1500_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZH1500_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZH1500_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZH3000_p80_m30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZH3000_m80_p30"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZH3000_p80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZH3000_m80_0"] = bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZllH250_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZllH250_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZllH250_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZllH250_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZllH350_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZllH350_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZllH350_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZllH350_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZllH365_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZllH365_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZllH365_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZllH365_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZllH380_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZllH380_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZllH380_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZllH380_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZllH500_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZllH500_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZllH500_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZllH500_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1000_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZllH1000_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZllH1000_p80_m20"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeZllH1000_m80_p20"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeZllH1000_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZllH1000_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1400_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZllH1400_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZllH1400_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZllH1400_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1500_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZllH1500_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZllH1500_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZllH1500_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZllH3000_p80_m30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZllH3000_m80_p30"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZllH3000_p80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZllH3000_m80_0"] = bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH250_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZqqH250_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZqqH250_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZqqH250_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH350_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZqqH350_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZqqH350_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZqqH350_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH365_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZqqH365_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZqqH365_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZqqH365_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH380_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZqqH380_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZqqH380_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZqqH380_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH500_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZqqH500_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZqqH500_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZqqH500_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1000_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZqqH1000_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZqqH1000_p80_m20"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeZqqH1000_m80_p20"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeZqqH1000_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZqqH1000_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1400_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZqqH1400_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZqqH1400_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZqqH1400_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1500_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZqqH1500_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZqqH1500_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZqqH1500_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH3000_p80_m30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZqqH3000_m80_p30"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZqqH3000_p80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZqqH3000_m80_0"] = bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["aPsk250_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["aPsk250_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    //
    obsThFactory["aPsk350_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["aPsk350_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    //
    obsThFactory["aPsk500_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["aPsk500_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    //
    obsThFactory["bPsk250_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["bPsk250_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    //
    obsThFactory["bPsk350_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["bPsk350_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);    
    //
    obsThFactory["bPsk500_p80_m30"] = bind(boost::factory<aPsk*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["bPsk500_m80_p30"] = bind(boost::factory<bPsk*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    //
    obsThFactory["eeWBF240"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeWBF250"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeWBF350"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeWBF365"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeWBF380"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeWBF500"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeWBF1000"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeWBF1400"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeWBF1500"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeWBF3000"] = bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeWBF250_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeWBF250_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeWBF250_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeWBF250_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeWBF350_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeWBF350_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeWBF350_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeWBF350_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeWBF365_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeWBF365_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeWBF365_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeWBF365_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeWBF380_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeWBF380_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeWBF380_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeWBF380_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeWBF500_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeWBF500_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeWBF500_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeWBF500_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1000_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeWBF1000_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeWBF1000_p80_m20"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeWBF1000_m80_p20"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeWBF1000_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeWBF1000_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1400_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeWBF1400_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeWBF1400_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeWBF1400_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1500_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeWBF1500_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeWBF1500_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeWBF1500_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeWBF3000_p80_m30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeWBF3000_m80_p30"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeWBF3000_p80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeWBF3000_m80_0"] = bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeHvv240"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeHvv250"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeHvv350"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeHvv365"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeHvv380"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeHvv500"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeHvv1000"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeHvv1400"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeHvv1500"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeHvv3000"] = bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mumuHvv3000"] = bind(boost::factory<mummHvv*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumuHvv10000"] = bind(boost::factory<mummHvv*>(), _1, sqrt_s_leptcoll_10000);
     //
    obsThFactory["eeHvv250_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeHvv250_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeHvv250_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeHvv250_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeHvv350_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeHvv350_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeHvv350_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeHvv350_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeHvv365_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeHvv365_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeHvv365_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeHvv365_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeHvv380_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeHvv380_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeHvv380_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeHvv380_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeHvv500_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeHvv500_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeHvv500_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeHvv500_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1000_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeHvv1000_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeHvv1000_p80_m20"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeHvv1000_m80_p20"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeHvv1000_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeHvv1000_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1400_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeHvv1400_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeHvv1400_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeHvv1400_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1500_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeHvv1500_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeHvv1500_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeHvv1500_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeHvv3000_p80_m30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeHvv3000_m80_p30"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeHvv3000_p80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeHvv3000_m80_0"] = bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZBF240"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZBF250"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZBF350"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZBF365"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZBF380"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZBF500"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZBF1000"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZBF1400"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZBF1500"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZBF3000"] = bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZBF250_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZBF250_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZBF250_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZBF250_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZBF350_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZBF350_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZBF350_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZBF350_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZBF365_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZBF365_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZBF365_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZBF365_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZBF380_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZBF380_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZBF380_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZBF380_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZBF500_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZBF500_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZBF500_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZBF500_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1000_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZBF1000_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZBF1000_p80_m20"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eeZBF1000_m80_p20"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eeZBF1000_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZBF1000_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1400_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZBF1400_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZBF1400_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZBF1400_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1500_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZBF1500_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZBF1500_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZBF1500_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZBF3000_p80_m30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZBF3000_m80_p30"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZBF3000_p80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZBF3000_m80_0"] = bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eettH500"] = bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eettH1000"] = bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eettH1400"] = bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eettH1500"] = bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eettH3000"] = bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mumuttH1500"] = bind(boost::factory<mummttH*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumuttH3000"] = bind(boost::factory<mummttH*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["eettH500_p80_m30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eettH500_m80_p30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eettH500_p80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eettH500_m80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eettH1000_p80_m30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eettH1000_m80_p30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eettH1000_p80_m20"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["eettH1000_m80_p20"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["eettH1000_p80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eettH1000_m80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eettH1400_p80_m30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eettH1400_m80_p30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eettH1400_p80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eettH1400_m80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eettH1500_p80_m30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eettH1500_m80_p30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eettH1500_p80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eettH1500_m80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eettH3000_p80_m30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eettH3000_m80_p30"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eettH3000_p80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eettH3000_m80_0"] = bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mumuH125"] = bind(boost::factory<mummH*>(), _1, sqrt_s_leptcoll_125);
    //
    obsThFactory["epWBF1300"] = bind(boost::factory<muepWBF*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBF1800"] = bind(boost::factory<muepWBF*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["epWBF3500"] = bind(boost::factory<muepWBF*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBF5000"] = bind(boost::factory<muepWBF*>(), _1, sqrt_s_FCCep_5);
    //
    obsThFactory["epZBF1300"] = bind(boost::factory<muepZBF*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBF1800"] = bind(boost::factory<muepZBF*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["epZBF3500"] = bind(boost::factory<muepZBF*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBF5000"] = bind(boost::factory<muepZBF*>(), _1, sqrt_s_FCCep_5);
    //-----  Decay width and Branching ratios (ratios with SM)  ----------
    obsThFactory["GammaHRatio"] = boost::factory<GammaHRatio*>();
    obsThFactory["BrHinvisible"] = boost::factory<BrHinvisible*>();
    obsThFactory["BrHinvisibleNP"] = boost::factory<BrHinvisibleNP*>();
    obsThFactory["BrHexotic"] = boost::factory<BrHexotic*>();
    obsThFactory["BrHvisRatio"] = boost::factory<BrHtovisRatio*>();
    obsThFactory["BrHtoinvRatio"] = boost::factory<BrHtoinvRatio*>();
    obsThFactory["BrHggRatio"] = boost::factory<BrHtoggRatio*>();
    obsThFactory["BrHWWRatio"] = boost::factory<BrHtoWWRatio*>();
    obsThFactory["BrHWW2l2vRatio"] = boost::factory<BrHtoWW2l2vRatio*>();
    obsThFactory["BrHZZRatio"] = boost::factory<BrHtoZZRatio*>();
    obsThFactory["BrHZZ4lRatio"] = boost::factory<BrHtoZZ4lRatio*>();
    obsThFactory["BrHZZ4eRatio"] = boost::factory<BrHtoZZ4eRatio*>();
    obsThFactory["BrHZZ2e2muRatio"] = boost::factory<BrHtoZZ2e2muRatio*>();
    obsThFactory["BrHZZ4muRatio"] = boost::factory<BrHtoZZ4muRatio*>();
    obsThFactory["BrHZgaRatio"] = boost::factory<BrHtoZgaRatio*>();
    obsThFactory["BrHZgallRatio"] = boost::factory<BrHtoZgallRatio*>();
    obsThFactory["BrHZgaeeRatio"] = boost::factory<BrHtoZgaeeRatio*>();
    obsThFactory["BrHZgamumuRatio"] = boost::factory<BrHtoZgamumuRatio*>();
    obsThFactory["BrHgagaRatio"] = boost::factory<BrHtogagaRatio*>();
    obsThFactory["BrHmumuRatio"] = boost::factory<BrHtomumuRatio*>();
    obsThFactory["BrHtautauRatio"] = boost::factory<BrHtotautauRatio*>();
    obsThFactory["BrHccRatio"] = boost::factory<BrHtoccRatio*>();
    obsThFactory["BrHbbRatio"] = boost::factory<BrHtobbRatio*>();
    // Dedicated 4 lepton decays
    obsThFactory["BrHto2l2vRatio"] = boost::factory<BrHto2l2vRatio*>();
    obsThFactory["BrHtoevmuvRatio"] = boost::factory<BrHtoevmuvRatio*>();
    obsThFactory["BrHto2e2vRatio"] = boost::factory<BrHto2e2vRatio*>();
    obsThFactory["BrHto2mu2vRatio"] = boost::factory<BrHto2mu2vRatio*>();
    obsThFactory["BrHto4lRatio"] = boost::factory<BrHto4lRatio*>();
    obsThFactory["BrHto4eRatio"] = boost::factory<BrHto4eRatio*>();
    obsThFactory["BrHto4muRatio"] = boost::factory<BrHto4muRatio*>();
    obsThFactory["BrHto2e2muRatio"] = boost::factory<BrHto2e2muRatio*>();
    // Other dedicated (semi-)leptonic 4 fermion decays
    obsThFactory["BrHtolljjRatio"] = boost::factory<BrHtolljjRatio*>();
    obsThFactory["BrHtolvjjRatio"] = boost::factory<BrHtolvjjRatio*>();
    obsThFactory["BrHtolv_lvorjjRatio"] = boost::factory<BrHtolv_lvorjjRatio*>();
    obsThFactory["BrHtoll_vvorjjRatio"] = boost::factory<BrHtoll_vvorjjRatio*>();
    //-----  Ratios of BR (ratios with SM)  ----------
    obsThFactory["BrHtogaga_over_mumu_Ratio"] = boost::factory<BrHtogaga_over_mumu_Ratio*>();
    obsThFactory["BrHtoZga_over_mumu_Ratio"] = boost::factory<BrHtoZga_over_mumu_Ratio*>();
    obsThFactory["BrHtoZmumuga_over_mumu_Ratio"] = boost::factory<BrHtoZmumuga_over_mumu_Ratio*>();
    obsThFactory["BrHtogaga_over_4l_Ratio"] = boost::factory<BrHtogaga_over_4l_Ratio*>(); 
    obsThFactory["BrHtobb_over_4l_Ratio"] = boost::factory<BrHtobb_over_4l_Ratio*>();
    obsThFactory["BrHto2l2v_over_4l_Ratio"] = boost::factory<BrHto2l2v_over_4l_Ratio*>();
    obsThFactory["BrHtotautau_over_4l_Ratio"] = boost::factory<BrHtotautau_over_4l_Ratio*>();
    obsThFactory["BrHtogaga_over_2e2mu_Ratio"] = boost::factory<BrHtogaga_over_2e2mu_Ratio*>();
    obsThFactory["BrHtoZga_over_4l_Ratio"] = boost::factory<BrHtoZga_over_4l_Ratio*>();
    obsThFactory["BrHtomumu_over_4l_Ratio"] = boost::factory<BrHtomumu_over_4l_Ratio*>();
    obsThFactory["BrHtomumu_over_4mu_Ratio"] = boost::factory<BrHtomumu_over_4mu_Ratio*>();
    obsThFactory["BrHto4l_over_gaga_Ratio"] = boost::factory<BrHto4l_over_gaga_Ratio*>();
    obsThFactory["BrHtoZga_over_gaga_Ratio"] = boost::factory<BrHtoZga_over_gaga_Ratio*>();
    obsThFactory["BrHtomumu_over_gaga_Ratio"] = boost::factory<BrHtomumu_over_gaga_Ratio*>();
    obsThFactory["BrHto2l2v_over_gaga_Ratio"] = boost::factory<BrHto2l2v_over_gaga_Ratio*>();
    //-----  Special observables --------
    obsThFactory["muttHZbb_boost100"] = bind(boost::factory<muttHZbbboost*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggHgagaInt14"] = bind(boost::factory<muggHgagaInt*>(), _1, sqrt_s_LHC14);
    //-----  Full Signal strengths per prod and decay: Hadron colliders  ----------
    obsThFactory["muggHgaga"] = bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHgaga"] = bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHgaga"] = bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHgaga"] = bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZZ"] = bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ"] = bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ"] = bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ"] = bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHWW"] = bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW"] = bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW"] = bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW"] = bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHtautau"] = bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHtautau"] = bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHtautau"] = bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHtautau"] = bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHbb"] = bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHbb"] = bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHbb"] = bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHbb"] = bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC13);
    // Indicating the energy explicit in the observable name
    obsThFactory["muggHgaga13"] = bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHgaga13"] = bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHgaga13"] = bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHgaga13"] = bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHgaga13"] = bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHgaga13"] = bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZga13"] = bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZga13"] = bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZga13"] = bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZga13"] = bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZga13"] = bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZga13"] = bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZZ13"] = bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ13"] = bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZZ13"] = bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZZ13"] = bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ13"] = bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ13"] = bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHZZ4l13"] = bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ4l13"] = bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZZ4l13"] = bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZZ4l13"] = bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ4l13"] = bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ4l13"] = bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHWW13"] = bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW13"] = bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHWW13"] = bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHWW13"] = bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW13"] = bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW13"] = bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHWW2l2v13"] = bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW2l2v13"] = bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHWW2l2v13"] = bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHWW2l2v13"] = bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW2l2v13"] = bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW2l2v13"] = bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muttHVV13"] = bind(boost::factory<muttHVV*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHmumu13"] = bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHmumu13"] = bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHmumu13"] = bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHmumu13"] = bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHmumu13"] = bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHmumu13"] = bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHtautau13"] = bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHtautau13"] = bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFpVHtautau13"] = bind(boost::factory<muVBFpVHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHtautau13"] = bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHtautau13"] = bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHtautau13"] = bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHtautau13"] = bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHbb13"] = bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHbb13"] = bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHbb13"] = bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHbb13"] = bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHbb13"] = bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHbb13"] = bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muVBFBRinv13"] = bind(boost::factory<muVBFBRinv*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHBRinv13"] = bind(boost::factory<muVHBRinv*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muVBFHinv13"] = bind(boost::factory<muVBFHinv*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHinv13"] = bind(boost::factory<muVHinv*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHgaga14"] = bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHgaga14"] = bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHgaga14"] = bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHgaga14"] = bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHgaga14"] = bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHgaga14"] = bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHZga14"] = bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZga14"] = bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZga14"] = bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZga14"] = bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZga14"] = bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZga14"] = bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHZZ14"] = bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZZ14"] = bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZZ14"] = bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZZ14"] = bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZZ14"] = bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZZ14"] = bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHZZ4l14"] = bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZZ4l14"] = bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZZ4l14"] = bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZZ4l14"] = bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZZ4l14"] = bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZZ4l14"] = bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHWW14"] = bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHWW14"] = bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHWW14"] = bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHWW14"] = bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHWW14"] = bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHWW14"] = bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHWW2l2v14"] = bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHWW2l2v14"] = bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHWW2l2v14"] = bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHWW2l2v14"] = bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHWW2l2v14"] = bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHWW2l2v14"] = bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHmumu14"] = bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHmumu14"] = bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHmumu14"] = bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHmumu14"] = bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHmumu14"] = bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHmumu14"] = bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHtautau14"] = bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHtautau14"] = bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHtautau14"] = bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHtautau14"] = bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHtautau14"] = bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHtautau14"] = bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHbb14"] = bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHbb14"] = bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHbb14"] = bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHbb14"] = bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHbb14"] = bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHbb14"] = bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muVBFBRinv14"] = bind(boost::factory<muVBFBRinv*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHBRinv14"] = bind(boost::factory<muVHBRinv*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muVBFHinv14"] = bind(boost::factory<muVBFHinv*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHinv14"] = bind(boost::factory<muVHinv*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHgaga27"] = bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHgaga27"] = bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHgaga27"] = bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHgaga27"] = bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHgaga27"] = bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHgaga27"] = bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHZga27"] = bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZga27"] = bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZga27"] = bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZga27"] = bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZga27"] = bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZga27"] = bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHZZ27"] = bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZZ27"] = bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZZ27"] = bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZZ27"] = bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZZ27"] = bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZZ27"] = bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHZZ4l27"] = bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZZ4l27"] = bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZZ4l27"] = bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZZ4l27"] = bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZZ4l27"] = bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZZ4l27"] = bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHWW27"] = bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHWW27"] = bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHWW27"] = bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHWW27"] = bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHWW27"] = bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHWW27"] = bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHWW2l2v27"] = bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHWW2l2v27"] = bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHWW2l2v27"] = bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHWW2l2v27"] = bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHWW2l2v27"] = bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHWW2l2v27"] = bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHmumu27"] = bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHmumu27"] = bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHmumu27"] = bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHmumu27"] = bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHmumu27"] = bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHmumu27"] = bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHtautau27"] = bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHtautau27"] = bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHtautau27"] = bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHtautau27"] = bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHtautau27"] = bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHtautau27"] = bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHbb27"] = bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHbb27"] = bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHbb27"] = bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHbb27"] = bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHbb27"] = bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHbb27"] = bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muVBFBRinv27"] = bind(boost::factory<muVBFBRinv*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHBRinv27"] = bind(boost::factory<muVHBRinv*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muVBFHinv27"] = bind(boost::factory<muVBFHinv*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHinv27"] = bind(boost::factory<muVHinv*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHgaga100"] = bind(boost::factory<muggHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHgaga100"] = bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHgaga100"] = bind(boost::factory<muZHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHgaga100"] = bind(boost::factory<muWHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHgaga100"] = bind(boost::factory<muVHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHgaga100"] = bind(boost::factory<muttHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZga100"] = bind(boost::factory<muggHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZgamumu100"] = bind(boost::factory<muggHZgamumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZga100"] = bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZga100"] = bind(boost::factory<muZHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZga100"] = bind(boost::factory<muWHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZga100"] = bind(boost::factory<muVHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZga100"] = bind(boost::factory<muttHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZZ100"] = bind(boost::factory<muggHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZZ100"] = bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZZ100"] = bind(boost::factory<muZHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZZ100"] = bind(boost::factory<muWHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZZ100"] = bind(boost::factory<muVHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZZ100"] = bind(boost::factory<muttHZZ*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHZZ4l100"] = bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZZ4mu100"] = bind(boost::factory<muggHZZ4mu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZZ4l100"] = bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZZ4l100"] = bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZZ4l100"] = bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZZ4l100"] = bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZZ4l100"] = bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHWW100"] = bind(boost::factory<muggHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHWW100"] = bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHWW100"] = bind(boost::factory<muZHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHWW100"] = bind(boost::factory<muWHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHWW100"] = bind(boost::factory<muVHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHWW100"] = bind(boost::factory<muttHWW*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHWW2l2v100"] = bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHWW2l2v100"] = bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHWW2l2v100"] = bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHWW2l2v100"] = bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHWW2l2v100"] = bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHWW2l2v100"] = bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHmumu100"] = bind(boost::factory<muggHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHmumu100"] = bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHmumu100"] = bind(boost::factory<muZHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHmumu100"] = bind(boost::factory<muWHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHmumu100"] = bind(boost::factory<muVHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHmumu100"] = bind(boost::factory<muttHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHtautau100"] = bind(boost::factory<muggHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHtautau100"] = bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHtautau100"] = bind(boost::factory<muZHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHtautau100"] = bind(boost::factory<muWHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHtautau100"] = bind(boost::factory<muVHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHtautau100"] = bind(boost::factory<muttHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHbb100"] = bind(boost::factory<muggHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHbb100"] = bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHbb100"] = bind(boost::factory<muZHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHbb100"] = bind(boost::factory<muWHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHbb100"] = bind(boost::factory<muVHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHbb100"] = bind(boost::factory<muttHbb*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muVBFBRinv100"] = bind(boost::factory<muVBFBRinv*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHBRinv100"] = bind(boost::factory<muVHBRinv*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muVBFHinv100"] = bind(boost::factory<muVBFHinv*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHinv100"] = bind(boost::factory<muVHinv*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muppHmumu8"] = bind(boost::factory<muppHmumu*>(), _1, sqrt_s_LHC8);
    obsThFactory["muppHmumu13"] = bind(boost::factory<muppHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muppHZga8"] = bind(boost::factory<muppHZga*>(), _1, sqrt_s_LHC8);
    obsThFactory["muppHZga13"] = bind(boost::factory<muppHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHH2ga2b14"] = bind(boost::factory<muggHH2ga2b*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHH2ga2b100"] = bind(boost::factory<muggHH2ga2b*>(), _1, sqrt_s_FCC100);
    //
    // Special version of the H signal strength at Hadron collider with separate theory uncertainty per prod x decay channel
    // Only for 14 and 27 TeV  
    //
    obsThFactory["muTHUggHgaga14"] = bind(boost::factory<muTHUggHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHgaga14"] = bind(boost::factory<muTHUVBFHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHgaga14"] = bind(boost::factory<muTHUZHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHgaga14"] = bind(boost::factory<muTHUWHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHgaga14"] = bind(boost::factory<muTHUVHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHgaga14"] = bind(boost::factory<muTHUttHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUggHZga14"] = bind(boost::factory<muTHUggHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHZga14"] = bind(boost::factory<muTHUVBFHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHZga14"] = bind(boost::factory<muTHUZHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHZga14"] = bind(boost::factory<muTHUWHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHZga14"] = bind(boost::factory<muTHUVHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHZga14"] = bind(boost::factory<muTHUttHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUggHZZ14"] = bind(boost::factory<muTHUggHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHZZ14"] = bind(boost::factory<muTHUVBFHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHZZ14"] = bind(boost::factory<muTHUZHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHZZ14"] = bind(boost::factory<muTHUWHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHZZ14"] = bind(boost::factory<muTHUVHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHZZ14"] = bind(boost::factory<muTHUttHZZ*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUggHZZ4l14"] = bind(boost::factory<muTHUggHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHZZ4l14"] = bind(boost::factory<muTHUVBFHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHZZ4l14"] = bind(boost::factory<muTHUZHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHZZ4l14"] = bind(boost::factory<muTHUWHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHZZ4l14"] = bind(boost::factory<muTHUVHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHZZ4l14"] = bind(boost::factory<muTHUttHZZ4l*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUggHWW14"] = bind(boost::factory<muTHUggHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHWW14"] = bind(boost::factory<muTHUVBFHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHWW14"] = bind(boost::factory<muTHUZHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHWW14"] = bind(boost::factory<muTHUWHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHWW14"] = bind(boost::factory<muTHUVHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHWW14"] = bind(boost::factory<muTHUttHWW*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUggHWW2l2v14"] = bind(boost::factory<muTHUggHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHWW2l2v14"] = bind(boost::factory<muTHUVBFHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHWW2l2v14"] = bind(boost::factory<muTHUZHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHWW2l2v14"] = bind(boost::factory<muTHUWHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHWW2l2v14"] = bind(boost::factory<muTHUVHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHWW2l2v14"] = bind(boost::factory<muTHUttHWW2l2v*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUggHmumu14"] = bind(boost::factory<muTHUggHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHmumu14"] = bind(boost::factory<muTHUVBFHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHmumu14"] = bind(boost::factory<muTHUZHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHmumu14"] = bind(boost::factory<muTHUWHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHmumu14"] = bind(boost::factory<muTHUVHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHmumu14"] = bind(boost::factory<muTHUttHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUggHtautau14"] = bind(boost::factory<muTHUggHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHtautau14"] = bind(boost::factory<muTHUVBFHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHtautau14"] = bind(boost::factory<muTHUZHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHtautau14"] = bind(boost::factory<muTHUWHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHtautau14"] = bind(boost::factory<muTHUVHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHtautau14"] = bind(boost::factory<muTHUttHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUggHbb14"] = bind(boost::factory<muTHUggHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVBFHbb14"] = bind(boost::factory<muTHUVBFHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUZHbb14"] = bind(boost::factory<muTHUZHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUWHbb14"] = bind(boost::factory<muTHUWHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHbb14"] = bind(boost::factory<muTHUVHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUttHbb14"] = bind(boost::factory<muTHUttHbb*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUVBFBRinv14"] = bind(boost::factory<muTHUVBFBRinv*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHBRinv14"] = bind(boost::factory<muTHUVHBRinv*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUVBFHinv14"] = bind(boost::factory<muTHUVBFHinv*>(), _1, sqrt_s_LHC14);
    obsThFactory["muTHUVHinv14"] = bind(boost::factory<muTHUVHinv*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muTHUggHgaga27"] = bind(boost::factory<muTHUggHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHgaga27"] = bind(boost::factory<muTHUVBFHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHgaga27"] = bind(boost::factory<muTHUZHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHgaga27"] = bind(boost::factory<muTHUWHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHgaga27"] = bind(boost::factory<muTHUVHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHgaga27"] = bind(boost::factory<muTHUttHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUggHZga27"] = bind(boost::factory<muTHUggHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHZga27"] = bind(boost::factory<muTHUVBFHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHZga27"] = bind(boost::factory<muTHUZHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHZga27"] = bind(boost::factory<muTHUWHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHZga27"] = bind(boost::factory<muTHUVHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHZga27"] = bind(boost::factory<muTHUttHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUggHZZ27"] = bind(boost::factory<muTHUggHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHZZ27"] = bind(boost::factory<muTHUVBFHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHZZ27"] = bind(boost::factory<muTHUZHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHZZ27"] = bind(boost::factory<muTHUWHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHZZ27"] = bind(boost::factory<muTHUVHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHZZ27"] = bind(boost::factory<muTHUttHZZ*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUggHZZ4l27"] = bind(boost::factory<muTHUggHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHZZ4l27"] = bind(boost::factory<muTHUVBFHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHZZ4l27"] = bind(boost::factory<muTHUZHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHZZ4l27"] = bind(boost::factory<muTHUWHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHZZ4l27"] = bind(boost::factory<muTHUVHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHZZ4l27"] = bind(boost::factory<muTHUttHZZ4l*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUggHWW27"] = bind(boost::factory<muTHUggHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHWW27"] = bind(boost::factory<muTHUVBFHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHWW27"] = bind(boost::factory<muTHUZHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHWW27"] = bind(boost::factory<muTHUWHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHWW27"] = bind(boost::factory<muTHUVHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHWW27"] = bind(boost::factory<muTHUttHWW*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUggHWW2l2v27"] = bind(boost::factory<muTHUggHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHWW2l2v27"] = bind(boost::factory<muTHUVBFHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHWW2l2v27"] = bind(boost::factory<muTHUZHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHWW2l2v27"] = bind(boost::factory<muTHUWHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHWW2l2v27"] = bind(boost::factory<muTHUVHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHWW2l2v27"] = bind(boost::factory<muTHUttHWW2l2v*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUggHmumu27"] = bind(boost::factory<muTHUggHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHmumu27"] = bind(boost::factory<muTHUVBFHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHmumu27"] = bind(boost::factory<muTHUZHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHmumu27"] = bind(boost::factory<muTHUWHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHmumu27"] = bind(boost::factory<muTHUVHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHmumu27"] = bind(boost::factory<muTHUttHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUggHtautau27"] = bind(boost::factory<muTHUggHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHtautau27"] = bind(boost::factory<muTHUVBFHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHtautau27"] = bind(boost::factory<muTHUZHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHtautau27"] = bind(boost::factory<muTHUWHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHtautau27"] = bind(boost::factory<muTHUVHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHtautau27"] = bind(boost::factory<muTHUttHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUggHbb27"] = bind(boost::factory<muTHUggHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVBFHbb27"] = bind(boost::factory<muTHUVBFHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUZHbb27"] = bind(boost::factory<muTHUZHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUWHbb27"] = bind(boost::factory<muTHUWHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHbb27"] = bind(boost::factory<muTHUVHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUttHbb27"] = bind(boost::factory<muTHUttHbb*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUVBFBRinv27"] = bind(boost::factory<muTHUVBFBRinv*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHBRinv27"] = bind(boost::factory<muTHUVHBRinv*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muTHUVBFHinv27"] = bind(boost::factory<muTHUVBFHinv*>(), _1, sqrt_s_LHC27);
    obsThFactory["muTHUVHinv27"] = bind(boost::factory<muTHUVHinv*>(), _1, sqrt_s_LHC27);
    //
    //
    //----- STXS bins at hadron colliders
    // Stage 0
    obsThFactory["STXS_0_qqH_13"] = bind(boost::factory<STXS_0_qqH*>(), _1, sqrt_s_LHC13);    
    // Stage 1: 4l final state
    obsThFactory["STXSggH_VBFtopo_j3v_4l_13"] = bind(boost::factory<STXSggH_VBFtopo_j3v_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSggH_VBFtopo_j3_4l_13"] = bind(boost::factory<STXSggH_VBFtopo_j3_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSggH0j4l_13"] = bind(boost::factory<STXSggH0j4l*>(), _1, sqrt_s_LHC13); 
    obsThFactory["STXSggH1j_pTH_0_60_4l_13"] = bind(boost::factory<STXSggH1j_pTH_0_60_4l*>(), _1, sqrt_s_LHC13); 
    obsThFactory["STXSggH1j_pTH_60_120_4l_13"] = bind(boost::factory<STXSggH1j_pTH_60_120_4l*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSggH1j_pTH_120_200_4l_13"] = bind(boost::factory<STXSggH1j_pTH_120_200_4l*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSggH1j_pTH_200_4l_13"] = bind(boost::factory<STXSggH1j_pTH_200_4l*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSggH2j_pTH_0_200_4l_13"] = bind(boost::factory<STXSggH2j_pTH_0_200_4l*>(), _1, sqrt_s_LHC13);  
    obsThFactory["STXSggH2j_pTH_0_60_4l_13"] = bind(boost::factory<STXSggH2j_pTH_0_60_4l*>(), _1, sqrt_s_LHC13);  
    obsThFactory["STXSggH2j_pTH_60_120_4l_13"] = bind(boost::factory<STXSggH2j_pTH_60_120_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSggH2j_pTH_120_200_4l_13"] = bind(boost::factory<STXSggH2j_pTH_120_200_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSggH2j_pTH_200_4l_13"] = bind(boost::factory<STXSggH2j_pTH_200_4l*>(), _1, sqrt_s_LHC13);    
    //
    obsThFactory["STXSqqHqq_VBFtopo_Rest_4l_13"] = bind(boost::factory<STXSqqHqq_VBFtopo_Rest_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHqq_VBFtopo_j3v_4l_13"] = bind(boost::factory<STXSqqHqq_VBFtopo_j3v_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHqq_VBFtopo_j3_4l_13"] = bind(boost::factory<STXSqqHqq_VBFtopo_j3_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHqq_nonVHtopo_4l_13"] = bind(boost::factory<STXSqqHqq_nonVHtopo_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHqq_VHtopo_4l_13"] = bind(boost::factory<STXSqqHqq_VHtopo_4l*>(), _1, sqrt_s_LHC13); 
    obsThFactory["STXSqqHqq_Rest_4l_13"] = bind(boost::factory<STXSqqHqq_Rest_4l*>(), _1, sqrt_s_LHC13);      
    obsThFactory["STXSqqHqq_pTj_200_4l_13"] = bind(boost::factory<STXSqqHqq_pTj_200_4l*>(), _1, sqrt_s_LHC13);    
    //
    obsThFactory["STXSqqHlv_pTV_0_250_4l_13"] = bind(boost::factory<STXSqqHlv_pTV_0_250_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_0_150_4l_13"] = bind(boost::factory<STXSqqHlv_pTV_0_150_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_150_250_0j_4l_13"] = bind(boost::factory<STXSqqHlv_pTV_150_250_0j_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_150_250_1j_4l_13"] = bind(boost::factory<STXSqqHlv_pTV_150_250_1j_4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_250_4l_13"] = bind(boost::factory<STXSqqHlv_pTV_250_4l*>(), _1, sqrt_s_LHC13);     
    //
    obsThFactory["STXSqqHll_pTV_0_150_4l_13"] = bind(boost::factory<STXSqqHll_pTV_0_150_4l*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSqqHll_pTV_150_250_4l_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_4l*>(), _1, sqrt_s_LHC13);    
    obsThFactory["STXSqqHll_pTV_150_250_0j_4l_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_0j_4l*>(), _1, sqrt_s_LHC13);    
    obsThFactory["STXSqqHll_pTV_150_250_1j_4l_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_1j_4l*>(), _1, sqrt_s_LHC13);         
    obsThFactory["STXSqqHll_pTV_250_4l_13"] = bind(boost::factory<STXSqqHll_pTV_250_4l*>(), _1, sqrt_s_LHC13);     
    //
    obsThFactory["STXSttHtH4l_13"] = bind(boost::factory<STXSttHtH4l*>(), _1, sqrt_s_LHC13);   
    // bb
    obsThFactory["STXSqqHlv_pTV_0_250_bb_13"] = bind(boost::factory<STXSqqHlv_pTV_0_250_bb*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_0_150_bb_13"] = bind(boost::factory<STXSqqHlv_pTV_0_150_bb*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_150_250_0j_bb_13"] = bind(boost::factory<STXSqqHlv_pTV_150_250_0j_bb*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_150_250_1j_bb_13"] = bind(boost::factory<STXSqqHlv_pTV_150_250_1j_bb*>(), _1, sqrt_s_LHC13);
    obsThFactory["STXSqqHlv_pTV_250_bb_13"] = bind(boost::factory<STXSqqHlv_pTV_250_bb*>(), _1, sqrt_s_LHC13);     
    //
    obsThFactory["STXSqqHll_pTV_0_150_bb_13"] = bind(boost::factory<STXSqqHll_pTV_0_150_bb*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSqqHll_pTV_150_250_bb_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_bb*>(), _1, sqrt_s_LHC13);    
    obsThFactory["STXSqqHll_pTV_150_250_0j_bb_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_0j_bb*>(), _1, sqrt_s_LHC13);    
    obsThFactory["STXSqqHll_pTV_150_250_1j_bb_13"] = bind(boost::factory<STXSqqHll_pTV_150_250_1j_bb*>(), _1, sqrt_s_LHC13);         
    obsThFactory["STXSqqHll_pTV_250_bb_13"] = bind(boost::factory<STXSqqHll_pTV_250_bb*>(), _1, sqrt_s_LHC13);    
    //
    obsThFactory["STXSWHqqHqq_VBFtopo_j3v_2b"] = bind(boost::factory<STXSWHqqHqq_VBFtopo_j3v_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSWHqqHqq_VBFtopo_j3_2b"] = bind(boost::factory<STXSWHqqHqq_VBFtopo_j3_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSWHqqHqq_VH2j_2b"] = bind(boost::factory<STXSWHqqHqq_VH2j_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSWHqqHqq_Rest_2b"] = bind(boost::factory<STXSWHqqHqq_Rest_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSWHqqHqq_pTj1_200_2b"] = bind(boost::factory<STXSWHqqHqq_pTj1_200_2b*>(), _1, sqrt_s_LHC13);     
    //
    obsThFactory["STXSZHqqHqq_VBFtopo_j3v_2b"] = bind(boost::factory<STXSZHqqHqq_VBFtopo_j3v_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSZHqqHqq_VBFtopo_j3_2b"] = bind(boost::factory<STXSZHqqHqq_VBFtopo_j3_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSZHqqHqq_VH2j_2b"] = bind(boost::factory<STXSZHqqHqq_VH2j_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSZHqqHqq_Rest_2b"] = bind(boost::factory<STXSZHqqHqq_Rest_2b*>(), _1, sqrt_s_LHC13);     
    obsThFactory["STXSZHqqHqq_pTj1_200_2b"] = bind(boost::factory<STXSZHqqHqq_pTj1_200_2b*>(), _1, sqrt_s_LHC13);     
    //
    // Stage 1.2: 4l final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_4l"] = bind(boost::factory<STXS12_ggH_pTH200_300_Nj01*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_4l"] = bind(boost::factory<STXS12_ggH_pTH300_450_Nj01*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_4l"] = bind(boost::factory<STXS12_ggH_pTH450_650_Nj01*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_4l"] = bind(boost::factory<STXS12_ggH_pTH650_Inf_Nj01*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_4l"] = bind(boost::factory<STXS12_ggH_pTH0_10_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_4l"] = bind(boost::factory<STXS12_ggH_pTH10_Inf_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_4l"] = bind(boost::factory<STXS12_ggH_pTH0_60_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_4l"] = bind(boost::factory<STXS12_ggH_pTH60_120_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_4l"] = bind(boost::factory<STXS12_ggH_pTH120_200_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH0_60_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH60_120_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH120_200_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_4l"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggHll_pTV0_75_4l"] = bind(boost::factory<STXS12_ggHll_pTV0_75*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggHll_pTV75_150_4l"] = bind(boost::factory<STXS12_ggHll_pTV75_150*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_4l"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_4l"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ggHll_pTV250_Inf_4l"] = bind(boost::factory<STXS12_ggHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_Nj0_4l"] = bind(boost::factory<STXS12_qqHqq_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_Nj1_4l"] = bind(boost::factory<STXS12_qqHqq_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj0_60_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj60_120_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj120_350_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_4l"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHlv_pTV0_75_4l"] = bind(boost::factory<STXS12_qqHlv_pTV0_75*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHlv_pTV75_150_4l"] = bind(boost::factory<STXS12_qqHlv_pTV75_150*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_4l"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_4l"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHlv_pTV250_Inf_4l"] = bind(boost::factory<STXS12_qqHlv_pTV250_Inf*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHll_pTV0_75_4l"] = bind(boost::factory<STXS12_qqHll_pTV0_75*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHll_pTV75_150_4l"] = bind(boost::factory<STXS12_qqHll_pTV75_150*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_4l"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_4l"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_qqHll_pTV250_Inf_4l"] = bind(boost::factory<STXS12_qqHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ttH_pTH0_60_4l"] = bind(boost::factory<STXS12_ttH_pTH0_60*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ttH_pTH60_120_4l"] = bind(boost::factory<STXS12_ttH_pTH60_120*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ttH_pTH120_200_4l"] = bind(boost::factory<STXS12_ttH_pTH120_200*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ttH_pTH200_300_4l"] = bind(boost::factory<STXS12_ttH_pTH200_300*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_ttH_pTH300_Inf_4l"] = bind(boost::factory<STXS12_ttH_pTH300_Inf*>(), _1, sqrt_s_LHC13, 1); 
    obsThFactory["STXS12_tH_4l"] = bind(boost::factory<STXS12_tH*>(), _1, sqrt_s_LHC13, 1); 
    //
    // Stage 1.2: ga ga final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_gaga"] = bind(boost::factory<STXS12_ggH_pTH200_300_Nj01*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_gaga"] = bind(boost::factory<STXS12_ggH_pTH300_450_Nj01*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_gaga"] = bind(boost::factory<STXS12_ggH_pTH450_650_Nj01*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_gaga"] = bind(boost::factory<STXS12_ggH_pTH650_Inf_Nj01*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_gaga"] = bind(boost::factory<STXS12_ggH_pTH0_10_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_gaga"] = bind(boost::factory<STXS12_ggH_pTH10_Inf_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_gaga"] = bind(boost::factory<STXS12_ggH_pTH0_60_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_gaga"] = bind(boost::factory<STXS12_ggH_pTH60_120_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_gaga"] = bind(boost::factory<STXS12_ggH_pTH120_200_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH0_60_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH60_120_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH120_200_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_gaga"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggHll_pTV0_75_gaga"] = bind(boost::factory<STXS12_ggHll_pTV0_75*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggHll_pTV75_150_gaga"] = bind(boost::factory<STXS12_ggHll_pTV75_150*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_gaga"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_gaga"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ggHll_pTV250_Inf_gaga"] = bind(boost::factory<STXS12_ggHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_Nj0_gaga"] = bind(boost::factory<STXS12_qqHqq_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_Nj1_gaga"] = bind(boost::factory<STXS12_qqHqq_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj0_60_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj60_120_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj120_350_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_gaga"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHlv_pTV0_75_gaga"] = bind(boost::factory<STXS12_qqHlv_pTV0_75*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHlv_pTV75_150_gaga"] = bind(boost::factory<STXS12_qqHlv_pTV75_150*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_gaga"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_gaga"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHlv_pTV250_Inf_gaga"] = bind(boost::factory<STXS12_qqHlv_pTV250_Inf*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHll_pTV0_75_gaga"] = bind(boost::factory<STXS12_qqHll_pTV0_75*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHll_pTV75_150_gaga"] = bind(boost::factory<STXS12_qqHll_pTV75_150*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_gaga"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_gaga"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_qqHll_pTV250_Inf_gaga"] = bind(boost::factory<STXS12_qqHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ttH_pTH0_60_gaga"] = bind(boost::factory<STXS12_ttH_pTH0_60*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ttH_pTH60_120_gaga"] = bind(boost::factory<STXS12_ttH_pTH60_120*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ttH_pTH120_200_gaga"] = bind(boost::factory<STXS12_ttH_pTH120_200*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ttH_pTH200_300_gaga"] = bind(boost::factory<STXS12_ttH_pTH200_300*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_ttH_pTH300_Inf_gaga"] = bind(boost::factory<STXS12_ttH_pTH300_Inf*>(), _1, sqrt_s_LHC13, 2); 
    obsThFactory["STXS12_tH_gaga"] = bind(boost::factory<STXS12_tH*>(), _1, sqrt_s_LHC13, 2); 
    //
    // Stage 1.2: bb final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_bb"] = bind(boost::factory<STXS12_ggH_pTH200_300_Nj01*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_bb"] = bind(boost::factory<STXS12_ggH_pTH300_450_Nj01*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_bb"] = bind(boost::factory<STXS12_ggH_pTH450_650_Nj01*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_bb"] = bind(boost::factory<STXS12_ggH_pTH650_Inf_Nj01*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_bb"] = bind(boost::factory<STXS12_ggH_pTH0_10_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_bb"] = bind(boost::factory<STXS12_ggH_pTH10_Inf_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_bb"] = bind(boost::factory<STXS12_ggH_pTH0_60_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_bb"] = bind(boost::factory<STXS12_ggH_pTH60_120_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_bb"] = bind(boost::factory<STXS12_ggH_pTH120_200_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH0_60_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH60_120_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH120_200_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_bb"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggHll_pTV0_75_bb"] = bind(boost::factory<STXS12_ggHll_pTV0_75*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggHll_pTV75_150_bb"] = bind(boost::factory<STXS12_ggHll_pTV75_150*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_bb"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_bb"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ggHll_pTV250_Inf_bb"] = bind(boost::factory<STXS12_ggHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_Nj0_bb"] = bind(boost::factory<STXS12_qqHqq_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_Nj1_bb"] = bind(boost::factory<STXS12_qqHqq_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj0_60_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj60_120_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj120_350_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_bb"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHlv_pTV0_75_bb"] = bind(boost::factory<STXS12_qqHlv_pTV0_75*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHlv_pTV75_150_bb"] = bind(boost::factory<STXS12_qqHlv_pTV75_150*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_bb"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_bb"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHlv_pTV250_Inf_bb"] = bind(boost::factory<STXS12_qqHlv_pTV250_Inf*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHll_pTV0_75_bb"] = bind(boost::factory<STXS12_qqHll_pTV0_75*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHll_pTV75_150_bb"] = bind(boost::factory<STXS12_qqHll_pTV75_150*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_bb"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_bb"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_qqHll_pTV250_Inf_bb"] = bind(boost::factory<STXS12_qqHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ttH_pTH0_60_bb"] = bind(boost::factory<STXS12_ttH_pTH0_60*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ttH_pTH60_120_bb"] = bind(boost::factory<STXS12_ttH_pTH60_120*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ttH_pTH120_200_bb"] = bind(boost::factory<STXS12_ttH_pTH120_200*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ttH_pTH200_300_bb"] = bind(boost::factory<STXS12_ttH_pTH200_300*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_ttH_pTH300_Inf_bb"] = bind(boost::factory<STXS12_ttH_pTH300_Inf*>(), _1, sqrt_s_LHC13, 3); 
    obsThFactory["STXS12_tH_bb"] = bind(boost::factory<STXS12_tH*>(), _1, sqrt_s_LHC13, 3); 
    //
    // Stage 1.2: evmuv final state
    obsThFactory["STXS12_ggH_pTH200_300_Nj01_evmuv"] = bind(boost::factory<STXS12_ggH_pTH200_300_Nj01*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH300_450_Nj01_evmuv"] = bind(boost::factory<STXS12_ggH_pTH300_450_Nj01*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH450_650_Nj01_evmuv"] = bind(boost::factory<STXS12_ggH_pTH450_650_Nj01*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH650_Inf_Nj01_evmuv"] = bind(boost::factory<STXS12_ggH_pTH650_Inf_Nj01*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH0_10_Nj0_evmuv"] = bind(boost::factory<STXS12_ggH_pTH0_10_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH10_Inf_Nj0_evmuv"] = bind(boost::factory<STXS12_ggH_pTH10_Inf_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH0_60_Nj1_evmuv"] = bind(boost::factory<STXS12_ggH_pTH0_60_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH60_120_Nj1_evmuv"] = bind(boost::factory<STXS12_ggH_pTH60_120_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_pTH120_200_Nj1_evmuv"] = bind(boost::factory<STXS12_ggH_pTH120_200_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH0_60_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH0_60_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH60_120_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH60_120_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj0_350_pTH120_200_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj0_350_pTH120_200_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2_evmuv"] = bind(boost::factory<STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggHll_pTV0_75_evmuv"] = bind(boost::factory<STXS12_ggHll_pTV0_75*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggHll_pTV75_150_evmuv"] = bind(boost::factory<STXS12_ggHll_pTV75_150*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj0_evmuv"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggHll_pTV150_250_Nj1_evmuv"] = bind(boost::factory<STXS12_ggHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ggHll_pTV250_Inf_evmuv"] = bind(boost::factory<STXS12_ggHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_Nj0_evmuv"] = bind(boost::factory<STXS12_qqHqq_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_Nj1_evmuv"] = bind(boost::factory<STXS12_qqHqq_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj0_60_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj0_60_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj60_120_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj60_120_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj120_350_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj120_350_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2_evmuv"] = bind(boost::factory<STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHlv_pTV0_75_evmuv"] = bind(boost::factory<STXS12_qqHlv_pTV0_75*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHlv_pTV75_150_evmuv"] = bind(boost::factory<STXS12_qqHlv_pTV75_150*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj0_evmuv"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHlv_pTV150_250_Nj1_evmuv"] = bind(boost::factory<STXS12_qqHlv_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHlv_pTV250_Inf_evmuv"] = bind(boost::factory<STXS12_qqHlv_pTV250_Inf*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHll_pTV0_75_evmuv"] = bind(boost::factory<STXS12_qqHll_pTV0_75*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHll_pTV75_150_evmuv"] = bind(boost::factory<STXS12_qqHll_pTV75_150*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj0_evmuv"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj0*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHll_pTV150_250_Nj1_evmuv"] = bind(boost::factory<STXS12_qqHll_pTV150_250_Nj1*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_qqHll_pTV250_Inf_evmuv"] = bind(boost::factory<STXS12_qqHll_pTV250_Inf*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ttH_pTH0_60_evmuv"] = bind(boost::factory<STXS12_ttH_pTH0_60*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ttH_pTH60_120_evmuv"] = bind(boost::factory<STXS12_ttH_pTH60_120*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ttH_pTH120_200_evmuv"] = bind(boost::factory<STXS12_ttH_pTH120_200*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ttH_pTH200_300_evmuv"] = bind(boost::factory<STXS12_ttH_pTH200_300*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_ttH_pTH300_Inf_evmuv"] = bind(boost::factory<STXS12_ttH_pTH300_Inf*>(), _1, sqrt_s_LHC13, 4); 
    obsThFactory["STXS12_tH_evmuv"] = bind(boost::factory<STXS12_tH*>(), _1, sqrt_s_LHC13, 4); 
    //
    //-----  Full Signal strengths per prod and decay: e+ e- colliders  ----------
    //
    // Pure WBF
    obsThFactory["mueeWBFbb240"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFbb250"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFbb350"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFbb365"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFbb380"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFbb500"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFbb1000"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFbb1400"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFbb1500"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFbb3000"] = bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFbb250_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeWBFbb250_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeWBFbb250_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeWBFbb250_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb350_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeWBFbb350_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeWBFbb350_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeWBFbb350_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb380_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeWBFbb380_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeWBFbb380_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeWBFbb380_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb500_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeWBFbb500_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeWBFbb500_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeWBFbb500_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb1000_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeWBFbb1000_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeWBFbb1000_p80_m20"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeWBFbb1000_m80_p20"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeWBFbb1000_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeWBFbb1000_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb1400_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeWBFbb1400_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeWBFbb1400_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeWBFbb1400_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb1500_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeWBFbb1500_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeWBFbb1500_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeWBFbb1500_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb3000_p80_m30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeWBFbb3000_m80_p30"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeWBFbb3000_p80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeWBFbb3000_m80_0"] = bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFcc240"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFcc250"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFcc350"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFcc365"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFcc380"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFcc500"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFcc1000"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFcc1400"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFcc1500"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFcc3000"] = bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFcc250_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc250_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc250_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFcc250_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc350_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc350_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc350_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFcc350_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc380_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc380_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc380_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFcc380_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc500_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc500_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc500_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFcc500_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1000_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1000_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1000_p80_m20"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFcc1000_m80_p20"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFcc1000_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1000_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1400_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1400_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1400_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1400_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc1500_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc1500_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc1500_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFcc1500_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFcc3000_p80_m30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFcc3000_m80_p30"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFcc3000_p80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFcc3000_m80_0"] = bind(boost::factory<mueeWBFccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFgg240"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFgg250"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFgg350"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFgg365"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFgg380"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFgg500"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFgg1000"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFgg1400"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFgg1500"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFgg3000"] = bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFgg250_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg250_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg250_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFgg250_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg350_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg350_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg350_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFgg350_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg380_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg380_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg380_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFgg380_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg500_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg500_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg500_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFgg500_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1000_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1000_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1000_p80_m20"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFgg1000_m80_p20"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFgg1000_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1000_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1400_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1400_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1400_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1400_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg1500_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg1500_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg1500_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFgg1500_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgg3000_p80_m30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgg3000_m80_p30"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgg3000_p80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFgg3000_m80_0"] = bind(boost::factory<mueeWBFggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFWW240"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFWW250"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFWW350"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFWW365"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFWW380"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFWW500"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFWW1000"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFWW1400"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFWW1500"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFWW3000"] = bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFWW250_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW250_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW250_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFWW250_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW350_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW350_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW350_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFWW350_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW380_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW380_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW380_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFWW380_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW500_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW500_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW500_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFWW500_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1000_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1000_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1000_p80_m20"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFWW1000_m80_p20"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFWW1000_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1000_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1400_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1400_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1400_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1400_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW1500_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW1500_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW1500_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFWW1500_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFWW3000_p80_m30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFWW3000_m80_p30"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFWW3000_p80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFWW3000_m80_0"] = bind(boost::factory<mueeWBFWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFtautau240"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFtautau250"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFtautau350"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFtautau365"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFtautau380"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFtautau500"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFtautau1000"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFtautau1400"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFtautau1500"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFtautau3000"] = bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFtautau250_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau250_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau250_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau250_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau350_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau350_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau350_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau350_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau380_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau380_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau380_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau380_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau500_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau500_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau500_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau500_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1000_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1000_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1000_p80_m20"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFtautau1000_m80_p20"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFtautau1000_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1000_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1400_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1400_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1400_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1400_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau1500_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau1500_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau1500_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau1500_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFtautau3000_p80_m30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFtautau3000_m80_p30"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFtautau3000_p80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFtautau3000_m80_0"] = bind(boost::factory<mueeWBFtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFZZ240"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFZZ250"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFZZ350"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFZZ365"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFZZ380"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFZZ500"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFZZ1000"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFZZ1400"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFZZ1500"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFZZ3000"] = bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFZZ250_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ250_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ250_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ250_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ350_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ350_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ350_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ350_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ380_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ380_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ380_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ380_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ500_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ500_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ500_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ500_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1000_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1000_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1000_p80_m20"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFZZ1000_m80_p20"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFZZ1000_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1000_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1400_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1400_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1400_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1400_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ1500_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ1500_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ1500_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ1500_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZZ3000_p80_m30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZZ3000_m80_p30"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZZ3000_p80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFZZ3000_m80_0"] = bind(boost::factory<mueeWBFZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFZga240"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFZga250"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFZga350"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFZga365"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFZga380"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFZga500"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFZga1000"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFZga1400"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFZga1500"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFZga3000"] = bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFZga250_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga250_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga250_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFZga250_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga350_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga350_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga350_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFZga350_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga380_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga380_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga380_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFZga380_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga500_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga500_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga500_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFZga500_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1000_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1000_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1000_p80_m20"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFZga1000_m80_p20"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFZga1000_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1000_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1400_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1400_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1400_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1400_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga1500_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga1500_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga1500_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFZga1500_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFZga3000_p80_m30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFZga3000_m80_p30"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFZga3000_p80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFZga3000_m80_0"] = bind(boost::factory<mueeWBFZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFgaga240"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFgaga250"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFgaga350"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFgaga365"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFgaga380"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFgaga500"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFgaga1000"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFgaga1400"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFgaga1500"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFgaga3000"] = bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFgaga250_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga250_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga250_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga250_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga350_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga350_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga350_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga350_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga380_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga380_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga380_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga380_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga500_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga500_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga500_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga500_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1000_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1000_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1000_p80_m20"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFgaga1000_m80_p20"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFgaga1000_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1000_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1400_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1400_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1400_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1400_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga1500_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga1500_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga1500_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga1500_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFgaga3000_p80_m30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFgaga3000_m80_p30"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFgaga3000_p80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFgaga3000_m80_0"] = bind(boost::factory<mueeWBFgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFmumu240"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFmumu250"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFmumu350"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFmumu365"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFmumu380"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFmumu500"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFmumu1000"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFmumu1400"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFmumu1500"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFmumu3000"] = bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
//    obsThFactory["mueeWBFmumu250_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu250_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu250_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu250_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu350_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu350_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu350_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu350_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu380_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu380_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu380_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu380_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu500_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu500_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu500_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu500_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1000_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1000_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1000_p80_m20"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
//    obsThFactory["mueeWBFmumu1000_m80_p20"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
//    obsThFactory["mueeWBFmumu1000_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1000_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1400_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1400_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1400_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1400_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu1500_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu1500_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu1500_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu1500_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
//    obsThFactory["mueeWBFmumu3000_p80_m30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
//    obsThFactory["mueeWBFmumu3000_m80_p30"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
//    obsThFactory["mueeWBFmumu3000_p80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
//    obsThFactory["mueeWBFmumu3000_m80_0"] = bind(boost::factory<mueeWBFmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // vvH
    obsThFactory["mueeHvvbb240"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvbb250"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvbb350"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvbb365"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvbb380"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvbb500"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvbb1000"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvbb1400"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvbb1500"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvbb3000"] = bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvbb250_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvbb250_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvbb250_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvbb250_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb350_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvbb350_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvbb350_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvbb350_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb380_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvbb380_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvbb380_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvbb380_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb500_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvbb500_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvbb500_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvbb500_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb1000_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvbb1000_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvbb1000_p80_m20"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvbb1000_m80_p20"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvbb1000_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvbb1000_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb1400_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvbb1400_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvbb1400_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvbb1400_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb1500_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvbb1500_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvbb1500_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvbb1500_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb3000_p80_m30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvbb3000_m80_p30"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvbb3000_p80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvbb3000_m80_0"] = bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc240"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvcc250"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvcc350"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvcc365"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvcc380"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvcc500"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvcc1000"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvcc1400"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvcc1500"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvcc3000"] = bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvcc250_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvcc250_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvcc250_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvcc250_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc350_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvcc350_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvcc350_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvcc350_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc380_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvcc380_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvcc380_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvcc380_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc500_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvcc500_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvcc500_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvcc500_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc1000_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvcc1000_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvcc1000_p80_m20"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvcc1000_m80_p20"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvcc1000_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvcc1000_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc1400_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvcc1400_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvcc1400_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvcc1400_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc1500_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvcc1500_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvcc1500_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvcc1500_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc3000_p80_m30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvcc3000_m80_p30"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvcc3000_p80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvcc3000_m80_0"] = bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg240"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvgg250"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvgg350"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvgg365"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvgg380"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvgg500"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvgg1000"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvgg1400"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvgg1500"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvgg3000"] = bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvgg250_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvgg250_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvgg250_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvgg250_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg350_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvgg350_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvgg350_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvgg350_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg380_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvgg380_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvgg380_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvgg380_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg500_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvgg500_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvgg500_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvgg500_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg1000_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvgg1000_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvgg1000_p80_m20"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvgg1000_m80_p20"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvgg1000_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvgg1000_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg1400_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvgg1400_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvgg1400_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvgg1400_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg1500_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvgg1500_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvgg1500_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvgg1500_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg3000_p80_m30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvgg3000_m80_p30"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvgg3000_p80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvgg3000_m80_0"] = bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW240"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvWW250"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvWW350"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvWW365"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvWW380"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvWW500"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvWW1000"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvWW1400"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvWW1500"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvWW3000"] = bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvWW250_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvWW250_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvWW250_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvWW250_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW350_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvWW350_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvWW350_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvWW350_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW380_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvWW380_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvWW380_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvWW380_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW500_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvWW500_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvWW500_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvWW500_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW1000_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvWW1000_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvWW1000_p80_m20"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvWW1000_m80_p20"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvWW1000_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvWW1000_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW1400_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvWW1400_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvWW1400_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvWW1400_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW1500_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvWW1500_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvWW1500_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvWW1500_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW3000_p80_m30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvWW3000_m80_p30"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvWW3000_p80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvWW3000_m80_0"] = bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau240"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvtautau250"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvtautau350"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvtautau365"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvtautau380"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvtautau500"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvtautau1000"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvtautau1400"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvtautau1500"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvtautau3000"] = bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvtautau250_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau250_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau250_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvtautau250_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau350_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau350_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau350_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvtautau350_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau380_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau380_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau380_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvtautau380_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau500_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau500_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau500_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvtautau500_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau1000_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau1000_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau1000_p80_m20"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvtautau1000_m80_p20"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvtautau1000_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvtautau1000_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau1400_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau1400_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau1400_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvtautau1400_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau1500_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau1500_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau1500_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvtautau1500_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau3000_p80_m30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau3000_m80_p30"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau3000_p80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvtautau3000_m80_0"] = bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ240"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvZZ250"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvZZ350"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvZZ365"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvZZ380"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvZZ500"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvZZ1000"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvZZ1400"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvZZ1500"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvZZ3000"] = bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvZZ250_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ250_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ250_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvZZ250_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ350_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ350_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ350_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvZZ350_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ380_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ380_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ380_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvZZ380_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ500_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ500_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ500_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvZZ500_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ1000_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ1000_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ1000_p80_m20"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvZZ1000_m80_p20"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvZZ1000_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvZZ1000_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ1400_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ1400_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ1400_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvZZ1400_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ1500_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ1500_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ1500_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvZZ1500_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ3000_p80_m30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ3000_m80_p30"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ3000_p80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvZZ3000_m80_0"] = bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga240"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvZga250"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvZga350"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvZga365"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvZga380"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvZga500"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvZga1000"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvZga1400"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvZga1500"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvZga3000"] = bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvZga250_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvZga250_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvZga250_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvZga250_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga350_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvZga350_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvZga350_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvZga350_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga380_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvZga380_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvZga380_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvZga380_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga500_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvZga500_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvZga500_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvZga500_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga1000_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvZga1000_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvZga1000_p80_m20"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvZga1000_m80_p20"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvZga1000_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvZga1000_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga1400_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvZga1400_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvZga1400_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvZga1400_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga1500_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvZga1500_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvZga1500_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvZga1500_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga3000_p80_m30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvZga3000_m80_p30"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvZga3000_p80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvZga3000_m80_0"] = bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga240"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvgaga250"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvgaga350"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvgaga365"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvgaga380"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvgaga500"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvgaga1000"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvgaga1400"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvgaga1500"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvgaga3000"] = bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvgaga250_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga250_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga250_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvgaga250_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga350_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga350_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga350_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvgaga350_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga380_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga380_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga380_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvgaga380_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga500_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga500_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga500_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvgaga500_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga1000_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga1000_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga1000_p80_m20"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvgaga1000_m80_p20"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvgaga1000_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvgaga1000_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga1400_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga1400_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga1400_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvgaga1400_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga1500_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga1500_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga1500_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvgaga1500_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga3000_p80_m30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga3000_m80_p30"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga3000_p80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvgaga3000_m80_0"] = bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu240"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvmumu250"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvmumu350"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvmumu365"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvmumu380"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvmumu500"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvmumu1000"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvmumu1400"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvmumu1500"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvmumu3000"] = bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvmumu250_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu250_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu250_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvmumu250_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu350_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu350_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu350_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvmumu350_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu380_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu380_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu380_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvmumu380_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu500_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu500_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu500_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvmumu500_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu1000_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu1000_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu1000_p80_m20"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeHvvmumu1000_m80_p20"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeHvvmumu1000_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeHvvmumu1000_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu1400_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu1400_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu1400_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvmumu1400_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu1500_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu1500_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu1500_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvmumu1500_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu3000_p80_m30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu3000_m80_p30"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu3000_p80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvmumu3000_m80_0"] = bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // ZH
    obsThFactory["mueeZHbb240"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHbb250"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHbb350"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHbb365"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHbb380"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHbb500"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHbb1000"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHbb1400"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHbb1500"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHbb3000"] = bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHbb250_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHbb250_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHbb250_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHbb250_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb350_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHbb350_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHbb350_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHbb350_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb380_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHbb380_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHbb380_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHbb380_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb500_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHbb500_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHbb500_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHbb500_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1000_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHbb1000_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHbb1000_p80_m20"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHbb1000_m80_p20"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHbb1000_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHbb1000_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1400_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHbb1400_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHbb1400_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHbb1400_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1500_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHbb1500_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHbb1500_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHbb1500_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb3000_p80_m30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHbb3000_m80_p30"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHbb3000_p80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHbb3000_m80_0"] = bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc240"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHcc250"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHcc350"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHcc365"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHcc380"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHcc500"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHcc1000"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHcc1400"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHcc1500"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHcc3000"] = bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHcc250_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHcc250_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHcc250_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHcc250_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc350_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHcc350_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHcc350_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHcc350_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc380_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHcc380_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHcc380_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHcc380_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc500_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHcc500_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHcc500_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHcc500_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc1000_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHcc1000_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHcc1000_p80_m20"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHcc1000_m80_p20"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHcc1000_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHcc1000_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc1400_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHcc1400_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHcc1400_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHcc1400_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc1500_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHcc1500_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHcc1500_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHcc1500_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc3000_p80_m30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHcc3000_m80_p30"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHcc3000_p80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHcc3000_m80_0"] = bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg240"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHgg250"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHgg350"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHgg365"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHgg380"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHgg500"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHgg1000"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHgg1400"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHgg1500"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHgg3000"] = bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHgg250_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHgg250_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHgg250_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHgg250_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg350_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHgg350_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHgg350_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHgg350_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg380_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHgg380_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHgg380_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHgg380_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg500_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHgg500_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHgg500_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHgg500_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg1000_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHgg1000_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHgg1000_p80_m20"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHgg1000_m80_p20"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHgg1000_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHgg1000_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg1400_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHgg1400_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHgg1400_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHgg1400_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg1500_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHgg1500_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHgg1500_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHgg1500_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg3000_p80_m30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHgg3000_m80_p30"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHgg3000_p80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHgg3000_m80_0"] = bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW240"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHWW250"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHWW350"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHWW365"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHWW380"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHWW500"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHWW1000"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHWW1400"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHWW1500"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHWW3000"] = bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHWW250_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHWW250_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHWW250_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHWW250_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW350_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHWW350_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHWW350_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHWW350_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW380_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHWW380_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHWW380_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHWW380_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW500_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHWW500_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHWW500_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHWW500_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW1000_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHWW1000_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHWW1000_p80_m20"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHWW1000_m80_p20"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHWW1000_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHWW1000_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW1400_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHWW1400_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHWW1400_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHWW1400_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW1500_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHWW1500_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHWW1500_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHWW1500_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW3000_p80_m30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHWW3000_m80_p30"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHWW3000_p80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHWW3000_m80_0"] = bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau240"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHtautau250"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHtautau350"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHtautau365"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHtautau380"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHtautau500"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHtautau1000"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHtautau1400"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHtautau1500"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHtautau3000"] = bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHtautau250_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHtautau250_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHtautau250_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHtautau250_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau350_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHtautau350_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHtautau350_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHtautau350_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau380_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHtautau380_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHtautau380_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHtautau380_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau500_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHtautau500_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHtautau500_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHtautau500_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau1000_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHtautau1000_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHtautau1000_p80_m20"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHtautau1000_m80_p20"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHtautau1000_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHtautau1000_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau1400_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHtautau1400_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHtautau1400_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHtautau1400_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau1500_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHtautau1500_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHtautau1500_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHtautau1500_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau3000_p80_m30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHtautau3000_m80_p30"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHtautau3000_p80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHtautau3000_m80_0"] = bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ240"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHZZ250"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHZZ350"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHZZ365"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHZZ380"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHZZ500"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHZZ1000"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHZZ1400"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHZZ1500"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHZZ3000"] = bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHZZ250_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHZZ250_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHZZ250_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHZZ250_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ350_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHZZ350_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHZZ350_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHZZ350_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ380_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHZZ380_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHZZ380_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHZZ380_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ500_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHZZ500_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHZZ500_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHZZ500_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ1000_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHZZ1000_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHZZ1000_p80_m20"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHZZ1000_m80_p20"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHZZ1000_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHZZ1000_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ1400_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHZZ1400_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHZZ1400_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHZZ1400_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ1500_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHZZ1500_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHZZ1500_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHZZ1500_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ3000_p80_m30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHZZ3000_m80_p30"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHZZ3000_p80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHZZ3000_m80_0"] = bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga240"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHZga250"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHZga350"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHZga365"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHZga380"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHZga500"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHZga1000"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHZga1400"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHZga1500"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHZga3000"] = bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHZga250_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHZga250_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHZga250_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHZga250_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga350_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHZga350_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHZga350_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHZga350_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga380_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHZga380_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHZga380_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHZga380_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga500_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHZga500_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHZga500_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHZga500_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga1000_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHZga1000_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHZga1000_p80_m20"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHZga1000_m80_p20"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHZga1000_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHZga1000_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga1400_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHZga1400_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHZga1400_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHZga1400_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga1500_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHZga1500_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHZga1500_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHZga1500_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga3000_p80_m30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHZga3000_m80_p30"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHZga3000_p80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHZga3000_m80_0"] = bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga240"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHgaga250"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHgaga350"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHgaga365"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHgaga380"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHgaga500"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHgaga1000"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHgaga1400"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHgaga1500"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHgaga3000"] = bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHgaga250_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHgaga250_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHgaga250_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHgaga250_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga350_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHgaga350_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHgaga350_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHgaga350_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga380_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHgaga380_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHgaga380_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHgaga380_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga500_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHgaga500_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHgaga500_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHgaga500_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga1000_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHgaga1000_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHgaga1000_p80_m20"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHgaga1000_m80_p20"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHgaga1000_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHgaga1000_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga1400_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHgaga1400_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHgaga1400_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHgaga1400_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga1500_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHgaga1500_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHgaga1500_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHgaga1500_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga3000_p80_m30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHgaga3000_m80_p30"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHgaga3000_p80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHgaga3000_m80_0"] = bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu240"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHmumu250"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHmumu350"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHmumu365"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHmumu380"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHmumu500"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHmumu1000"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHmumu1400"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHmumu1500"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHmumu3000"] = bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHmumu250_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHmumu250_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHmumu250_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHmumu250_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu350_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHmumu350_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHmumu350_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHmumu350_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu380_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHmumu380_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHmumu380_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHmumu380_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu500_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHmumu500_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHmumu500_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHmumu500_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu1000_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHmumu1000_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHmumu1000_p80_m20"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHmumu1000_m80_p20"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHmumu1000_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHmumu1000_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu1400_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHmumu1400_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHmumu1400_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHmumu1400_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu1500_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHmumu1500_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHmumu1500_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHmumu1500_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu3000_p80_m30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHmumu3000_m80_p30"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHmumu3000_p80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHmumu3000_m80_0"] = bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv240"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHinv250"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHinv350"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHinv365"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHinv380"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHinv500"] = bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_500);
    //
    obsThFactory["mueeZHinv250_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHinv250_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHinv250_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHinv250_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv350_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHinv350_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHinv350_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHinv350_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv380_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHinv380_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHinv380_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHinv380_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv500_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHinv500_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHinv500_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHinv500_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv1000_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHinv1000_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHinv1000_p80_m20"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHinv1000_m80_p20"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHinv1000_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHinv1000_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv1400_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHinv1400_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHinv1400_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHinv1400_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv1500_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHinv1500_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHinv1500_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHinv1500_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv3000_p80_m30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHinv3000_m80_p30"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHinv3000_p80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHinv3000_m80_0"] = bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv240"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHBRinv250"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHBRinv350"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHBRinv365"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHBRinv380"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHBRinv500"] = bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_500);
    //
    obsThFactory["mueeZHBRinv250_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv250_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv250_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHBRinv250_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv350_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv350_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv350_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHBRinv350_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv380_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv380_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv380_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHBRinv380_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv500_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv500_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv500_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHBRinv500_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv1000_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv1000_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv1000_p80_m20"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZHBRinv1000_m80_p20"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZHBRinv1000_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHBRinv1000_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv1400_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv1400_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv1400_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHBRinv1400_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv1500_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv1500_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv1500_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHBRinv1500_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv3000_p80_m30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv3000_m80_p30"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv3000_p80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHBRinv3000_m80_0"] = bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // ZBF
    obsThFactory["mueeZBFbb240"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZBFbb250"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZBFbb350"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZBFbb365"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZBFbb380"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZBFbb500"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZBFbb1000"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZBFbb1400"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZBFbb1500"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZBFbb3000"] = bind(boost::factory<mueeZBFbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZBFbb250_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZBFbb250_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZBFbb250_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZBFbb250_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb350_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZBFbb350_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZBFbb350_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZBFbb350_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb380_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZBFbb380_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZBFbb380_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZBFbb380_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //    
    obsThFactory["mueeZBFbb500_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZBFbb500_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZBFbb500_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZBFbb500_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //    
    obsThFactory["mueeZBFbb1000_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZBFbb1000_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZBFbb1000_p80_m20"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueeZBFbb1000_m80_p20"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueeZBFbb1000_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZBFbb1000_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //    
    obsThFactory["mueeZBFbb1400_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZBFbb1400_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZBFbb1400_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZBFbb1400_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb1500_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZBFbb1500_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZBFbb1500_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZBFbb1500_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb3000_p80_m30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZBFbb3000_m80_p30"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZBFbb3000_p80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZBFbb3000_m80_0"] = bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // eettH
    obsThFactory["mueettHbb500"] = bind(boost::factory<mueettHbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueettHbb1000"] = bind(boost::factory<mueettHbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueettHbb1400"] = bind(boost::factory<mueettHbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueettHbb1500"] = bind(boost::factory<mueettHbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueettHbb3000"] = bind(boost::factory<mueettHbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueettHbb500_p80_m30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueettHbb500_m80_p30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueettHbb500_p80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueettHbb500_m80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb1000_p80_m30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueettHbb1000_m80_p30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueettHbb1000_p80_m20"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_20);
    obsThFactory["mueettHbb1000_m80_p20"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_20);
    obsThFactory["mueettHbb1000_p80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueettHbb1000_m80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb1400_p80_m30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueettHbb1400_m80_p30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueettHbb1400_p80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueettHbb1400_m80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb1500_p80_m30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueettHbb1500_m80_p30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueettHbb1500_p80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueettHbb1500_m80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb3000_p80_m30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueettHbb3000_m80_p30"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueettHbb3000_p80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueettHbb3000_m80_0"] = bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    //-----  Full Signal strengths per prod and decay: mu+ mu- colliders  ----------
    //
    obsThFactory["mumumuHbb125"] = bind(boost::factory<mummHbb*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHcc125"] = bind(boost::factory<mummHcc*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHgg125"] = bind(boost::factory<mummHgg*>(), _1, sqrt_s_leptcoll_125);    
    obsThFactory["mumumuHWW125"] = bind(boost::factory<mummHWW*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHtautau125"] = bind(boost::factory<mummHtautau*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHZZ125"] = bind(boost::factory<mummHZZ*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHZga125"] = bind(boost::factory<mummHZga*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHgaga125"] = bind(boost::factory<mummHgaga*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHmumu125"] = bind(boost::factory<mummHmumu*>(), _1, sqrt_s_leptcoll_125);   
    // The same in the narrow width approximation
    obsThFactory["mumumuHbbNWA125"] = bind(boost::factory<mummHbbNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHccNWA125"] = bind(boost::factory<mummHccNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHggNWA125"] = bind(boost::factory<mummHggNWA*>(), _1, sqrt_s_leptcoll_125);    
    obsThFactory["mumumuHWWNWA125"] = bind(boost::factory<mummHWWNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHtautauNWA125"] = bind(boost::factory<mummHtautauNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHZZNWA125"] = bind(boost::factory<mummHZZNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHZgaNWA125"] = bind(boost::factory<mummHZgaNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHgagaNWA125"] = bind(boost::factory<mummHgagaNWA*>(), _1, sqrt_s_leptcoll_125);
    obsThFactory["mumumuHmumuNWA125"] = bind(boost::factory<mummHmumuNWA*>(), _1, sqrt_s_leptcoll_125);
    //
    // Signal strengths above the pole
    //
    // Hvv
    obsThFactory["mumumuHvvbb3000"] = bind(boost::factory<mummHvvbb*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvbb10000"] = bind(boost::factory<mummHvvbb*>(), _1, sqrt_s_leptcoll_10000);    
    //
    obsThFactory["mumumuHvvcc3000"] = bind(boost::factory<mummHvvcc*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvcc10000"] = bind(boost::factory<mummHvvcc*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvgg3000"] = bind(boost::factory<mummHvvgg*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvgg10000"] = bind(boost::factory<mummHvvgg*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvWW3000"] = bind(boost::factory<mummHvvWW*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvWW10000"] = bind(boost::factory<mummHvvWW*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvtautau3000"] = bind(boost::factory<mummHvvtautau*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvtautau10000"] = bind(boost::factory<mummHvvtautau*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvZZ3000"] = bind(boost::factory<mummHvvZZ*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvZZ10000"] = bind(boost::factory<mummHvvZZ*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvZga3000"] = bind(boost::factory<mummHvvZga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvZga10000"] = bind(boost::factory<mummHvvZga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvgaga3000"] = bind(boost::factory<mummHvvgaga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvgaga10000"] = bind(boost::factory<mummHvvgaga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHvvmumu3000"] = bind(boost::factory<mummHvvmumu*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHvvmumu10000"] = bind(boost::factory<mummHvvmumu*>(), _1, sqrt_s_leptcoll_10000);
    //
    // Hmumu
    obsThFactory["mumumuHmumubb3000"] = bind(boost::factory<mummHmmbb*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumubb10000"] = bind(boost::factory<mummHmmbb*>(), _1, sqrt_s_leptcoll_10000);    
    //
    obsThFactory["mumumuHmumucc3000"] = bind(boost::factory<mummHmmcc*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumucc10000"] = bind(boost::factory<mummHmmcc*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumugg3000"] = bind(boost::factory<mummHmmgg*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumugg10000"] = bind(boost::factory<mummHmmgg*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumuWW3000"] = bind(boost::factory<mummHmmWW*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumuWW10000"] = bind(boost::factory<mummHmmWW*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumutautau3000"] = bind(boost::factory<mummHmmtautau*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumutautau10000"] = bind(boost::factory<mummHmmtautau*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumuZZ3000"] = bind(boost::factory<mummHmmZZ*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumuZZ10000"] = bind(boost::factory<mummHmmZZ*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumuZga3000"] = bind(boost::factory<mummHmmZga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumuZga10000"] = bind(boost::factory<mummHmmZga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumugaga3000"] = bind(boost::factory<mummHmmgaga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumugaga10000"] = bind(boost::factory<mummHmmgaga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuHmumumumu3000"] = bind(boost::factory<mummHmmmumu*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuHmumumumu10000"] = bind(boost::factory<mummHmmmumu*>(), _1, sqrt_s_leptcoll_10000);
    //
    // ZH
    obsThFactory["mumumuZHbb3000"] = bind(boost::factory<mummZHbb*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHbb10000"] = bind(boost::factory<mummZHbb*>(), _1, sqrt_s_leptcoll_10000);    
    //
    obsThFactory["mumumuZHcc3000"] = bind(boost::factory<mummZHcc*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHcc10000"] = bind(boost::factory<mummZHcc*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHgg3000"] = bind(boost::factory<mummZHgg*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHgg10000"] = bind(boost::factory<mummZHgg*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHWW3000"] = bind(boost::factory<mummZHWW*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHWW10000"] = bind(boost::factory<mummZHWW*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHtautau3000"] = bind(boost::factory<mummZHtautau*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHtautau10000"] = bind(boost::factory<mummZHtautau*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHZZ3000"] = bind(boost::factory<mummZHZZ*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHZZ10000"] = bind(boost::factory<mummZHZZ*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHZga3000"] = bind(boost::factory<mummZHZga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHZga10000"] = bind(boost::factory<mummZHZga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHgaga3000"] = bind(boost::factory<mummZHgaga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHgaga10000"] = bind(boost::factory<mummZHgaga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuZHmumu3000"] = bind(boost::factory<mummZHmumu*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuZHmumu10000"] = bind(boost::factory<mummZHmumu*>(), _1, sqrt_s_leptcoll_10000);
    //
    // mumuttH 
    obsThFactory["mumumuttHbb3000"] = bind(boost::factory<mummttHbb*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHbb10000"] = bind(boost::factory<mummttHbb*>(), _1, sqrt_s_leptcoll_10000);    
    //
    obsThFactory["mumumuttHcc3000"] = bind(boost::factory<mummttHcc*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHcc10000"] = bind(boost::factory<mummttHcc*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHgg3000"] = bind(boost::factory<mummttHgg*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHgg10000"] = bind(boost::factory<mummttHgg*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHWW3000"] = bind(boost::factory<mummttHWW*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHWW10000"] = bind(boost::factory<mummttHWW*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHtautau3000"] = bind(boost::factory<mummttHtautau*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHtautau10000"] = bind(boost::factory<mummttHtautau*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHZZ3000"] = bind(boost::factory<mummttHZZ*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHZZ10000"] = bind(boost::factory<mummttHZZ*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHZga3000"] = bind(boost::factory<mummttHZga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHZga10000"] = bind(boost::factory<mummttHZga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHgaga3000"] = bind(boost::factory<mummttHgaga*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHgaga10000"] = bind(boost::factory<mummttHgaga*>(), _1, sqrt_s_leptcoll_10000);
    //
    obsThFactory["mumumuttHmumu3000"] = bind(boost::factory<mummttHmumu*>(), _1, sqrt_s_leptcoll_3000);
    obsThFactory["mumumuttHmumu10000"] = bind(boost::factory<mummttHmumu*>(), _1, sqrt_s_leptcoll_10000);    
    //
    //-----  Full Signal strengths per prod and decay: Lepton-Hadron colliders  ----------
    //
    obsThFactory["muepWBFbb1300"] = bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFcc1300"] = bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFgg1300"] = bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFWW2l2v1300"] = bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFZZ4l1300"] = bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFgaga1300"] = bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepWBFtautau1300"] = bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_LHeC_1_3);
    //
    obsThFactory["muepWBFbb1800"] = bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFcc1800"] = bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFgg1800"] = bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFWW2l2v1800"] = bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFZZ4l1800"] = bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFgaga1800"] = bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepWBFtautau1800"] = bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_LHeC_1_8);
    //
    obsThFactory["muepWBFbb3500"] = bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFcc3500"] = bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFgg3500"] = bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFWW2l2v3500"] = bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFZZ4l3500"] = bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFgaga3500"] = bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepWBFtautau3500"] = bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_FCCep_3_5);
    //
    obsThFactory["muepWBFbb5000"] = bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFcc5000"] = bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFgg5000"] = bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFWW2l2v5000"] = bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFZZ4l5000"] = bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFgaga5000"] = bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepWBFtautau5000"] = bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_FCCep_5);
    //
    obsThFactory["muepZBFbb1300"] = bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFcc1300"] = bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFgg1300"] = bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFWW2l2v1300"] = bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFZZ4l1300"] = bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFgaga1300"] = bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["muepZBFtautau1300"] = bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_LHeC_1_3);
    //
    obsThFactory["muepZBFbb1800"] = bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFcc1800"] = bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFgg1800"] = bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFWW2l2v1800"] = bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFZZ4l1800"] = bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFgaga1800"] = bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["muepZBFtautau1800"] = bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_LHeC_1_8);
    //
    obsThFactory["muepZBFbb3500"] = bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFcc3500"] = bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFgg3500"] = bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFWW2l2v3500"] = bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFZZ4l3500"] = bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFgaga3500"] = bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["muepZBFtautau3500"] = bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_FCCep_3_5);
    //
    obsThFactory["muepZBFbb5000"] = bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFcc5000"] = bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFgg5000"] = bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFWW2l2v5000"] = bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFZZ4l5000"] = bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFgaga5000"] = bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_FCCep_5);
    obsThFactory["muepZBFtautau5000"] = bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_FCCep_5);
    //
    //-----  Limits  ----------
    obsThFactory["UpperLimit_ppHZgammaA13"] = bind(boost::factory<UpperLimit_ppHZgammaA13*>(), _1, sqrt_s_LHC13);
    obsThFactory["UpperLimit_ppHZgammaC13"] = bind(boost::factory<UpperLimit_ppHZgammaC13*>(), _1, sqrt_s_LHC13);
    obsThFactory["UpperLimit_ppHZgammaA"] = bind(boost::factory<UpperLimit_ppHZgammaA*>(), _1, sqrt_s_LHC8);
    obsThFactory["UpperLimit_ppHZgammaC"] = bind(boost::factory<UpperLimit_ppHZgammaC*>(), _1, sqrt_s_LHC8);
    //-----  Others  ----------
    obsThFactory["cg_plus_ct"] = boost::factory<cg_plus_ct*>();
    obsThFactory["cga_plus_ct"] = boost::factory<cga_plus_ct*>();
    obsThFactory["cg_minus_cga"] = boost::factory<cg_minus_cga*>();
    obsThFactory["cV_plus_cb"] = boost::factory<cV_plus_cb*>();
    obsThFactory["cV_plus_ctau"] = boost::factory<cV_plus_ctau*>();
    obsThFactory["cb_minus_cc"] = boost::factory<cb_minus_cc*>();
    obsThFactory["cb_minus_ctau"] = boost::factory<cb_minus_ctau*>();
    obsThFactory["cc_minus_ctau"] = boost::factory<cc_minus_ctau*>();

    //-----  e+e- -> W+ W- Optimized Observables  -----
    obsThFactory["eeWWOO"] = boost::factory<eeWW*>();

    //-----  Epsilon parameters  -----
    obsThFactory["epsilon1"] = boost::factory<Epsilon1*>();
    obsThFactory["epsilon2"] = boost::factory<Epsilon2*>();
    obsThFactory["epsilon3"] = boost::factory<Epsilon3*>();
    obsThFactory["epsilonb"] = boost::factory<Epsilonb*>();
    
    
    //----- NPSMEFT6dtopquark  -----
    
    obsThFactory["C_phit"] = boost::factory<C_phit*>();
    obsThFactory["C_phiQ3"] = boost::factory<C_phiQ3*>();
    obsThFactory["C_phiQ1"] = boost::factory<C_phiQ1*>();
    obsThFactory["C_phiQm"] = boost::factory<C_phiQm*>();
    obsThFactory["C_tW"] = boost::factory<C_tW*>();
    obsThFactory["C_tZ"] = boost::factory<C_tZ*>();
    obsThFactory["C_tB"] = boost::factory<C_tB*>();
    obsThFactory["C_tphi"] = boost::factory<C_tphi*>();
    obsThFactory["C_phib"] = boost::factory<C_phib*>();
    obsThFactory["C_bW"] = boost::factory<C_bW*>();
    obsThFactory["C_bB"] = boost::factory<C_bB*>();
    obsThFactory["C_bZ"] = boost::factory<C_bZ*>();
    
    
    obsThFactory["C_ed"] = boost::factory<C_ed*>();
    obsThFactory["C_eq"] = boost::factory<C_eq*>();
    obsThFactory["C_ld"] = boost::factory<C_ld*>();
    obsThFactory["C_lqP"] = boost::factory<C_lqP*>();
    obsThFactory["C_eu"] = boost::factory<C_eu*>();
    obsThFactory["C_lu"] = boost::factory<C_lu*>();
    obsThFactory["C_lqM"] = boost::factory<C_lqM*>();
    
    
    obsThFactory["C_phitb"] = boost::factory<C_phitb*>();
    obsThFactory["C_tG"] = boost::factory<C_tG*>();
    obsThFactory["C_tu8"] = boost::factory<C_tu8*>();
    obsThFactory["C_td8"] = boost::factory<C_td8*>();
    obsThFactory["C_Qq18"] = boost::factory<C_Qq18*>();
    obsThFactory["C_tq8"] = boost::factory<C_tq8*>();
    obsThFactory["C_Qq38"] = boost::factory<C_Qq38*>();
    obsThFactory["C_Qu8"] = boost::factory<C_Qu8*>();
    obsThFactory["C_Qd8"] = boost::factory<C_Qd8*>();
   
    obsThFactory["Rb_NPSMEFT6dtopquark"] = boost::factory<Rb_NPSMEFT6dtopquark*>();
    obsThFactory["AFBLR"] = boost::factory<AFBLR*>();
    obsThFactory["SigmattZ"] = boost::factory<sigmattZ*>();
    obsThFactory["SigmattA"] = boost::factory<sigmattA*>();
    obsThFactory["SigmattH"] = boost::factory<sigmattH*>();
    obsThFactory["SigmattW"] = boost::factory<sigmattW*>();
    obsThFactory["SigmattWqEM"] = boost::factory<ttWqEM*>();
    obsThFactory["SigmattWqSUM"] = boost::factory<ttWqSUM*>();
    obsThFactory["Sigmatchannel13"] = boost::factory<sigmatchannel13*>();
    obsThFactory["Sigmatchannel8"] = boost::factory<sigmatchannel8*>();
    obsThFactory["SigmaschannelTev"] = boost::factory<sigmaschannelTev*>();
    obsThFactory["Sigmaschannel8"] = boost::factory<sigmaschannel8*>();
    obsThFactory["SigmatW"] = boost::factory<sigmatW*>();
     obsThFactory["SigmatW_8TeV"] = boost::factory<sigmatW_8TeV*>();
    obsThFactory["SigmatqZ"] = boost::factory<sigmatqZ*>();
    obsThFactory["SigmatAq"] = boost::factory<sigmatAq*>();
    obsThFactory["tH_Theo_Exp"] = boost::factory<tH_tchan*>(); 
    obsThFactory["ttHSUM"] = boost::factory<ttHSUM*>();
    obsThFactory["F0"] = boost::factory<F0*>();
    obsThFactory["FL"] = boost::factory<FL*>();
    
    obsThFactory["SigmattbarLHC13"] = boost::factory<sigmattbarLHC13*>();
    obsThFactory["SigmattbarLHC8"] = boost::factory<sigmattbarLHC8*>();
    obsThFactory["SigmattbarTev"] = boost::factory<sigmattbarTev*>();        
            
            
    obsThFactory["sigma_250_bb_eLpR"] = boost::factory<sigma_250_bb_eLpR*>();
    obsThFactory["a_250_bb_eLpR"] = boost::factory<a_250_bb_eLpR*>();
    obsThFactory["sigma_250_bb_eRpL"] = boost::factory<sigma_250_bb_eRpL*>();
    obsThFactory["a_250_bb_eRpL"] = boost::factory<a_250_bb_eRpL*>();
    obsThFactory["sigma_500_bb_eLpR"] = boost::factory<sigma_500_bb_eLpR*>();
    obsThFactory["a_500_bb_eLpR"] = boost::factory<a_500_bb_eLpR*>();
    obsThFactory["sigma_500_bb_eRpL"] = boost::factory<sigma_500_bb_eRpL*>();
    obsThFactory["a_500_bb_eRpL"] = boost::factory<a_500_bb_eRpL*>();
    obsThFactory["sigma_500_tt_eLpR"] = boost::factory<sigma_500_tt_eLpR*>();
    obsThFactory["a_500_tt_eLpR"] = boost::factory<a_500_tt_eLpR*>();
    obsThFactory["sigma_500_tt_eRpL"] = boost::factory<sigma_500_tt_eRpL*>();
    obsThFactory["a_500_tt_eRpL"] = boost::factory<a_500_tt_eRpL*>();
    obsThFactory["pt_500_tt_eLpR"] = boost::factory<pt_500_tt_eLpR*>();
    obsThFactory["pt_500_tt_eRpL"] = boost::factory<pt_500_tt_eRpL*>();
    
    obsThFactory["sigma_1000_bb_eLpR"] = boost::factory<sigma_1000_bb_eLpR*>();
    obsThFactory["a_1000_bb_eLpR"] = boost::factory<a_1000_bb_eLpR*>();
    obsThFactory["sigma_1000_bb_eRpL"] = boost::factory<sigma_1000_bb_eRpL*>();
    obsThFactory["a_1000_bb_eRpL"] = boost::factory<a_1000_bb_eRpL*>();   
  
    
    //SM ttZ bins
    
    obsThFactory["ttZ_bin_0_40"] = boost::factory<ttZ_bin_0_40*>();
    obsThFactory["ttZ_bin_40_70"] = boost::factory<ttZ_bin_40_70*>();
    obsThFactory["ttZ_bin_70_110"] = boost::factory<ttZ_bin_70_110*>();
    obsThFactory["ttZ_bin_110_160"] = boost::factory<ttZ_bin_110_160*>();
    obsThFactory["ttZ_bin_160_220"] = boost::factory<ttZ_bin_160_220*>();
    obsThFactory["ttZ_bin_220_290"] = boost::factory<ttZ_bin_220_290*>();
    obsThFactory["ttZ_bin_290_400"] = boost::factory<ttZ_bin_290_400*>();
    
   
    
    //ttA bins
    obsThFactory["ttA_bin_20_25"] = boost::factory<ttA_bin_20_25*>();
    obsThFactory["ttA_bin_25_30"] = boost::factory<ttA_bin_25_30*>();
    obsThFactory["ttA_bin_30_35"] = boost::factory<ttA_bin_30_35*>();
    obsThFactory["ttA_bin_35_40"] = boost::factory<ttA_bin_35_40*>();
    obsThFactory["ttA_bin_40_47"] = boost::factory<ttA_bin_40_47*>();
    obsThFactory["ttA_bin_47_55"] = boost::factory<ttA_bin_47_55*>();
    obsThFactory["ttA_bin_55_70"] = boost::factory<ttA_bin_55_70*>();
    obsThFactory["ttA_bin_70_85"] = boost::factory<ttA_bin_70_85*>();
    obsThFactory["ttA_bin_85_132"] = boost::factory<ttA_bin_85_132*>();
    obsThFactory["ttA_bin_132_180"] = boost::factory<ttA_bin_132_180*>();
    obsThFactory["ttA_bin_180_300"] = boost::factory<ttA_bin_180_300*>();
    
    //OPTIMIZED OBSERVABLES
    
    obsThFactory["op1"] = boost::factory<op1*>();
    obsThFactory["op2"] = boost::factory<op2*>();
    obsThFactory["op3"] = boost::factory<op3*>();
    obsThFactory["op4"] = boost::factory<op4*>();
    
    //OPTIMIZED OBSERVABLES 1000 GeV
    
    obsThFactory["op_1000_1"] = boost::factory<op_1000_1*>();
    obsThFactory["op_1000_2"] = boost::factory<op_1000_2*>();
    obsThFactory["op_1000_3"] = boost::factory<op_1000_3*>();
    obsThFactory["op_1000_4"] = boost::factory<op_1000_4*>();
    obsThFactory["op_1000_5"] = boost::factory<op_1000_5*>();
    obsThFactory["op_1000_6"] = boost::factory<op_1000_6*>();
    obsThFactory["op_1000_7"] = boost::factory<op_1000_7*>();
    obsThFactory["op_1000_8"] = boost::factory<op_1000_8*>();

   obsThFactory["gLt"] = boost::factory<gLt*>();
   obsThFactory["gLb"] = boost::factory<gLb*>();
   obsThFactory["gRt"] = boost::factory<gRt*>();
   obsThFactory["gRb"] = boost::factory<gRb*>();



    /* BEGIN: REMOVE FROM THE PACKAGE */
    //-----  LEP-II two-fermion processes  -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189.,
        192., 196., 200., 202., 205., 207.};
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192.,
        196., 200., 202., 205., 207.};
    for (int i = 0; i < 12; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
        obsThFactory["sigmaqLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2sigmaHadron*>(), _1, sqrt_s[i]);
        obsThFactory["sigmamuLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2sigmaMu*>(), _1, sqrt_s[i]);
        obsThFactory["sigmatauLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2sigmaTau*>(), _1, sqrt_s[i]);
        obsThFactory["AFBmuLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2AFBmu*>(), _1, sqrt_s[i]);
        obsThFactory["AFBtauLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2AFBtau*>(), _1, sqrt_s[i]);
    }
    for (int i = 0; i < 10; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
        obsThFactory["AFBbottomLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2AFBbottom*>(), _1, sqrt_s_HF[i]);
        obsThFactory["AFBcharmLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2AFBcharm*>(), _1, sqrt_s_HF[i]);
        obsThFactory["RbottomLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2Rbottom*>(), _1, sqrt_s_HF[i]);
        obsThFactory["RcharmLEP2_" + sqrt_s_str] = bind(boost::factory<LEP2Rcharm*>(), _1, sqrt_s_HF[i]);
    }
    /* END: REMOVE FROM THE PACKAGE */

    //-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["DmBd"] = boost::factory<DmBd*>();
    obsThFactory["DmBs"] = boost::factory<DmBs*>();
    obsThFactory["RmBs"] = boost::factory<RmBs*>();
    obsThFactory["alpha"] = boost::factory<Alpha*>();
    obsThFactory["alpha_2a"] = boost::factory<Alpha_2a*>();
    obsThFactory["SJPsiK"] = boost::factory<SJPsiK*>();
    obsThFactory["C2beta"] = boost::factory<C2beta*>();
    obsThFactory["Phis_JPsiPhi"] = boost::factory<Phis_JPsiPhi*>();
    obsThFactory["EpsilonK"] = boost::factory<EpsilonK*>();
    obsThFactory["DmK"] = boost::factory<DmK*>();
    obsThFactory["ImADC2"] = boost::factory<ImADC2*>();
    /* BEGIN: REMOVE FROM THE PACKAGE */
    obsThFactory["M12D"] = boost::factory<M12D*>();
    obsThFactory["ArgD"] = boost::factory<ArgD*>();
    //----- eps'/eps  -----
    obsThFactory["EpsilonP_O_Epsilon"] = boost::factory<EpsilonP_O_Epsilon*>();
    /* END: REMOVE FROM THE PACKAGE */
    //----- CKM  -----
    obsThFactory["Vud"] = bind(boost::factory<VCKM*>(), _1, 1, 1);
    obsThFactory["Vus"] = bind(boost::factory<VCKM*>(), _1, 1, 2);
    obsThFactory["Vub"] = bind(boost::factory<VCKM*>(), _1, 1, 3);
    obsThFactory["Vcd"] = bind(boost::factory<VCKM*>(), _1, 2, 1);
    obsThFactory["Vcs"] = bind(boost::factory<VCKM*>(), _1, 2, 2);
    obsThFactory["Vcb"] = bind(boost::factory<VCKM*>(), _1, 2, 3);
    obsThFactory["Vtd"] = bind(boost::factory<VCKM*>(), _1, 3, 1);
    obsThFactory["Vts"] = bind(boost::factory<VCKM*>(), _1, 3, 2);
    obsThFactory["Vtb"] = bind(boost::factory<VCKM*>(), _1, 3, 3);
    obsThFactory["CKM_alpha"] = boost::factory<CKM_Alpha*>();
    obsThFactory["CKM_gamma"] = boost::factory<CKM_Gamma*>();
    obsThFactory["CKM_beta"] = boost::factory<CKM_Beta*>();
    obsThFactory["CKM_betas"] = boost::factory<CKM_Betas*>();
    obsThFactory["CKM_2betapgamma"] = boost::factory<CKM_2BpG*>();
    obsThFactory["CKM_s2beta"] = boost::factory<CKM_S2Beta*>();
    obsThFactory["CKM_c2beta"] = boost::factory<CKM_C2Beta*>();
    obsThFactory["CKM_rho"] = boost::factory<CKM_rho*>();
    obsThFactory["CKM_eta"] = boost::factory<CKM_eta*>();
    obsThFactory["CKM_sintheta12"] = boost::factory<CKM_SinTheta12*>();
    obsThFactory["CKM_sintheta13"] = boost::factory<CKM_SinTheta13*>();
    obsThFactory["CKM_sintheta23"] = boost::factory<CKM_SinTheta23*>();
    obsThFactory["CKM_delta"] = boost::factory<CKM_Delta*>();
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
    obsThFactory["BR_Bdmumu"] = bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_D, StandardModel::MU);
    obsThFactory["BRbar_Bdmumu"] = bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_D, StandardModel::MU);
    obsThFactory["Amumu_Bd"] = bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_D, StandardModel::MU);
    obsThFactory["Smumu_Bd"] = bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_D, StandardModel::MU);

    obsThFactory["BR_Bsmumu"] = bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_S, StandardModel::MU);
    obsThFactory["BRbar_Bsmumu"] = bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_S, StandardModel::MU);
    obsThFactory["Amumu_Bs"] = bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_S, StandardModel::MU);
    obsThFactory["Smumu_Bs"] = bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_S, StandardModel::MU);
    obsThFactory["BR_Bsee"] = bind(boost::factory<Mll*>(), _1, 1, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["BRbar_Bsee"] = bind(boost::factory<Mll*>(), _1, 2, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["Aee_Bs"] = bind(boost::factory<Mll*>(), _1, 3, StandardModel::B_S, StandardModel::ELECTRON);
    obsThFactory["See_Bs"] = bind(boost::factory<Mll*>(), _1, 4, StandardModel::B_S, StandardModel::ELECTRON);

    obsThFactory["BR_BdmumuOBR_Bsmumu"] = boost::factory<BdmumuOBsmumu*>();
   //----- b to q gamma  -----
    obsThFactory["BR_bsgamma"] = bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 1);
    obsThFactory["ACP_bsgamma"] = bind(boost::factory<Bsgamma*>(), _1, StandardModel::STRANGE, 2);
    obsThFactory["BR_bdgamma"] = bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 1);
    obsThFactory["ACP_bdgamma"] = bind(boost::factory<Bsgamma*>(), _1, StandardModel::DOWN, 2);
    obsThFactory["BR_bqgamma"] = bind(boost::factory<Bsgamma*>(), _1, 1);
    obsThFactory["ACP_bqgamma"] = bind(boost::factory<Bsgamma*>(), _1, 2);

    //----- Wilson coefficients  -----
    obsThFactory["WC_real_C7g"] = bind(boost::factory<WC_C7g*>(), _1, 0);
    obsThFactory["WC_imag_C7g"] = bind(boost::factory<WC_C7g*>(), _1, 1);
    obsThFactory["WC_abs_C7g"] = bind(boost::factory<WC_C7g*>(), _1, 2);
    obsThFactory["WC_arg_C7g"] = bind(boost::factory<WC_C7g*>(), _1, 3);

    obsThFactory["WC_real_C9_mu"] = bind(boost::factory<WC_C9*>(), _1, 0, StandardModel::MU);
    obsThFactory["WC_imag_C9_mu"] = bind(boost::factory<WC_C9*>(), _1, 1, StandardModel::MU);
    obsThFactory["WC_abs_C9_mu"] = bind(boost::factory<WC_C9*>(), _1, 2, StandardModel::MU);
    obsThFactory["WC_arg_C9_mu"] = bind(boost::factory<WC_C9*>(), _1, 3, StandardModel::MU);

    obsThFactory["WC_real_C9_el"] = bind(boost::factory<WC_C9*>(), _1, 0, StandardModel::ELECTRON);
    obsThFactory["WC_imag_C9_el"] = bind(boost::factory<WC_C9*>(), _1, 1, StandardModel::ELECTRON);
    obsThFactory["WC_abs_C9_el"] = bind(boost::factory<WC_C9*>(), _1, 2, StandardModel::ELECTRON);
    obsThFactory["WC_arg_C9_el"] = bind(boost::factory<WC_C9*>(), _1, 3, StandardModel::ELECTRON);

    obsThFactory["WC_real_C10_mu"] = bind(boost::factory<WC_C10*>(), _1, 0, StandardModel::MU);
    obsThFactory["WC_imag_C10_mu"] = bind(boost::factory<WC_C10*>(), _1, 1, StandardModel::MU);
    obsThFactory["WC_abs_C10_mu"] = bind(boost::factory<WC_C10*>(), _1, 2, StandardModel::MU);
    obsThFactory["WC_arg_C10_mu"] = bind(boost::factory<WC_C10*>(), _1, 3, StandardModel::MU);

    obsThFactory["WC_real_C10_el"] = bind(boost::factory<WC_C10*>(), _1, 0, StandardModel::ELECTRON);
    obsThFactory["WC_imag_C10_el"] = bind(boost::factory<WC_C10*>(), _1, 1, StandardModel::ELECTRON);
    obsThFactory["WC_abs_C10_el"] = bind(boost::factory<WC_C10*>(), _1, 2, StandardModel::ELECTRON);
    obsThFactory["WC_arg_C10_el"] = bind(boost::factory<WC_C10*>(), _1, 3, StandardModel::ELECTRON);

    obsThFactory["WC_arg_C10_el"] = bind(boost::factory<WC_C10*>(), _1, 3, StandardModel::ELECTRON);
    
    //----- B to K* ll  -----
    obsThFactory["P_1_BdKstmu"] = bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_1_BdKste"] = bind(boost::factory<P_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_2_BdKstmu"] = bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2_BdKste"] = bind(boost::factory<P_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_3_BdKstmu"] = bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3_BdKste"] = bind(boost::factory<P_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_4p_BdKstmu"] = bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4p_BdKste"] = bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_5p_BdKstmu"] = bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5p_BdKste"] = bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_6p_BdKstmu"] = bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6p_BdKste"] = bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["P_8p_BdKstmu"] = bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_8p_BdKste"] = bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["Gammap_BdKstmu"] = bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_FB_BdKstmu"] = bind(boost::factory<A_FB*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BR_BdKstmu"] = bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BR_BdKste"] = bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["RKst_BdKstll"] = bind(boost::factory<R_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RKstL_BdKstll"] = bind(boost::factory<RL_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RKstT_BdKstll"] = bind(boost::factory<RT_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["R6_BdKstll"] = bind(boost::factory<R_6*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["ACP_BdKstmu"] = bind(boost::factory<ACP_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P3CP_BdKstmu"] = bind(boost::factory<P3CP*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_L_BdKstmu"] = bind(boost::factory<F_L*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_L_BdKste"] = bind(boost::factory<F_L*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
    obsThFactory["M_1p_BdKstmu"] = bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["M_2p_BdKstmu"] = bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_3_BdKstmu"] = bind(boost::factory<S_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_4_BdKstmu"] = bind(boost::factory<S_4*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_5_BdKstmu"] = bind(boost::factory<S_5*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_7_BdKstmu"] = bind(boost::factory<S_7*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_8_BdKstmu"] = bind(boost::factory<S_8*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_9_BdKstmu"] = bind(boost::factory<S_9*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_5_BdKstmu"] = bind(boost::factory<A_5*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_6_BdKstmu"] = bind(boost::factory<A_6*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_6c_BdKstmu"] = bind(boost::factory<A_6c*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_8_BdKstmu"] = bind(boost::factory<A_8*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_9_BdKstmu"] = bind(boost::factory<A_9*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["P_1f_BdKstmu"] = bind(boost::factory<P_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_2f_BdKstmu"] = bind(boost::factory<P_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_3f_BdKstmu"] = bind(boost::factory<P_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_4pf_BdKstmu"] = bind(boost::factory<P_4Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_5pf_BdKstmu"] = bind(boost::factory<P_5Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_6pf_BdKstmu"] = bind(boost::factory<P_6Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_8pf_BdKstmu"] = bind(boost::factory<P_8Primef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Gammapf_BdKstmu"] = bind(boost::factory<GammaPrimef*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["BRf_BdKstmu"] = bind(boost::factory<BRf_MVll*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["A_FBf_BdKstmu"] = bind(boost::factory<A_FBf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["F_Lf_BdKstmu"] = bind(boost::factory<F_Lf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_3f_BdKstmu"] = bind(boost::factory<S_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_4f_BdKstmu"] = bind(boost::factory<S_4f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_5f_BdKstmu"] = bind(boost::factory<S_5f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_7f_BdKstmu"] = bind(boost::factory<S_7f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_8f_BdKstmu"] = bind(boost::factory<S_8f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_9f_BdKstmu"] = bind(boost::factory<S_9f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_relationf"] = bind(boost::factory<P_relationf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["P_relation_exactf"] = bind(boost::factory<P_relation_exactf*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    
    obsThFactory["F_L_BpKstmu"] = bind(boost::factory<F_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_3_BpKstmu"] = bind(boost::factory<S_3*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_4_BpKstmu"] = bind(boost::factory<S_4*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_5_BpKstmu"] = bind(boost::factory<S_5*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["A_FB_BpKstmu"] = bind(boost::factory<A_FB*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_7_BpKstmu"] = bind(boost::factory<S_7*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_8_BpKstmu"] = bind(boost::factory<S_8*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["S_9_BpKstmu"] = bind(boost::factory<S_9*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    
    obsThFactory["P_1_BpKstmu"] = bind(boost::factory<P_1*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_2_BpKstmu"] = bind(boost::factory<P_2*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_3_BpKstmu"] = bind(boost::factory<P_3*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_4p_BpKstmu"] = bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_5p_BpKstmu"] = bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_6p_BpKstmu"] = bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["P_8p_BpKstmu"] = bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    
    obsThFactory["V0_BdKstmu"] = bind(boost::factory<V0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Vp_BdKstmu"] = bind(boost::factory<Vp*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Vm_BdKstmu"] = bind(boost::factory<Vm*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["T0_BdKstmu"] = bind(boost::factory<T0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Tp_BdKstmu"] = bind(boost::factory<Tp*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["Tm_BdKstmu"] = bind(boost::factory<Tm*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["S_BdKstmu"] = bind(boost::factory<S*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["QCDfC9_1f_BdKstmu"] = bind(boost::factory<QCDfC9_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9_2f_BdKstmu"] = bind(boost::factory<QCDfC9_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9_3f_BdKstmu"] = bind(boost::factory<QCDfC9_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["QCDfC9p_1f_BdKstmu"] = bind(boost::factory<QCDfC9p_1f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9p_2f_BdKstmu"] = bind(boost::factory<QCDfC9p_2f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
    obsThFactory["QCDfC9p_3f_BdKstmu"] = bind(boost::factory<QCDfC9p_3f*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);

    obsThFactory["Regtilde_1_BdKstmu"] = bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Regtilde_2_BdKstmu"] = bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Regtilde_3_BdKstmu"] = bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);

    obsThFactory["Imgtilde_1_BdKstmu"] = bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imgtilde_2_BdKstmu"] = bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imgtilde_3_BdKstmu"] = bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);

    obsThFactory["Absgtilde_1_BdKstmu"] = bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absgtilde_2_BdKstmu"] = bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absgtilde_3_BdKstmu"] = bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);

    obsThFactory["Arggtilde_1_BdKstmu"] = bind(boost::factory<gtilde_1*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Arggtilde_2_BdKstmu"] = bind(boost::factory<gtilde_2*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Arggtilde_3_BdKstmu"] = bind(boost::factory<gtilde_3*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);

    obsThFactory["Reh_0_BdKstmu"] = bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Reh_p_BdKstmu"] = bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);
    obsThFactory["Reh_m_BdKstmu"] = bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 1);

    obsThFactory["Imh_0_BdKstmu"] = bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imh_p_BdKstmu"] = bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);
    obsThFactory["Imh_m_BdKstmu"] = bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 2);

    obsThFactory["Absh_0_BdKstmu"] = bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absh_p_BdKstmu"] = bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);
    obsThFactory["Absh_m_BdKstmu"] = bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 3);

    obsThFactory["Argh_0_BdKstmu"] = bind(boost::factory<h_0*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Argh_p_BdKstmu"] = bind(boost::factory<h_p*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);
    obsThFactory["Argh_m_BdKstmu"] = bind(boost::factory<h_m*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::MU, 4);

    //----- B+ to K*+ ll  -----
    obsThFactory["A_FB_BpKstmu"] = bind(boost::factory<A_FB*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["F_L_BpKstmu"] = bind(boost::factory<F_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);
    obsThFactory["BR_BpKstmu"] = bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU);

/* BEGIN: REMOVE FROM THE PACKAGE */
    //----- B to X_q ll -----
    obsThFactory["R_BXsee"] = bind(boost::factory<R_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["HT_BXsee"] = bind(boost::factory<HT_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["HL_BXsee"] = bind(boost::factory<HL_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["HA_BXsee"] = bind(boost::factory<HA_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["BR_BXsee"] = bind(boost::factory<BR_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["AFB_BXsee"] = bind(boost::factory<AFB_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
/* END: REMOVE FROM THE PACKAGE */

    //----- B to K* gamma  -----
    obsThFactory["BR_BKstgamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ACP_BKstgamma"] = bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["C_BKstgamma"] = bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["S_BKstgamma"] = bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ADG_BKstgamma"] = bind(boost::factory<ADG_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_1"] = bind(boost::factory<DC7_1*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["DC7_2"] = bind(boost::factory<DC7_2*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_L"] = bind(boost::factory<AbsDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_R"] = bind(boost::factory<AbsDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_L_Bd"] = bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_R_Bd"] = bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_L_Bd"] = bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_R_Bd"] = bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["hp0_hm0"] = bind(boost::factory<hp0_hm0*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_QCDF_Bd"] = bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["AbsDC7_QCDF_Bd_bar"] = bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_QCDF_Bd"] = bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ReDC7_QCDF_Bd_bar"] = bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_QCDF_Bd"] = bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["ImDC7_QCDF_Bd_bar"] = bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_D, StandardModel::K_star);


    //----- B+ to K*+ gamma  -----
    obsThFactory["BR_BpKstgamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ACP_BpKstgamma"] = bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_L_Bp"] = bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_R_Bp"] = bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_L_Bp"] = bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_R_Bp"] = bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["AbsDC7_QCDF_Bp"] = bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["AbsDC7_QCDF_Bp_bar"] = bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_QCDF_Bp"] = bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ReDC7_QCDF_Bp_bar"] = bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_QCDF_Bp"] = bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["ImDC7_QCDF_Bp_bar"] = bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_P, StandardModel::K_star_P);

    //----- B to PHI gamma  -----
    obsThFactory["BR_Bsphigamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["C_Bsphigamma"] = bind(boost::factory<C_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["S_Bsphigamma"] = bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ADG_Bsphigamma"] = bind(boost::factory<ADG_MVgamma*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_L_Bs"] = bind(boost::factory<ReDC7_L*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_R_Bs"] = bind(boost::factory<ReDC7_R*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_L_Bs"] = bind(boost::factory<ImDC7_L*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_R_Bs"] = bind(boost::factory<ImDC7_R*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["AbsDC7_QCDF_Bs"] = bind(boost::factory<AbsDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["AbsDC7_QCDF_Bs_bar"] = bind(boost::factory<AbsDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_QCDF_Bs"] = bind(boost::factory<ReDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ReDC7_QCDF_Bs_bar"] = bind(boost::factory<ReDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_QCDF_Bs"] = bind(boost::factory<ImDC7_QCDF*>(), _1, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["ImDC7_QCDF_Bs_bar"] = bind(boost::factory<ImDC7_QCDF_bar*>(), _1, StandardModel::B_S, StandardModel::PHI);

    //----- B to RHO gamma  -----
    obsThFactory["BR_Brhogamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::RHO);
    obsThFactory["ACP_Brhogamma"] = bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_D, StandardModel::RHO);
    obsThFactory["S_Brhogamma"] = bind(boost::factory<S_MVgamma*>(), _1, StandardModel::B_D, StandardModel::RHO);
    
    //----- B+ to RHO+ gamma  -----
    obsThFactory["BR_Bprhogamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_P, StandardModel::RHO_P);
    obsThFactory["ACP_Bprhogamma"] = bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_P, StandardModel::RHO_P);
    
    //----- B to OMEGA gamma  -----
    obsThFactory["BR_Bomegagamma"] = bind(boost::factory<BR_MVgamma*>(), _1, StandardModel::B_D, StandardModel::OMEGA);
    obsThFactory["ACP_Bomegagamma"] = bind(boost::factory<ACP_MVgamma*>(), _1, StandardModel::B_D, StandardModel::OMEGA);
    
    //----- B to V gamma  -----
    obsThFactory["R_BVgamma"] = bind(boost::factory<R_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::B_S, StandardModel::PHI);
    obsThFactory["D0p_BKstgamma"] = bind(boost::factory<D0p_MVgamma*>(), _1, StandardModel::B_D, StandardModel::K_star, StandardModel::B_P, StandardModel::K_star_P);
    obsThFactory["DACP_BKstgamma"] = bind(boost::factory<DACP_MVgamma*>(), _1, StandardModel::B_P, StandardModel::K_star_P, StandardModel::B_D, StandardModel::K_star);
    obsThFactory["D0p_Brhogamma"] = bind(boost::factory<D0p_MVgamma*>(), _1, StandardModel::B_D, StandardModel::RHO, StandardModel::B_P, StandardModel::RHO_P);
    obsThFactory["DACP_Brhogamma"] = bind(boost::factory<DACP_MVgamma*>(), _1, StandardModel::B_P, StandardModel::RHO_P, StandardModel::B_D, StandardModel::RHO);
    
    //----- B to phi ll  -----
    obsThFactory["P_1_Bsphimu"] = bind(boost::factory<P_1*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_2_Bsphimu"] = bind(boost::factory<P_2*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_3_Bsphimu"] = bind(boost::factory<P_3*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_4p_Bsphimu"] = bind(boost::factory<P_4Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_5p_Bsphimu"] = bind(boost::factory<P_5Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_6p_Bsphimu"] = bind(boost::factory<P_6Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P_8p_Bsphimu"] = bind(boost::factory<P_8Prime*>(), _1, StandardModel::B_D, StandardModel::PHI, StandardModel::MU);
    obsThFactory["Gammap_Bsphimu"] = bind(boost::factory<GammaPrime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_FB_Bsphimu"] = bind(boost::factory<A_FB*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["BR_Bsphimu"] = bind(boost::factory<BR_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["Rphi_Bsphill"] = bind(boost::factory<R_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RphiL_Bsphill"] = bind(boost::factory<RL_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RphiT_Bsphill"] = bind(boost::factory<RT_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["R6_Bsphill"] = bind(boost::factory<R_6*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["ACP_Bsphimu"] = bind(boost::factory<ACP_MVll*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["P3CP_Bsphimu"] = bind(boost::factory<P3CP*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["F_L_Bsphimu"] = bind(boost::factory<F_L*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_1p_Bsphimu"] = bind(boost::factory<M_1Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["M_2p_Bsphimu"] = bind(boost::factory<M_2Prime*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_3_Bsphimu"] = bind(boost::factory<S_3*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_4_Bsphimu"] = bind(boost::factory<S_4*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_5_Bsphimu"] = bind(boost::factory<S_5*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_7_Bsphimu"] = bind(boost::factory<S_7*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_8_Bsphimu"] = bind(boost::factory<S_8*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["S_9_Bsphimu"] = bind(boost::factory<S_9*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_5_Bsphimu"] = bind(boost::factory<A_5*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_6_Bsphimu"] = bind(boost::factory<A_6*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_6c_Bsphimu"] = bind(boost::factory<A_6c*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_8_Bsphimu"] = bind(boost::factory<A_8*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
    obsThFactory["A_9_Bsphimu"] = bind(boost::factory<A_9*>(), _1, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);

    //----- B+ to K+ ll  -----
    obsThFactory["BR_BpKmu"] = bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["BR_BpKe"] = bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["dBR_BpKmu"] = bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
    obsThFactory["dBR_BpKe"] = bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
    obsThFactory["RK_BpKll"] = bind(boost::factory<R_MPll*>(), _1, StandardModel::B_P, StandardModel::K_P, StandardModel::MU, StandardModel::ELECTRON);

    //----- B0 to K0 ll  -----
    obsThFactory["BR_B0Kmu"] = bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);
    obsThFactory["BR_B0Ke"] = bind(boost::factory<BR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON);
    obsThFactory["dBR_B0Kmu"] = bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);
    obsThFactory["dBR_B0Ke"] = bind(boost::factory<dBR_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON);
    obsThFactory["RK_B0Kll"] = bind(boost::factory<R_MPll*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU, StandardModel::ELECTRON);
    
    obsThFactory["DC9_hlambda"] = bind(boost::factory<DC9_hlambda*>(), _1, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);

    //----- B to D*lnu -----
    obsThFactory["Gammaw_MVlnu"] = bind(boost::factory<Gammaw_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RDstar_MVlnu"] = bind(boost::factory<RDstar_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["Gammacl_MVlnu"] = bind(boost::factory<Gammacl_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["GammacV_MVlnu"] = bind(boost::factory<GammacV_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["Gammachi_MVlnu"] = bind(boost::factory<Gammachi_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_MVlnu"] = bind(boost::factory<UnitarityV_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityA_MVlnu"] = bind(boost::factory<UnitarityA_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_D_Dst"] = bind(boost::factory<UnitarityV_D_Dst*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["hA1_at_w1"] = bind(boost::factory<FF_hA1atw1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hV_w"] = bind(boost::factory<FF_hV*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA1_w"] = bind(boost::factory<FF_hA1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA2_w"] = bind(boost::factory<FF_hA2*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["hA3_w"] = bind(boost::factory<FF_hA3*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R1_w"] = bind(boost::factory<FF_R1*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R2_w"] = bind(boost::factory<FF_R2*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["R0_w"] = bind(boost::factory<FF_R0*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::MU);
    obsThFactory["FL_MVtaunu"] = bind(boost::factory<FL_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU);
    obsThFactory["Ptau_MVtaunu"] = bind(boost::factory<Plep_MVlnu*>(), _1, StandardModel::B_D, StandardModel::D_star_P, StandardModel::TAU);
    //----- B to Dlnu -----
    obsThFactory["Gammaw_MPlnu"] = bind(boost::factory<Gammaw_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["RD_MPlnu"] = bind(boost::factory<RD_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::TAU, StandardModel::MU, StandardModel::ELECTRON);
    obsThFactory["UnitarityV_MPlnu"] = bind(boost::factory<UnitarityV_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["UnitarityA_MPlnu"] = bind(boost::factory<UnitarityA_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["Unitarity_Strong_MPlnu"] = bind(boost::factory<Unitarity_Strong_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["af0_0"] = bind(boost::factory<af0_0*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["FF0_MPlnu"] = bind(boost::factory<FF0_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);
    obsThFactory["FFplus_MPlnu"] = bind(boost::factory<FFplus_MPlnu*>(), _1, StandardModel::B_D, StandardModel::D_P, StandardModel::ELECTRON);

    //----- B to tau nu  -----
    obsThFactory["btaunu"] = bind(boost::factory<Btaunu*>(), _1, StandardModel::B_P);
    obsThFactory["bctaunu"] = bind(boost::factory<Btaunu*>(), _1, StandardModel::B_C);
    
    /* BEGIN: REMOVE FROM THE PACKAGE */
    obsThFactory["Deltaamu"] = boost::factory<Deltaamu*>();
    /* END: REMOVE FROM THE PACKAGE */

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
    obsThFactory["mueconversion_Ti"] = boost::factory<mueconversion_Ti*>();

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
/* BEGIN: REMOVE FROM THE PACKAGE */
#if FEYNHIGGS
    obsThFactory["OutputSLHAfromFH"] = boost::factory<OutputSLHAfromFH*>(); // for debug
#endif
/* END: REMOVE FROM THE PACKAGE */
    obsThFactory["MHl"] = bind(boost::factory<Mhiggs*>(), _1, 0);
    obsThFactory["MHh"] = bind(boost::factory<Mhiggs*>(), _1, 1);
    obsThFactory["MHa"] = bind(boost::factory<Mhiggs*>(), _1, 2);
    obsThFactory["MHp"] = bind(boost::factory<Mhiggs*>(), _1, 3);
    obsThFactory["Msu1"] = bind(boost::factory<Msup*>(), _1, 0);
    obsThFactory["Msu2"] = bind(boost::factory<Msup*>(), _1, 1);
    obsThFactory["Msu3"] = bind(boost::factory<Msup*>(), _1, 2);
    obsThFactory["Msu4"] = bind(boost::factory<Msup*>(), _1, 3);
    obsThFactory["Msu5"] = bind(boost::factory<Msup*>(), _1, 4);
    obsThFactory["Msu6"] = bind(boost::factory<Msup*>(), _1, 5);
    obsThFactory["Msd1"] = bind(boost::factory<Msdown*>(), _1, 0);
    obsThFactory["Msd2"] = bind(boost::factory<Msdown*>(), _1, 1);
    obsThFactory["Msd3"] = bind(boost::factory<Msdown*>(), _1, 2);
    obsThFactory["Msd4"] = bind(boost::factory<Msdown*>(), _1, 3);
    obsThFactory["Msd5"] = bind(boost::factory<Msdown*>(), _1, 4);
    obsThFactory["Msd6"] = bind(boost::factory<Msdown*>(), _1, 5);
    obsThFactory["Msl1"] = bind(boost::factory<Mslepton*>(), _1, 0);
    obsThFactory["Msl2"] = bind(boost::factory<Mslepton*>(), _1, 1);
    obsThFactory["Msl3"] = bind(boost::factory<Mslepton*>(), _1, 2);
    obsThFactory["Msl4"] = bind(boost::factory<Mslepton*>(), _1, 3);
    obsThFactory["Msl5"] = bind(boost::factory<Mslepton*>(), _1, 4);
    obsThFactory["Msl6"] = bind(boost::factory<Mslepton*>(), _1, 5);
    obsThFactory["Msnu1"] = bind(boost::factory<Msneutrino*>(), _1, 0);
    obsThFactory["Msnu2"] = bind(boost::factory<Msneutrino*>(), _1, 1);
    obsThFactory["Msnu3"] = bind(boost::factory<Msneutrino*>(), _1, 2);
    obsThFactory["Mch1"] = bind(boost::factory<Mchargino*>(), _1, 0);
    obsThFactory["Mch2"] = bind(boost::factory<Mchargino*>(), _1, 1);
    obsThFactory["Mneu1"] = bind(boost::factory<Mneutralino*>(), _1, 0);
    obsThFactory["Mneu2"] = bind(boost::factory<Mneutralino*>(), _1, 1);
    obsThFactory["Mneu3"] = bind(boost::factory<Mneutralino*>(), _1, 2);
    obsThFactory["Mneu4"] = bind(boost::factory<Mneutralino*>(), _1, 3);
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
    obsThFactory["mu_VBF_Vh_htogg"] = boost::factory<VBF_Vh_htogg*>();
    obsThFactory["mu_VBF_Vh_htocc"] = boost::factory<VBF_Vh_htocc*>();
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
    obsThFactory["mu_htobb"] = boost::factory<mu_htobb*>();
    obsThFactory["mu_htoWW"] = boost::factory<mu_htoWW*>();
    obsThFactory["mu_htotautau"] = boost::factory<mu_htotautau*>();
    obsThFactory["mu_htoZga"] = boost::factory<mu_htoZga*>();
    obsThFactory["Gamma_h_THDM"] = boost::factory<Gamma_h_THDM*>();
    obsThFactory["rh_gaga_THDM"] = boost::factory<rh_gaga_THDM*>();
    obsThFactory["rh_Zga_THDM"] = boost::factory<rh_Zga_THDM*>();
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
    obsThFactory["Hobs_pp_H_AZ_bbll_CMS8"] = boost::factory<Hobs_pp_H_AZ_bbll_CMS8*>();
    obsThFactory["Hobs_pp_H_AZ_tautaull_CMS8"] = boost::factory<Hobs_pp_H_AZ_tautaull_CMS8*>();
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
    obsThFactory["Robs_pp_H_AZ_bbll_CMS8"] = boost::factory<Robs_pp_H_AZ_bbll_CMS8*>();
    obsThFactory["Robs_pp_H_AZ_tautaull_CMS8"] = boost::factory<Robs_pp_H_AZ_tautaull_CMS8*>();
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
    obsThFactory["log10_pp_H_AZ_bbll_TH8"] = boost::factory<log10_pp_H_AZ_bbll_TH8*>();
    obsThFactory["log10_pp_H_AZ_tautaull_TH8"] = boost::factory<log10_pp_H_AZ_tautaull_TH8*>();
    obsThFactory["Hobs_ttF_H_tt_ATLAS13"] = boost::factory<Hobs_ttF_H_tt_ATLAS13*>();
    obsThFactory["Hobs_bbF_H_tt_ATLAS13"] = boost::factory<Hobs_bbF_H_tt_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_tautau_ATLAS13"] = boost::factory<Hobs_ggF_H_tautau_ATLAS13*>();
    obsThFactory["Hobs_bbF_H_tautau_ATLAS13"] = boost::factory<Hobs_bbF_H_tautau_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_tautau_CMS13"] = boost::factory<Hobs_ggF_H_tautau_CMS13*>();
    obsThFactory["Hobs_bbF_H_tautau_CMS13"] = boost::factory<Hobs_bbF_H_tautau_CMS13*>();
    obsThFactory["Hobs_pp_H_gaga_ATLAS13"] = boost::factory<Hobs_pp_H_gaga_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_gaga_CMS13"] = boost::factory<Hobs_ggF_H_gaga_CMS13*>();
    obsThFactory["Hobs_pp_H_Zga_llga_ATLAS13"] = boost::factory<Hobs_pp_H_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_Zga_llga_ATLAS13"] = boost::factory<Hobs_ggF_H_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_pp_H_Zga_llga_CMS13"] = boost::factory<Hobs_pp_H_Zga_llga_CMS13*>();
    obsThFactory["Hobs_pp_H_Zga_qqga_CMS13"] = boost::factory<Hobs_pp_H_Zga_qqga_CMS13*>();
    obsThFactory["Hobs_ggF_H_Zga_CMS13"] = boost::factory<Hobs_ggF_H_Zga_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llnunu_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llnunu_ATLAS13*>();
    obsThFactory["Hobs_pp_H_ZZ_llnunu_CMS13"] = boost::factory<Hobs_pp_H_ZZ_llnunu_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llnunu_CMS13"] = boost::factory<Hobs_ggF_H_ZZ_llnunu_CMS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llnunu_CMS13"] = boost::factory<Hobs_VBF_H_ZZ_llnunu_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llll_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llll_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_llll_ATLAS13*>();
    obsThFactory["Hobs_pp_H_ZZ_llll_CMS13"] = boost::factory<Hobs_pp_H_ZZ_llll_CMS13*>();
    obsThFactory["Hobs_VBF_VH_H_ZZ_llll_CMS13"] = boost::factory<Hobs_VBF_VH_H_ZZ_llll_CMS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_qqllnunu_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_qqllnunu_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_qqllnunu_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_qqllnunu_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_llqq_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_ZZ_llqq_ATLAS13"] = boost::factory<Hobs_VBF_H_ZZ_llqq_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_ZZ_nunuqq_ATLAS13"] = boost::factory<Hobs_ggF_H_ZZ_nunuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_H_ZZ_llqq_CMS13"] = boost::factory<Hobs_pp_H_ZZ_llqq_CMS13*>();
    obsThFactory["Hobs_ggF_H_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_ggF_H_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_VBF_H_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_ggF_H_WW_enumunu_ATLAS13"] = boost::factory<Hobs_ggF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VBF_H_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VBF_H_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_ggF_VBF_H_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggF_VBF_H_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_pp_H_VV_qqqq_ATLAS13"] = boost::factory<Hobs_pp_H_VV_qqqq_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbgaga_ATLAS13"] = boost::factory<Hobs_pp_H_hh_bbgaga_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbgaga_CMS13"] = boost::factory<Hobs_pp_H_hh_bbgaga_CMS13*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_ATLAS13"] = boost::factory<Hobs_pp_H_hh_bbbb_ATLAS13*>();
    obsThFactory["Hobs_pp_H_hh_bbbb_CMS13"] = boost::factory<Hobs_pp_H_hh_bbbb_CMS13*>();
    obsThFactory["Hobs_ggF_H_hh_bbbb_CMS13"] = boost::factory<Hobs_ggF_H_hh_bbbb_CMS13*>();
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
    obsThFactory["Robs_ggF_H_Zga_CMS13"] = boost::factory<Robs_ggF_H_Zga_CMS13*>();
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
    obsThFactory["Robs_ggF_H_hh_bbbb_CMS13"] = boost::factory<Robs_ggF_H_hh_bbbb_CMS13*>();
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
    obsThFactory["log10_ggF_H_Zga_TH13"] = boost::factory<log10_ggF_H_Zga_TH13*>();
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
    obsThFactory["log10_ggF_H_hh_bbbb_TH13"] = boost::factory<log10_ggF_H_hh_bbbb_TH13*>();
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
    obsThFactory["Hobs_pp_A_HZ_bbll_CMS8"] = boost::factory<Hobs_pp_A_HZ_bbll_CMS8*>();
    obsThFactory["Hobs_pp_A_HZ_tautaull_CMS8"] = boost::factory<Hobs_pp_A_HZ_tautaull_CMS8*>();
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
    obsThFactory["Robs_pp_A_HZ_bbll_CMS8"] = boost::factory<Robs_pp_A_HZ_bbll_CMS8*>();
    obsThFactory["Robs_pp_A_HZ_tautaull_CMS8"] = boost::factory<Robs_pp_A_HZ_tautaull_CMS8*>();
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
    obsThFactory["log10_pp_A_HZ_bbll_TH8"] = boost::factory<log10_pp_A_HZ_bbll_TH8*>();
    obsThFactory["log10_pp_A_HZ_tautaull_TH8"] = boost::factory<log10_pp_A_HZ_tautaull_TH8*>();
    obsThFactory["Hobs_ttF_A_tt_ATLAS13"] = boost::factory<Hobs_ttF_A_tt_ATLAS13*>();
    obsThFactory["Hobs_bbF_A_tt_ATLAS13"] = boost::factory<Hobs_bbF_A_tt_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_tautau_ATLAS13"] = boost::factory<Hobs_ggF_A_tautau_ATLAS13*>();
    obsThFactory["Hobs_bbF_A_tautau_ATLAS13"] = boost::factory<Hobs_bbF_A_tautau_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_tautau_CMS13"] = boost::factory<Hobs_ggF_A_tautau_CMS13*>();
    obsThFactory["Hobs_bbF_A_tautau_CMS13"] = boost::factory<Hobs_bbF_A_tautau_CMS13*>();
    obsThFactory["Hobs_pp_A_gaga_ATLAS13"] = boost::factory<Hobs_pp_A_gaga_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_gaga_CMS13"] = boost::factory<Hobs_ggF_A_gaga_CMS13*>();
    obsThFactory["Hobs_pp_A_Zga_llga_ATLAS13"] = boost::factory<Hobs_pp_A_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_ggF_A_Zga_llga_ATLAS13"] = boost::factory<Hobs_ggF_A_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_pp_A_Zga_llga_CMS13"] = boost::factory<Hobs_pp_A_Zga_llga_CMS13*>();
    obsThFactory["Hobs_pp_A_Zga_qqga_CMS13"] = boost::factory<Hobs_pp_A_Zga_qqga_CMS13*>();
    obsThFactory["Hobs_ggF_A_Zga_CMS13"] = boost::factory<Hobs_ggF_A_Zga_CMS13*>();
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
    obsThFactory["Robs_ggF_A_Zga_CMS13"] = boost::factory<Robs_ggF_A_Zga_CMS13*>();
    obsThFactory["Robs_ggF_A_hZ_bbZ_ATLAS13"] = boost::factory<Robs_ggF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Robs_bbF_A_hZ_bbZ_ATLAS13"] = boost::factory<Robs_bbF_A_hZ_bbZ_ATLAS13*>();
    obsThFactory["Robs_pp_A_bb_CMS13"] = boost::factory<Robs_pp_A_bb_CMS13*>();
    obsThFactory["log10_ggF_A_tautau_TH13"] = boost::factory<log10_ggF_A_tautau_TH13*>();
    obsThFactory["log10_bbF_A_tautau_TH13"] = boost::factory<log10_bbF_A_tautau_TH13*>();
    obsThFactory["log10_pp_A_gaga_TH13"] = boost::factory<log10_pp_A_gaga_TH13*>();
    obsThFactory["log10_ggF_A_gaga_TH13"] = boost::factory<log10_ggF_A_gaga_TH13*>();
    obsThFactory["log10_pp_A_Zga_TH13"] = boost::factory<log10_pp_A_Zga_TH13*>();
    obsThFactory["log10_ggF_A_Zga_TH13"] = boost::factory<log10_ggF_A_Zga_TH13*>();
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
    obsThFactory["Hobs_pp_Hpm_tb_ATLAS8"] = boost::factory<Hobs_pp_Hpm_tb_ATLAS8*>();
    obsThFactory["Hobs_pp_Hp_tb_CMS8"] = boost::factory<Hobs_pp_Hp_tb_CMS8*>();
    obsThFactory["Hobs_pp_Hpm_taunu_ATLAS13"] = boost::factory<Hobs_pp_Hpm_taunu_ATLAS13*>();
    obsThFactory["Hobs_pp_Hpm_taunu_CMS13"] = boost::factory<Hobs_pp_Hpm_taunu_CMS13*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13_1"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13_1*>();
    obsThFactory["Hobs_pp_Hp_tb_ATLAS13_2"] = boost::factory<Hobs_pp_Hp_tb_ATLAS13_2*>();
    obsThFactory["Robs_pp_Hpm_taunu_ATLAS8"] = boost::factory<Robs_pp_Hpm_taunu_ATLAS8*>();
    obsThFactory["Robs_pp_Hp_taunu_CMS8"] = boost::factory<Robs_pp_Hp_taunu_CMS8*>();
    obsThFactory["Robs_pp_Hpm_tb_ATLAS8"] = boost::factory<Robs_pp_Hpm_tb_ATLAS8*>();
    obsThFactory["Robs_pp_Hp_tb_CMS8"] = boost::factory<Robs_pp_Hp_tb_CMS8*>();
    obsThFactory["Robs_pp_Hpm_taunu_ATLAS13"] = boost::factory<Robs_pp_Hpm_taunu_ATLAS13*>();
    obsThFactory["Robs_pp_Hpm_taunu_CMS13"] = boost::factory<Robs_pp_Hpm_taunu_CMS13*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13"] = boost::factory<Robs_pp_Hp_tb_ATLAS13*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13_1"] = boost::factory<Robs_pp_Hp_tb_ATLAS13_1*>();
    obsThFactory["Robs_pp_Hp_tb_ATLAS13_2"] = boost::factory<Robs_pp_Hp_tb_ATLAS13_2*>();
    obsThFactory["log10_pp_Hpm_taunu_TH8"] = boost::factory<log10_pp_Hpm_taunu_TH8*>();
    obsThFactory["log10_pp_Hp_tb_TH8"] = boost::factory<log10_pp_Hp_tb_TH8*>();
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

    /* BEGIN: REMOVE FROM THE PACKAGE */
    //-----  GeneralTHDM observables  -----
    obsThFactory["mH1"] = boost::factory<mH1_GTHDM*>();
    obsThFactory["mH2"] = boost::factory<mH2_GTHDM*>();
    obsThFactory["mH3"] = boost::factory<mH3_GTHDM*>();
    obsThFactory["m1_2"] = boost::factory<m1_2*>();
    obsThFactory["m2_2"] = boost::factory<m2_2*>();
    obsThFactory["m3_2"] = boost::factory<m3_2*>();
    obsThFactory["mHlight"] = boost::factory<mHlight_GTHDM*>();
    obsThFactory["mHmedium"] = boost::factory<mHmedium_GTHDM*>();
    obsThFactory["mHheavy"] = boost::factory<mHheavy_GTHDM*>();
    obsThFactory["mHp_GTHDM"] = boost::factory<mHp_GTHDM*>();
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
    obsThFactory["Rem12_2_GTHDM"] = boost::factory<Rem12_2_GTHDM*>();
    obsThFactory["Imm12_2_GTHDM"] = boost::factory<Imm12_2_GTHDM*>();
    obsThFactory["lambda1_GTHDM"] = boost::factory<lambda1_GTHDM*>();
    obsThFactory["lambda2_GTHDM"] = boost::factory<lambda2_GTHDM*>();
    obsThFactory["lambda3_GTHDM"] = boost::factory<lambda3_GTHDM*>();
    obsThFactory["lambda4_GTHDM"] = boost::factory<lambda4_GTHDM*>();
    obsThFactory["v1_GTHDM"] = boost::factory<v1_GTHDM*>();
    obsThFactory["v2_GTHDM"] = boost::factory<v2_GTHDM*>();
    obsThFactory["tanbeta"] = boost::factory<tanbeta_GTHDM*>();

    obsThFactory["lambda1H_GTHDM"] = boost::factory<lambda1H_GTHDM*>();
    obsThFactory["lambda2H_GTHDM"] = boost::factory<lambda2H_GTHDM*>();
    obsThFactory["lambda3H_GTHDM"] = boost::factory<lambda3H_GTHDM*>();
    obsThFactory["lambda4H_GTHDM"] = boost::factory<lambda4H_GTHDM*>();
    obsThFactory["Relambda5H_GTHDM"] = boost::factory<Relambda5H_GTHDM*>();
    obsThFactory["Imlambda5H_GTHDM"] = boost::factory<Imlambda5H_GTHDM*>();
    obsThFactory["Relambda6H_GTHDM"] = boost::factory<Relambda6H_GTHDM*>();
    obsThFactory["Imlambda6H_GTHDM"] = boost::factory<Imlambda6H_GTHDM*>();
    obsThFactory["Relambda7H_GTHDM"] = boost::factory<Relambda7H_GTHDM*>();
    obsThFactory["Imlambda7H_GTHDM"] = boost::factory<Imlambda7H_GTHDM*>();


    obsThFactory["R11"]= boost::factory<R11_GTHDM*>();
    obsThFactory["R12"]= boost::factory<R12_GTHDM*>();
    obsThFactory["R13"]= boost::factory<R13_GTHDM*>();
    obsThFactory["R21"]= boost::factory<R21_GTHDM*>();
    obsThFactory["R22"]= boost::factory<R22_GTHDM*>();
    obsThFactory["R23"]= boost::factory<R23_GTHDM*>();
    obsThFactory["R31"]= boost::factory<R31_GTHDM*>();
    obsThFactory["R32"]= boost::factory<R32_GTHDM*>();
    obsThFactory["R33"]= boost::factory<R33_GTHDM*>();

    obsThFactory["cosalpha1"]= boost::factory<cosalpha1_GTHDM*>();

    obsThFactory["Q_stGTHDM"]= boost::factory<Q_stGTHDM*>();
    obsThFactory["DeltaQ_GTHDM"]= boost::factory<DeltaQ_GTHDM*>();
    obsThFactory["g1atQGTHDM"]= boost::factory<g1atQGTHDM*>();
    obsThFactory["g2atQGTHDM"]= boost::factory<g2atQGTHDM*>();
    obsThFactory["g3atQGTHDM"]= boost::factory<g3atQGTHDM*>();
    obsThFactory["etaU1atQGTHDM"]= boost::factory<etaU1atQGTHDM*>();
    obsThFactory["etaU2atQGTHDM"]= boost::factory<etaU2atQGTHDM*>();
    obsThFactory["etaD1atQGTHDM"]= boost::factory<etaD1atQGTHDM*>();
    obsThFactory["etaD2atQGTHDM"]= boost::factory<etaD2atQGTHDM*>();
    obsThFactory["etaL1atQGTHDM"]= boost::factory<etaL1atQGTHDM*>();
    obsThFactory["etaL2atQGTHDM"]= boost::factory<etaL2atQGTHDM*>();
    obsThFactory["lambda1atQGTHDM"]= boost::factory<lambda1atQGTHDM*>();
    obsThFactory["lambda2atQGTHDM"]= boost::factory<lambda2atQGTHDM*>();
    obsThFactory["lambda3atQGTHDM"]= boost::factory<lambda3atQGTHDM*>();
    obsThFactory["lambda4atQGTHDM"]= boost::factory<lambda4atQGTHDM*>();
    obsThFactory["Relambda5atQGTHDM"]= boost::factory<Relambda5atQGTHDM*>();
    obsThFactory["Relambda6atQGTHDM"]= boost::factory<Relambda6atQGTHDM*>();
    obsThFactory["Relambda7atQGTHDM"]= boost::factory<Relambda7atQGTHDM*>();

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

    obsThFactory["EffectivePotMin1_GTHDM"] = boost::factory<EffectivePotMin1_GTHDM*>();
    obsThFactory["EffectivePotMin2_GTHDM"] = boost::factory<EffectivePotMin2_GTHDM*>();


    obsThFactory["DeltaS_GTHDM"] = boost::factory<GTHDMDeltaS*>();
    obsThFactory["DeltaT_GTHDM"] = boost::factory<GTHDMDeltaT*>();
    obsThFactory["DeltaU_GTHDM"] = boost::factory<GTHDMDeltaU*>();

    obsThFactory["Rb0_GTHDM"] = boost::factory<Rb0GTHDM*>();

    obsThFactory["GTHDMgminus2_mu"] = boost::factory<GeneralTHDMgminus2_mu*>();

    obsThFactory["GTHDM_mu_ggF_tth_htobb8"] = boost::factory<GTHDM_ggF_tth_htobb8*>();
    obsThFactory["GTHDM_mu_ggF_tth_htoWW8"] = boost::factory<GTHDM_ggF_tth_htoWW8*>();
    obsThFactory["GTHDM_mu_ggF_tth_htotautau8"] = boost::factory<GTHDM_ggF_tth_htotautau8*>();
    obsThFactory["GTHDM_mu_ggF_tth_htoZZ8"] = boost::factory<GTHDM_ggF_tth_htoZZ8*>();
    obsThFactory["GTHDM_mu_ggF_tth_htogaga8"] = boost::factory<GTHDM_ggF_tth_htogaga8*>();
    obsThFactory["GTHDM_mu_ggF_tth_htobb13"] = boost::factory<GTHDM_ggF_tth_htobb13*>();
    obsThFactory["GTHDM_mu_ggF_tth_htoWW13"] = boost::factory<GTHDM_ggF_tth_htoWW13*>();
    obsThFactory["GTHDM_mu_ggF_tth_htotautau13"] = boost::factory<GTHDM_ggF_tth_htotautau13*>();
    obsThFactory["GTHDM_mu_ggF_tth_htoZZ13"] = boost::factory<GTHDM_ggF_tth_htoZZ13*>();
    obsThFactory["GTHDM_mu_ggF_tth_htogaga13"] = boost::factory<GTHDM_ggF_tth_htogaga13*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htobb"] = boost::factory<GTHDM_VBF_Vh_htobb*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htoWW"] = boost::factory<GTHDM_VBF_Vh_htoWW*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htotautau"] = boost::factory<GTHDM_VBF_Vh_htotautau*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htoZZ"] = boost::factory<GTHDM_VBF_Vh_htoZZ*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htogaga"] = boost::factory<GTHDM_VBF_Vh_htogaga*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htogg"] = boost::factory<GTHDM_VBF_Vh_htogg*>();
    obsThFactory["GTHDM_mu_VBF_Vh_htocc"] = boost::factory<GTHDM_VBF_Vh_htocc*>();
    obsThFactory["GTHDM_mu_ggF_htobb"] = boost::factory<GTHDM_ggF_htobb*>();
    obsThFactory["GTHDM_mu_ggF_htoWW"] = boost::factory<GTHDM_ggF_htoWW*>();
    obsThFactory["GTHDM_mu_ggF_htotautau"] = boost::factory<GTHDM_ggF_htotautau*>();
    obsThFactory["GTHDM_mu_ggF_htoZZ"] = boost::factory<GTHDM_ggF_htoZZ*>();
    obsThFactory["GTHDM_mu_ggF_htogaga"] = boost::factory<GTHDM_ggF_htogaga*>();
    obsThFactory["GTHDM_mu_tth_htobb"] = boost::factory<GTHDM_tth_htobb*>();
    obsThFactory["GTHDM_mu_tth_htoWW"] = boost::factory<GTHDM_tth_htoWW*>();
    obsThFactory["GTHDM_mu_tth_htotautau"] = boost::factory<GTHDM_tth_htotautau*>();
    obsThFactory["GTHDM_mu_tth_htoZZ"] = boost::factory<GTHDM_tth_htoZZ*>();
    obsThFactory["GTHDM_mu_tth_htogaga"] = boost::factory<GTHDM_tth_htogaga*>();
    obsThFactory["GTHDM_mu_htobb"] = boost::factory<GTHDM_mu_htobb*>();
    obsThFactory["GTHDM_mu_htoWW"] = boost::factory<GTHDM_mu_htoWW*>();
    obsThFactory["GTHDM_mu_htotautau"] = boost::factory<GTHDM_mu_htotautau*>();
    obsThFactory["GTHDM_mu_htoZga"] = boost::factory<GTHDM_mu_htoZga*>();

   obsThFactory["yu1r_GTHDM"] = boost::factory<yu1r_GTHDM*>();
   obsThFactory["yd1r_GTHDM"] = boost::factory<yd1r_GTHDM*>();
   obsThFactory["yl1r_GTHDM"] = boost::factory<yl1r_GTHDM*>();
   
   obsThFactory["yu2r_GTHDM"] = boost::factory<yu2r_GTHDM*>();
   obsThFactory["yd2r_GTHDM"] = boost::factory<yd2r_GTHDM*>();
   obsThFactory["yl2r_GTHDM"] = boost::factory<yl2r_GTHDM*>();
   
   obsThFactory["yu3r_GTHDM"] = boost::factory<yu3r_GTHDM*>();
   obsThFactory["yd3r_GTHDM"] = boost::factory<yd3r_GTHDM*>();
   obsThFactory["yl3r_GTHDM"] = boost::factory<yl3r_GTHDM*>();

   obsThFactory["suR_GTHDM"] = boost::factory<suR_GTHDM*>();
   obsThFactory["sdR_GTHDM"] = boost::factory<sdR_GTHDM*>();
   obsThFactory["slR_GTHDM"] = boost::factory<slR_GTHDM*>();
   

   obsThFactory["rh_gg_GTHDM"] = boost::factory<rh_gg_GTHDM*>();
   obsThFactory["rh_gaga_GTHDM"] = boost::factory<rh_gaga_GTHDM*>();
   obsThFactory["rh_Zga_GTHDM"] = boost::factory<rh_Zga_GTHDM*>();



   /* obsThFactory["Gamma_HH_THDM"] = boost::factory<Gamma_HH_THDM*>();
    obsThFactory["rHH_gg_THDM"] = boost::factory<rHH_gg_THDM*>();
    obsThFactory["BR_HH_hh_THDM"] = boost::factory<BR_HH_hh_THDM*>();
    obsThFactory["BR_HH_AA_THDM"] = boost::factory<BR_HH_AA_THDM*>();
    obsThFactory["BR_HH_HpHm_THDM"] = boost::factory<BR_HH_HpHm_THDM*>();
    obsThFactory["BR_HH_AZ_THDM"] = boost::factory<BR_HH_AZ_THDM*>();
    obsThFactory["BR_HH_HpW_THDM"] = boost::factory<BR_HH_HpW_THDM*>(); */

    obsThFactory["Hobs_tt_phi2_tt_ATLAS13"] = boost::factory<Hobs_tt_phi2_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_phi2_tt_ATLAS13"] = boost::factory<Hobs_bb_phi2_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_phi2_bb_CMS8"] = boost::factory<Hobs_bb_phi2_bb_CMS8*>();
    obsThFactory["Hobs_gg_phi2_bb_CMS8"] = boost::factory<Hobs_gg_phi2_bb_CMS8*>();
    obsThFactory["Hobs_pp_phi2_bb_CMS13"] = boost::factory<Hobs_pp_phi2_bb_CMS13*>();
    obsThFactory["Hobs_bb_phi2_bb_CMS13"] = boost::factory<Hobs_bb_phi2_bb_CMS13*>();
    obsThFactory["Hobs_gg_phi2_tautau_ATLAS8"] = boost::factory<Hobs_gg_phi2_tautau_ATLAS8*>();
    obsThFactory["Hobs_bb_phi2_tautau_ATLAS8"] = boost::factory<Hobs_bb_phi2_tautau_ATLAS8*>();
    obsThFactory["Hobs_gg_phi2_tautau_CMS8"] = boost::factory<Hobs_gg_phi2_tautau_CMS8*>();
    obsThFactory["Hobs_bb_phi2_tautau_CMS8"] = boost::factory<Hobs_bb_phi2_tautau_CMS8*>();
    obsThFactory["Hobs_gg_phi2_tautau_ATLAS13"] = boost::factory<Hobs_gg_phi2_tautau_ATLAS13*>();
    obsThFactory["Hobs_bb_phi2_tautau_ATLAS13"] = boost::factory<Hobs_bb_phi2_tautau_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_tautau_CMS13"] = boost::factory<Hobs_gg_phi2_tautau_CMS13*>();
    obsThFactory["Hobs_bb_phi2_tautau_CMS13"] = boost::factory<Hobs_bb_phi2_tautau_CMS13*>();
    obsThFactory["Hobs_gg_phi2_gaga_ATLAS8"] = boost::factory<Hobs_gg_phi2_gaga_ATLAS8*>();
    obsThFactory["Hobs_pp_phi2_gaga_ATLAS13"] = boost::factory<Hobs_pp_phi2_gaga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_gaga_CMS13"] = boost::factory<Hobs_gg_phi2_gaga_CMS13*>();
    obsThFactory["Hobs_pp_phi2_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_phi2_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_phi2_Zga_llga_CMS8"] = boost::factory<Hobs_pp_phi2_Zga_llga_CMS8*>();
    obsThFactory["Hobs_gg_phi2_Zga_llga_ATLAS13"] = boost::factory<Hobs_gg_phi2_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_Zga_qqga_ATLAS13"] = boost::factory<Hobs_gg_phi2_Zga_qqga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_Zga_CMS13"] = boost::factory<Hobs_gg_phi2_Zga_CMS13*>();
    obsThFactory["Hobs_gg_phi2_ZZ_ATLAS8"] = boost::factory<Hobs_gg_phi2_ZZ_ATLAS8*>();
    obsThFactory["Hobs_VV_phi2_ZZ_ATLAS8"] = boost::factory<Hobs_VV_phi2_ZZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi2_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_gg_phi2_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_VV_phi2_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_VV_phi2_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13"] = boost::factory<Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13*>();
    obsThFactory["Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13"] = boost::factory<Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_phi2_ZZ_llqqnunull_CMS13"] = boost::factory<Hobs_pp_phi2_ZZ_llqqnunull_CMS13*>();
    obsThFactory["Hobs_pp_phi2_ZZ_qqnunu_CMS13"] = boost::factory<Hobs_pp_phi2_ZZ_qqnunu_CMS13*>();
    obsThFactory["Hobs_gg_phi2_WW_ATLAS8"] = boost::factory<Hobs_gg_phi2_WW_ATLAS8*>();
    obsThFactory["Hobs_VV_phi2_WW_ATLAS8"] = boost::factory<Hobs_VV_phi2_WW_ATLAS8*>();
    obsThFactory["Hobs_gg_phi2_WW_enumunu_ATLAS13"] = boost::factory<Hobs_gg_phi2_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VV_phi2_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VV_phi2_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_ggVV_phi2_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggVV_phi2_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_gg_phi2_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_gg_phi2_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_VV_phi2_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_VV_phi2_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_phi2_WW_lnuqq_CMS13"] = boost::factory<Hobs_pp_phi2_WW_lnuqq_CMS13*>();
    obsThFactory["Hobs_mu_pp_phi2_VV_CMS8"] = boost::factory<Hobs_mu_pp_phi2_VV_CMS8*>();
    obsThFactory["Hobs_pp_phi2_VV_qqqq_ATLAS13"] = boost::factory<Hobs_pp_phi2_VV_qqqq_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_phi1phi1_ATLAS8"] = boost::factory<Hobs_gg_phi2_phi1phi1_ATLAS8*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbbb_CMS8"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbbb_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbgaga_CMS8"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbgaga_CMS8*>();
    obsThFactory["Hobs_gg_phi2_phi1phi1_bbtautau_CMS8"] = boost::factory<Hobs_gg_phi2_phi1phi1_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbtautau_CMS8"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbgaga_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbgaga_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbVV_CMS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbVV_CMS13*>();
    obsThFactory["Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13"] = boost::factory<Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13"] = boost::factory<Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_bbZ_ATLAS8"] = boost::factory<Hobs_gg_phi2_phi1Z_bbZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_bbll_CMS8"] = boost::factory<Hobs_gg_phi2_phi1Z_bbll_CMS8*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8"] = boost::factory<Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_tautaull_CMS8"] = boost::factory<Hobs_gg_phi2_phi1Z_tautaull_CMS8*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_bbZ_ATLAS13"] = boost::factory<Hobs_gg_phi2_phi1Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_bbZ_1_CMS13"] = boost::factory<Hobs_gg_phi2_phi1Z_bbZ_1_CMS13*>();
    obsThFactory["Hobs_gg_phi2_phi1Z_bbZ_2_CMS13"] = boost::factory<Hobs_gg_phi2_phi1Z_bbZ_2_CMS13*>();
    obsThFactory["Hobs_bb_phi2_phi1Z_bbZ_ATLAS13"] = boost::factory<Hobs_bb_phi2_phi1Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bb_phi2_phi1Z_bbZ_1_CMS13"] = boost::factory<Hobs_bb_phi2_phi1Z_bbZ_1_CMS13*>();
    obsThFactory["Hobs_bb_phi2_phi1Z_bbZ_2_CMS13"] = boost::factory<Hobs_bb_phi2_phi1Z_bbZ_2_CMS13*>();

    obsThFactory["Hobs_tt_phi3_tt_ATLAS13"] = boost::factory<Hobs_tt_phi3_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_phi3_tt_ATLAS13"] = boost::factory<Hobs_bb_phi3_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_phi3_bb_CMS8"] = boost::factory<Hobs_bb_phi3_bb_CMS8*>();
    obsThFactory["Hobs_gg_phi3_bb_CMS8"] = boost::factory<Hobs_gg_phi3_bb_CMS8*>();
    obsThFactory["Hobs_pp_phi3_bb_CMS13"] = boost::factory<Hobs_pp_phi3_bb_CMS13*>();
    obsThFactory["Hobs_bb_phi3_bb_CMS13"] = boost::factory<Hobs_bb_phi3_bb_CMS13*>();
    obsThFactory["Hobs_gg_phi3_tautau_ATLAS8"] = boost::factory<Hobs_gg_phi3_tautau_ATLAS8*>();
    obsThFactory["Hobs_bb_phi3_tautau_ATLAS8"] = boost::factory<Hobs_bb_phi3_tautau_ATLAS8*>();
    obsThFactory["Hobs_gg_phi3_tautau_CMS8"] = boost::factory<Hobs_gg_phi3_tautau_CMS8*>();
    obsThFactory["Hobs_bb_phi3_tautau_CMS8"] = boost::factory<Hobs_bb_phi3_tautau_CMS8*>();
    obsThFactory["Hobs_gg_phi3_tautau_ATLAS13"] = boost::factory<Hobs_gg_phi3_tautau_ATLAS13*>();
    obsThFactory["Hobs_bb_phi3_tautau_ATLAS13"] = boost::factory<Hobs_bb_phi3_tautau_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_tautau_CMS13"] = boost::factory<Hobs_gg_phi3_tautau_CMS13*>();
    obsThFactory["Hobs_bb_phi3_tautau_CMS13"] = boost::factory<Hobs_bb_phi3_tautau_CMS13*>();
    obsThFactory["Hobs_gg_phi3_gaga_ATLAS8"] = boost::factory<Hobs_gg_phi3_gaga_ATLAS8*>();
    obsThFactory["Hobs_pp_phi3_gaga_ATLAS13"] = boost::factory<Hobs_pp_phi3_gaga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_gaga_CMS13"] = boost::factory<Hobs_gg_phi3_gaga_CMS13*>();
    obsThFactory["Hobs_pp_phi3_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_phi3_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_phi3_Zga_llga_CMS8"] = boost::factory<Hobs_pp_phi3_Zga_llga_CMS8*>();
    obsThFactory["Hobs_gg_phi3_Zga_llga_ATLAS13"] = boost::factory<Hobs_gg_phi3_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_Zga_qqga_ATLAS13"] = boost::factory<Hobs_gg_phi3_Zga_qqga_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_Zga_CMS13"] = boost::factory<Hobs_gg_phi3_Zga_CMS13*>();
    obsThFactory["Hobs_gg_phi3_ZZ_ATLAS8"] = boost::factory<Hobs_gg_phi3_ZZ_ATLAS8*>();
    obsThFactory["Hobs_VV_phi3_ZZ_ATLAS8"] = boost::factory<Hobs_VV_phi3_ZZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi3_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_gg_phi3_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_VV_phi3_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_VV_phi3_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13"] = boost::factory<Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13*>();
    obsThFactory["Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13"] = boost::factory<Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_phi3_ZZ_llqqnunull_CMS13"] = boost::factory<Hobs_pp_phi3_ZZ_llqqnunull_CMS13*>();
    obsThFactory["Hobs_pp_phi3_ZZ_qqnunu_CMS13"] = boost::factory<Hobs_pp_phi3_ZZ_qqnunu_CMS13*>();
    obsThFactory["Hobs_gg_phi3_WW_ATLAS8"] = boost::factory<Hobs_gg_phi3_WW_ATLAS8*>();
    obsThFactory["Hobs_VV_phi3_WW_ATLAS8"] = boost::factory<Hobs_VV_phi3_WW_ATLAS8*>();
    obsThFactory["Hobs_gg_phi3_WW_enumunu_ATLAS13"] = boost::factory<Hobs_gg_phi3_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VV_phi3_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VV_phi3_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_ggVV_phi3_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggVV_phi3_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_gg_phi3_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_gg_phi3_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_VV_phi3_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_VV_phi3_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_pp_phi3_WW_lnuqq_CMS13"] = boost::factory<Hobs_pp_phi3_WW_lnuqq_CMS13*>();
    obsThFactory["Hobs_mu_pp_phi3_VV_CMS8"] = boost::factory<Hobs_mu_pp_phi3_VV_CMS8*>();
    obsThFactory["Hobs_pp_phi3_VV_qqqq_ATLAS13"] = boost::factory<Hobs_pp_phi3_VV_qqqq_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_phi1phi1_ATLAS8"] = boost::factory<Hobs_gg_phi3_phi1phi1_ATLAS8*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbbb_CMS8"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbbb_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbgaga_CMS8"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbgaga_CMS8*>();
    obsThFactory["Hobs_gg_phi3_phi1phi1_bbtautau_CMS8"] = boost::factory<Hobs_gg_phi3_phi1phi1_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbtautau_CMS8"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbgaga_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbgaga_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbVV_CMS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbVV_CMS13*>();
    obsThFactory["Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13"] = boost::factory<Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13"] = boost::factory<Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13*>();
   
    
    
    obsThFactory["Hobs_gg_phi3_phi1Z_bbZ_ATLAS8"] = boost::factory<Hobs_gg_phi3_phi1Z_bbZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_bbll_CMS8"] = boost::factory<Hobs_gg_phi3_phi1Z_bbll_CMS8*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8"] = boost::factory<Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_tautaull_CMS8"] = boost::factory<Hobs_gg_phi3_phi1Z_tautaull_CMS8*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_bbZ_ATLAS13"] = boost::factory<Hobs_gg_phi3_phi1Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_bbZ_1_CMS13"] = boost::factory<Hobs_gg_phi3_phi1Z_bbZ_1_CMS13*>();
    obsThFactory["Hobs_gg_phi3_phi1Z_bbZ_2_CMS13"] = boost::factory<Hobs_gg_phi3_phi1Z_bbZ_2_CMS13*>();
    obsThFactory["Hobs_bb_phi3_phi1Z_bbZ_ATLAS13"] = boost::factory<Hobs_bb_phi3_phi1Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bb_phi3_phi1Z_bbZ_1_CMS13"] = boost::factory<Hobs_bb_phi3_phi1Z_bbZ_1_CMS13*>();
    obsThFactory["Hobs_bb_phi3_phi1Z_bbZ_2_CMS13"] = boost::factory<Hobs_bb_phi3_phi1Z_bbZ_2_CMS13*>();

    obsThFactory["Hobs_pp_phi3_phi2Z_bbll_1_CMS8"] = boost::factory<Hobs_pp_phi3_phi2Z_bbll_1_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi3Z_bbll_1_CMS8"] = boost::factory<Hobs_pp_phi2_phi3Z_bbll_1_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi2Z_bbll_2_CMS8"] = boost::factory<Hobs_pp_phi3_phi2Z_bbll_2_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi3Z_bbll_2_CMS8"] = boost::factory<Hobs_pp_phi2_phi3Z_bbll_2_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi2Z_tautaull_1_CMS8"] = boost::factory<Hobs_pp_phi3_phi2Z_tautaull_1_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi3Z_tautaull_1_CMS8"] = boost::factory<Hobs_pp_phi2_phi3Z_tautaull_1_CMS8*>();
    obsThFactory["Hobs_pp_phi3_phi2Z_tautaull_2_CMS8"] = boost::factory<Hobs_pp_phi3_phi2Z_tautaull_2_CMS8*>();
    obsThFactory["Hobs_pp_phi2_phi3Z_tautaull_2_CMS8"] = boost::factory<Hobs_pp_phi2_phi3Z_tautaull_2_CMS8*>();
    obsThFactory["Hobs_gg_phi3_phi2Z_bbZ_ATLAS13"] = boost::factory<Hobs_gg_phi3_phi2Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_gg_phi2_phi3Z_bbZ_ATLAS13"] = boost::factory<Hobs_gg_phi2_phi3Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bb_phi3_phi2Z_bbZ_ATLAS13"] = boost::factory<Hobs_bb_phi3_phi2Z_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bb_phi2_phi3Z_bbZ_ATLAS13"] = boost::factory<Hobs_bb_phi2_phi3Z_bbZ_ATLAS13*>();

    obsThFactory["Hobs_pp_Hpm_taunu_ATLAS8"] = boost::factory<Hobs_pp_Hpm_taunu_ATLAS8_GTHDM*>();
    obsThFactory["Hobs_pp_Hp_taunu_CMS8"] = boost::factory<Hobs_pp_Hp_taunu_CMS8_GTHDM*>();
    obsThFactory["Hobs_pp_Hpm_taunu_ATLAS13"] = boost::factory<Hobs_pp_Hpm_taunu_ATLAS13_GTHDM*>();
    obsThFactory["Hobs_pp_Hpm_taunu_CMS13"] = boost::factory<Hobs_pp_Hpm_taunu_CMS13_GTHDM*>();
    obsThFactory["Hobs_pp_Hpm_tb_ATLAS8"] = boost::factory<Hobs_pp_Hpm_tb_ATLAS8_GTHDM*>();
    obsThFactory["Hobs_pp_Hp_tb_CMS8"] = boost::factory<Hobs_pp_Hp_tb_CMS8_GTHDM*>();
    obsThFactory["Hobs_pp_Hpm_tb_ATLAS13"] = boost::factory<Hobs_pp_Hpm_tb_ATLAS13*>();

    obsThFactory["log10_tt_phi2_tt_TH13"] = boost::factory<log10_tt_phi2_tt_TH13*>();
    obsThFactory["log10_tt_phi3_tt_TH13"] = boost::factory<log10_tt_phi3_tt_TH13*>();
    obsThFactory["log10_bb_phi2_tt_TH13"] = boost::factory<log10_bb_phi2_tt_TH13*>();
    obsThFactory["log10_bb_phi3_tt_TH13"] = boost::factory<log10_bb_phi3_tt_TH13*>();

    obsThFactory["log10_bb_phi2_bb_TH8"] = boost::factory<log10_bb_phi2_bb_TH8*>();
    obsThFactory["log10_bb_phi3_bb_TH8"] = boost::factory<log10_bb_phi3_bb_TH8*>();
    obsThFactory["log10_gg_phi2_bb_TH8"] = boost::factory<log10_gg_phi2_bb_TH8*>();
    obsThFactory["log10_gg_phi3_bb_TH8"] = boost::factory<log10_gg_phi3_bb_TH8*>();
    obsThFactory["log10_pp_phi2_bb_TH13"] = boost::factory<log10_pp_phi2_bb_TH13*>();
    obsThFactory["log10_pp_phi3_bb_TH13"] = boost::factory<log10_pp_phi3_bb_TH13*>();
    obsThFactory["log10_bb_phi2_bb_TH13"] = boost::factory<log10_bb_phi2_bb_TH13*>();
    obsThFactory["log10_bb_phi3_bb_TH13"] = boost::factory<log10_bb_phi3_bb_TH13*>();

    obsThFactory["log10_gg_phi2_tautau_TH8"] = boost::factory<log10_gg_phi2_tautau_TH8*>();
    obsThFactory["log10_gg_phi3_tautau_TH8"] = boost::factory<log10_gg_phi3_tautau_TH8*>();
    obsThFactory["log10_bb_phi2_tautau_TH8"] = boost::factory<log10_bb_phi2_tautau_TH8*>();
    obsThFactory["log10_bb_phi3_tautau_TH8"] = boost::factory<log10_bb_phi3_tautau_TH8*>();
    obsThFactory["log10_gg_phi2_tautau_TH13"] = boost::factory<log10_gg_phi2_tautau_TH13*>();
    obsThFactory["log10_gg_phi3_tautau_TH13"] = boost::factory<log10_gg_phi3_tautau_TH13*>();
    obsThFactory["log10_bb_phi2_tautau_TH13"] = boost::factory<log10_bb_phi2_tautau_TH13*>();
    obsThFactory["log10_bb_phi3_tautau_TH13"] = boost::factory<log10_bb_phi3_tautau_TH13*>();

    obsThFactory["log10_gg_phi2_gaga_TH8"] = boost::factory<log10_gg_phi2_gaga_TH8*>();
    obsThFactory["log10_gg_phi3_gaga_TH8"] = boost::factory<log10_gg_phi3_gaga_TH8*>();
    obsThFactory["log10_pp_phi2_gaga_TH13"] = boost::factory<log10_pp_phi2_gaga_TH13*>();
    obsThFactory["log10_pp_phi3_gaga_TH13"] = boost::factory<log10_pp_phi3_gaga_TH13*>();
    obsThFactory["log10_gg_phi2_gaga_TH13"] = boost::factory<log10_gg_phi2_gaga_TH13*>();
    obsThFactory["log10_gg_phi3_gaga_TH13"] = boost::factory<log10_gg_phi3_gaga_TH13*>();


    obsThFactory["log10_pp_phi2_Zga_llga_TH8"] = boost::factory<log10_pp_phi2_Zga_llga_TH8*>();
    obsThFactory["log10_pp_phi3_Zga_llga_TH8"] = boost::factory<log10_pp_phi3_Zga_llga_TH8*>();
    obsThFactory["log10_gg_phi2_Zga_TH13"] = boost::factory<log10_gg_phi2_Zga_TH13*>();
    obsThFactory["log10_gg_phi3_Zga_TH13"] = boost::factory<log10_gg_phi3_Zga_TH13*>();

    obsThFactory["log10_gg_phi2_ZZ_TH8"] = boost::factory<log10_gg_phi2_ZZ_TH8*>();
    obsThFactory["log10_gg_phi3_ZZ_TH8"] = boost::factory<log10_gg_phi3_ZZ_TH8*>();
    obsThFactory["log10_VV_phi2_ZZ_TH8"] = boost::factory<log10_VV_phi2_ZZ_TH8*>();
    obsThFactory["log10_VV_phi3_ZZ_TH8"] = boost::factory<log10_VV_phi3_ZZ_TH8*>();
    obsThFactory["log10_gg_phi2_ZZ_TH13"] = boost::factory<log10_gg_phi2_ZZ_TH13*>();
    obsThFactory["log10_gg_phi3_ZZ_TH13"] = boost::factory<log10_gg_phi3_ZZ_TH13*>();
    obsThFactory["log10_VV_phi2_ZZ_TH13"] = boost::factory<log10_VV_phi2_ZZ_TH13*>();
    obsThFactory["log10_VV_phi3_ZZ_TH13"] = boost::factory<log10_VV_phi3_ZZ_TH13*>();
    obsThFactory["log10_pp_phi2_ZZ_TH13"] = boost::factory<log10_pp_phi2_ZZ_TH13*>();
    obsThFactory["log10_pp_phi3_ZZ_TH13"] = boost::factory<log10_pp_phi3_ZZ_TH13*>();

    obsThFactory["log10_gg_phi2_WW_TH8"] = boost::factory<log10_gg_phi2_WW_TH8*>();
    obsThFactory["log10_gg_phi3_WW_TH8"] = boost::factory<log10_gg_phi3_WW_TH8*>();
    obsThFactory["log10_VV_phi2_WW_TH8"] = boost::factory<log10_VV_phi2_WW_TH8*>();
    obsThFactory["log10_VV_phi3_WW_TH8"] = boost::factory<log10_VV_phi3_WW_TH8*>();
    obsThFactory["log10_gg_phi2_WW_TH13"] = boost::factory<log10_gg_phi2_WW_TH13*>();
    obsThFactory["log10_gg_phi3_WW_TH13"] = boost::factory<log10_gg_phi3_WW_TH13*>();
    obsThFactory["log10_VV_phi2_WW_TH13"] = boost::factory<log10_VV_phi2_WW_TH13*>();
    obsThFactory["log10_VV_phi3_WW_TH13"] = boost::factory<log10_VV_phi3_WW_TH13*>();
    obsThFactory["log10_ggVV_phi2_WW_lnulnu_TH13"] = boost::factory<log10_ggVV_phi2_WW_lnulnu_TH13*>();
    obsThFactory["log10_ggVV_phi3_WW_lnulnu_TH13"] = boost::factory<log10_ggVV_phi3_WW_lnulnu_TH13*>();
    obsThFactory["log10_pp_phi2_WW_TH13"] = boost::factory<log10_pp_phi2_WW_TH13*>();
    obsThFactory["log10_pp_phi3_WW_TH13"] = boost::factory<log10_pp_phi3_WW_TH13*>();

    obsThFactory["log10_mu_pp_phi2_VV_TH8"] = boost::factory<log10_mu_pp_phi2_VV_TH8*>();
    obsThFactory["log10_mu_pp_phi3_VV_TH8"] = boost::factory<log10_mu_pp_phi3_VV_TH8*>();
    obsThFactory["log10_pp_phi2_VV_TH13"] = boost::factory<log10_pp_phi2_VV_TH13*>();
    obsThFactory["log10_pp_phi3_VV_TH13"] = boost::factory<log10_pp_phi3_VV_TH13*>();

    obsThFactory["log10_gg_phi2_phi1phi1_TH8"] = boost::factory<log10_gg_phi2_phi1phi1_TH8*>();
    obsThFactory["log10_gg_phi3_phi1phi1_TH8"] = boost::factory<log10_gg_phi3_phi1phi1_TH8*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbbb_TH8"] = boost::factory<log10_pp_phi2_phi1phi1_bbbb_TH8*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbbb_TH8"] = boost::factory<log10_pp_phi3_phi1phi1_bbbb_TH8*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbgaga_TH8"] = boost::factory<log10_pp_phi2_phi1phi1_bbgaga_TH8*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbgaga_TH8"] = boost::factory<log10_pp_phi3_phi1phi1_bbgaga_TH8*>();
    obsThFactory["log10_gg_phi2_phi1phi1_bbtautau_TH8"] = boost::factory<log10_gg_phi2_phi1phi1_bbtautau_TH8*>();
    obsThFactory["log10_gg_phi3_phi1phi1_bbtautau_TH8"] = boost::factory<log10_gg_phi3_phi1phi1_bbtautau_TH8*>();
    obsThFactory["log10_pp_phi2_phi1phi1_TH8"] = boost::factory<log10_pp_phi2_phi1phi1_TH8*>();
    obsThFactory["log10_pp_phi3_phi1phi1_TH8"] = boost::factory<log10_pp_phi3_phi1phi1_TH8*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbbb_TH13"] = boost::factory<log10_pp_phi2_phi1phi1_bbbb_TH13*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbbb_TH13"] = boost::factory<log10_pp_phi3_phi1phi1_bbbb_TH13*>();
    obsThFactory["log10_pp_phi2_phi1phi1_TH13"] = boost::factory<log10_pp_phi2_phi1phi1_TH13*>();
    obsThFactory["log10_pp_phi3_phi1phi1_TH13"] = boost::factory<log10_pp_phi3_phi1phi1_TH13*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbgaga_TH13"] = boost::factory<log10_pp_phi2_phi1phi1_bbgaga_TH13*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbgaga_TH13"] = boost::factory<log10_pp_phi3_phi1phi1_bbgaga_TH13*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbtautau_TH13"] = boost::factory<log10_pp_phi2_phi1phi1_bbtautau_TH13*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbtautau_TH13"] = boost::factory<log10_pp_phi3_phi1phi1_bbtautau_TH13*>();
    obsThFactory["log10_pp_phi2_phi1phi1_bbVV_TH13"] = boost::factory<log10_pp_phi2_phi1phi1_bbVV_TH13*>();
    obsThFactory["log10_pp_phi3_phi1phi1_bbVV_TH13"] = boost::factory<log10_pp_phi3_phi1phi1_bbVV_TH13*>();
    obsThFactory["log10_gg_phi2_phi1phi1_gagaWW_TH13"] = boost::factory<log10_gg_phi2_phi1phi1_gagaWW_TH13*>();
    obsThFactory["log10_gg_phi3_phi1phi1_gagaWW_TH13"] = boost::factory<log10_gg_phi3_phi1phi1_gagaWW_TH13*>();

    obsThFactory["log10_gg_phi2_phi1Z_bbZ_TH8"] = boost::factory<log10_gg_phi2_phi1Z_bbZ_TH8*>();
    obsThFactory["log10_gg_phi3_phi1Z_bbZ_TH8"] = boost::factory<log10_gg_phi3_phi1Z_bbZ_TH8*>();
    obsThFactory["log10_gg_phi2_phi1Z_bbll_TH8"] = boost::factory<log10_gg_phi2_phi1Z_bbll_TH8*>();
    obsThFactory["log10_gg_phi3_phi1Z_bbll_TH8"] = boost::factory<log10_gg_phi3_phi1Z_bbll_TH8*>();
    obsThFactory["log10_gg_phi2_phi1Z_tautauZ_TH8"] = boost::factory<log10_gg_phi2_phi1Z_tautauZ_TH8*>();
    obsThFactory["log10_gg_phi3_phi1Z_tautauZ_TH8"] = boost::factory<log10_gg_phi3_phi1Z_tautauZ_TH8*>();
    obsThFactory["log10_gg_phi2_phi1Z_tautaull_TH8"] = boost::factory<log10_gg_phi2_phi1Z_tautaull_TH8*>();
    obsThFactory["log10_gg_phi3_phi1Z_tautaull_TH8"] = boost::factory<log10_gg_phi3_phi1Z_tautaull_TH8*>();
    obsThFactory["log10_gg_phi2_phi1Z_bbZ_TH13"] = boost::factory<log10_gg_phi2_phi1Z_bbZ_TH13*>();
    obsThFactory["log10_gg_phi3_phi1Z_bbZ_TH13"] = boost::factory<log10_gg_phi3_phi1Z_bbZ_TH13*>();
    obsThFactory["log10_bb_phi2_phi1Z_bbZ_TH13"] = boost::factory<log10_bb_phi2_phi1Z_bbZ_TH13*>();
    obsThFactory["log10_bb_phi3_phi1Z_bbZ_TH13"] = boost::factory<log10_bb_phi3_phi1Z_bbZ_TH13*>();

    obsThFactory["log10_pp_phi3_phi2Z_bbll_TH8"] = boost::factory<log10_pp_phi3_phi2Z_bbll_TH8*>();
    obsThFactory["log10_pp_phi3_phi2Z_tautaull_TH8"] = boost::factory<log10_pp_phi3_phi2Z_tautaull_TH8*>();
    obsThFactory["log10_gg_phi3_phi2Z_bbZ_TH13"] = boost::factory<log10_gg_phi3_phi2Z_bbZ_TH13*>();
    obsThFactory["log10_bb_phi3_phi2Z_bbZ_TH13"] = boost::factory<log10_bb_phi3_phi2Z_bbZ_TH13*>();

    obsThFactory["log10_pp_Hpm_taunu_TH8_GTHDM"] = boost::factory<log10_pp_Hpm_taunu_TH8_GTHDM*>();
    obsThFactory["log10_pp_Hp_taunu_TH8_GTHDM"] = boost::factory<log10_pp_Hp_taunu_TH8_GTHDM*>();
    obsThFactory["log10_pp_Hpm_taunu_TH13_GTHDM"] = boost::factory<log10_pp_Hpm_taunu_TH13_GTHDM*>();
    obsThFactory["log10_pp_Hpm_tb_TH8_GTHDM"] = boost::factory<log10_pp_Hpm_tb_TH8_GTHDM*>();
    obsThFactory["log10_pp_Hp_tb_TH8_GTHDM"] = boost::factory<log10_pp_Hp_tb_TH8_GTHDM*>();
    obsThFactory["log10_pp_Hpm_tb_TH13_GTHDM"] = boost::factory<log10_pp_Hpm_tb_TH13_GTHDM*>();

        /* END: REMOVE FROM THE PACKAGE */

    //-----  LeftRightSymmetric model observables  -----
    obsThFactory["mu1_2_LRSM"] = boost::factory<mu1_2_LRSM*>();
    obsThFactory["mu2_2_LRSM"] = boost::factory<mu2_2_LRSM*>();
    obsThFactory["mu3_2_LRSM"] = boost::factory<mu3_2_LRSM*>();
    obsThFactory["rho2_LRSM"] = boost::factory<rho2_LRSM*>();
    obsThFactory["rho3_LRSM"] = boost::factory<rho3_LRSM*>();
    obsThFactory["alpha3_LRSM"] = boost::factory<alpha3_LRSM*>();
    obsThFactory["mH00_LRSM"] = bind(boost::factory<MH0_LRSM*>(), _1, 0);
    obsThFactory["mH01_LRSM"] = bind(boost::factory<MH0_LRSM*>(), _1, 1);
    obsThFactory["mH02_LRSM"] = bind(boost::factory<MH0_LRSM*>(), _1, 2);
    obsThFactory["mH03_LRSM"] = bind(boost::factory<MH0_LRSM*>(), _1, 3);
    obsThFactory["mH04_LRSM"] = bind(boost::factory<MH0_LRSM*>(), _1, 4);
    obsThFactory["mH05_LRSM"] = boost::factory<MH05_LRSM*>();
    obsThFactory["mH06_LRSM"] = boost::factory<MH06_LRSM*>();
    obsThFactory["MH01_app1"] = boost::factory<MH01_app1*>();
    obsThFactory["MH01_app"] = bind(boost::factory<MH0_app*>(), _1, 0);
    obsThFactory["MH02_app"] = bind(boost::factory<MH0_app*>(), _1, 1);
    obsThFactory["MH03_app"] = bind(boost::factory<MH0_app*>(), _1, 2);
    obsThFactory["MH04_app"] = bind(boost::factory<MH0_app*>(), _1, 3);

    /* BEGIN: REMOVE FROM THE PACKAGE */
    //-----  THDMW model observables  -----
    obsThFactory["Q_stTHDMW"] = boost::factory<Q_stTHDMW*>();
    obsThFactory["DeltaQ_THDMW"] = boost::factory<DeltaQ_THDMW*>();
    obsThFactory["lambda1atQTHDMW"] = boost::factory<lambda1atQTHDMW*>();
    obsThFactory["lambda2atQTHDMW"] = boost::factory<lambda2atQTHDMW*>();
    obsThFactory["lambda3atQTHDMW"] = boost::factory<lambda3atQTHDMW*>();
    obsThFactory["lambda4atQTHDMW"] = boost::factory<lambda4atQTHDMW*>();
    obsThFactory["mu1atQTHDMW"] = boost::factory<mu1atQTHDMW*>();
    obsThFactory["mu2atQTHDMW"] = boost::factory<mu2atQTHDMW*>();
    obsThFactory["mu3atQTHDMW"] = boost::factory<mu3atQTHDMW*>();
    obsThFactory["mu4atQTHDMW"] = boost::factory<mu4atQTHDMW*>();
    obsThFactory["mu5atQTHDMW"] = boost::factory<mu5atQTHDMW*>();
    obsThFactory["mu6atQTHDMW"] = boost::factory<mu6atQTHDMW*>();
    obsThFactory["nu1atQTHDMW"] = boost::factory<nu1atQTHDMW*>();
    obsThFactory["omega1atQTHDMW"] = boost::factory<omega1atQTHDMW*>();
    obsThFactory["kappa1atQTHDMW"] = boost::factory<kappa1atQTHDMW*>();
    obsThFactory["nu2atQTHDMW"] = boost::factory<nu2atQTHDMW*>();
    obsThFactory["omega2atQTHDMW"] = boost::factory<omega2atQTHDMW*>();
    obsThFactory["kappa2atQTHDMW"] = boost::factory<kappa2atQTHDMW*>();
    obsThFactory["nu4atQTHDMW"] = boost::factory<nu4atQTHDMW*>();
    obsThFactory["omega4atQTHDMW"] = boost::factory<omega4atQTHDMW*>();
    obsThFactory["nu3atQTHDMW"] = boost::factory<nu3atQTHDMW*>();
    obsThFactory["nu5atQTHDMW"] = boost::factory<nu5atQTHDMW*>();
    //-----  Positivity constraints  -----
    obsThFactory["THDMWpositivity1"] = boost::factory<THDMWpositivity1*>();
    obsThFactory["THDMWpositivity2"] = boost::factory<THDMWpositivity2*>();
    obsThFactory["THDMWpositivity3"] = boost::factory<THDMWpositivity3*>();
    obsThFactory["THDMWpositivity4"] = boost::factory<THDMWpositivity4*>();
    obsThFactory["THDMWpositivity5"] = boost::factory<THDMWpositivity5*>();
    obsThFactory["THDMWpositivity6"] = boost::factory<THDMWpositivity6*>();
    obsThFactory["THDMWpositivity7"] = boost::factory<THDMWpositivity7*>();
    obsThFactory["THDMWpositivity8"] = boost::factory<THDMWpositivity8*>();
    obsThFactory["THDMWpositivity9"] = boost::factory<THDMWpositivity9*>();
    obsThFactory["THDMWpositivity10"] = boost::factory<THDMWpositivity10*>();
    obsThFactory["THDMWpositivity11"] = boost::factory<THDMWpositivity11*>();
    obsThFactory["THDMWpositivity12"] = boost::factory<THDMWpositivity12*>();
    obsThFactory["THDMWpositiveMassSquares"] = boost::factory<THDMWpositiveMassSquares*>();
    //-----  Tree-level unitarity constraints  -----
    obsThFactory["THDMWunitarity1"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 0);
    obsThFactory["THDMWunitarity2"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 1);
    obsThFactory["THDMWunitarity3"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 2);
    obsThFactory["THDMWunitarity4"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 3);
    obsThFactory["THDMWunitarity5"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 4);
    obsThFactory["THDMWunitarity6"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 5);
    obsThFactory["THDMWunitarity7"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 6);
    obsThFactory["THDMWunitarity8"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 7);
    obsThFactory["THDMWunitarity9"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 8);
    obsThFactory["THDMWunitarity10"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 9);
    obsThFactory["THDMWunitarity11"] = bind(boost::factory<THDMWunitarityLO*>(), _1, 10);
    //-----  One-loop unitarity constraints  -----
    obsThFactory["THDMWNLOunitarity1"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 0);
    obsThFactory["THDMWNLOunitarity2"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 1);
    obsThFactory["THDMWNLOunitarity3"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 2);
    obsThFactory["THDMWNLOunitarity4"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 3);
    obsThFactory["THDMWNLOunitarity5"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 4);
    obsThFactory["THDMWNLOunitarity6"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 5);
    obsThFactory["THDMWNLOunitarity7"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 6);
    obsThFactory["THDMWNLOunitarity8"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 7);
    obsThFactory["THDMWNLOunitarity9"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 8);
    obsThFactory["THDMWNLOunitarity10"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 9);
    obsThFactory["THDMWNLOunitarity11"] = bind(boost::factory<THDMWunitarityNLO*>(), _1, 10);
    //-----  One-loop "plus" unitarity constraints  -----
    obsThFactory["THDMWNLOpunitarity1"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 0);
    obsThFactory["THDMWNLOpunitarity2"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 1);
    obsThFactory["THDMWNLOpunitarity3"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 2);
    obsThFactory["THDMWNLOpunitarity4"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 3);
    obsThFactory["THDMWNLOpunitarity5"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 4);
    obsThFactory["THDMWNLOpunitarity6"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 5);
    obsThFactory["THDMWNLOpunitarity7"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 6);
    obsThFactory["THDMWNLOpunitarity8"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 7);
    obsThFactory["THDMWNLOpunitarity9"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 8);
    obsThFactory["THDMWNLOpunitarity10"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 9);
    obsThFactory["THDMWNLOpunitarity11"] = bind(boost::factory<THDMWunitarityNLOp*>(), _1, 10);
    //-----   R' criteria for perturbative unitarity  -----
    obsThFactory["THDMWunitarityRp1"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 0);
    obsThFactory["THDMWunitarityRp2"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 1);
    obsThFactory["THDMWunitarityRp3"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 2);
    obsThFactory["THDMWunitarityRp4"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 3);
    obsThFactory["THDMWunitarityRp5"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 4);
    obsThFactory["THDMWunitarityRp6"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 5);
    obsThFactory["THDMWunitarityRp7"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 6);
    obsThFactory["THDMWunitarityRp8"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 7);
    obsThFactory["THDMWunitarityRp9"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 8);
    obsThFactory["THDMWunitarityRp10"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 9);
    obsThFactory["THDMWunitarityRp11"] = bind(boost::factory<THDMWunitarityRp*>(), _1, 10);
    //-----   Physical parameters  -----
    obsThFactory["m12sqTHDMW"] = boost::factory<m12sqTHDMW*>();
    obsThFactory["m11sqTHDMW"] = boost::factory<m11sqTHDMW*>();
    obsThFactory["m22sqTHDMW"] = boost::factory<m22sqTHDMW*>();
    obsThFactory["mhsqTHDMW"] = boost::factory<mhsqTHDMW*>();
    obsThFactory["mhTHDMW"] = boost::factory<mhTHDMW*>();
    obsThFactory["mHHsqTHDMW"] = boost::factory<mHHsqTHDMW*>();
    obsThFactory["mHHTHDMW"] = boost::factory<mHHTHDMW*>();
    obsThFactory["mAsqTHDMW"] = boost::factory<mAsqTHDMW*>();
    obsThFactory["mATHDMW"] = boost::factory<mATHDMW*>();
    obsThFactory["mSpsqTHDMW"] = boost::factory<mSpsqTHDMW*>();
    obsThFactory["mSpTHDMW"] = boost::factory<mSpTHDMW*>();
    obsThFactory["mSRsqTHDMW"] = boost::factory<mSRsqTHDMW*>();
    obsThFactory["mSRTHDMW"] = boost::factory<mSRTHDMW*>();
    obsThFactory["mSIsqTHDMW"] = boost::factory<mSIsqTHDMW*>();
    obsThFactory["mSITHDMW"] = boost::factory<mSITHDMW*>();
    obsThFactory["mAmmHH_THDMW"] = boost::factory<mAmmHH_THDMW*>();
    obsThFactory["mHHmmA_THDMW"] = boost::factory<mHHmmA_THDMW*>();
    obsThFactory["mAmmSR_THDMW"] = boost::factory<mAmmSR_THDMW*>();
    obsThFactory["mSRmmA_THDMW"] = boost::factory<mSRmmA_THDMW*>();
    obsThFactory["mAmmSI_THDMW"] = boost::factory<mAmmSI_THDMW*>();
    obsThFactory["mSImmA_THDMW"] = boost::factory<mSImmA_THDMW*>();
    obsThFactory["mHHmmSR_THDMW"] = boost::factory<mHHmmSR_THDMW*>();
    obsThFactory["mSRmmHH_THDMW"] = boost::factory<mSRmmHH_THDMW*>();
    obsThFactory["mHHmmSI_THDMW"] = boost::factory<mHHmmSI_THDMW*>();
    obsThFactory["mSImmHH_THDMW"] = boost::factory<mSImmHH_THDMW*>();
    obsThFactory["mSRmmSI_THDMW"] = boost::factory<mSRmmSI_THDMW*>();
    obsThFactory["mSImmSR_THDMW"] = boost::factory<mSImmSR_THDMW*>();
    obsThFactory["mSpmmSI_THDMW"] = boost::factory<mSpmmSI_THDMW*>();
    obsThFactory["mSpmmSR_THDMW"] = boost::factory<mSpmmSR_THDMW*>();
    obsThFactory["mSRmmSp_THDMW"] = boost::factory<mSRmmSp_THDMW*>();
    obsThFactory["mSImmSp_THDMW"] = boost::factory<mSImmSp_THDMW*>();
    //-----   Higgs observables -----
    obsThFactory["rh_gg_THDMW"] = boost::factory<rh_gg_THDMW*>();
    obsThFactory["rh_gaga_THDMW"] = boost::factory<rh_gaga_THDMW*>();
    obsThFactory["rh_Zga_THDMW"] = boost::factory<rh_Zga_THDMW*>();
    //-----   Direct Searches -----
    obsThFactory["Hobs_pp_Sr_tt_ATLAS13"] = boost::factory<Hobs_pp_Sr_tt_ATLAS13*>();
    obsThFactory["log10_pp_Sr_tt_TH13"] = boost::factory<log10_pp_Sr_tt_TH13*>();
    obsThFactory["Hobs_pp_Srtt_tttt_ATLAS13"] = boost::factory<Hobs_pp_Srtt_tttt_ATLAS13*>();
    obsThFactory["log10_pp_Srtt_tttt_TH13"] = boost::factory<log10_pp_Srtt_tttt_TH13*>();
    obsThFactory["Hobs_pp_Sr_jj_CMS13"] = boost::factory<Hobs_pp_Sr_jj_CMS13*>();
    obsThFactory["log10_pp_Sr_jj_TH13"] = boost::factory<log10_pp_Sr_jj_TH13*>();
    obsThFactory["Hobs_pp_SrSr_jjjj_ATLAS13"] = boost::factory<Hobs_pp_SrSr_jjjj_ATLAS13*>();
    obsThFactory["log10_pp_SrSr_jjjj_TH13"] = boost::factory<log10_pp_SrSr_jjjj_TH13*>();
    obsThFactory["Hobs_pp_Stb_tbtb_ATLAS13"] = boost::factory<Hobs_pp_Stb_tbtb_ATLAS13*>();
    obsThFactory["log10_pp_Stb_tbtb_TH13"] = boost::factory<log10_pp_Stb_tbtb_TH13*>();
    obsThFactory["Hobs_pp_Sitt_tttt_ATLAS13"] = boost::factory<Hobs_pp_Sitt_tttt_ATLAS13*>();
    obsThFactory["log10_pp_Sitt_tttt_TH13"] = boost::factory<log10_pp_Sitt_tttt_TH13*>();
    obsThFactory["Hobs_pp_Srbb_bbbb_CMS13"] = boost::factory<Hobs_pp_Srbb_bbbb_CMS13*>();
    obsThFactory["log10_pp_Srbb_bbbb_TH13"] = boost::factory<log10_pp_Srbb_bbbb_TH13*>();
    obsThFactory["Hobs_pp_Srbb_bbbb_CMS8"] = boost::factory<Hobs_pp_Srbb_bbbb_CMS8*>();
    obsThFactory["log10_pp_Srbb_bbbb_TH8"] = boost::factory<log10_pp_Srbb_bbbb_TH8*>();
    obsThFactory["Hobs_pp_Sibb_bbbb_CMS13"] = boost::factory<Hobs_pp_Sibb_bbbb_CMS13*>();
    obsThFactory["log10_pp_Sibb_bbbb_TH13"] = boost::factory<log10_pp_Sibb_bbbb_TH13*>();
    obsThFactory["Hobs_pp_Sibb_bbbb_CMS8"] = boost::factory<Hobs_pp_Sibb_bbbb_CMS8*>();
    obsThFactory["log10_pp_Sibb_bbbb_TH8"] = boost::factory<log10_pp_Sibb_bbbb_TH8*>();
    obsThFactory["Hobs_pp_Sr_bb_CMS13"] = boost::factory<Hobs_pp_Sr_bb_CMS13*>();
    obsThFactory["log10_pp_Sr_bb_TH13"] = boost::factory<log10_pp_Sr_bb_TH13*>();
    obsThFactory["Hobs_pp_Sr_bb_CMS8"] = boost::factory<Hobs_pp_Sr_bb_CMS8*>();
    obsThFactory["log10_pp_Sr_bb_TH8"] = boost::factory<log10_pp_Sr_bb_TH8*>();
    obsThFactory["Hobs_pp_Si_bb_CMS13"] = boost::factory<Hobs_pp_Si_bb_CMS13*>();
    obsThFactory["log10_pp_Si_bb_TH13"] = boost::factory<log10_pp_Si_bb_TH13*>();
    obsThFactory["log10_pp_Si_bb_TH8"] = boost::factory<log10_pp_Si_bb_TH8*>();
    obsThFactory["Hobs_pp_Si_bb_CMS8"] = boost::factory<Hobs_pp_Si_bb_CMS8*>();

    //obsThFactory["logpp_SrSr_jjjj_TH13"] = boost::factory<logpp_SrSr_jjjj_TH13*>();
    //-----        EWPO       ------
    obsThFactory["Rb0_THDMW"] = boost::factory<Rb0THDMW*>();
    //-----        STU       ------
    obsThFactory["DeltaS_THDMW"] = boost::factory<THDMWDeltaS*>();
    obsThFactory["DeltaT_THDMW"] = boost::factory<THDMWDeltaT*>();
    obsThFactory["DeltaU_THDMW"] = boost::factory<THDMWDeltaU*>();

    /* END: REMOVE FROM THE PACKAGE */

    //-----  GeorgiMachacek observables  -----
    //-----  GeorgiMachacek quantities -----
    obsThFactory["tanbetaGM"] = boost::factory<tanbetaGM*>();
    obsThFactory["m1sqGM"] = boost::factory<m1sqGM*>();
    obsThFactory["m2sqGM"] = boost::factory<m2sqGM*>();
    obsThFactory["lambda1GM"] = boost::factory<lambda1GM*>();
    obsThFactory["lambda2GM"] = boost::factory<lambda2GM*>();
    obsThFactory["lambda3GM"] = boost::factory<lambda3GM*>();
    obsThFactory["lambda4GM"] = boost::factory<lambda4GM*>();
    obsThFactory["lambda5GM"] = boost::factory<lambda5GM*>();
    obsThFactory["vPhiGM"] = boost::factory<vPhiGM*>();
    obsThFactory["GMmHh"] = boost::factory<GMmass_mHh*>();
    obsThFactory["GMmA"] = boost::factory<GMmass_mA*>();
    obsThFactory["GMmH5"] = boost::factory<GMmass_mH5*>();
    obsThFactory["GMmHlmmHh"] = boost::factory<GMmassdifference_mHlmmHh*>();
    obsThFactory["GMmHhmmHl"] = boost::factory<GMmassdifference_mHhmmHl*>();
    obsThFactory["GMmHlmmA"] = boost::factory<GMmassdifference_mHlmmA*>();
    obsThFactory["GMmAmmHl"] = boost::factory<GMmassdifference_mAmmHl*>();
    obsThFactory["GMmHlmmH5"] = boost::factory<GMmassdifference_mHlmmH5*>();
    obsThFactory["GMmH5mmHl"] = boost::factory<GMmassdifference_mH5mmHl*>();
    obsThFactory["GMmHhmmA"] = boost::factory<GMmassdifference_mHhmmA*>();
    obsThFactory["GMmAmmHh"] = boost::factory<GMmassdifference_mAmmHh*>();
    obsThFactory["GMmHhmmH5"] = boost::factory<GMmassdifference_mHhmmH5*>();
    obsThFactory["GMmH5mmHh"] = boost::factory<GMmassdifference_mH5mmHh*>();
    obsThFactory["GMmAmmH5"] = boost::factory<GMmassdifference_mAmmH5*>();
    obsThFactory["GMmH5mmA"] = boost::factory<GMmassdifference_mH5mmA*>();
    obsThFactory["GMGammah"] = boost::factory<GMGammah*>();
    obsThFactory["GMGammaH1"] = boost::factory<GMGammaH1*>();
    obsThFactory["GMGammaH3"] = boost::factory<GMGammaH3*>();
    obsThFactory["GMGammaH3p"] = boost::factory<GMGammaH3p*>();
    obsThFactory["GMGammaH5"] = boost::factory<GMGammaH5*>();
    obsThFactory["GMGammaH5p"] = boost::factory<GMGammaH5p*>();
    obsThFactory["GMGammaH5pp"] = boost::factory<GMGammaH5pp*>();
    obsThFactory["GMghhh"] = boost::factory<GMghhh*>();
    //-----  Tree-level unitarity constraints  -----
    obsThFactory["GMunitarity1"] = bind(boost::factory<GMunitarityLO*>(), _1, 0);
    obsThFactory["GMunitarity2"] = bind(boost::factory<GMunitarityLO*>(), _1, 1);
    obsThFactory["GMunitarity3"] = bind(boost::factory<GMunitarityLO*>(), _1, 2);
    obsThFactory["GMunitarity4"] = bind(boost::factory<GMunitarityLO*>(), _1, 3);
    obsThFactory["GMunitarity5"] = bind(boost::factory<GMunitarityLO*>(), _1, 4);
    obsThFactory["GMunitarity6"] = bind(boost::factory<GMunitarityLO*>(), _1, 5);
    obsThFactory["GMunitarity7"] = bind(boost::factory<GMunitarityLO*>(), _1, 6);
    obsThFactory["GMunitarity8"] = bind(boost::factory<GMunitarityLO*>(), _1, 7);
    obsThFactory["GMunitarity9"] = bind(boost::factory<GMunitarityLO*>(), _1, 8);
    obsThFactory["GMunitarity10"] = bind(boost::factory<GMunitarityLO*>(), _1, 9);
    obsThFactory["GMunitarity11"] = bind(boost::factory<GMunitarityLO*>(), _1, 10);
    obsThFactory["GMunitarity12"] = bind(boost::factory<GMunitarityLO*>(), _1, 11);
    obsThFactory["GMunitarity13"] = bind(boost::factory<GMunitarityLO*>(), _1, 12);
    obsThFactory["GMunitarity14"] = bind(boost::factory<GMunitarityLO*>(), _1, 13);
    obsThFactory["GMunitarity15"] = bind(boost::factory<GMunitarityLO*>(), _1, 14);
    obsThFactory["GMunitarity16"] = bind(boost::factory<GMunitarityLO*>(), _1, 15);
    obsThFactory["GMunitarity17"] = bind(boost::factory<GMunitarityLO*>(), _1, 16);
    //-----  Positivity constraints  -----
    obsThFactory["GMpositivity1"] = boost::factory<GMpositivity1*>();
    obsThFactory["GMpositivity2"] = boost::factory<GMpositivity2*>();
    obsThFactory["GMpositivity3"] = boost::factory<GMpositivity3*>();
    obsThFactory["GMpositivity4"] = boost::factory<GMpositivity4*>();
    //-----  Higgs observables  -----
    obsThFactory["rh_gaga_GM"] = boost::factory<rh_gaga_GM*>();
    obsThFactory["rh_Zga_GM"] = boost::factory<rh_Zga_GM*>();
    //-----  Direct Higgs searches -----
    obsThFactory["BR_H1_tt_GM"] = boost::factory<BR_H1_tt_GM*>();
    obsThFactory["BR_H1_bb_GM"] = boost::factory<BR_H1_bb_GM*>();
    obsThFactory["BR_H1_tautau_GM"] = boost::factory<BR_H1_tautau_GM*>();
    obsThFactory["BR_H1_WW_GM"] = boost::factory<BR_H1_WW_GM*>();
    obsThFactory["BR_H1_ZZ_GM"] = boost::factory<BR_H1_ZZ_GM*>();
    obsThFactory["BR_H1_gaga_GM"] = boost::factory<BR_H1_gaga_GM*>();
    obsThFactory["BR_H1_Zga_GM"] = boost::factory<BR_H1_Zga_GM*>();
    obsThFactory["BR_H1_H3Z_GM"] = boost::factory<BR_H1_H3Z_GM*>();
    obsThFactory["BR_H1_H3pW_GM"] = boost::factory<BR_H1_H3pW_GM*>();
    obsThFactory["BR_H1_hh_GM"] = boost::factory<BR_H1_hh_GM*>();
    obsThFactory["BR_H1_H3H3_GM"] = boost::factory<BR_H1_H3H3_GM*>();
    obsThFactory["BR_H1_H3pH3m_GM"] = boost::factory<BR_H1_H3pH3m_GM*>();
    obsThFactory["BR_H1_H5H5_GM"] = boost::factory<BR_H1_H5H5_GM*>();
    obsThFactory["BR_H1_H5pH5m_GM"] = boost::factory<BR_H1_H5pH5m_GM*>();
    obsThFactory["BR_H1_H5ppH5mm_GM"] = boost::factory<BR_H1_H5ppH5mm_GM*>();
    obsThFactory["BR_H3_tt_GM"] = boost::factory<BR_H3_tt_GM*>();
    obsThFactory["BR_H3_bb_GM"] = boost::factory<BR_H3_bb_GM*>();
    obsThFactory["BR_H3_tautau_GM"] = boost::factory<BR_H3_tautau_GM*>();
    obsThFactory["BR_H3_gaga_GM"] = boost::factory<BR_H3_gaga_GM*>();
    obsThFactory["BR_H3_Zga_GM"] = boost::factory<BR_H3_Zga_GM*>();
    obsThFactory["BR_H3_hZ_GM"] = boost::factory<BR_H3_hZ_GM*>();
    obsThFactory["BR_H3_H1Z_GM"] = boost::factory<BR_H3_H1Z_GM*>();
    obsThFactory["BR_H3_H5Z_GM"] = boost::factory<BR_H3_H5Z_GM*>();
    obsThFactory["BR_H3_H5pW_GM"] = boost::factory<BR_H3_H5pW_GM*>();
    obsThFactory["BR_H3p_taunu_GM"] = boost::factory<BR_H3p_taunu_GM*>();
    obsThFactory["BR_H3p_tb_GM"] = boost::factory<BR_H3p_tb_GM*>();
    obsThFactory["BR_H3p_hW_GM"] = boost::factory<BR_H3p_hW_GM*>();
    obsThFactory["BR_H3p_H1W_GM"] = boost::factory<BR_H3p_H1W_GM*>();
    obsThFactory["BR_H3p_H5W_GM"] = boost::factory<BR_H3p_H5W_GM*>();
    obsThFactory["BR_H3p_H5pZ_GM"] = boost::factory<BR_H3p_H5pZ_GM*>();
    obsThFactory["BR_H3p_H5ppW_GM"] = boost::factory<BR_H3p_H5ppW_GM*>();
    obsThFactory["BR_H5_WW_GM"] = boost::factory<BR_H5_WW_GM*>();
    obsThFactory["BR_H5_ZZ_GM"] = boost::factory<BR_H5_ZZ_GM*>();
    obsThFactory["BR_H5_gaga_GM"] = boost::factory<BR_H5_gaga_GM*>();
    obsThFactory["BR_H5_Zga_GM"] = boost::factory<BR_H5_Zga_GM*>();
    obsThFactory["BR_H5_H3Z_GM"] = boost::factory<BR_H5_H3Z_GM*>();
    obsThFactory["BR_H5_H3pW_GM"] = boost::factory<BR_H5_H3pW_GM*>();
    obsThFactory["BR_H5_H3H3_GM"] = boost::factory<BR_H5_H3H3_GM*>();
    obsThFactory["BR_H5_H3pH3m_GM"] = boost::factory<BR_H5_H3pH3m_GM*>();
    obsThFactory["BR_H5p_WZ_GM"] = boost::factory<BR_H5p_WZ_GM*>();
    obsThFactory["BR_H5p_H3W_GM"] = boost::factory<BR_H5p_H3W_GM*>();
    obsThFactory["BR_H5p_H3pZ_GM"] = boost::factory<BR_H5p_H3pZ_GM*>();
    obsThFactory["BR_H5p_H3pH3_GM"] = boost::factory<BR_H5p_H3pH3_GM*>();
    obsThFactory["BR_H5pp_WW_GM"] = boost::factory<BR_H5pp_WW_GM*>();
    obsThFactory["BR_H5pp_H3pW_GM"] = boost::factory<BR_H5pp_H3pW_GM*>();
    obsThFactory["BR_H5pp_H3pH3p_GM"] = boost::factory<BR_H5pp_H3pH3p_GM*>();
    obsThFactory["Hobs_tt_H1_tt_ATLAS13"] = boost::factory<Hobs_tt_H1_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_H1_tt_ATLAS13"] = boost::factory<Hobs_bb_H1_tt_ATLAS13*>();
    obsThFactory["Hobs_tt_H3_tt_ATLAS13"] = boost::factory<Hobs_tt_H3_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_H3_tt_ATLAS13"] = boost::factory<Hobs_bb_H3_tt_ATLAS13*>();
    obsThFactory["Hobs_bb_H1_bb_CMS8"] = boost::factory<Hobs_bb_H1_bb_CMS8*>();
    obsThFactory["Hobs_gg_H1_bb_CMS8"] = boost::factory<Hobs_gg_H1_bb_CMS8*>();
    obsThFactory["Hobs_pp_H1_bb_CMS13"] = boost::factory<Hobs_pp_H1_bb_CMS13*>();
    obsThFactory["Hobs_bb_H1_bb_CMS13"] = boost::factory<Hobs_bb_H1_bb_CMS13*>();
    obsThFactory["Hobs_bb_H3_bb_CMS8"] = boost::factory<Hobs_bb_H3_bb_CMS8*>();
    obsThFactory["Hobs_gg_H3_bb_CMS8"] = boost::factory<Hobs_gg_H3_bb_CMS8*>();
    obsThFactory["Hobs_pp_H3_bb_CMS13"] = boost::factory<Hobs_pp_H3_bb_CMS13*>();
    obsThFactory["Hobs_bb_H3_bb_CMS13"] = boost::factory<Hobs_bb_H3_bb_CMS13*>();
    obsThFactory["Hobs_gg_H1_tautau_CMS8"] = boost::factory<Hobs_gg_H1_tautau_CMS8*>();
    obsThFactory["Hobs_bb_H1_tautau_CMS8"] = boost::factory<Hobs_bb_H1_tautau_CMS8*>();
    obsThFactory["Hobs_gg_H1_tautau_ATLAS13"] = boost::factory<Hobs_gg_H1_tautau_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_tautau_CMS13"] = boost::factory<Hobs_gg_H1_tautau_CMS13*>();
    obsThFactory["Hobs_bb_H1_tautau_ATLAS13"] = boost::factory<Hobs_bb_H1_tautau_ATLAS13*>();
    obsThFactory["Hobs_bb_H1_tautau_CMS13"] = boost::factory<Hobs_bb_H1_tautau_CMS13*>();
    obsThFactory["Hobs_gg_H1_tautau_ATLAS8"] = boost::factory<Hobs_gg_H1_tautau_ATLAS8*>();
    obsThFactory["Hobs_bb_H1_tautau_ATLAS8"] = boost::factory<Hobs_bb_H1_tautau_ATLAS8*>();
    obsThFactory["Hobs_gg_H3_tautau_ATLAS8"] = boost::factory<Hobs_gg_H3_tautau_ATLAS8*>();
    obsThFactory["Hobs_gg_H3_tautau_CMS8"] = boost::factory<Hobs_gg_H3_tautau_CMS8*>();
    obsThFactory["Hobs_bb_H3_tautau_ATLAS8"] = boost::factory<Hobs_bb_H3_tautau_ATLAS8*>();
    obsThFactory["Hobs_bb_H3_tautau_CMS8"] = boost::factory<Hobs_bb_H3_tautau_CMS8*>();
    obsThFactory["Hobs_gg_H3_tautau_ATLAS13"] = boost::factory<Hobs_gg_H3_tautau_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_tautau_CMS13"] = boost::factory<Hobs_gg_H3_tautau_CMS13*>();
    obsThFactory["Hobs_bb_H3_tautau_ATLAS13"] = boost::factory<Hobs_bb_H3_tautau_ATLAS13*>();
    obsThFactory["Hobs_bb_H3_tautau_CMS13"] = boost::factory<Hobs_bb_H3_tautau_CMS13*>();
    obsThFactory["Hobs_gg_H1_gaga_ATLAS8"] = boost::factory<Hobs_gg_H1_gaga_ATLAS8*>();
    obsThFactory["Hobs_pp_H1_gaga_ATLAS13"] = boost::factory<Hobs_pp_H1_gaga_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_gaga_CMS13"] = boost::factory<Hobs_gg_H1_gaga_CMS13*>();
    obsThFactory["Hobs_gg_H3_gaga_ATLAS8"] = boost::factory<Hobs_gg_H3_gaga_ATLAS8*>();
    obsThFactory["Hobs_pp_H3_gaga_ATLAS13"] = boost::factory<Hobs_pp_H3_gaga_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_gaga_CMS13"] = boost::factory<Hobs_gg_H3_gaga_CMS13*>();
    obsThFactory["Hobs_pp_H5_gaga_ATLAS13"] = boost::factory<Hobs_pp_H5_gaga_ATLAS13*>();
    obsThFactory["Hobs_pp_H1_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_H1_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_H1_Zga_llga_CMS8"] = boost::factory<Hobs_pp_H1_Zga_llga_CMS8*>();
    obsThFactory["Hobs_gg_H1_Zga_llga_ATLAS13"] = boost::factory<Hobs_gg_H1_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_Zga_qqga_ATLAS13"] = boost::factory<Hobs_gg_H1_Zga_qqga_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_Zga_CMS13"] = boost::factory<Hobs_gg_H1_Zga_CMS13*>();
    obsThFactory["Hobs_pp_H3_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_H3_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_H3_Zga_llga_CMS8"] = boost::factory<Hobs_pp_H3_Zga_llga_CMS8*>();
    obsThFactory["Hobs_gg_H3_Zga_llga_ATLAS13"] = boost::factory<Hobs_gg_H3_Zga_llga_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_Zga_qqga_ATLAS13"] = boost::factory<Hobs_gg_H3_Zga_qqga_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_Zga_CMS13"] = boost::factory<Hobs_gg_H3_Zga_CMS13*>();
    obsThFactory["Hobs_pp_H5_Zga_llga_ATLAS8"] = boost::factory<Hobs_pp_H5_Zga_llga_ATLAS8*>();
    obsThFactory["Hobs_pp_H5_Zga_llga_CMS8"] = boost::factory<Hobs_pp_H5_Zga_llga_CMS8*>();
    obsThFactory["Hobs_gg_H1_ZZ_ATLAS8"] = boost::factory<Hobs_gg_H1_ZZ_ATLAS8*>();
    obsThFactory["Hobs_VV_H1_ZZ_ATLAS8"] = boost::factory<Hobs_VV_H1_ZZ_ATLAS8*>();
    obsThFactory["Hobs_gg_H1_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_gg_H1_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_VV_H1_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_VV_H1_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_ZZ_qqllnunu_ATLAS13"] = boost::factory<Hobs_gg_H1_ZZ_qqllnunu_ATLAS13*>();
    obsThFactory["Hobs_VV_H1_ZZ_qqllnunu_ATLAS13"] = boost::factory<Hobs_VV_H1_ZZ_qqllnunu_ATLAS13*>();
    obsThFactory["Hobs_pp_H1_ZZ_llqqnunull_CMS13"] = boost::factory<Hobs_pp_H1_ZZ_llqqnunull_CMS13*>();
    obsThFactory["Hobs_pp_H1_ZZ_qqnunu_CMS13"] = boost::factory<Hobs_pp_H1_ZZ_qqnunu_CMS13*>();
    obsThFactory["Hobs_VV_H5_ZZ_ATLAS8"] = boost::factory<Hobs_VV_H5_ZZ_ATLAS8*>();
    obsThFactory["Hobs_VV_H5_ZZ_llllnunu_ATLAS13"] = boost::factory<Hobs_VV_H5_ZZ_llllnunu_ATLAS13*>();
    obsThFactory["Hobs_VV_H5_ZZ_qqllnunu_ATLAS13"] = boost::factory<Hobs_VV_H5_ZZ_qqllnunu_ATLAS13*>();
    obsThFactory["Hobs_pp_H5_ZZ_llqqnunull_CMS13"] = boost::factory<Hobs_pp_H5_ZZ_llqqnunull_CMS13*>();
    obsThFactory["Hobs_pp_H5_ZZ_qqnunu_CMS13"] = boost::factory<Hobs_pp_H5_ZZ_qqnunu_CMS13*>();
    obsThFactory["Hobs_gg_H1_WW_ATLAS8"] = boost::factory<Hobs_gg_H1_WW_ATLAS8*>();
    obsThFactory["Hobs_VV_H1_WW_ATLAS8"] = boost::factory<Hobs_VV_H1_WW_ATLAS8*>();
    obsThFactory["Hobs_gg_H1_WW_enumunu_ATLAS13"] = boost::factory<Hobs_gg_H1_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VV_H1_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VV_H1_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_gg_H1_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_VV_H1_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_VV_H1_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_ggVV_H1_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggVV_H1_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_pp_H1_WW_lnuqq_CMS13"] = boost::factory<Hobs_pp_H1_WW_lnuqq_CMS13*>();
    obsThFactory["Hobs_VV_H5_WW_ATLAS8"] = boost::factory<Hobs_VV_H5_WW_ATLAS8*>();
    obsThFactory["Hobs_VV_H5_WW_enumunu_ATLAS13"] = boost::factory<Hobs_VV_H5_WW_enumunu_ATLAS13*>();
    obsThFactory["Hobs_VV_H5_WW_lnuqq_ATLAS13"] = boost::factory<Hobs_VV_H5_WW_lnuqq_ATLAS13*>();
    obsThFactory["Hobs_ggVV_H5_WW_lnulnu_CMS13"] = boost::factory<Hobs_ggVV_H5_WW_lnulnu_CMS13*>();
    obsThFactory["Hobs_pp_H5_WW_lnuqq_CMS13"] = boost::factory<Hobs_pp_H5_WW_lnuqq_CMS13*>();
    obsThFactory["Hobs_mu_pp_H1_VV_CMS8"] = boost::factory<Hobs_mu_pp_H1_VV_CMS8*>();
    obsThFactory["Hobs_pp_H1_VV_qqqq_ATLAS13"] = boost::factory<Hobs_pp_H1_VV_qqqq_ATLAS13*>();
    obsThFactory["Hobs_mu_pp_H5_VV_CMS8"] = boost::factory<Hobs_mu_pp_H5_VV_CMS8*>();
    obsThFactory["Hobs_pp_H5_VV_qqqq_ATLAS13"] = boost::factory<Hobs_pp_H5_VV_qqqq_ATLAS13*>();
    obsThFactory["Hobs_gg_H1_hh_ATLAS8"] = boost::factory<Hobs_gg_H1_hh_ATLAS8*>();
    obsThFactory["Hobs_pp_H1_hh_bbbb_CMS8"] = boost::factory<Hobs_pp_H1_hh_bbbb_CMS8*>();
    obsThFactory["Hobs_pp_H1_hh_gagabb_CMS8"] = boost::factory<Hobs_pp_H1_hh_gagabb_CMS8*>();
    obsThFactory["Hobs_gg_H1_hh_bbtautau_CMS8"] = boost::factory<Hobs_gg_H1_hh_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_H1_hh_bbtautau_CMS8"] = boost::factory<Hobs_pp_H1_hh_bbtautau_CMS8*>();
    obsThFactory["Hobs_pp_H1_hh_bbbb_ATLAS13"] = boost::factory<Hobs_pp_H1_hh_bbbb_ATLAS13*>();
    obsThFactory["Hobs_pp_H1_hh_bbbb_1_CMS13"] = boost::factory<Hobs_pp_H1_hh_bbbb_1_CMS13*>();
    obsThFactory["Hobs_pp_H1_hh_bbbb_2_CMS13"] = boost::factory<Hobs_pp_H1_hh_bbbb_2_CMS13*>();
    obsThFactory["Hobs_gg_H1_hh_bbbb_CMS13"] = boost::factory<Hobs_gg_H1_hh_bbbb_CMS13*>();
    obsThFactory["Hobs_pp_H1_hh_gagabb_ATLAS13"] = boost::factory<Hobs_pp_H1_hh_gagabb_ATLAS13*>();
    obsThFactory["Hobs_pp_H1_hh_gagabb_CMS13"] = boost::factory<Hobs_pp_H1_hh_gagabb_CMS13*>();
    obsThFactory["Hobs_pp_H1_hh_bbtautau_ATLAS13"] = boost::factory<Hobs_pp_H1_hh_bbtautau_ATLAS13*>();
    obsThFactory["Hobs_pp_H1_hh_bbtautau_1_CMS13"] = boost::factory<Hobs_pp_H1_hh_bbtautau_1_CMS13*>();
    obsThFactory["Hobs_pp_H1_hh_bbtautau_2_CMS13"] = boost::factory<Hobs_pp_H1_hh_bbtautau_2_CMS13*>();
    obsThFactory["Hobs_pp_H1_hh_bblnulnu_CMS13"] = boost::factory<Hobs_pp_H1_hh_bblnulnu_CMS13*>();
    obsThFactory["Hobs_gg_H1_hh_gagaWW_ATLAS13"] = boost::factory<Hobs_gg_H1_hh_gagaWW_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_hZ_bbZ_ATLAS8"] = boost::factory<Hobs_gg_H3_hZ_bbZ_ATLAS8*>();
    obsThFactory["Hobs_gg_H3_hZ_bbll_CMS8"] = boost::factory<Hobs_gg_H3_hZ_bbll_CMS8*>();
    obsThFactory["Hobs_gg_H3_hZ_tautauZ_ATLAS8"] = boost::factory<Hobs_gg_H3_hZ_tautauZ_ATLAS8*>();
    obsThFactory["Hobs_gg_H3_hZ_tautaull_CMS8"] = boost::factory<Hobs_gg_H3_hZ_tautaull_CMS8*>();
    obsThFactory["Hobs_gg_H3_hZ_bbZ_ATLAS13"] = boost::factory<Hobs_gg_H3_hZ_bbZ_ATLAS13*>();
    obsThFactory["Hobs_bb_H3_hZ_bbZ_ATLAS13"] = boost::factory<Hobs_bb_H3_hZ_bbZ_ATLAS13*>();
    obsThFactory["Hobs_gg_H3_hZ_bbZ_1_CMS13"] = boost::factory<Hobs_gg_H3_hZ_bbZ_1_CMS13*>();
    obsThFactory["Hobs_bb_H3_hZ_bbZ_1_CMS13"] = boost::factory<Hobs_bb_H3_hZ_bbZ_1_CMS13*>();
    obsThFactory["Hobs_gg_H3_hZ_bbZ_2_CMS13"] = boost::factory<Hobs_gg_H3_hZ_bbZ_2_CMS13*>();
    obsThFactory["Hobs_bb_H3_hZ_bbZ_2_CMS13"] = boost::factory<Hobs_bb_H3_hZ_bbZ_2_CMS13*>();
    obsThFactory["Hobs_pp_H3_H1Z_bbll_CMS8"] = boost::factory<Hobs_pp_H3_H1Z_bbll_CMS8*>();
    obsThFactory["Hobs_pp_H1_H3Z_bbll_CMS8"] = boost::factory<Hobs_pp_H1_H3Z_bbll_CMS8*>();
    obsThFactory["Hobs_pp_H5_H3Z_bbll_CMS8"] = boost::factory<Hobs_pp_H5_H3Z_bbll_CMS8*>();
    obsThFactory["Hobs_gg_H3_H1Z_bbll_ATLAS13"] = boost::factory<Hobs_gg_H3_H1Z_bbll_ATLAS13*>();
    obsThFactory["Hobs_bb_H3_H1Z_bbll_ATLAS13"] = boost::factory<Hobs_bb_H3_H1Z_bbll_ATLAS13*>();
    obsThFactory["Hobs_pp_H3pm_taunu_ATLAS8"] = boost::factory<Hobs_pp_H3pm_taunu_ATLAS8*>();
    obsThFactory["Hobs_pp_H3p_taunu_CMS8"] = boost::factory<Hobs_pp_H3p_taunu_CMS8*>();
    obsThFactory["Hobs_pp_H3pm_taunu_ATLAS13"] = boost::factory<Hobs_pp_H3pm_taunu_ATLAS13*>();
    obsThFactory["Hobs_pp_H3pm_taunu_CMS13"] = boost::factory<Hobs_pp_H3pm_taunu_CMS13*>();
    obsThFactory["Hobs_pp_H3pm_tb_ATLAS8"] = boost::factory<Hobs_pp_H3pm_tb_ATLAS8*>();
    obsThFactory["Hobs_pp_H3p_tb_CMS8"] = boost::factory<Hobs_pp_H3p_tb_CMS8*>();
    obsThFactory["Hobs_pp_H3pm_tb_ATLAS13"] = boost::factory<Hobs_pp_H3pm_tb_ATLAS13*>();
    obsThFactory["Hobs_WZ_H5pm_WZ_qqll_ATLAS8"] = boost::factory<Hobs_WZ_H5pm_WZ_qqll_ATLAS8*>();
    obsThFactory["Hobs_WZ_H5pm_WZ_lnull_ATLAS13"] = boost::factory<Hobs_WZ_H5pm_WZ_lnull_ATLAS13*>();
    obsThFactory["Robs_WZ_H5pm_WZ_lnull_ATLAS13"] = boost::factory<Robs_WZ_H5pm_WZ_lnull_ATLAS13*>();
    obsThFactory["Hobs_WZ_H5pm_WZ_lnull_1_CMS13"] = boost::factory<Hobs_WZ_H5pm_WZ_lnull_1_CMS13*>();
    obsThFactory["Hobs_WZ_H5pm_WZ_lnull_2_CMS13"] = boost::factory<Hobs_WZ_H5pm_WZ_lnull_2_CMS13*>();
    obsThFactory["Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8"] = boost::factory<Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8*>();
    obsThFactory["Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8"] = boost::factory<Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8*>();
    obsThFactory["Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8"] = boost::factory<Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8*>();
    obsThFactory["Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13"] = boost::factory<Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13*>();
    obsThFactory["Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13"] = boost::factory<Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13*>();
    obsThFactory["Hobs_VV_H5ppmm_WW_jjll_CMS8"] = boost::factory<Hobs_VV_H5ppmm_WW_jjll_CMS8*>();
    obsThFactory["Hobs_VV_H5ppmm_WW_jjll_CMS13"] = boost::factory<Hobs_VV_H5ppmm_WW_jjll_CMS13*>();
    obsThFactory["log10_tt_H1_tt_TH13"] = boost::factory<log10_tt_H1_tt_TH13*>();
    obsThFactory["log10_bb_H1_tt_TH13"] = boost::factory<log10_bb_H1_tt_TH13*>();
    obsThFactory["log10_tt_H3_tt_TH13"] = boost::factory<log10_tt_H3_tt_TH13*>();
    obsThFactory["log10_bb_H3_tt_TH13"] = boost::factory<log10_bb_H3_tt_TH13*>();
    obsThFactory["log10_bb_H1_bb_TH8"] = boost::factory<log10_bb_H1_bb_TH8*>();
    obsThFactory["log10_gg_H1_bb_TH8"] = boost::factory<log10_gg_H1_bb_TH8*>();
    obsThFactory["log10_pp_H1_bb_TH13"] = boost::factory<log10_pp_H1_bb_TH13*>();
    obsThFactory["log10_bb_H1_bb_TH13"] = boost::factory<log10_bb_H1_bb_TH13*>();
    obsThFactory["log10_bb_H3_bb_TH8"] = boost::factory<log10_bb_H3_bb_TH8*>();
    obsThFactory["log10_gg_H3_bb_TH8"] = boost::factory<log10_gg_H3_bb_TH8*>();
    obsThFactory["log10_pp_H3_bb_TH13"] = boost::factory<log10_pp_H3_bb_TH13*>();
    obsThFactory["log10_bb_H3_bb_TH13"] = boost::factory<log10_bb_H3_bb_TH13*>();
    obsThFactory["log10_gg_H1_tautau_TH8"] = boost::factory<log10_gg_H1_tautau_TH8*>();
    obsThFactory["log10_bb_H1_tautau_TH8"] = boost::factory<log10_bb_H1_tautau_TH8*>();
    obsThFactory["log10_gg_H1_tautau_TH13"] = boost::factory<log10_gg_H1_tautau_TH13*>();
    obsThFactory["log10_bb_H1_tautau_TH13"] = boost::factory<log10_bb_H1_tautau_TH13*>();
    obsThFactory["log10_gg_H3_tautau_TH8"] = boost::factory<log10_gg_H3_tautau_TH8*>();
    obsThFactory["log10_bb_H3_tautau_TH8"] = boost::factory<log10_bb_H3_tautau_TH8*>();
    obsThFactory["log10_gg_H3_tautau_TH13"] = boost::factory<log10_gg_H3_tautau_TH13*>();
    obsThFactory["log10_bb_H3_tautau_TH13"] = boost::factory<log10_bb_H3_tautau_TH13*>();
    obsThFactory["log10_gg_H1_gaga_TH8"] = boost::factory<log10_gg_H1_gaga_TH8*>();
    obsThFactory["log10_pp_H1_gaga_TH13"] = boost::factory<log10_pp_H1_gaga_TH13*>();
    obsThFactory["log10_gg_H1_gaga_TH13"] = boost::factory<log10_gg_H1_gaga_TH13*>();
    obsThFactory["log10_gg_H3_gaga_TH8"] = boost::factory<log10_gg_H3_gaga_TH8*>();
    obsThFactory["log10_pp_H3_gaga_TH13"] = boost::factory<log10_pp_H3_gaga_TH13*>();
    obsThFactory["log10_gg_H3_gaga_TH13"] = boost::factory<log10_gg_H3_gaga_TH13*>();
    obsThFactory["log10_pp_H5_gaga_TH13"] = boost::factory<log10_pp_H5_gaga_TH13*>();
    obsThFactory["log10_pp_H1_Zga_llga_TH8"] = boost::factory<log10_pp_H1_Zga_llga_TH8*>();
    obsThFactory["log10_gg_H1_Zga_TH13"] = boost::factory<log10_gg_H1_Zga_TH13*>();
    obsThFactory["log10_pp_H3_Zga_llga_TH8"] = boost::factory<log10_pp_H3_Zga_llga_TH8*>();
    obsThFactory["log10_gg_H3_Zga_TH13"] = boost::factory<log10_gg_H3_Zga_TH13*>();
    obsThFactory["log10_pp_H5_Zga_llga_TH8"] = boost::factory<log10_pp_H5_Zga_llga_TH8*>();
    obsThFactory["log10_gg_H1_ZZ_TH8"] = boost::factory<log10_gg_H1_ZZ_TH8*>();
    obsThFactory["log10_VV_H1_ZZ_TH8"] = boost::factory<log10_VV_H1_ZZ_TH8*>();
    obsThFactory["log10_gg_H1_ZZ_TH13"] = boost::factory<log10_gg_H1_ZZ_TH13*>();
    obsThFactory["log10_VV_H1_ZZ_TH13"] = boost::factory<log10_VV_H1_ZZ_TH13*>();
    obsThFactory["log10_pp_H1_ZZ_TH13"] = boost::factory<log10_pp_H1_ZZ_TH13*>();
    obsThFactory["log10_VV_H5_ZZ_TH8"] = boost::factory<log10_VV_H5_ZZ_TH8*>();
    obsThFactory["log10_VV_H5_ZZ_TH13"] = boost::factory<log10_VV_H5_ZZ_TH13*>();
    obsThFactory["log10_pp_H5_ZZ_TH13"] = boost::factory<log10_pp_H5_ZZ_TH13*>();
    obsThFactory["log10_gg_H1_WW_TH8"] = boost::factory<log10_gg_H1_WW_TH8*>();
    obsThFactory["log10_VV_H1_WW_TH8"] = boost::factory<log10_VV_H1_WW_TH8*>();
    obsThFactory["log10_gg_H1_WW_TH13"] = boost::factory<log10_gg_H1_WW_TH13*>();
    obsThFactory["log10_VV_H1_WW_TH13"] = boost::factory<log10_VV_H1_WW_TH13*>();
    obsThFactory["log10_ggVV_H1_WW_lnulnu_TH13"] = boost::factory<log10_ggVV_H1_WW_lnulnu_TH13*>();
    obsThFactory["log10_pp_H1_WW_TH13"] = boost::factory<log10_pp_H1_WW_TH13*>();
    obsThFactory["log10_VV_H5_WW_TH8"] = boost::factory<log10_VV_H5_WW_TH8*>();
    obsThFactory["log10_VV_H5_WW_TH13"] = boost::factory<log10_VV_H5_WW_TH13*>();
    obsThFactory["log10_ggVV_H5_WW_lnulnu_TH13"] = boost::factory<log10_ggVV_H5_WW_lnulnu_TH13*>();
    obsThFactory["log10_pp_H5_WW_TH13"] = boost::factory<log10_pp_H5_WW_TH13*>();
    obsThFactory["log10_pp_H1_VV_TH8"] = boost::factory<log10_pp_H1_VV_TH8*>();
    obsThFactory["log10_mu_pp_H1_VV_TH8"] = boost::factory<log10_mu_pp_H1_VV_TH8*>();
    obsThFactory["log10_pp_H1_VV_TH13"] = boost::factory<log10_pp_H1_VV_TH13*>();
    obsThFactory["log10_pp_H5_VV_TH8"] = boost::factory<log10_pp_H5_VV_TH8*>();
    obsThFactory["log10_mu_pp_H5_VV_TH8"] = boost::factory<log10_mu_pp_H5_VV_TH8*>();
    obsThFactory["log10_pp_H5_VV_TH13"] = boost::factory<log10_pp_H5_VV_TH13*>();
    obsThFactory["log10_gg_H1_hh_TH8"] = boost::factory<log10_gg_H1_hh_TH8*>();
    obsThFactory["log10_pp_H1_hh_TH8"] = boost::factory<log10_pp_H1_hh_TH8*>();
    obsThFactory["log10_pp_H1_hh_bbbb_TH8"] = boost::factory<log10_pp_H1_hh_bbbb_TH8*>();
    obsThFactory["log10_pp_H1_hh_gagabb_TH8"] = boost::factory<log10_pp_H1_hh_gagabb_TH8*>();
    obsThFactory["log10_gg_H1_hh_bbtautau_TH8"] = boost::factory<log10_gg_H1_hh_bbtautau_TH8*>();
    obsThFactory["log10_pp_H1_hh_TH13"] = boost::factory<log10_pp_H1_hh_TH13*>();
    obsThFactory["log10_gg_H1_hh_TH13"] = boost::factory<log10_gg_H1_hh_TH13*>();
    obsThFactory["log10_pp_H1_hh_bbbb_TH13"] = boost::factory<log10_pp_H1_hh_bbbb_TH13*>();
    obsThFactory["log10_gg_H1_hh_bbbb_TH13"] = boost::factory<log10_gg_H1_hh_bbbb_TH13*>();
    obsThFactory["log10_pp_H1_hh_gagabb_TH13"] = boost::factory<log10_pp_H1_hh_gagabb_TH13*>();
    obsThFactory["log10_pp_H1_hh_bbtautau_TH13"] = boost::factory<log10_pp_H1_hh_bbtautau_TH13*>();
    obsThFactory["log10_pp_H1_hh_bblnulnu_TH13"] = boost::factory<log10_pp_H1_hh_bblnulnu_TH13*>();
    obsThFactory["log10_gg_H1_hh_gagaWW_TH13"] = boost::factory<log10_gg_H1_hh_gagaWW_TH13*>();
    obsThFactory["log10_gg_H3_hZ_bbZ_TH8"] = boost::factory<log10_gg_H3_hZ_bbZ_TH8*>();
    obsThFactory["log10_gg_H3_hZ_bbll_TH8"] = boost::factory<log10_gg_H3_hZ_bbll_TH8*>();
    obsThFactory["log10_gg_H3_hZ_tautauZ_TH8"] = boost::factory<log10_gg_H3_hZ_tautauZ_TH8*>();
    obsThFactory["log10_gg_H3_hZ_tautaull_TH8"] = boost::factory<log10_gg_H3_hZ_tautaull_TH8*>();
    obsThFactory["log10_gg_H3_hZ_bbZ_TH13"] = boost::factory<log10_gg_H3_hZ_bbZ_TH13*>();
    obsThFactory["log10_bb_H3_hZ_bbZ_TH13"] = boost::factory<log10_bb_H3_hZ_bbZ_TH13*>();
    obsThFactory["log10_pp_H3_H1Z_bbll_TH8"] = boost::factory<log10_pp_H3_H1Z_bbll_TH8*>();
    obsThFactory["log10_pp_H3_H5Z_bbll_TH8"] = boost::factory<log10_pp_H3_H5Z_bbll_TH8*>();
    obsThFactory["log10_pp_H1_H3Z_bbll_TH8"] = boost::factory<log10_pp_H1_H3Z_bbll_TH8*>();
    obsThFactory["log10_pp_H5_H3Z_bbll_TH8"] = boost::factory<log10_pp_H5_H3Z_bbll_TH8*>();
    obsThFactory["log10_pp_H3pm_taunu_TH8"] = boost::factory<log10_pp_H3pm_taunu_TH8*>();
    obsThFactory["log10_pp_H3p_taunu_TH8"] = boost::factory<log10_pp_H3p_taunu_TH8*>();
    obsThFactory["log10_pp_H3pm_taunu_TH13"] = boost::factory<log10_pp_H3pm_taunu_TH13*>();
    obsThFactory["log10_pp_H3pm_tb_TH8"] = boost::factory<log10_pp_H3pm_tb_TH8*>();
    obsThFactory["log10_pp_H3p_tb_TH8"] = boost::factory<log10_pp_H3p_tb_TH8*>();
    obsThFactory["log10_pp_H3pm_tb_TH13"] = boost::factory<log10_pp_H3pm_tb_TH13*>();
    obsThFactory["log10_WZ_H5pm_WZ_TH8"] = boost::factory<log10_WZ_H5pm_WZ_TH8*>();
    obsThFactory["log10_WZ_H5pm_WZ_TH13"] = boost::factory<log10_WZ_H5pm_WZ_TH13*>();
    obsThFactory["log10_pp_H5ppmmH5mmpp_TH8"] = boost::factory<log10_pp_H5ppmmH5mmpp_TH8*>();
    obsThFactory["log10_pp_H5ppmmH5mmpp_TH13"] = boost::factory<log10_pp_H5ppmmH5mmpp_TH13*>();
    obsThFactory["log10_pp_H5ppmmH5mmpp_WWWW_TH13"] = boost::factory<log10_pp_H5ppmmH5mmpp_WWWW_TH13*>();
    obsThFactory["log10_VV_H5ppmm_WW_TH8"] = boost::factory<log10_VV_H5ppmm_WW_TH8*>();
    obsThFactory["log10_VV_H5ppmm_WW_TH13"] = boost::factory<log10_VV_H5ppmm_WW_TH13*>();

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
