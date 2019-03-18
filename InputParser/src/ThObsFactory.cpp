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
#include "THDMObservables.h"
#include "LRSMObservables.h"
/* BEGIN: REMOVE FROM THE PACKAGE */
#include "GeneralTHDMObservables.h"
#include "THDMWObservables.h"
/* END: REMOVE FROM THE PACKAGE */
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

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
    const double sqrt_s_LEP2_207 = 0.2066; ///< the center-of-mass energy in TeV
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
    //
    const double sqrt_s_LHeC_1_3 = 1.3; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHeC_1_8 = 1.8; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_FCCep_3_5 = 3.5; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCCep_5 = 5.0; ///< the center-of-mass energy in TeV
    // Polarizations at lepton colliders
    const double pol_0 = 0.0;
    const double pol_30 = 30.0;
    const double pol_80 = 80.0;
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
    //-----  Triple gauge coupling observables  -----
    obsThFactory["deltag1Z"] = boost::factory<deltag1Z*>();
    obsThFactory["deltaKgamma"] = boost::factory<deltaKgamma*>();
    obsThFactory["lambdaZ"] = boost::factory<lambdaZ*>();
    //-----  ee -> WW observables: LEP2 total cross section  -----
    obsThFactory["eeWW_LEP2_161"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_161);
    obsThFactory["eeWW_LEP2_172"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_172);
    obsThFactory["eeWW_LEP2_183"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_183);
    obsThFactory["eeWW_LEP2_189"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_189);
    obsThFactory["eeWW_LEP2_192"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_192);
    obsThFactory["eeWW_LEP2_196"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_196);
    obsThFactory["eeWW_LEP2_200"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_200);
    obsThFactory["eeWW_LEP2_202"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_202);
    obsThFactory["eeWW_LEP2_205"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_205);
    obsThFactory["eeWW_LEP2_207"] = boost::bind(boost::factory<xseeWW*>(), _1, sqrt_s_LEP2_207);
    //-----  ee -> WW observables: LEP2 differential cross section  -----
    obsThFactory["deeWWdcos_LEP2_183_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_183_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos1, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_189_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_189_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos2, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_198_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_198_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos3, cos10_LEP2_WW, cos11_LEP2_WW);
    //
    obsThFactory["deeWWdcos_LEP2_206_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos1_LEP2_WW, cos2_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos2_LEP2_WW, cos3_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos3_LEP2_WW, cos4_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos4_LEP2_WW, cos5_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos5_LEP2_WW, cos6_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos6_LEP2_WW, cos7_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos7_LEP2_WW, cos8_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos8_LEP2_WW, cos9_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos9_LEP2_WW, cos10_LEP2_WW);
    obsThFactory["deeWWdcos_LEP2_206_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_LEP2_WWcos4, cos10_LEP2_WW, cos11_LEP2_WW);
    //-----  ee -> WW observables: Future colliders differential cross section  -----
    obsThFactory["deeWWdcos_161_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_161_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_161_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_161_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_161_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_161_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_161_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_161_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_161_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_161_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_161, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_240_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_240_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_240_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_240_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_240_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_240_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_240_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_240_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_240_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_240_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_240, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_250_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_250_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_250_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_250_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_250_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_250_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_250_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_250_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_250_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_250_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_250, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_350_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_350_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_350_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_350_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_350_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_350_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_350_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_350_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_350_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_350_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_350, cos10_ee_WW, cos11_ee_WW);
    //
    obsThFactory["deeWWdcos_365_Bin1"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos1_ee_WW, cos2_ee_WW);
    obsThFactory["deeWWdcos_365_Bin2"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos2_ee_WW, cos3_ee_WW);
    obsThFactory["deeWWdcos_365_Bin3"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos3_ee_WW, cos4_ee_WW);
    obsThFactory["deeWWdcos_365_Bin4"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos4_ee_WW, cos5_ee_WW);
    obsThFactory["deeWWdcos_365_Bin5"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos5_ee_WW, cos6_ee_WW);
    obsThFactory["deeWWdcos_365_Bin6"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos6_ee_WW, cos7_ee_WW);
    obsThFactory["deeWWdcos_365_Bin7"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos7_ee_WW, cos8_ee_WW);
    obsThFactory["deeWWdcos_365_Bin8"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos8_ee_WW, cos9_ee_WW);
    obsThFactory["deeWWdcos_365_Bin9"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos9_ee_WW, cos10_ee_WW);
    obsThFactory["deeWWdcos_365_Bin10"] = boost::bind(boost::factory<dxseeWWdcosBin*>(), _1, sqrt_s_leptcoll_365, cos10_ee_WW, cos11_ee_WW);
    //-----  ee -> WW observables: total rates (ratio with the SM)  -----
    obsThFactory["eeWW161"] = boost::bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_161);
    //
    obsThFactory["eeWW240"] = boost::bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeWW240_p80_m30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_240, pol_80, -pol_30);
    obsThFactory["eeWW240_m80_p30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_240, -pol_80, pol_30);
    //
    obsThFactory["eeWW250"] = boost::bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeWW250_p80_m30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeWW250_m80_p30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    //
    obsThFactory["eeWW350"] = boost::bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeWW350_p80_m30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeWW350_m80_p30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    //
    obsThFactory["eeWW380_p80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeWW380_m80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeWW365"] = boost::bind(boost::factory<mueeWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeWW365_p80_m30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeWW365_m80_p30"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    //
    obsThFactory["eeWW1500_p80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeWW1500_m80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeWW3000_p80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeWW3000_m80_0"] = boost::bind(boost::factory<mueeWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //----- High Energy diboson observables at hadron colliders
    obsThFactory["ppZHprobe14"] = boost::bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_LHC14);
    obsThFactory["ppZHprobe27"] = boost::bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_LHC27);
    obsThFactory["ppZHprobe100"] = boost::bind(boost::factory<ppZHprobe*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["mupTVppWZ_14_Bin1"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 100., 150.);
    obsThFactory["mupTVppWZ_14_Bin2"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 150., 220.);
    obsThFactory["mupTVppWZ_14_Bin3"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 220., 300.);
    obsThFactory["mupTVppWZ_14_Bin4"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 300., 500.);
    obsThFactory["mupTVppWZ_14_Bin5"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 500., 750.);
    obsThFactory["mupTVppWZ_14_Bin6"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC14, 750., 1200.);
    //
    obsThFactory["mupTVppWZ_27_Bin1"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 150., 220.);
    obsThFactory["mupTVppWZ_27_Bin2"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 220., 300.);
    obsThFactory["mupTVppWZ_27_Bin3"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 300., 500.);
    obsThFactory["mupTVppWZ_27_Bin4"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 500., 750.);
    obsThFactory["mupTVppWZ_27_Bin5"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 750., 1200.);
    obsThFactory["mupTVppWZ_27_Bin6"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_LHC27, 1200., 1800.);
    //
    obsThFactory["mupTVppWZ_100_Bin1"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 220., 300.);
    obsThFactory["mupTVppWZ_100_Bin2"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 300., 500.);
    obsThFactory["mupTVppWZ_100_Bin3"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 500., 750.);
    obsThFactory["mupTVppWZ_100_Bin4"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 750., 1200.);
    obsThFactory["mupTVppWZ_100_Bin5"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 1200., 1800.);
    obsThFactory["mupTVppWZ_100_Bin6"] = boost::bind(boost::factory<mupTVppWZ*>(), _1, sqrt_s_FCC100, 1800., 2400.);
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
    obsThFactory["lambz_HB"] = boost::factory<lambzHB*>();
    //-----  Other useful observables to work with new physics  ----------
    //-----  Z couplings with leptons ---------
    obsThFactory["delgZeL"] = boost::bind(boost::factory<delgZlL*>(), _1, StandardModel::ELECTRON);
    obsThFactory["delgZeR"] = boost::bind(boost::factory<delgZlR*>(), _1, StandardModel::ELECTRON);
    obsThFactory["delgZmuL"] = boost::bind(boost::factory<delgZlL*>(), _1, StandardModel::MU);
    obsThFactory["delgZmuR"] = boost::bind(boost::factory<delgZlR*>(), _1, StandardModel::MU);
    obsThFactory["delgZtaL"] = boost::bind(boost::factory<delgZlL*>(), _1, StandardModel::TAU);
    obsThFactory["delgZtaR"] = boost::bind(boost::factory<delgZlR*>(), _1, StandardModel::TAU);
    //-----  Z couplings with up sector quarks ---------
    obsThFactory["delgZuL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::UP);
    obsThFactory["delgZuR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::UP);
    obsThFactory["delgZcL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::CHARM);
    obsThFactory["delgZcR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::CHARM);
    obsThFactory["delgZtL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::TOP);
    obsThFactory["delgZtR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::TOP);
    //-----  Z couplings with down sector quarks ---------
    obsThFactory["delgZdL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::DOWN);
    obsThFactory["delgZdR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::DOWN);
    obsThFactory["delgZsL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::STRANGE);
    obsThFactory["delgZsR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::STRANGE);
    obsThFactory["delgZbL"] = boost::bind(boost::factory<delgZqL*>(), _1, StandardModel::BOTTOM);
    obsThFactory["delgZbR"] = boost::bind(boost::factory<delgZqR*>(), _1, StandardModel::BOTTOM);
    obsThFactory["deltaMW"] = boost::factory<deltaMW*>();
    obsThFactory["oblWpar"] = boost::factory<oblW*>();
    obsThFactory["oblYpar"] = boost::factory<oblY*>();
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

    //-----  Higgs Extension observables  ----------

    //-----  Production cross sections (ratios with SM)  ----------
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
    obsThFactory["VBFgamma8"] = boost::bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC8);
    obsThFactory["VH8"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC8);
    obsThFactory["WH8"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ZH8"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ttH8"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC8);
    obsThFactory["ggH13"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH+ttH13"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF13"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBF+VH13"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["VBFgamma13"] = boost::bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC13);
    obsThFactory["VH13"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC13);
    obsThFactory["WH13"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ZH13"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ttH13"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC13);
    obsThFactory["ggH14"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH+ttH14"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF14"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBF+VH14"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["VBFgamma14"] = boost::bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC14);
    obsThFactory["VH14"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC14);
    obsThFactory["WH14"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ZH14"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ttH14"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC14);
    obsThFactory["ggH27"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ggH+ttH27"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBF27"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBF+VH27"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_LHC27);
    obsThFactory["VBFgamma27"] = boost::bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_LHC27);
    obsThFactory["VH27"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_LHC27);
    obsThFactory["WH27"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ZH27"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ttH27"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_LHC27);
    obsThFactory["ggH100"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH+ttH100"] = boost::bind(boost::factory<muggHpttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF100"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBF+VH100"] = boost::bind(boost::factory<muVBFpVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["VBFgamma100"] = boost::bind(boost::factory<muVBFgamma*>(), _1, sqrt_s_FCC100);
    obsThFactory["VH100"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_FCC100);
    obsThFactory["WH100"] = boost::bind(boost::factory<muWH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ZH100"] = boost::bind(boost::factory<muZH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ttH100"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggH196"] = boost::bind(boost::factory<muggH*>(), _1, sqrt_s_TeV);
    obsThFactory["VBF196"] = boost::bind(boost::factory<muVBF*>(), _1, sqrt_s_TeV);
    obsThFactory["VH196"] = boost::bind(boost::factory<muVH*>(), _1, sqrt_s_TeV);
    obsThFactory["ttH196"] = boost::bind(boost::factory<muttH*>(), _1, sqrt_s_TeV);
    //
    obsThFactory["eeZH240"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZH250"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZH350"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZH365"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZH380"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZH500"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZH1000"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZH1400"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZH1500"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZH3000"] = boost::bind(boost::factory<mueeZH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZllH240"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZllH250"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZllH350"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZllH365"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZllH380"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZllH500"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZllH1000"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZllH1400"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZllH1500"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZllH3000"] = boost::bind(boost::factory<mueeZllH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZqqH240"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZqqH250"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZqqH350"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZqqH365"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZqqH380"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZqqH500"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZqqH1000"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZqqH1400"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZqqH1500"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZqqH3000"] = boost::bind(boost::factory<mueeZqqH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZH250_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZH250_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZH250_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZH250_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZH350_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZH350_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZH350_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZH350_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZH365_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZH365_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZH365_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZH365_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZH380_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZH380_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZH380_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZH380_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZH500_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZH500_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZH500_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZH500_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZH1000_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZH1000_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZH1000_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZH1000_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZH1400_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZH1400_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZH1400_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZH1400_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZH1500_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZH1500_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZH1500_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZH1500_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZH3000_p80_m30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZH3000_m80_p30"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZH3000_p80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZH3000_m80_0"] = boost::bind(boost::factory<mueeZHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZllH250_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZllH250_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZllH250_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZllH250_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZllH350_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZllH350_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZllH350_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZllH350_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZllH365_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZllH365_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZllH365_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZllH365_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZllH380_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZllH380_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZllH380_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZllH380_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZllH500_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZllH500_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZllH500_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZllH500_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1000_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZllH1000_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZllH1000_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZllH1000_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1400_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZllH1400_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZllH1400_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZllH1400_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZllH1500_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZllH1500_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZllH1500_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZllH1500_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZllH3000_p80_m30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZllH3000_m80_p30"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZllH3000_p80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZllH3000_m80_0"] = boost::bind(boost::factory<mueeZllHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH250_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZqqH250_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZqqH250_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZqqH250_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH350_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZqqH350_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZqqH350_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZqqH350_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH365_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZqqH365_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZqqH365_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZqqH365_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH380_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZqqH380_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZqqH380_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZqqH380_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH500_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZqqH500_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZqqH500_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZqqH500_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1000_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZqqH1000_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZqqH1000_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZqqH1000_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1400_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZqqH1400_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZqqH1400_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZqqH1400_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH1500_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZqqH1500_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZqqH1500_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZqqH1500_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZqqH3000_p80_m30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZqqH3000_m80_p30"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZqqH3000_p80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZqqH3000_m80_0"] = boost::bind(boost::factory<mueeZqqHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeWBF240"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeWBF250"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeWBF350"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeWBF365"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeWBF380"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeWBF500"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeWBF1000"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeWBF1400"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeWBF1500"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeWBF3000"] = boost::bind(boost::factory<mueeWBF*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeWBF250_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeWBF250_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeWBF250_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeWBF250_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeWBF350_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeWBF350_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeWBF350_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeWBF350_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeWBF365_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeWBF365_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeWBF365_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeWBF365_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeWBF380_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeWBF380_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeWBF380_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeWBF380_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeWBF500_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeWBF500_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeWBF500_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeWBF500_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1000_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeWBF1000_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeWBF1000_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeWBF1000_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1400_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeWBF1400_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeWBF1400_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeWBF1400_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeWBF1500_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeWBF1500_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeWBF1500_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeWBF1500_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeWBF3000_p80_m30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeWBF3000_m80_p30"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeWBF3000_p80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeWBF3000_m80_0"] = boost::bind(boost::factory<mueeWBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeHvv240"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeHvv250"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeHvv350"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeHvv365"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeHvv380"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeHvv500"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeHvv1000"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeHvv1400"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeHvv1500"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeHvv3000"] = boost::bind(boost::factory<mueeHvv*>(), _1, sqrt_s_leptcoll_3000);
     //
    obsThFactory["eeHvv250_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeHvv250_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeHvv250_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeHvv250_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeHvv350_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeHvv350_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeHvv350_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeHvv350_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeHvv365_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeHvv365_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeHvv365_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeHvv365_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeHvv380_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeHvv380_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeHvv380_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeHvv380_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeHvv500_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeHvv500_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeHvv500_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeHvv500_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1000_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeHvv1000_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeHvv1000_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeHvv1000_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1400_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeHvv1400_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeHvv1400_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeHvv1400_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeHvv1500_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeHvv1500_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeHvv1500_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeHvv1500_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeHvv3000_p80_m30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeHvv3000_m80_p30"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeHvv3000_p80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeHvv3000_m80_0"] = boost::bind(boost::factory<mueeHvvPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eeZBF240"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["eeZBF250"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["eeZBF350"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["eeZBF365"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["eeZBF380"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["eeZBF500"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eeZBF1000"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eeZBF1400"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eeZBF1500"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eeZBF3000"] = boost::bind(boost::factory<mueeZBF*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eeZBF250_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["eeZBF250_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["eeZBF250_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["eeZBF250_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["eeZBF350_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["eeZBF350_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["eeZBF350_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["eeZBF350_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["eeZBF365_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, -pol_30);
    obsThFactory["eeZBF365_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_30);
    obsThFactory["eeZBF365_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, pol_80, pol_0);
    obsThFactory["eeZBF365_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_365, -pol_80, pol_0);
    //
    obsThFactory["eeZBF380_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["eeZBF380_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["eeZBF380_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["eeZBF380_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["eeZBF500_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eeZBF500_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eeZBF500_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eeZBF500_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1000_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eeZBF1000_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eeZBF1000_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eeZBF1000_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1400_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eeZBF1400_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eeZBF1400_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eeZBF1400_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eeZBF1500_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eeZBF1500_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eeZBF1500_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eeZBF1500_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eeZBF3000_p80_m30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eeZBF3000_m80_p30"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eeZBF3000_p80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eeZBF3000_m80_0"] = boost::bind(boost::factory<mueeZBFPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["eettH500"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["eettH1000"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["eettH1400"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["eettH1500"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["eettH3000"] = boost::bind(boost::factory<mueettH*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["eettH500_p80_m30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["eettH500_m80_p30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["eettH500_p80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["eettH500_m80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["eettH1000_p80_m30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["eettH1000_m80_p30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["eettH1000_p80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["eettH1000_m80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["eettH1400_p80_m30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["eettH1400_m80_p30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["eettH1400_p80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["eettH1400_m80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["eettH1500_p80_m30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["eettH1500_m80_p30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["eettH1500_p80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["eettH1500_m80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["eettH3000_p80_m30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["eettH3000_m80_p30"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["eettH3000_p80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["eettH3000_m80_0"] = boost::bind(boost::factory<mueettHPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mumuH125"] = boost::bind(boost::factory<mummH*>(), _1, sqrt_s_leptcoll_125);
    //
    obsThFactory["epWBF1300"] = boost::bind(boost::factory<muepWBF*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBF1800"] = boost::bind(boost::factory<muepWBF*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["epWBF3500"] = boost::bind(boost::factory<muepWBF*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBF5000"] = boost::bind(boost::factory<muepWBF*>(), _1, sqrt_s_FCCep_5);
    //
    obsThFactory["epZBF1300"] = boost::bind(boost::factory<muepZBF*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBF1800"] = boost::bind(boost::factory<muepZBF*>(), _1, sqrt_s_LHeC_1_8);
    obsThFactory["epZBF3500"] = boost::bind(boost::factory<muepZBF*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBF5000"] = boost::bind(boost::factory<muepZBF*>(), _1, sqrt_s_FCCep_5);
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
    //-----  Ratios of BR (ratios with SM)  ----------
    obsThFactory["BrHtogaga_over_mumu_Ratio"] = boost::factory<BrHtogaga_over_mumu_Ratio*>();
    obsThFactory["BrHtoZga_over_mumu_Ratio"] = boost::factory<BrHtoZga_over_mumu_Ratio*>();
    obsThFactory["BrHtoZmumuga_over_mumu_Ratio"] = boost::factory<BrHtoZmumuga_over_mumu_Ratio*>();
    obsThFactory["BrHtogaga_over_4l_Ratio"] = boost::factory<BrHtogaga_over_4l_Ratio*>();
    obsThFactory["BrHtogaga_over_2e2mu_Ratio"] = boost::factory<BrHtogaga_over_2e2mu_Ratio*>();
    obsThFactory["BrHtoZga_over_4l_Ratio"] = boost::factory<BrHtoZga_over_4l_Ratio*>();
    obsThFactory["BrHtomumu_over_4l_Ratio"] = boost::factory<BrHtomumu_over_4l_Ratio*>();
    obsThFactory["BrHtomumu_over_4mu_Ratio"] = boost::factory<BrHtomumu_over_4mu_Ratio*>();
    obsThFactory["BrHto4l_over_gaga_Ratio"] = boost::factory<BrHto4l_over_gaga_Ratio*>();
    obsThFactory["BrHtoZga_over_gaga_Ratio"] = boost::factory<BrHtoZga_over_gaga_Ratio*>();
    obsThFactory["BrHtomumu_over_gaga_Ratio"] = boost::factory<BrHtomumu_over_gaga_Ratio*>();
    //-----  Special observables --------
    obsThFactory["muttHZbb_boost100"] = boost::bind(boost::factory<muttHZbbboost*>(), _1, sqrt_s_FCC100);
    obsThFactory["ggHgagaInt14"] = boost::bind(boost::factory<muggHgagaInt*>(), _1, sqrt_s_LHC14);
    //-----  Full Signal strengths per prod and decay: Hadron colliders  ----------
    obsThFactory["muggHgaga"] = boost::bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHgaga"] = boost::bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHgaga"] = boost::bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHgaga"] = boost::bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZZ"] = boost::bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ"] = boost::bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ"] = boost::bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ"] = boost::bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHWW"] = boost::bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW"] = boost::bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW"] = boost::bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW"] = boost::bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHtautau"] = boost::bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHtautau"] = boost::bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHtautau"] = boost::bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHtautau"] = boost::bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHbb"] = boost::bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHbb"] = boost::bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHbb"] = boost::bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHbb"] = boost::bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC13);
    // Indicating the energy explicit in the observable name
    obsThFactory["muggHgaga13"] = boost::bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHgaga13"] = boost::bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHgaga13"] = boost::bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHgaga13"] = boost::bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHgaga13"] = boost::bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHgaga13"] = boost::bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZga13"] = boost::bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZga13"] = boost::bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZga13"] = boost::bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZga13"] = boost::bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZga13"] = boost::bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZga13"] = boost::bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHZZ13"] = boost::bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ13"] = boost::bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZZ13"] = boost::bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZZ13"] = boost::bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ13"] = boost::bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ13"] = boost::bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHZZ4l13"] = boost::bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHZZ4l13"] = boost::bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHZZ4l13"] = boost::bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHZZ4l13"] = boost::bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHZZ4l13"] = boost::bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHZZ4l13"] = boost::bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHWW13"] = boost::bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW13"] = boost::bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHWW13"] = boost::bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHWW13"] = boost::bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW13"] = boost::bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW13"] = boost::bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHWW2l2v13"] = boost::bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHWW2l2v13"] = boost::bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHWW2l2v13"] = boost::bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHWW2l2v13"] = boost::bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHWW2l2v13"] = boost::bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHWW2l2v13"] = boost::bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHmumu13"] = boost::bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHmumu13"] = boost::bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHmumu13"] = boost::bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHmumu13"] = boost::bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHmumu13"] = boost::bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHmumu13"] = boost::bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHtautau13"] = boost::bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHtautau13"] = boost::bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHtautau13"] = boost::bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHtautau13"] = boost::bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHtautau13"] = boost::bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHtautau13"] = boost::bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHbb13"] = boost::bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVBFHbb13"] = boost::bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muZHbb13"] = boost::bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muWHbb13"] = boost::bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muVHbb13"] = boost::bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC13);
    obsThFactory["muttHbb13"] = boost::bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC13);
    //
    obsThFactory["muggHgaga14"] = boost::bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHgaga14"] = boost::bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHgaga14"] = boost::bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHgaga14"] = boost::bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHgaga14"] = boost::bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHgaga14"] = boost::bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHZga14"] = boost::bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZga14"] = boost::bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZga14"] = boost::bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZga14"] = boost::bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZga14"] = boost::bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZga14"] = boost::bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHZZ14"] = boost::bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZZ14"] = boost::bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZZ14"] = boost::bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZZ14"] = boost::bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZZ14"] = boost::bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZZ14"] = boost::bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHZZ4l14"] = boost::bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHZZ4l14"] = boost::bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHZZ4l14"] = boost::bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHZZ4l14"] = boost::bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHZZ4l14"] = boost::bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHZZ4l14"] = boost::bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHWW14"] = boost::bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHWW14"] = boost::bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHWW14"] = boost::bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHWW14"] = boost::bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHWW14"] = boost::bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHWW14"] = boost::bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHWW2l2v14"] = boost::bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHWW2l2v14"] = boost::bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHWW2l2v14"] = boost::bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHWW2l2v14"] = boost::bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHWW2l2v14"] = boost::bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHWW2l2v14"] = boost::bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHmumu14"] = boost::bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHmumu14"] = boost::bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHmumu14"] = boost::bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHmumu14"] = boost::bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHmumu14"] = boost::bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHmumu14"] = boost::bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHtautau14"] = boost::bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHtautau14"] = boost::bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHtautau14"] = boost::bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHtautau14"] = boost::bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHtautau14"] = boost::bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHtautau14"] = boost::bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHbb14"] = boost::bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVBFHbb14"] = boost::bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muZHbb14"] = boost::bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muWHbb14"] = boost::bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muVHbb14"] = boost::bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC14);
    obsThFactory["muttHbb14"] = boost::bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC14);
    //
    obsThFactory["muggHgaga27"] = boost::bind(boost::factory<muggHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHgaga27"] = boost::bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHgaga27"] = boost::bind(boost::factory<muZHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHgaga27"] = boost::bind(boost::factory<muWHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHgaga27"] = boost::bind(boost::factory<muVHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHgaga27"] = boost::bind(boost::factory<muttHgaga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHZga27"] = boost::bind(boost::factory<muggHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZga27"] = boost::bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZga27"] = boost::bind(boost::factory<muZHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZga27"] = boost::bind(boost::factory<muWHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZga27"] = boost::bind(boost::factory<muVHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZga27"] = boost::bind(boost::factory<muttHZga*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHZZ27"] = boost::bind(boost::factory<muggHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZZ27"] = boost::bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZZ27"] = boost::bind(boost::factory<muZHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZZ27"] = boost::bind(boost::factory<muWHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZZ27"] = boost::bind(boost::factory<muVHZZ*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZZ27"] = boost::bind(boost::factory<muttHZZ*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHZZ4l27"] = boost::bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHZZ4l27"] = boost::bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHZZ4l27"] = boost::bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHZZ4l27"] = boost::bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHZZ4l27"] = boost::bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHZZ4l27"] = boost::bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHWW27"] = boost::bind(boost::factory<muggHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHWW27"] = boost::bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHWW27"] = boost::bind(boost::factory<muZHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHWW27"] = boost::bind(boost::factory<muWHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHWW27"] = boost::bind(boost::factory<muVHWW*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHWW27"] = boost::bind(boost::factory<muttHWW*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHWW2l2v27"] = boost::bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHWW2l2v27"] = boost::bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHWW2l2v27"] = boost::bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHWW2l2v27"] = boost::bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHWW2l2v27"] = boost::bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHWW2l2v27"] = boost::bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHmumu27"] = boost::bind(boost::factory<muggHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHmumu27"] = boost::bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHmumu27"] = boost::bind(boost::factory<muZHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHmumu27"] = boost::bind(boost::factory<muWHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHmumu27"] = boost::bind(boost::factory<muVHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHmumu27"] = boost::bind(boost::factory<muttHmumu*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHtautau27"] = boost::bind(boost::factory<muggHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHtautau27"] = boost::bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHtautau27"] = boost::bind(boost::factory<muZHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHtautau27"] = boost::bind(boost::factory<muWHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHtautau27"] = boost::bind(boost::factory<muVHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHtautau27"] = boost::bind(boost::factory<muttHtautau*>(), _1, sqrt_s_LHC27);
    obsThFactory["muggHbb27"] = boost::bind(boost::factory<muggHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVBFHbb27"] = boost::bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muZHbb27"] = boost::bind(boost::factory<muZHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muWHbb27"] = boost::bind(boost::factory<muWHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muVHbb27"] = boost::bind(boost::factory<muVHbb*>(), _1, sqrt_s_LHC27);
    obsThFactory["muttHbb27"] = boost::bind(boost::factory<muttHbb*>(), _1, sqrt_s_LHC27);
    //
    obsThFactory["muggHgaga100"] = boost::bind(boost::factory<muggHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHgaga100"] = boost::bind(boost::factory<muVBFHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHgaga100"] = boost::bind(boost::factory<muZHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHgaga100"] = boost::bind(boost::factory<muWHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHgaga100"] = boost::bind(boost::factory<muVHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHgaga100"] = boost::bind(boost::factory<muttHgaga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZga100"] = boost::bind(boost::factory<muggHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZgamumu100"] = boost::bind(boost::factory<muggHZgamumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZga100"] = boost::bind(boost::factory<muVBFHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZga100"] = boost::bind(boost::factory<muZHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZga100"] = boost::bind(boost::factory<muWHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZga100"] = boost::bind(boost::factory<muVHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZga100"] = boost::bind(boost::factory<muttHZga*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZZ100"] = boost::bind(boost::factory<muggHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZZ100"] = boost::bind(boost::factory<muVBFHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZZ100"] = boost::bind(boost::factory<muZHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZZ100"] = boost::bind(boost::factory<muWHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZZ100"] = boost::bind(boost::factory<muVHZZ*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZZ100"] = boost::bind(boost::factory<muttHZZ*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHZZ4l100"] = boost::bind(boost::factory<muggHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHZZ4mu100"] = boost::bind(boost::factory<muggHZZ4mu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHZZ4l100"] = boost::bind(boost::factory<muVBFHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHZZ4l100"] = boost::bind(boost::factory<muZHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHZZ4l100"] = boost::bind(boost::factory<muWHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHZZ4l100"] = boost::bind(boost::factory<muVHZZ4l*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHZZ4l100"] = boost::bind(boost::factory<muttHZZ4l*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHWW100"] = boost::bind(boost::factory<muggHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHWW100"] = boost::bind(boost::factory<muVBFHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHWW100"] = boost::bind(boost::factory<muZHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHWW100"] = boost::bind(boost::factory<muWHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHWW100"] = boost::bind(boost::factory<muVHWW*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHWW100"] = boost::bind(boost::factory<muttHWW*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHWW2l2v100"] = boost::bind(boost::factory<muggHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHWW2l2v100"] = boost::bind(boost::factory<muVBFHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHWW2l2v100"] = boost::bind(boost::factory<muZHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHWW2l2v100"] = boost::bind(boost::factory<muWHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHWW2l2v100"] = boost::bind(boost::factory<muVHWW2l2v*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHWW2l2v100"] = boost::bind(boost::factory<muttHWW2l2v*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muggHmumu100"] = boost::bind(boost::factory<muggHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHmumu100"] = boost::bind(boost::factory<muVBFHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHmumu100"] = boost::bind(boost::factory<muZHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHmumu100"] = boost::bind(boost::factory<muWHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHmumu100"] = boost::bind(boost::factory<muVHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHmumu100"] = boost::bind(boost::factory<muttHmumu*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHtautau100"] = boost::bind(boost::factory<muggHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHtautau100"] = boost::bind(boost::factory<muVBFHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHtautau100"] = boost::bind(boost::factory<muZHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHtautau100"] = boost::bind(boost::factory<muWHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHtautau100"] = boost::bind(boost::factory<muVHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHtautau100"] = boost::bind(boost::factory<muttHtautau*>(), _1, sqrt_s_FCC100);
    obsThFactory["muggHbb100"] = boost::bind(boost::factory<muggHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVBFHbb100"] = boost::bind(boost::factory<muVBFHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muZHbb100"] = boost::bind(boost::factory<muZHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muWHbb100"] = boost::bind(boost::factory<muWHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muVHbb100"] = boost::bind(boost::factory<muVHbb*>(), _1, sqrt_s_FCC100);
    obsThFactory["muttHbb100"] = boost::bind(boost::factory<muttHbb*>(), _1, sqrt_s_FCC100);
    //
    obsThFactory["muppHmumu8"] = boost::bind(boost::factory<muppHmumu*>(), _1, sqrt_s_LHC8);
    obsThFactory["muppHmumu13"] = boost::bind(boost::factory<muppHmumu*>(), _1, sqrt_s_LHC13);
    obsThFactory["muppHZga8"] = boost::bind(boost::factory<muppHZga*>(), _1, sqrt_s_LHC8);
    obsThFactory["muppHZga13"] = boost::bind(boost::factory<muppHZga*>(), _1, sqrt_s_LHC13);
    obsThFactory["muggHH2ga2b14"] = boost::bind(boost::factory<muggHH2ga2b*>(), _1, sqrt_s_LHC14);
    obsThFactory["muggHH2ga2b100"] = boost::bind(boost::factory<muggHH2ga2b*>(), _1, sqrt_s_FCC100);
    //-----  Full Signal strengths per prod and decay: Lepton colliders  ----------
    //
    // Pure WBF
    obsThFactory["mueeWBFbb240"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFbb250"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFbb350"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFbb365"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFbb380"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFbb500"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFbb1000"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFbb1400"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFbb1500"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFbb3000"] = boost::bind(boost::factory<mueeWBFbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFbb250_p80_m30"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeWBFbb250_m80_p30"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeWBFbb250_p80_0"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeWBFbb250_m80_0"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFbb350_p80_m30"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeWBFbb350_m80_p30"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeWBFbb350_p80_0"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeWBFbb350_m80_0"] = boost::bind(boost::factory<mueeWBFbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeWBFcc240"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFcc250"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFcc350"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFcc365"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFcc380"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFcc500"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFcc1000"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFcc1400"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFcc1500"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFcc3000"] = boost::bind(boost::factory<mueeWBFcc*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFgg240"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFgg250"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFgg350"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFgg365"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFgg380"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFgg500"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFgg1000"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFgg1400"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFgg1500"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFgg3000"] = boost::bind(boost::factory<mueeWBFgg*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFWW240"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFWW250"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFWW350"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFWW365"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFWW380"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFWW500"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFWW1000"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFWW1400"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFWW1500"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFWW3000"] = boost::bind(boost::factory<mueeWBFWW*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFtautau240"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFtautau250"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFtautau350"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFtautau365"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFtautau380"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFtautau500"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFtautau1000"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFtautau1400"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFtautau1500"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFtautau3000"] = boost::bind(boost::factory<mueeWBFtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFZZ240"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFZZ250"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFZZ350"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFZZ365"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFZZ380"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFZZ500"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFZZ1000"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFZZ1400"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFZZ1500"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFZZ3000"] = boost::bind(boost::factory<mueeWBFZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFZga240"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFZga250"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFZga350"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFZga365"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFZga380"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFZga500"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFZga1000"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFZga1400"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFZga1500"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFZga3000"] = boost::bind(boost::factory<mueeWBFZga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFgaga240"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFgaga250"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFgaga350"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFgaga365"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFgaga380"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFgaga500"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFgaga1000"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFgaga1400"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFgaga1500"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFgaga3000"] = boost::bind(boost::factory<mueeWBFgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeWBFmumu240"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeWBFmumu250"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeWBFmumu350"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeWBFmumu365"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeWBFmumu380"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeWBFmumu500"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeWBFmumu1000"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeWBFmumu1400"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeWBFmumu1500"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeWBFmumu3000"] = boost::bind(boost::factory<mueeWBFmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
    // vvH
    obsThFactory["mueeHvvbb240"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvbb250"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvbb350"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvbb365"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvbb380"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvbb500"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvbb1000"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvbb1400"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvbb1500"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvbb3000"] = boost::bind(boost::factory<mueeHvvbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvbb250_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvbb250_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvbb250_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvbb250_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb350_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvbb350_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvbb350_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvbb350_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb380_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvbb380_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvbb380_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvbb380_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb500_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvbb500_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvbb500_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvbb500_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb1400_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvbb1400_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvbb1400_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvbb1400_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb1500_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvbb1500_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvbb1500_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvbb1500_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvbb3000_p80_m30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvbb3000_m80_p30"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvbb3000_p80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvbb3000_m80_0"] = boost::bind(boost::factory<mueeHvvbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc240"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvcc250"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvcc350"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvcc365"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvcc380"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvcc500"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvcc1000"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvcc1400"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvcc1500"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvcc3000"] = boost::bind(boost::factory<mueeHvvcc*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvcc250_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvcc250_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvcc250_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvcc250_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc350_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvcc350_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvcc350_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvcc350_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc380_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvcc380_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvcc380_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvcc380_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc500_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvcc500_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvcc500_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvcc500_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc1400_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvcc1400_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvcc1400_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvcc1400_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc1500_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvcc1500_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvcc1500_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvcc1500_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvcc3000_p80_m30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvcc3000_m80_p30"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvcc3000_p80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvcc3000_m80_0"] = boost::bind(boost::factory<mueeHvvccPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg240"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvgg250"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvgg350"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvgg365"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvgg380"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvgg500"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvgg1000"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvgg1400"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvgg1500"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvgg3000"] = boost::bind(boost::factory<mueeHvvgg*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvgg250_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvgg250_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvgg250_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvgg250_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg350_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvgg350_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvgg350_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvgg350_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg380_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvgg380_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvgg380_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvgg380_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg500_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvgg500_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvgg500_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvgg500_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg1400_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvgg1400_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvgg1400_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvgg1400_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg1500_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvgg1500_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvgg1500_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvgg1500_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgg3000_p80_m30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvgg3000_m80_p30"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvgg3000_p80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvgg3000_m80_0"] = boost::bind(boost::factory<mueeHvvggPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW240"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvWW250"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvWW350"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvWW365"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvWW380"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvWW500"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvWW1000"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvWW1400"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvWW1500"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvWW3000"] = boost::bind(boost::factory<mueeHvvWW*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvWW250_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvWW250_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvWW250_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvWW250_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW350_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvWW350_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvWW350_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvWW350_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW380_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvWW380_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvWW380_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvWW380_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW500_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvWW500_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvWW500_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvWW500_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW1400_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvWW1400_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvWW1400_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvWW1400_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW1500_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvWW1500_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvWW1500_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvWW1500_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvWW3000_p80_m30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvWW3000_m80_p30"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvWW3000_p80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvWW3000_m80_0"] = boost::bind(boost::factory<mueeHvvWWPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau240"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvtautau250"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvtautau350"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvtautau365"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvtautau380"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvtautau500"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvtautau1000"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvtautau1400"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvtautau1500"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvtautau3000"] = boost::bind(boost::factory<mueeHvvtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvtautau250_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau250_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau250_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvtautau250_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau350_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau350_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau350_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvtautau350_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau380_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau380_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau380_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvtautau380_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau500_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau500_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau500_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvtautau500_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau1400_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau1400_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau1400_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvtautau1400_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau1500_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau1500_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau1500_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvtautau1500_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvtautau3000_p80_m30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvtautau3000_m80_p30"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvtautau3000_p80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvtautau3000_m80_0"] = boost::bind(boost::factory<mueeHvvtautauPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ240"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvZZ250"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvZZ350"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvZZ365"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvZZ380"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvZZ500"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvZZ1000"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvZZ1400"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvZZ1500"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvZZ3000"] = boost::bind(boost::factory<mueeHvvZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvZZ250_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ250_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ250_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvZZ250_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ350_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ350_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ350_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvZZ350_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ380_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ380_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ380_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvZZ380_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ500_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ500_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ500_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvZZ500_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ1400_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ1400_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ1400_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvZZ1400_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ1500_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ1500_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ1500_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvZZ1500_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZZ3000_p80_m30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvZZ3000_m80_p30"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvZZ3000_p80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvZZ3000_m80_0"] = boost::bind(boost::factory<mueeHvvZZPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga240"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvZga250"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvZga350"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvZga365"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvZga380"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvZga500"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvZga1000"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvZga1400"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvZga1500"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvZga3000"] = boost::bind(boost::factory<mueeHvvZga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvZga250_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvZga250_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvZga250_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvZga250_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga350_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvZga350_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvZga350_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvZga350_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga380_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvZga380_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvZga380_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvZga380_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga500_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvZga500_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvZga500_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvZga500_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga1400_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvZga1400_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvZga1400_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvZga1400_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga1500_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvZga1500_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvZga1500_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvZga1500_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvZga3000_p80_m30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvZga3000_m80_p30"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvZga3000_p80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvZga3000_m80_0"] = boost::bind(boost::factory<mueeHvvZgaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga240"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvgaga250"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvgaga350"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvgaga365"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvgaga380"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvgaga500"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvgaga1000"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvgaga1400"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvgaga1500"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvgaga3000"] = boost::bind(boost::factory<mueeHvvgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvgaga250_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga250_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga250_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvgaga250_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga350_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga350_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga350_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvgaga350_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga380_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga380_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga380_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvgaga380_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga500_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga500_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga500_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvgaga500_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga1400_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga1400_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga1400_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvgaga1400_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga1500_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga1500_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga1500_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvgaga1500_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvgaga3000_p80_m30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvgaga3000_m80_p30"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvgaga3000_p80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvgaga3000_m80_0"] = boost::bind(boost::factory<mueeHvvgagaPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu240"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeHvvmumu250"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeHvvmumu350"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeHvvmumu365"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeHvvmumu380"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeHvvmumu500"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeHvvmumu1000"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeHvvmumu1400"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeHvvmumu1500"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeHvvmumu3000"] = boost::bind(boost::factory<mueeHvvmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeHvvmumu250_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu250_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu250_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeHvvmumu250_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu350_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu350_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu350_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeHvvmumu350_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu380_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu380_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu380_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeHvvmumu380_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu500_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu500_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu500_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeHvvmumu500_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu1400_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu1400_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu1400_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeHvvmumu1400_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu1500_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu1500_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu1500_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeHvvmumu1500_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeHvvmumu3000_p80_m30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeHvvmumu3000_m80_p30"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeHvvmumu3000_p80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeHvvmumu3000_m80_0"] = boost::bind(boost::factory<mueeHvvmumuPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // ZH
    obsThFactory["mueeZHbb240"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHbb250"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHbb350"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHbb365"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHbb380"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHbb500"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHbb1000"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHbb1400"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHbb1500"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHbb3000"] = boost::bind(boost::factory<mueeZHbb*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHbb250_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHbb250_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHbb250_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHbb250_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb350_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHbb350_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHbb350_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHbb350_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb380_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHbb380_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHbb380_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHbb380_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb500_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHbb500_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHbb500_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHbb500_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1000_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, -pol_30);
    obsThFactory["mueeZHbb1000_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_30);
    obsThFactory["mueeZHbb1000_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, pol_80, pol_0);
    obsThFactory["mueeZHbb1000_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1400_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZHbb1400_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZHbb1400_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZHbb1400_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb1500_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZHbb1500_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZHbb1500_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZHbb1500_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHbb3000_p80_m30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZHbb3000_m80_p30"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZHbb3000_p80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mueeZHbb3000_m80_0"] = boost::bind(boost::factory<mueeZHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc240"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHcc250"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHcc350"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHcc365"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHcc380"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHcc500"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHcc1000"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHcc1400"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHcc1500"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHcc3000"] = boost::bind(boost::factory<mueeZHcc*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHcc250_p80_m30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHcc250_m80_p30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHcc250_p80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHcc250_m80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc350_p80_m30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHcc350_m80_p30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHcc350_p80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHcc350_m80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc380_p80_m30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHcc380_m80_p30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHcc380_p80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHcc380_m80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHcc500_p80_m30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHcc500_m80_p30"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHcc500_p80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHcc500_m80_0"] = boost::bind(boost::factory<mueeZHccPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg240"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHgg250"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHgg350"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHgg365"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHgg380"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHgg500"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHgg1000"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHgg1400"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHgg1500"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHgg3000"] = boost::bind(boost::factory<mueeZHgg*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHgg250_p80_m30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHgg250_m80_p30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHgg250_p80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHgg250_m80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg350_p80_m30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHgg350_m80_p30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHgg350_p80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHgg350_m80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg380_p80_m30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHgg380_m80_p30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHgg380_p80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHgg380_m80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgg500_p80_m30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHgg500_m80_p30"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHgg500_p80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHgg500_m80_0"] = boost::bind(boost::factory<mueeZHggPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW240"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHWW250"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHWW350"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHWW365"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHWW380"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHWW500"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHWW1000"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHWW1400"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHWW1500"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHWW3000"] = boost::bind(boost::factory<mueeZHWW*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHWW250_p80_m30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHWW250_m80_p30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHWW250_p80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHWW250_m80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW350_p80_m30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHWW350_m80_p30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHWW350_p80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHWW350_m80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW380_p80_m30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHWW380_m80_p30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHWW380_p80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHWW380_m80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHWW500_p80_m30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHWW500_m80_p30"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHWW500_p80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHWW500_m80_0"] = boost::bind(boost::factory<mueeZHWWPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau240"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHtautau250"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHtautau350"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHtautau365"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHtautau380"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHtautau500"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHtautau1000"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHtautau1400"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHtautau1500"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHtautau3000"] = boost::bind(boost::factory<mueeZHtautau*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHtautau250_p80_m30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHtautau250_m80_p30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHtautau250_p80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHtautau250_m80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau350_p80_m30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHtautau350_m80_p30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHtautau350_p80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHtautau350_m80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau380_p80_m30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHtautau380_m80_p30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHtautau380_p80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHtautau380_m80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHtautau500_p80_m30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHtautau500_m80_p30"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHtautau500_p80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHtautau500_m80_0"] = boost::bind(boost::factory<mueeZHtautauPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ240"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHZZ250"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHZZ350"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHZZ365"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHZZ380"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHZZ500"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHZZ1000"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHZZ1400"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHZZ1500"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHZZ3000"] = boost::bind(boost::factory<mueeZHZZ*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHZZ250_p80_m30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHZZ250_m80_p30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHZZ250_p80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHZZ250_m80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ350_p80_m30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHZZ350_m80_p30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHZZ350_p80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHZZ350_m80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ380_p80_m30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHZZ380_m80_p30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHZZ380_p80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHZZ380_m80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZZ500_p80_m30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHZZ500_m80_p30"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHZZ500_p80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHZZ500_m80_0"] = boost::bind(boost::factory<mueeZHZZPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga240"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHZga250"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHZga350"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHZga365"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHZga380"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHZga500"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHZga1000"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHZga1400"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHZga1500"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHZga3000"] = boost::bind(boost::factory<mueeZHZga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHZga250_p80_m30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHZga250_m80_p30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHZga250_p80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHZga250_m80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga350_p80_m30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHZga350_m80_p30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHZga350_p80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHZga350_m80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga380_p80_m30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHZga380_m80_p30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHZga380_p80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHZga380_m80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHZga500_p80_m30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHZga500_m80_p30"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHZga500_p80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHZga500_m80_0"] = boost::bind(boost::factory<mueeZHZgaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga240"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHgaga250"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHgaga350"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHgaga365"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHgaga380"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHgaga500"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHgaga1000"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHgaga1400"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHgaga1500"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHgaga3000"] = boost::bind(boost::factory<mueeZHgaga*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHgaga250_p80_m30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHgaga250_m80_p30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHgaga250_p80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHgaga250_m80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga350_p80_m30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHgaga350_m80_p30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHgaga350_p80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHgaga350_m80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga380_p80_m30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHgaga380_m80_p30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHgaga380_p80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHgaga380_m80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHgaga500_p80_m30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHgaga500_m80_p30"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHgaga500_p80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHgaga500_m80_0"] = boost::bind(boost::factory<mueeZHgagaPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu240"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHmumu250"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHmumu350"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHmumu365"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHmumu380"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHmumu500"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_500);
    obsThFactory["mueeZHmumu1000"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1000);
    obsThFactory["mueeZHmumu1400"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1400);
    obsThFactory["mueeZHmumu1500"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_1500);
    obsThFactory["mueeZHmumu3000"] = boost::bind(boost::factory<mueeZHmumu*>(), _1, sqrt_s_leptcoll_3000);
    //
    obsThFactory["mueeZHmumu250_p80_m30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHmumu250_m80_p30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHmumu250_p80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHmumu250_m80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu350_p80_m30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHmumu350_m80_p30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHmumu350_p80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHmumu350_m80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu380_p80_m30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHmumu380_m80_p30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHmumu380_p80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHmumu380_m80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHmumu500_p80_m30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHmumu500_m80_p30"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHmumu500_p80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHmumu500_m80_0"] = boost::bind(boost::factory<mueeZHmumuPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv240"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHinv250"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHinv350"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHinv365"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHinv380"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHinv500"] = boost::bind(boost::factory<mueeZHinv*>(), _1, sqrt_s_leptcoll_500);
    //
    obsThFactory["mueeZHinv250_p80_m30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHinv250_m80_p30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHinv250_p80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHinv250_m80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv350_p80_m30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHinv350_m80_p30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHinv350_p80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHinv350_m80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv380_p80_m30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHinv380_m80_p30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHinv380_p80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHinv380_m80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHinv500_p80_m30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHinv500_m80_p30"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHinv500_p80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHinv500_m80_0"] = boost::bind(boost::factory<mueeZHinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv240"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_240);
    obsThFactory["mueeZHBRinv250"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_250);
    obsThFactory["mueeZHBRinv350"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_350);
    obsThFactory["mueeZHBRinv365"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_365);
    obsThFactory["mueeZHBRinv380"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_380);
    obsThFactory["mueeZHBRinv500"] = boost::bind(boost::factory<mueeZHBRinv*>(), _1, sqrt_s_leptcoll_500);
    //
    obsThFactory["mueeZHBRinv250_p80_m30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv250_m80_p30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv250_p80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, pol_80, pol_0);
    obsThFactory["mueeZHBRinv250_m80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_250, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv350_p80_m30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv350_m80_p30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv350_p80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, pol_80, pol_0);
    obsThFactory["mueeZHBRinv350_m80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_350, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv380_p80_m30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv380_m80_p30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv380_p80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, pol_80, pol_0);
    obsThFactory["mueeZHBRinv380_m80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_380, -pol_80, pol_0);
    //
    obsThFactory["mueeZHBRinv500_p80_m30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, -pol_30);
    obsThFactory["mueeZHBRinv500_m80_p30"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_30);
    obsThFactory["mueeZHBRinv500_p80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, pol_80, pol_0);
    obsThFactory["mueeZHBRinv500_m80_0"] = boost::bind(boost::factory<mueeZHBRinvPol*>(), _1, sqrt_s_leptcoll_500, -pol_80, pol_0);
    //
    // ZBF
    obsThFactory["mueeZBFbb1400_p80_m30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueeZBFbb1400_m80_p30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueeZBFbb1400_p80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueeZBFbb1400_m80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb1500_p80_m30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueeZBFbb1500_m80_p30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueeZBFbb1500_p80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueeZBFbb1500_m80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueeZBFbb3000_p80_m30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueeZBFbb3000_m80_p30"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueeZBFbb3000_p80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mmueeZBFbb3000_m80_0"] = boost::bind(boost::factory<mueeZBFbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    // eettH
    obsThFactory["mueettHbb1400_p80_m30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, -pol_30);
    obsThFactory["mueettHbb1400_m80_p30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_30);
    obsThFactory["mueettHbb1400_p80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, pol_80, pol_0);
    obsThFactory["mueettHbb1400_m80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1400, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb1500_p80_m30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, -pol_30);
    obsThFactory["mueettHbb1500_m80_p30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_30);
    obsThFactory["mueettHbb1500_p80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, pol_80, pol_0);
    obsThFactory["mueettHbb1500_m80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_1500, -pol_80, pol_0);
    //
    obsThFactory["mueettHbb3000_p80_m30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, -pol_30);
    obsThFactory["mueettHbb3000_m80_p30"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_30);
    obsThFactory["mueettHbb3000_p80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, pol_80, pol_0);
    obsThFactory["mmueettHbb3000_m80_0"] = boost::bind(boost::factory<mueettHbbPol*>(), _1, sqrt_s_leptcoll_3000, -pol_80, pol_0);
    //
    //-----  Full Signal strengths per prod and decay: Lepton-Hadron colliders  ----------
    //
    obsThFactory["epWBFbb1300"] = boost::bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFcc1300"] = boost::bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFgg1300"] = boost::bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFWW2l2v1300"] = boost::bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFZZ4l1300"] = boost::bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFgaga1300"] = boost::bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epWBFtautau1300"] = boost::bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_LHeC_1_3);
    //
    obsThFactory["epWBFbb3500"] = boost::bind(boost::factory<muepWBFbb*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFcc3500"] = boost::bind(boost::factory<muepWBFcc*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFgg3500"] = boost::bind(boost::factory<muepWBFgg*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFWW2l2v3500"] = boost::bind(boost::factory<muepWBFWW2l2v*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFZZ4l3500"] = boost::bind(boost::factory<muepWBFZZ4l*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFgaga3500"] = boost::bind(boost::factory<muepWBFgaga*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epWBFtautau3500"] = boost::bind(boost::factory<muepWBFtautau*>(), _1, sqrt_s_FCCep_3_5);
    //
    obsThFactory["epZBFbb1300"] = boost::bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFcc1300"] = boost::bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFgg1300"] = boost::bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFWW2l2v1300"] = boost::bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFZZ4l1300"] = boost::bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFgaga1300"] = boost::bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_LHeC_1_3);
    obsThFactory["epZBFtautau1300"] = boost::bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_LHeC_1_3);
    //
    obsThFactory["epZBFbb3500"] = boost::bind(boost::factory<muepZBFbb*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFcc3500"] = boost::bind(boost::factory<muepZBFcc*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFgg3500"] = boost::bind(boost::factory<muepZBFgg*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFWW2l2v3500"] = boost::bind(boost::factory<muepZBFWW2l2v*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFZZ4l3500"] = boost::bind(boost::factory<muepZBFZZ4l*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFgaga3500"] = boost::bind(boost::factory<muepZBFgaga*>(), _1, sqrt_s_FCCep_3_5);
    obsThFactory["epZBFtautau3500"] = boost::bind(boost::factory<muepZBFtautau*>(), _1, sqrt_s_FCCep_3_5);
    //
    //-----  Limits  ----------
    obsThFactory["UpperLimit_ppHZgammaA13"] = boost::bind(boost::factory<UpperLimit_ppHZgammaA13*>(), _1, sqrt_s_LHC13);
    obsThFactory["UpperLimit_ppHZgammaC13"] = boost::bind(boost::factory<UpperLimit_ppHZgammaC13*>(), _1, sqrt_s_LHC13);
    obsThFactory["UpperLimit_ppHZgammaA"] = boost::bind(boost::factory<UpperLimit_ppHZgammaA*>(), _1, sqrt_s_LHC8);
    obsThFactory["UpperLimit_ppHZgammaC"] = boost::bind(boost::factory<UpperLimit_ppHZgammaC*>(), _1, sqrt_s_LHC8);
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
    obsThFactory["C_tW"] = boost::factory<C_tW*>();
    obsThFactory["C_tB"] = boost::factory<C_tB*>();
    obsThFactory["C_tphi"] = boost::factory<C_tphi*>();
    obsThFactory["C_phib"] = boost::factory<C_phib*>();
    obsThFactory["C_bW"] = boost::factory<C_bW*>();
    obsThFactory["C_bB"] = boost::factory<C_bB*>();
    
    obsThFactory["Rb_NPSMEFT6dtopquark"] = boost::factory<Rb_NPSMEFT6dtopquark*>();
    obsThFactory["AFBLR"] = boost::factory<AFBLR*>();
    obsThFactory["SigmattZ"] = boost::factory<sigmattZ*>();
    obsThFactory["SigmattA_1"] = boost::factory<sigmattA_1*>();
    obsThFactory["SigmattA_2"] = boost::factory<sigmattA_2*>();
    obsThFactory["SigmattH"] = boost::factory<sigmattH*>();
    obsThFactory["SigmattW"] = boost::factory<sigmattW*>();
    obsThFactory["Sigmatq"] = boost::factory<sigmatq*>();
    obsThFactory["SigmatW"] = boost::factory<sigmatW*>();
    obsThFactory["SigmatqZ"] = boost::factory<sigmatqZ*>();
    obsThFactory["F0"] = boost::factory<F0*>();
    obsThFactory["FL"] = boost::factory<FL*>();
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

    /* BEGIN: REMOVE FROM THE PACKAGE */
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
    /* END: REMOVE FROM THE PACKAGE */

    //-----  Flavour observables  -----
    //----- DF = 2  -----
    obsThFactory["DmBd"] = boost::factory<DmBd*>();
    obsThFactory["DmBs"] = boost::factory<DmBs*>();
    obsThFactory["RmBs"] = boost::factory<RmBs*>();
    obsThFactory["SJPsiK"] = boost::factory<SJPsiK*>();
    obsThFactory["Betas_JPsiPhi"] = boost::factory<Betas_JPsiPhi*>();
    obsThFactory["EpsilonK"] = boost::factory<EpsilonK*>();
    obsThFactory["DmK"] = boost::factory<DmK*>();
    /* BEGIN: REMOVE FROM THE PACKAGE */
    obsThFactory["M12D"] = boost::factory<M12D*>();
    obsThFactory["ArgD"] = boost::factory<ArgD*>();
    //----- eps'/eps  -----
    obsThFactory["EpsiloP_o_Epsilon"] = boost::factory<EpsilonP_O_Epsilon*>();
    /* END: REMOVE FROM THE PACKAGE */
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

/* BEGIN: REMOVE FROM THE PACKAGE */
    //----- B to X_q ll -----
    obsThFactory["R_BXsee"] = boost::bind(boost::factory<R_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["Rlow_BXsee"] = boost::bind(boost::factory<Rlow_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
    obsThFactory["Rhigh_BXsee"] = boost::bind(boost::factory<Rhigh_BXqll*>(), _1, StandardModel::STRANGE, StandardModel::ELECTRON);
/* END: REMOVE FROM THE PACKAGE */

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
    obsThFactory["OutputSLHAfromFH"] = boost::factory<OutputSLHAfromFH*>(); // for debug
/* END: REMOVE FROM THE PACKAGE */
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

   obsThFactory["yu1R_GTHDM"] = boost::factory<yu1R_GTHDM*>();
   obsThFactory["yd1R_GTHDM"] = boost::factory<yd1R_GTHDM*>();
   obsThFactory["yl1R_GTHDM"] = boost::factory<yl1R_GTHDM*>();

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
    obsThFactory["THDMWunitarity1"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 0);
    obsThFactory["THDMWunitarity2"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 1);
    obsThFactory["THDMWunitarity3"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 2);
    obsThFactory["THDMWunitarity4"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 3);
    obsThFactory["THDMWunitarity5"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 4);
    obsThFactory["THDMWunitarity6"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 5);
    obsThFactory["THDMWunitarity7"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 6);
    obsThFactory["THDMWunitarity8"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 7);
    obsThFactory["THDMWunitarity9"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 8);
    obsThFactory["THDMWunitarity10"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 9);
    obsThFactory["THDMWunitarity11"] = boost::bind(boost::factory<THDMWunitarityLO*>(), _1, 10);
    //-----  One-loop unitarity constraints  -----
    obsThFactory["THDMWNLOunitarity1"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 0);
    obsThFactory["THDMWNLOunitarity2"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 1);
    obsThFactory["THDMWNLOunitarity3"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 2);
    obsThFactory["THDMWNLOunitarity4"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 3);
    obsThFactory["THDMWNLOunitarity5"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 4);
    obsThFactory["THDMWNLOunitarity6"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 5);
    obsThFactory["THDMWNLOunitarity7"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 6);
    obsThFactory["THDMWNLOunitarity8"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 7);
    obsThFactory["THDMWNLOunitarity9"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 8);
    obsThFactory["THDMWNLOunitarity10"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 9);
    obsThFactory["THDMWNLOunitarity11"] = boost::bind(boost::factory<THDMWunitarityNLO*>(), _1, 10);
    //-----  One-loop "plus" unitarity constraints  -----
    obsThFactory["THDMWNLOpunitarity1"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 0);
    obsThFactory["THDMWNLOpunitarity2"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 1);
    obsThFactory["THDMWNLOpunitarity3"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 2);
    obsThFactory["THDMWNLOpunitarity4"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 3);
    obsThFactory["THDMWNLOpunitarity5"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 4);
    obsThFactory["THDMWNLOpunitarity6"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 5);
    obsThFactory["THDMWNLOpunitarity7"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 6);
    obsThFactory["THDMWNLOpunitarity8"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 7);
    obsThFactory["THDMWNLOpunitarity9"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 8);
    obsThFactory["THDMWNLOpunitarity10"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 9);
    obsThFactory["THDMWNLOpunitarity11"] = boost::bind(boost::factory<THDMWunitarityNLOp*>(), _1, 10);
    //-----   R' criteria for perturbative unitarity  -----
    obsThFactory["THDMWunitarityRp1"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 0);
    obsThFactory["THDMWunitarityRp2"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 1);
    obsThFactory["THDMWunitarityRp3"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 2);
    obsThFactory["THDMWunitarityRp4"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 3);
    obsThFactory["THDMWunitarityRp5"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 4);
    obsThFactory["THDMWunitarityRp6"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 5);
    obsThFactory["THDMWunitarityRp7"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 6);
    obsThFactory["THDMWunitarityRp8"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 7);
    obsThFactory["THDMWunitarityRp9"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 8);
    obsThFactory["THDMWunitarityRp10"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 9);
    obsThFactory["THDMWunitarityRp11"] = boost::bind(boost::factory<THDMWunitarityRp*>(), _1, 10);
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
    obsThFactory["GMunitarity1"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 0);
    obsThFactory["GMunitarity2"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 1);
    obsThFactory["GMunitarity3"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 2);
    obsThFactory["GMunitarity4"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 3);
    obsThFactory["GMunitarity5"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 4);
    obsThFactory["GMunitarity6"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 5);
    obsThFactory["GMunitarity7"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 6);
    obsThFactory["GMunitarity8"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 7);
    obsThFactory["GMunitarity9"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 8);
    obsThFactory["GMunitarity10"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 9);
    obsThFactory["GMunitarity11"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 10);
    obsThFactory["GMunitarity12"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 11);
    obsThFactory["GMunitarity13"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 12);
    obsThFactory["GMunitarity14"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 13);
    obsThFactory["GMunitarity15"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 14);
    obsThFactory["GMunitarity16"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 15);
    obsThFactory["GMunitarity17"] = boost::bind(boost::factory<GMunitarityLO*>(), _1, 16);
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
