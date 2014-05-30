/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <vector>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include <EWObservables.h>
#include <HiggsThObservables.h>
#include <ParamObs.h>
#include "ThFactory.h"

ThFactory::ThFactory(const StandardModel& myModel)
{
    //-----  EW precision observables  -----
    thobs["Mw"] = new Mw(myModel);
    thobs["GammaW"] = new GammaW(myModel);
    thobs["GammaZ"] = new GammaZ(myModel);
    thobs["sigmaHadron"] = new sigmaHadron(myModel);
    thobs["sin2thetaEff"] = new sin2thetaEff(myModel);
    thobs["PtauPol"] = new PtauPol(myModel);
    thobs["Alepton"] = new Alepton(myModel);
    thobs["Acharm"] = new Acharm(myModel);
    thobs["Abottom"] = new Abottom(myModel);
    thobs["AFBlepton"] = new AFBlepton(myModel);
    thobs["AFBcharm"] = new AFBcharm(myModel);
    thobs["AFBbottom"] = new AFBbottom(myModel);
    thobs["Rlepton"] = new Rlepton(myModel);
    thobs["Rcharm"] = new Rcharm(myModel);
    thobs["Rbottom"] = new Rbottom(myModel);

    //----- Higgs Extension observables ----------

    if (myModel.ModelName().compare(0, 5, "Higgs") == 0) {
        thobs["ggH"] = new muggH(myModel);
        thobs["VBF"] = new muVBF(myModel);
        thobs["VH"] = new muWH(myModel);
        thobs["ttH"] = new muttH(myModel);
        thobs["BrHWW"] = new BrWW(myModel);
        thobs["BrHZZ"] = new BrZZ(myModel);
        thobs["BrHgaga"] = new Brgaga(myModel);
        thobs["BrHtautau"] = new Brtautau(myModel);
    }

    //-----  Epsilon parameters  -----
    thobs["epsilon1"] = new Epsilon1(myModel);
    thobs["epsilon2"] = new Epsilon2(myModel);
    thobs["epsilon3"] = new Epsilon3(myModel);
    thobs["epsilonb"] = new Epsilonb(myModel);

    //-----  LEP-II two-fermion processes  -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189.,
        192., 196., 200., 202., 205., 207.};
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192.,
        196., 200., 202., 205., 207.};
    LEP2sigmaHadron * myLEP2sigmaHadron[12];
    LEP2sigmaMu * myLEP2sigmaMu[12];
    LEP2sigmaTau * myLEP2sigmaTau[12];
    LEP2AFBmu * myLEP2AFBmu[12];
    LEP2AFBtau * myLEP2AFBtau[12];
    LEP2AFBbottom * myLEP2AFBbottom[10];
    LEP2AFBcharm * myLEP2AFBcharm[10];
    LEP2Rbottom * myLEP2Rbottom[10];
    LEP2Rcharm * myLEP2Rcharm[10];
    for (int i = 0; i < 12; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
        myLEP2sigmaHadron[i] = new LEP2sigmaHadron(myModel, sqrt_s[i]);
        thobs["sigmaqLEP2_" + sqrt_s_str] = myLEP2sigmaHadron[i];
        myLEP2sigmaMu[i] = new LEP2sigmaMu(myModel, sqrt_s[i]);
        thobs["sigmamuLEP2_" + sqrt_s_str] = myLEP2sigmaMu[i];
        myLEP2sigmaTau[i] = new LEP2sigmaTau(myModel, sqrt_s[i]);
        thobs["sigmatauLEP2_" + sqrt_s_str] = myLEP2sigmaTau[i];
        myLEP2AFBmu[i] = new LEP2AFBmu(myModel, sqrt_s[i]);
        thobs["AFBmuLEP2_" + sqrt_s_str] = myLEP2AFBmu[i];
        myLEP2AFBtau[i] = new LEP2AFBtau(myModel, sqrt_s[i]);
        thobs["AFBtauLEP2_" + sqrt_s_str] = myLEP2AFBtau[i];
    }
    for (int i = 0; i < 10; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
        myLEP2AFBbottom[i] = new LEP2AFBbottom(myModel, sqrt_s_HF[i]);
        thobs["AFBbottomLEP2_" + sqrt_s_str] = myLEP2AFBbottom[i];
        myLEP2AFBcharm[i] = new LEP2AFBcharm(myModel, sqrt_s_HF[i]);
        thobs["AFBcharmLEP2_" + sqrt_s_str] = myLEP2AFBcharm[i];
        myLEP2Rbottom[i] = new LEP2Rbottom(myModel, sqrt_s_HF[i]);
        thobs["RbottomLEP2_" + sqrt_s_str] = myLEP2Rbottom[i];
        myLEP2Rcharm[i] = new LEP2Rcharm(myModel, sqrt_s_HF[i]);
        thobs["RcharmLEP2_" + sqrt_s_str] = myLEP2Rcharm[i];
    }

//    //-----  Flavour observables  -----
//    thobs["Dmd1"] = new DmBd(myFlavour);
//    thobs["Dms1"] = new DmBs(myFlavour);
//    thobs["M12D"] = new M12D(myFlavour);
//    thobs["ArgD"] = new ArgD(myFlavour);
//    thobs["EpsilonK"] = new EpsilonK(myFlavour);
//    thobs["EpsiloP_o_Epsilon"] = new EpsilonP_O_Epsilon(myFlavour);
//    thobs["DmK"] = new DmK(myFlavour);
//    thobs["Vud"] = new Vud(myFlavour);
//    thobs["Vus"] = new Vus(myFlavour);
//    thobs["Vub"] = new Vub(myFlavour);
//    thobs["Vcb"] = new Vcb(myFlavour);
//    thobs["alpha"] = new Alpha(myFlavour);
//    thobs["alpha_2a"] = new Alpha_2a(myFlavour);
//    thobs["gamma"] = new CKMGamma(myFlavour);
//    thobs["SJPsiK"] = new SJPsiK(myFlavour);
//    thobs["SJPsiPhi"] = new SJPsiPhi(myFlavour);
//    thobs["BR_Bdmumu"] = new Bdmumu(myFlavour, 1);
//    thobs["BRbar_Bdmumu"] = new Bdmumu(myFlavour, 2);
//    thobs["Amumu_Bd"] = new Bdmumu(myFlavour, 3);
//    thobs["Smumu_Bd"] = new Bdmumu(myFlavour, 4);
//    thobs["BR_Bsmumu"] = new Bsmumu(myFlavour, 1);
//    thobs["BRbar_Bsmumu"] = new Bsmumu(myFlavour, 2);
//    thobs["Amumu_Bs"] = new Bsmumu(myFlavour, 3);
//    thobs["Smumu_Bs"] = new Bsmumu(myFlavour, 4);
//
//    //-----  Lepton Flavour observables  -----
//    thobs["li_lj_gamma"] = new li_lj_gamma(myLeptonFlavour);

    //-----  SM input parameters  -----
//    thobs["Mz"] = new StandardModelParams(myMO, "Mz");
//    thobs["mHl"] = new StandardModelParams(myMO, "mHl");

    //-----  NP input parameters, etc.   -----
//    if (myModel.ModelName().compare("NPEffective1") == 0
//            || myModel.ModelName().compare("NPEffective2") == 0) {
//        thobs["cHQ1pPLUScHQ2p_NP"] = new NewPhysicsParams(myMO, "cHQ1pPLUScHQ2p_NP");
//        thobs["cHQ2pMINUScHQ2_NP"] = new NewPhysicsParams(myMO, "cHQ2pMINUScHQ2_NP");
//        thobs["cHQ3pPLUScHQ3_NP"] = new NewPhysicsParams(myMO, "cHQ3pPLUScHQ3_NP");
//        thobs["c_Ae_NP"] = new NewPhysicsParams(myMO, "c_Ae_NP");
//        thobs["c_GammaZ_uds_NP"] = new NewPhysicsParams(myMO, "c_GammaZ_uds_NP");
//    }
//    if ((myModel.ModelName().compare("NPZbbbar") == 0)
//            || (myModel.ModelName().compare("NPZbbbarLR") == 0)) {
//        thobs["deltaGVb"] = new NewPhysicsParams(myMO, "deltaGVb");
//        thobs["deltaGAb"] = new NewPhysicsParams(myMO, "deltaGAb");
//        thobs["deltaGLb"] = new NewPhysicsParams(myMO, "deltaGLb");
//        thobs["deltaGRb"] = new NewPhysicsParams(myMO, "deltaGRb");
//        thobs["deltaRhoZb"] = new NewPhysicsParams(myMO, "deltaRhoZb");
//        thobs["deltaKappaZb"] = new NewPhysicsParams(myMO, "deltaKappaZb");
//    }
//
//    //-----  SUSY spectra and observables  -----
//    if (myModel.ModelName().compare("SUSY") == 0
//            || myModel.ModelName().compare("SUSYMassInsertion") == 0
//            || myModel.ModelName().compare("GeneralSUSY") == 0
//            || myModel.ModelName().compare("pMSSM") == 0
//            || myModel.ModelName().compare("MFV") == 0) {
//        thobs["OutputSLHAfromFH"] = new OutputSLHAfromFH(myMO); // for debug
//        thobs["MHl"] = new Mhiggs(myMO, 0);
//        thobs["MHh"] = new Mhiggs(myMO, 1);
//        thobs["MHa"] = new Mhiggs(myMO, 2);
//        thobs["MHp"] = new Mhiggs(myMO, 3);
//        thobs["Msu1"] = new Msup(myMO, 0);
//        thobs["Msu2"] = new Msup(myMO, 1);
//        thobs["Msu3"] = new Msup(myMO, 2);
//        thobs["Msu4"] = new Msup(myMO, 3);
//        thobs["Msu5"] = new Msup(myMO, 4);
//        thobs["Msu6"] = new Msup(myMO, 5);
//        thobs["Msd1"] = new Msdown(myMO, 0);
//        thobs["Msd2"] = new Msdown(myMO, 1);
//        thobs["Msd3"] = new Msdown(myMO, 2);
//        thobs["Msd4"] = new Msdown(myMO, 3);
//        thobs["Msd5"] = new Msdown(myMO, 4);
//        thobs["Msd6"] = new Msdown(myMO, 5);
//        thobs["Mch1"] = new Mchargino(myMO, 0);
//        thobs["Mch2"] = new Mchargino(myMO, 1);
//        thobs["Mneu1"] = new Mneutralino(myMO, 0);
//        thobs["Mneu2"] = new Mneutralino(myMO, 1);
//        thobs["Mneu3"] = new Mneutralino(myMO, 2);
//        thobs["Mneu4"] = new Mneutralino(myMO, 3);
//        thobs["Mw_dRho"] = new Mw_dRho(myMO);
//    }
}

ThFactory::~ThFactory()
{
    for (std::map<std::string, ThObservable *>::iterator it = thobs.begin();
            it != thobs.end(); it++)
        if (it->second != NULL)
            delete it->second;
}

ThObservable * ThFactory::getThMethod(const std::string& name)
{
    if (thobs.find(name) == thobs.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThFactory");
    return (thobs[name]);
}
