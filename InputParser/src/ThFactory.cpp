/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <vector>
#include <boost/lexical_cast.hpp>
#include <EWObservables.h>
#include <FlavourObservables.h>
#include <StandardModelParams.h>
#include <NewPhysicsParams.h>
#include <SUSYObservables.h>
#include "ThFactory.h"

ThFactory::ThFactory(const StandardModel& myModel) 
: myEW(myModel), myFlavour(myModel), myMO(myModel)
{    
    //-----  EW precision observables  -----
    thobs["Mw"] = new Mw(myEW);
    thobs["GammaW"] = new GammaW(myEW);
    thobs["GammaZ"] = new GammaZ(myEW);
    thobs["sigmaHadron"] = new sigmaHadron(myEW);
    thobs["sin2thetaEff"] = new sin2thetaEff(myEW);
    thobs["PtauPol"] = new PtauPol(myEW);
    thobs["Alepton"] = new Alepton(myEW);
    thobs["Acharm"] = new Acharm(myEW);
    thobs["Abottom"] = new Abottom(myEW);
    thobs["AFBlepton"] = new AFBlepton(myEW);
    thobs["AFBcharm"] = new AFBcharm(myEW);
    thobs["AFBbottom"] = new AFBbottom(myEW);
    thobs["Rlepton"] = new Rlepton(myEW);
    thobs["Rcharm"] = new Rcharm(myEW);
    thobs["Rbottom"] = new Rbottom(myEW);

    //-----  Epsilon parameters  -----
    thobs["epsilon1"] = new NewPhysicsParams(myMO, "epsilon1");
    thobs["epsilon2"] = new NewPhysicsParams(myMO, "epsilon2");
    thobs["epsilon3"] = new NewPhysicsParams(myMO, "epsilon3");
    thobs["epsilonb"] = new NewPhysicsParams(myMO, "epsilonb");

    //-----  LEP-II two-fermion processes  -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189.,
                               192., 196., 200., 202., 205., 207.};
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192.,
                                  196., 200., 202., 205., 207.};
    LEP2sigmaHadron* myLEP2sigmaHadron[12];
    LEP2sigmaMu* myLEP2sigmaMu[12];
    LEP2sigmaTau* myLEP2sigmaTau[12];
    LEP2AFBmu* myLEP2AFBmu[12];
    LEP2AFBtau* myLEP2AFBtau[12];
    LEP2AFBbottom* myLEP2AFBbottom[10];
    LEP2AFBcharm* myLEP2AFBcharm[10];
    LEP2Rbottom* myLEP2Rbottom[10];
    LEP2Rcharm* myLEP2Rcharm[10];
    for (int i=0; i<12; i++) { 
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
        myLEP2sigmaHadron[i] = new LEP2sigmaHadron(myEW, sqrt_s[i]);
        thobs["sigmaqLEP2_" + sqrt_s_str] = myLEP2sigmaHadron[i];
        myLEP2sigmaMu[i] = new LEP2sigmaMu(myEW, sqrt_s[i]);
        thobs["sigmamuLEP2_" + sqrt_s_str] = myLEP2sigmaMu[i];
        myLEP2sigmaTau[i] = new LEP2sigmaTau(myEW, sqrt_s[i]);
        thobs["sigmatauLEP2_" + sqrt_s_str] = myLEP2sigmaTau[i];
        myLEP2AFBmu[i] = new LEP2AFBmu(myEW, sqrt_s[i]);
        thobs["AFBmuLEP2_" + sqrt_s_str] = myLEP2AFBmu[i];
        myLEP2AFBtau[i] = new LEP2AFBtau(myEW, sqrt_s[i]);
        thobs["AFBtauLEP2_" + sqrt_s_str] = myLEP2AFBtau[i];
    }
    for (int i=0; i<10; i++) { 
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
        myLEP2AFBbottom[i] = new LEP2AFBbottom(myEW, sqrt_s_HF[i]);
        thobs["AFBbottomLEP2_" + sqrt_s_str] = myLEP2AFBbottom[i];
        myLEP2AFBcharm[i] = new LEP2AFBcharm(myEW, sqrt_s_HF[i]);
        thobs["AFBcharmLEP2_" + sqrt_s_str] = myLEP2AFBcharm[i];
        myLEP2Rbottom[i] = new LEP2Rbottom(myEW, sqrt_s_HF[i]);  
        thobs["RbottomLEP2_" + sqrt_s_str] = myLEP2Rbottom[i];
        myLEP2Rcharm[i] = new LEP2Rcharm(myEW, sqrt_s_HF[i]);
        thobs["RcharmLEP2_" + sqrt_s_str] = myLEP2Rcharm[i];
    }    

    //-----  Flavour observables  -----
    thobs["Dmd1"] = new DmBd(myFlavour);
    thobs["Dms1"] = new DmBs(myFlavour);
    thobs["M12D"] = new M12D(myFlavour);
    thobs["ArgD"] = new ArgD(myFlavour);
    thobs["EpsilonK"] = new EpsilonK(myFlavour);
    thobs["EpsiloP_o_Epsilon"] = new EpsilonP_O_Epsilon(myFlavour);
    thobs["DmK"] = new DmK(myFlavour);
    thobs["Vud"] = new Vud(myFlavour);
    thobs["Vus"] = new Vus(myFlavour);
    thobs["Vub"] = new Vub(myFlavour);
    thobs["Vcb"] = new Vcb(myFlavour);
    thobs["alpha"] = new Alpha(myFlavour);
    thobs["alpha_2a"] = new Alpha_2a(myFlavour);
    thobs["gamma"] = new CKMGamma(myFlavour);
    thobs["SJPsiK"] = new SJPsiK(myFlavour);
    thobs["SJPsiPhi"] = new SJPsiPhi(myFlavour);
    thobs["BR_Bdmumu"] = new Bdmumu(myFlavour, 1);
    thobs["BRbar_Bdmumu"] = new Bdmumu(myFlavour, 2);
    thobs["Amumu_Bd"] = new Bdmumu(myFlavour, 3);
    thobs["Smumu_Bd"] = new Bdmumu(myFlavour, 4);
    thobs["BR_Bsmumu"] = new Bsmumu(myFlavour, 1);
    thobs["BRbar_Bsmumu"] = new Bsmumu(myFlavour, 2);
    thobs["Amumu_Bs"] = new Bsmumu(myFlavour, 3);
    thobs["Smumu_Bs"] = new Bsmumu(myFlavour, 4);

    //-----  SM input parameters, etc.  -----
    thobs["AlsMz"] = new StandardModelParams(myMO, "AlsMz");
    thobs["dAle5Mz"] = new StandardModelParams(myMO, "dAle5Mz");
    thobs["Mz"] = new StandardModelParams(myMO, "Mz");
    thobs["mtop"] = new StandardModelParams(myMO, "mtop");
    thobs["mHl"] = new StandardModelParams(myMO, "mHl");
    thobs["delMw"] = new StandardModelParams(myMO, "delMw");
    thobs["delSin2th_l"] = new StandardModelParams(myMO, "delSin2th_l");
    thobs["delGammaZ"] = new StandardModelParams(myMO, "delGammaZ");
    thobs["delRhoZ_nu"] = new StandardModelParams(myMO, "delRhoZ_nu");
    thobs["delRhoZ_e"] = new StandardModelParams(myMO, "delRhoZ_e");
    if (myMO.getModel().IsFlagTestSubleadingTwoLoopEW()) {
        thobs["delRhoZ_u"] = new StandardModelParams(myMO, "delRhoZ_u");
        thobs["delRhoZ_d"] = new StandardModelParams(myMO, "delRhoZ_d");
    }
    thobs["delRhoZ_b"] = new StandardModelParams(myMO, "delRhoZ_b");

    //-----  NP input parameters, etc.   -----
    if(myModel.ModelName().compare("NPEffective1") == 0
            || myModel.ModelName().compare("NPEffective2") == 0 ) {
        thobs["cHLp_NP"] = new NewPhysicsParams(myMO, "cHLp_NP");
        thobs["cHQp_NP"] = new NewPhysicsParams(myMO, "cHQp_NP");
        thobs["cHQ_NP"] = new NewPhysicsParams(myMO, "cHQ_NP");
        thobs["cHL_NP"] = new NewPhysicsParams(myMO, "cHL_NP");
        thobs["cHE_NP"] = new NewPhysicsParams(myMO, "cHE_NP");
        thobs["c_Ae_NP"] = new NewPhysicsParams(myMO, "c_Ae_NP");
        thobs["c_GammaZ_uds_NP"] = new NewPhysicsParams(myMO, "c_GammaZ_uds_NP");
        thobs["cHU2_NP"] = new NewPhysicsParams(myMO, "cHU2_NP");
        thobs["cHD3_NP"] = new NewPhysicsParams(myMO, "cHD3_NP");
        thobs["cHQ1pPLUScHQ2p_NP"] = new NewPhysicsParams(myMO, "cHQ1pPLUScHQ2p_NP");
        thobs["cHQ2pMINUScHQ2_NP"] = new NewPhysicsParams(myMO, "cHQ2pMINUScHQ2_NP");
        thobs["cHQ3pPLUScHQ3_NP"] = new NewPhysicsParams(myMO, "cHQ3pPLUScHQ3_NP");
    }
    if(myModel.ModelName().compare("NPZbbbar") == 0) {
        thobs["deltaGVb"] = new NewPhysicsParams(myMO, "deltaGVb");
        thobs["deltaGAb"] = new NewPhysicsParams(myMO, "deltaGAb");
        thobs["deltaGLb"] = new NewPhysicsParams(myMO, "deltaGLb");
        thobs["deltaGRb"] = new NewPhysicsParams(myMO, "deltaGRb");
        thobs["deltaRhoZb"] = new NewPhysicsParams(myMO, "deltaRhoZb");
        thobs["deltaKappaZb"] = new NewPhysicsParams(myMO, "deltaKappaZb");
    }

    //-----  SUSY spectra and observables  -----
    if(myModel.ModelName().compare("SUSY") == 0
            || myModel.ModelName().compare("SUSYMassInsertion") == 0
            || myModel.ModelName().compare("GeneralSUSY") == 0
            || myModel.ModelName().compare("pMSSM") == 0
            || myModel.ModelName().compare("MFV") == 0) {
        thobs["OutputSLHAfromFH"] = new OutputSLHAfromFH(myMO); // for debug
        thobs["MHl"] = new Mhiggs(myMO, 0);
        thobs["MHh"] = new Mhiggs(myMO, 1);
        thobs["MHa"] = new Mhiggs(myMO, 2);
        thobs["MHp"] = new Mhiggs(myMO, 3);
        thobs["Msu1"] = new Msup(myMO, 0);
        thobs["Msu2"] = new Msup(myMO, 1);
        thobs["Msu3"] = new Msup(myMO, 2);
        thobs["Msu4"] = new Msup(myMO, 3);
        thobs["Msu5"] = new Msup(myMO, 4);
        thobs["Msu6"] = new Msup(myMO, 5);
        thobs["Msd1"] = new Msdown(myMO, 0);
        thobs["Msd2"] = new Msdown(myMO, 1);
        thobs["Msd3"] = new Msdown(myMO, 2);
        thobs["Msd4"] = new Msdown(myMO, 3);
        thobs["Msd5"] = new Msdown(myMO, 4);
        thobs["Msd6"] = new Msdown(myMO, 5);
        thobs["Mch1"] = new Mchargino(myMO, 0);
        thobs["Mch2"] = new Mchargino(myMO, 1);
        thobs["Mneu1"] = new Mneutralino(myMO, 0);
        thobs["Mneu2"] = new Mneutralino(myMO, 1);
        thobs["Mneu3"] = new Mneutralino(myMO, 2);
        thobs["Mneu4"] = new Mneutralino(myMO, 3);
        thobs["Mw_dRho"] = new Mw_dRho(myMO);
    }
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
    if (thobs.find(name) == thobs.end()) {
        std::cout << "wrong observable " << name << " in ThFactory" << std::endl;
        exit(EXIT_FAILURE);
    }
    return (thobs[name]);
}
