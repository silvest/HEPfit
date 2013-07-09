/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <vector>
#include <boost/lexical_cast.hpp>
#include "ThFactory.h"
#include <SMInputs.h>
#include <NPInputs.h>
#include <SUSYObservables.h>
#include <FlavourObservables.h>
#include <EWObservables.h>
//#include <ZFEWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel) 
: myFlavour(myModel), myEW(myModel), myMO(myModel)
//, myZFitter(myModel) 
{
    //-----   Flavour observables   -----
    thobs["Dmd0"] = new DmBd0(myFlavour);
    thobs["Dmd1"] = new DmBd(myFlavour);
    thobs["Dms0"] = new DmBs0(myFlavour);
    thobs["Dms1"] = new DmBs(myFlavour);
    thobs["M12D"] = new M12D(myFlavour);
    thobs["ArgD"] = new ArgD(myFlavour);
    thobs["EpsilonK"] = new EpsilonK(myFlavour) ;
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
    if(myModel.ModelName().compare("SUSY") || myModel.ModelName().compare("GeneralSUSY")
            || myModel.ModelName().compare("MFV")) {
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
    }
    
    //-----  SM input parameters, etc.  -----
    thobs["AlsMz"] = new alsMz(myMO);
    thobs["dAle5Mz"] = new dAle5Mz(myMO);
    thobs["mHl"] = new mHl(myMO);
    thobs["Mz"] = new mZ(myMO);
    thobs["mtop"] = new mtpole(myMO);
    thobs["delRhoZ_nu"] = new delRhoZ_nu(myMO);
    thobs["delRhoZ_e"] = new delRhoZ_e(myMO);
    thobs["delRhoZ_u"] = new delRhoZ_u(myMO);
    thobs["delRhoZ_d"] = new delRhoZ_d(myMO);
    thobs["delRhoZ_b"] = new delRhoZ_b(myMO);
    
    //-----  NP input parameters, etc.  -----
    thobs["cHLp_NP"] = new cHLp_NP(myMO);
    thobs["cHQp_NP"] = new cHQp_NP(myMO);
    thobs["cHQ_NP"] = new cHQ_NP(myMO);
    thobs["cHL_NP"] = new cHL_NP(myMO);
    thobs["cHE_NP"] = new cHE_NP(myMO);
    thobs["cHU2_NP"] = new cHU2_NP(myMO);
    thobs["cHD3_NP"] = new cHD3_NP(myMO);
    thobs["cHQ1pPLUScHQ2p_NP"] = new cHQ1pPLUScHQ2p_NP(myMO);
    thobs["cHQ2pMINUScHQ2_NP"] = new cHQ2pMINUScHQ2_NP(myMO);
    thobs["cHQ3pPLUScHQ3_NP"] = new cHQ3pPLUScHQ3_NP(myMO);
    thobs["c_Ae_NP"] = new c_Ae_NP(myMO);
    thobs["c_GammaZ_uds_NP"] = new c_GammaZ_uds_NP(myMO);
    thobs["deltaGVb"] = new deltaGVb(myMO);
    thobs["deltaGAb"] = new deltaGAb(myMO);
    thobs["deltaGLb"] = new deltaGLb(myMO);
    thobs["deltaGRb"] = new deltaGRb(myMO);
    thobs["deltaRhoZb"] = new deltaRhoZb(myMO);
    thobs["deltaKappaZb"] = new deltaKappaZb(myMO);
    
    //-----   Higgs mass   -----
    thobs["Mh0"] = new Mh0(myEW);

    //-----  Z-pole observables (with EW and StandardModel)  -----
    thobs["Mw"] = new Mw(myEW);
    thobs["sin2thetaEff"] = new sin2thetaEff(myEW);
    thobs["GammaW"] = new GammaW(myEW);
    thobs["GammaZ"] = new GammaZ(myEW);
    thobs["Alepton"] = new Alepton(myEW);
    thobs["Acharm"] = new Acharm(myEW);
    thobs["Abottom"] = new Abottom(myEW);
    thobs["PtauPol"] = new PtauPol(myEW);
    thobs["AFBlepton"] = new AFBlepton(myEW);
    thobs["AFBcharm"] = new AFBcharm(myEW);
    thobs["AFBbottom"] = new AFBbottom(myEW);
    thobs["Rlepton"] = new Rlepton(myEW);
    thobs["Rcharm"] = new Rcharm(myEW);
    thobs["Rbottom"] = new Rbottom(myEW);
    thobs["sigmaHadron"] = new sigmaHadron(myEW);

    //-----  epsilon parameters by Altarelli et al.  -----
    thobs["epsilon1"] = new epsilon1(myEW);
    thobs["epsilon2"] = new epsilon2(myEW);
    thobs["epsilon3"] = new epsilon3(myEW);   
    thobs["epsilonb"] = new epsilonb(myEW);   
    
    //-----   Z-pole observables (with ZFitter)   -----
    //thobs["Mw"] = new ZFMw(myZFitter);
    //thobs["sin2thetaEff"] = new ZFsin2thetaEff(myZFitter);
    //thobs["GammaW"] = new ZFGammaW(myZFitter);
    //thobs["GammaZ"] = new ZFGammaZ(myZFitter);
    //thobs["Alepton"] = new ZFAlepton(myZFitter);
    //thobs["Acharm"] = new ZFAcharm(myZFitter);
    //thobs["Abottom"] = new ZFAbottom(myZFitter);
    //thobs["PtauPol"] = new ZFPtauPol(myZFitter);
    //thobs["AFBlepton"] = new ZFAFBlepton(myZFitter);
    //thobs["AFBcharm"] = new ZFAFBcharm(myZFitter);
    //thobs["AFBbottom"] = new ZFAFBbottom(myZFitter);    
    //thobs["Rlepton"] = new ZFRlepton(myZFitter);
    //thobs["Rcharm"] = new ZFRcharm(myZFitter);
    //thobs["Rbottom"] = new ZFRbottom(myZFitter);
    //thobs["sigmaHadron"] = new ZFsigmaHadron(myZFitter);

    //-----   For LEP-II observables   -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189., 
                               192., 196., 200., 202., 205., 207.};
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192., 
                                  196., 200., 202., 205., 207.};

    //-----  LEP-II observables (with EWSMTwoFermionsLEP2 class in StandardModel)  -----
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

    //-----  LEP-II observables (with ZFitter)  -----
    //ZFsigmaQuarksLEP2* myZFsigmaQuarks[12];
    //ZFsigmaMuLEP2* myZFsigmaMu[12];
    //ZFsigmaTauLEP2* myZFsigmaTau[12];
    //ZFAFBmuLEP2* myZFAFBmu[12];
    //ZFAFBtauLEP2* myZFAFBtau[12];
    //ZFAFBbottomLEP2* myZFAFBbottomLEP2[10];
    //ZFAFBcharmLEP2* myZFAFBcharmLEP2[10];
    //ZFRbottomLEP2* myZFRbottomLEP2[10];
    //ZFRcharmLEP2* myZFRcharmLEP2[10];
    //for (int i=0; i<12; i++) { 
    //    std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
    //    myZFsigmaQuarks[i] = new ZFsigmaQuarksLEP2(myZFitter, sqrt_s[i]);
    //    thobs["sigmaqLEP2_" + sqrt_s_str] = myZFsigmaQuarks[i];
    //    myZFsigmaMu[i] = new ZFsigmaMuLEP2(myZFitter, sqrt_s[i]);
    //    thobs["sigmamuLEP2_" + sqrt_s_str] = myZFsigmaMu[i];
    //    myZFsigmaTau[i] = new ZFsigmaTauLEP2(myZFitter, sqrt_s[i]);
    //    thobs["sigmatauLEP2_" + sqrt_s_str] = myZFsigmaTau[i];
    //    myZFAFBmu[i] = new ZFAFBmuLEP2(myZFitter, sqrt_s[i]);
    //    thobs["AFBmuLEP2_" + sqrt_s_str] = myZFAFBmu[i];
    //    myZFAFBtau[i] = new ZFAFBtauLEP2(myZFitter, sqrt_s[i]);
    //    thobs["AFBtauLEP2_" + sqrt_s_str] = myZFAFBtau[i];
    //}
    //for (int i=0; i<10; i++) { 
    //    std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
    //    myZFAFBbottomLEP2[i] = new ZFAFBbottomLEP2(myZFitter, sqrt_s_HF[i]);
    //    thobs["AFBbottomLEP2_" + sqrt_s_str] = myZFAFBbottomLEP2[i];
    //    myZFAFBcharmLEP2[i] = new ZFAFBcharmLEP2(myZFitter, sqrt_s_HF[i]);
    //    thobs["AFBcharmLEP2_" + sqrt_s_str] = myZFAFBcharmLEP2[i];
    //    myZFRbottomLEP2[i] = new ZFRbottomLEP2(myZFitter, sqrt_s_HF[i]);  
    //    thobs["RbottomLEP2_" + sqrt_s_str] = myZFRbottomLEP2[i];
    //    myZFRcharmLEP2[i] = new ZFRcharmLEP2(myZFitter, sqrt_s_HF[i]);
    //    thobs["RcharmLEP2_" + sqrt_s_str] = myZFRcharmLEP2[i];  
    //}    
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
