/* 
 * File:   ThFactory.cpp
 * Author: silvest
 * 
 * Created on April 19, 2011, 12:23 PM
 */

#include <vector>
#include <boost/lexical_cast.hpp>
#include "ThFactory.h"
#include <FlavourObservables.h>
#include <EWObservables.h>
//#include <ZFEWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel) :
myFlavour(myModel), myEW(myModel)
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
    thobs["DmK"] = new DmK(myFlavour);
    thobs["Vud"] = new Vud(myFlavour);
    thobs["Vus"] = new Vus(myFlavour);
    thobs["Vub"] = new Vub(myFlavour);
    thobs["Vcb"] = new Vcb(myFlavour);
    thobs["alpha"] = new Alpha(myFlavour);
    thobs["alpha_2a"] = new Alpha_2a(myFlavour);
    thobs["gamma"] = new Gamma(myFlavour);
    thobs["SJPsiK"] = new SJPsiK(myFlavour);
    
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

    //-----  LEP-II observables (with EWSMOneLoopLEP2 class in StandardModel)  -----    
    //sigmaqLEP2* mySigmaqLEP2[12];
    //sigmamuLEP2* mySigmamuLEP2[12];
    //sigmatauLEP2* mySigmatauLEP2[12];
    //AFBmuLEP2* myAFBmuLEP2[12];
    //AFBtauLEP2* myAFBtauLEP2[12];
    //AFBbottomLEP2* myAFBbottomLEP2[10];
    //AFBcharmLEP2* myAFBcharmLEP2[10];
    //RbottomLEP2* myRbottomLEP2[10];
    //RcharmLEP2* myRcharmLEP2[10];
    //for (int i=0; i<12; i++) { 
    //    std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
    //    mySigmaqLEP2[i] = new sigmaqLEP2(myEW, sqrt_s[i]);
    //    thobs["sigmaqLEP2_" + sqrt_s_str] = mySigmaqLEP2[i];
    //    mySigmamuLEP2[i] = new sigmamuLEP2(myEW, sqrt_s[i]);
    //    thobs["sigmamuLEP2_" + sqrt_s_str] = mySigmamuLEP2[i];
    //    mySigmatauLEP2[i] = new sigmatauLEP2(myEW, sqrt_s[i]);
    //    thobs["sigmatauLEP2_" + sqrt_s_str] = mySigmatauLEP2[i];
    //    myAFBmuLEP2[i] = new AFBmuLEP2(myEW, sqrt_s[i]);
    //    thobs["AFBmuLEP2_" + sqrt_s_str] = myAFBmuLEP2[i];
    //    myAFBtauLEP2[i] = new AFBtauLEP2(myEW, sqrt_s[i]);
    //    thobs["AFBtauLEP2_" + sqrt_s_str] = myAFBtauLEP2[i];
    //    myAFBbottomLEP2[i] = new AFBbottomLEP2(myEW, sqrt_s[i]);
    //    thobs["AFBbottomLEP2_" + sqrt_s_str] = myAFBbottomLEP2[i];
    //    myAFBcharmLEP2[i] = new AFBcharmLEP2(myEW, sqrt_s[i]);
    //    thobs["AFBcharmLEP2_" + sqrt_s_str] = myAFBcharmLEP2[i];
    //    myRbottomLEP2[i] = new RbottomLEP2(myEW, sqrt_s[i]);  
    //    thobs["RbottomLEP2_" + sqrt_s_str] = myRbottomLEP2[i];
    //    myRcharmLEP2[i] = new RcharmLEP2(myEW, sqrt_s[i]);
    //    thobs["RcharmLEP2_" + sqrt_s_str] = myRcharmLEP2[i];        
    //}

}

ThFactory::~ThFactory() {
    for (std::map<std::string, ThObservable *>::iterator it = thobs.begin();
            it != thobs.end(); it++)
        if (it->second != NULL)
            delete it->second;
}

ThObservable * ThFactory::getThMethod(const std::string& name) {
    if (thobs.find(name) == thobs.end()) {
        std::cout << "wrong observable " << name << " in ThFactory" << std::endl;
        exit(EXIT_FAILURE);
    }
    return (thobs[name]);
}
