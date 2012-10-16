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
#include <SusyObservables.h>
#include <EWObservables.h>
//#include <ZFEWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel) :
myFlavour(myModel), myEW(myModel), myMO(myModel)
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
    thobs["gamma"] = new CKMGamma(myFlavour);
    thobs["SJPsiK"] = new SJPsiK(myFlavour);
    if(myModel.ModelName().compare("MFV")) {
        thobs["Msq"] = new Msq(myMO);
    }
    
    //-----   Higgs mass   -----
    thobs["Mh0"] = new Mh0(myEW);

    //-----  Z-pole observables (with EW and StandardModel)  -----
    //bool flag_CHMN = true; // TEST: use EW_CHMN class
    bool flag_CHMN = false;
    //bool flag_BURGESS = true; // TEST: use the formulae in hep-ph/9411257 by C.P. Burgess
    bool flag_BURGESS = false;
    thobs["Mw"] = new Mw(myEW, flag_CHMN, flag_BURGESS);
    thobs["sin2thetaEff"] = new sin2thetaEff(myEW, flag_CHMN, flag_BURGESS);
    thobs["GammaW"] = new GammaW(myEW, flag_CHMN, flag_BURGESS);
    thobs["GammaZ"] = new GammaZ(myEW, flag_CHMN, flag_BURGESS);
    thobs["Alepton"] = new Alepton(myEW, flag_CHMN, flag_BURGESS);
    thobs["Acharm"] = new Acharm(myEW, flag_CHMN, flag_BURGESS);
    thobs["Abottom"] = new Abottom(myEW, flag_CHMN, flag_BURGESS);
    thobs["PtauPol"] = new PtauPol(myEW, flag_CHMN, flag_BURGESS);
    thobs["AFBlepton"] = new AFBlepton(myEW, flag_CHMN, flag_BURGESS);
    thobs["AFBcharm"] = new AFBcharm(myEW, flag_CHMN, flag_BURGESS);
    thobs["AFBbottom"] = new AFBbottom(myEW, flag_CHMN, flag_BURGESS);
    thobs["Rlepton"] = new Rlepton(myEW, flag_CHMN, flag_BURGESS);
    thobs["Rcharm"] = new Rcharm(myEW, flag_CHMN, flag_BURGESS);
    thobs["Rbottom"] = new Rbottom(myEW, flag_CHMN, flag_BURGESS);
    thobs["sigmaHadron"] = new sigmaHadron(myEW, flag_CHMN, flag_BURGESS);

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
