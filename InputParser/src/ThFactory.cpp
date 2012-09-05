/* 
 * File:   ThFactory.cpp
 * Author: silvest
 * 
 * Created on April 19, 2011, 12:23 PM
 */

#include "ThFactory.h"
#include <FlavourObservables.h>
#include <EWObservables.h>
#include <ZFEWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel) : 
        myFlavour(myModel), myEW(myModel), myZFitter(myModel) {
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

//    thobs["Mw"] = new Mw(myEW);
//    thobs["sin2thetaEff"] = new sin2thetaEff(myEW);
//    thobs["GammaW"] = new GammaW(myEW);
//    thobs["GammaZ"] = new GammaZ(myEW);
//    thobs["Alepton"] = new Alepton(myEW);
//    thobs["Acharm"] = new Acharm(myEW);
//    thobs["Abottom"] = new Abottom(myEW);
//    thobs["PtauPol"] = new PtauPol(myEW);
//    thobs["AFBlepton"] = new AFBlepton(myEW);
//    thobs["AFBcharm"] = new AFBcharm(myEW);
//    thobs["AFBbottom"] = new AFBbottom(myEW);    
//    thobs["Rlepton"] = new Rlepton(myEW);
//    thobs["Rcharm"] = new Rcharm(myEW);
//    thobs["Rbottom"] = new Rbottom(myEW);
//    thobs["sigmaHadron"] = new sigmaHadron(myEW);
    
//    thobs["Mw"] = new ZFMw(myZFitter);
//    thobs["sin2thetaEff"] = new ZFsin2thetaEff(myZFitter);
//    thobs["GammaW"] = new ZFGammaW(myZFitter);
//    thobs["GammaZ"] = new ZFGammaZ(myZFitter);
//    thobs["Alepton"] = new ZFAlepton(myZFitter);
//    thobs["Acharm"] = new ZFAcharm(myZFitter);
//    thobs["Abottom"] = new ZFAbottom(myZFitter);
//    thobs["PtauPol"] = new ZFPtauPol(myZFitter);
//    thobs["AFBlepton"] = new ZFAFBlepton(myZFitter);
//    thobs["AFBcharm"] = new ZFAFBcharm(myZFitter);
//    thobs["AFBbottom"] = new ZFAFBbottom(myZFitter);    
//    thobs["Rlepton"] = new ZFRlepton(myZFitter);
//    thobs["Rcharm"] = new ZFRcharm(myZFitter);
//    thobs["Rbottom"] = new ZFRbottom(myZFitter);
//    thobs["sigmaHadron"] = new ZFsigmaHadron(myZFitter);    
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
