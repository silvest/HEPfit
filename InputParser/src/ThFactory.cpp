/* 
 * File:   ThFactory.cpp
 * Author: silvest
 * 
 * Created on April 19, 2011, 12:23 PM
 */

#include "ThFactory.h"
#include <Flavour.h>
#include <EW.h>
#include <EWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel, StandardModelMatching& myModelMatching) : 
myFlavour(myModel, myModelMatching), myEW(myModel) {
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
    thobs["Mw"] = new Mw(myEW);
    thobs["GammaW"] = new GammaW(myEW);
    thobs["GammaZ"] = new GammaZ(myEW);
    thobs["AFBlepton"] = new AFBlepton(myEW);
    thobs["AFBcharm"] = new AFBcharm(myEW);
    thobs["AFBbottom"] = new AFBbottom(myEW);
    thobs["PtauPol"] = new PtauPol(myEW);
    thobs["Alepton"] = new Alepton(myEW);
    thobs["Acharm"] = new Acharm(myEW);
    thobs["Abottom"] = new Abottom(myEW);
    thobs["Rlepton"] = new Rlepton(myEW);
    thobs["Rcharm"] = new Rcharm(myEW);
    thobs["Rbottom"] = new Rbottom(myEW);
    thobs["sigmaHadron"] = new sigmaHadron(myEW);
    thobs["sin2thetaEff"] = new sin2thetaEff(myEW);
    thobs["obliqueEpsilon1"] = new obliqueEpsilon1(myEW);
    thobs["obliqueEpsilon2"] = new obliqueEpsilon2(myEW);
    thobs["obliqueEpsilon3"] = new obliqueEpsilon3(myEW);
    thobs["obliqueS"] = new obliqueS(myEW);
    thobs["obliqueT"] = new obliqueT(myEW);
    thobs["obliqueU"] = new obliqueU(myEW);
    thobs["obliqueShat"] = new obliqueShat(myEW);
    thobs["ObliqueThat"] = new obliqueThat(myEW);
    thobs["obliqueW"] = new obliqueW(myEW);
    thobs["obliqueY"] = new obliqueY(myEW);
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
