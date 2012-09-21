/* 
 * File:   ThFactory.cpp
 * Author: silvest
 * 
 * Created on April 19, 2011, 12:23 PM
 */

#include "ThFactory.h"
#include <FlavourObservables.h>
#include <EWObservables.h>
//#include <ZFEWObservables.h>

ThFactory::ThFactory(const StandardModel& myModel) :
myFlavour(myModel), myEW(myModel)
//, myZFitter(myModel) 
{
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
    
    thobs["Mh0"] = new Mh0(myEW);

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

    //    I need to change s for LEP2 observables    
    //    thobs["sigmaqLEP2"] = new sigmaqLEP2(myEW); 

    //    thobs["sigmamuLEP2_130"] = new sigmamuLEP2(myEW,130.);
    //    thobs["sigmamuLEP2_136"] = new sigmamuLEP2(myEW,136.);
    //    thobs["sigmamuLEP2_161"] = new sigmamuLEP2(myEW,161.);
    //    thobs["sigmamuLEP2_172"] = new sigmamuLEP2(myEW,172.);
    //    thobs["sigmamuLEP2_182"] = new sigmamuLEP2(myEW,182.);
    //    thobs["sigmamuLEP2_189"] = new sigmamuLEP2(myEW,189.);
    //    thobs["sigmamuLEP2_192"] = new sigmamuLEP2(myEW,192.);
    //    thobs["sigmamuLEP2_196"] = new sigmamuLEP2(myEW,196.);
    //    thobs["sigmamuLEP2_200"] = new sigmamuLEP2(myEW,200.);
    //    thobs["sigmamuLEP2_202"] = new sigmamuLEP2(myEW,202.);
    //    thobs["sigmamuLEP2_205"] = new sigmamuLEP2(myEW,205.);    
    //    thobs["sigmamuLEP2_207"] = new sigmamuLEP2(myEW,207.);

    //    thobs["sigmatauLEP2_130"] = new sigmatauLEP2(myEW,130.);
    //    thobs["sigmatauLEP2_136"] = new sigmatauLEP2(myEW,136.);
    //    thobs["sigmatauLEP2_161"] = new sigmatauLEP2(myEW,161.);
    //    thobs["sigmatauLEP2_172"] = new sigmatauLEP2(myEW,172.);
    //    thobs["sigmatauLEP2_182"] = new sigmatauLEP2(myEW,182.);
    //    thobs["sigmatauLEP2_189"] = new sigmatauLEP2(myEW,189.);
    //    thobs["sigmatauLEP2_192"] = new sigmatauLEP2(myEW,192.);
    //    thobs["sigmatauLEP2_196"] = new sigmatauLEP2(myEW,196.);
    //    thobs["sigmatauLEP2_200"] = new sigmatauLEP2(myEW,200.);
    //    thobs["sigmatauLEP2_202"] = new sigmatauLEP2(myEW,202.);
    //    thobs["sigmatauLEP2_205"] = new sigmatauLEP2(myEW,205.);    
    //    thobs["sigmatauLEP2_207"] = new sigmatauLEP2(myEW,207.);

    //    thobs["sigmatauLEP2"] = new sigmatauLEP2(myEW);
    //    thobs["AFBmuLEP2"] = new AFBmuLEP2(myEW);
    //    thobs["AFBtauLEP2"] = new AFBtauLEP2(myEW);
    //    thobs["AFBbottomLEP2"] = new AFBbottomLEP2(myEW);
    //    thobs["AFBcharmLEP2"] = new AFBcharmLEP2(myEW);    
    //    thobs["RbottomLEP2"] = new RbottomLEP2(myEW);    
    //    thobs["RcharmLEP2"] = new RcharmLEP2(myEW); 



    //    thobs["obliqueS"] = new obliqueS(myEW);
    //    thobs["obliqueT"] = new obliqueT(myEW);
    //    thobs["obliqueU"] = new obliqueU(myEW);
    //    thobs["obliqueW"] = new obliqueW(myEW);
    //    thobs["obliqueX"] = new obliqueX(myEW);
    //    thobs["obliqueY"] = new obliqueY(myEW);  
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
