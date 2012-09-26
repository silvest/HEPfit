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

    //-----  LEP-II observables (with EWSMTwoFermionsLEP2 class in StandardModel)  -----
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189., 
                               192., 196., 200., 202., 205., 207.};
    LEP2sigmaHadron *myLEP2sigmaHadron[12];
    LEP2sigmaMu *myLEP2sigmaMu[12];
    LEP2sigmaTau *myLEP2sigmaTau[12];
    LEP2AFBmu *myLEP2AFBmu[12];
    LEP2AFBtau *myLEP2AFBtau[12];
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
    //thobs["sigmaqLEP2_130"] = new sigmaqLEP2(myEW,130.); 
    //thobs["sigmaqLEP2_136"] = new sigmaqLEP2(myEW,136.);
    //thobs["sigmaqLEP2_161"] = new sigmaqLEP2(myEW,161.);
    //thobs["sigmaqLEP2_172"] = new sigmaqLEP2(myEW,172.);
    //thobs["sigmaqLEP2_183"] = new sigmaqLEP2(myEW,183.);
    //thobs["sigmaqLEP2_189"] = new sigmaqLEP2(myEW,189.);
    //thobs["sigmaqLEP2_192"] = new sigmaqLEP2(myEW,192.);
    //thobs["sigmaqLEP2_196"] = new sigmaqLEP2(myEW,196.);
    //thobs["sigmaqLEP2_200"] = new sigmaqLEP2(myEW,200.);
    //thobs["sigmaqLEP2_202"] = new sigmaqLEP2(myEW,202.);
    //thobs["sigmaqLEP2_205"] = new sigmaqLEP2(myEW,205.);
    //thobs["sigmaqLEP2_207"] = new sigmaqLEP2(myEW,207.);
    
    //thobs["sigmamuLEP2_130"] = new sigmamuLEP2(myEW,130.);
    //thobs["sigmamuLEP2_136"] = new sigmamuLEP2(myEW,136.);
    //thobs["sigmamuLEP2_161"] = new sigmamuLEP2(myEW,161.);
    //thobs["sigmamuLEP2_172"] = new sigmamuLEP2(myEW,172.);
    //thobs["sigmamuLEP2_183"] = new sigmamuLEP2(myEW,183.);
    //thobs["sigmamuLEP2_189"] = new sigmamuLEP2(myEW,189.);
    //thobs["sigmamuLEP2_192"] = new sigmamuLEP2(myEW,192.);
    //thobs["sigmamuLEP2_196"] = new sigmamuLEP2(myEW,196.);
    //thobs["sigmamuLEP2_200"] = new sigmamuLEP2(myEW,200.);
    //thobs["sigmamuLEP2_202"] = new sigmamuLEP2(myEW,202.);
    //thobs["sigmamuLEP2_205"] = new sigmamuLEP2(myEW,205.);    
    //thobs["sigmamuLEP2_207"] = new sigmamuLEP2(myEW,207.);
    
    //thobs["sigmatauLEP2_130"] = new sigmatauLEP2(myEW,130.);
    //thobs["sigmatauLEP2_136"] = new sigmatauLEP2(myEW,136.);
    //thobs["sigmatauLEP2_161"] = new sigmatauLEP2(myEW,161.);
    //thobs["sigmatauLEP2_172"] = new sigmatauLEP2(myEW,172.);
    //thobs["sigmatauLEP2_183"] = new sigmatauLEP2(myEW,183.);
    //thobs["sigmatauLEP2_189"] = new sigmatauLEP2(myEW,189.);
    //thobs["sigmatauLEP2_192"] = new sigmatauLEP2(myEW,192.);
    //thobs["sigmatauLEP2_196"] = new sigmatauLEP2(myEW,196.);
    //thobs["sigmatauLEP2_200"] = new sigmatauLEP2(myEW,200.);
    //thobs["sigmatauLEP2_202"] = new sigmatauLEP2(myEW,202.);
    //thobs["sigmatauLEP2_205"] = new sigmatauLEP2(myEW,205.);    
    //thobs["sigmatauLEP2_207"] = new sigmatauLEP2(myEW,207.);

    //thobs["AFBmuLEP2_130"] = new AFBmuLEP2(myEW,130.);
    //thobs["AFBmuLEP2_136"] = new AFBmuLEP2(myEW,136.);
    //thobs["AFBmuLEP2_161"] = new AFBmuLEP2(myEW,161.);
    //thobs["AFBmuLEP2_172"] = new AFBmuLEP2(myEW,172.);
    //thobs["AFBmuLEP2_183"] = new AFBmuLEP2(myEW,183.);
    //thobs["AFBmuLEP2_189"] = new AFBmuLEP2(myEW,189.);
    //thobs["AFBmuLEP2_192"] = new AFBmuLEP2(myEW,192.);
    //thobs["AFBmuLEP2_196"] = new AFBmuLEP2(myEW,196.);
    //thobs["AFBmuLEP2_200"] = new AFBmuLEP2(myEW,200.);
    //thobs["AFBmuLEP2_202"] = new AFBmuLEP2(myEW,202.);
    //thobs["AFBmuLEP2_205"] = new AFBmuLEP2(myEW,205.);
    //thobs["AFBmuLEP2_207"] = new AFBmuLEP2(myEW,207.);
    
    //thobs["AFBtauLEP2_130"] = new AFBtauLEP2(myEW,130.);
    //thobs["AFBtauLEP2_136"] = new AFBtauLEP2(myEW,136.);
    //thobs["AFBtauLEP2_161"] = new AFBtauLEP2(myEW,161.);
    //thobs["AFBtauLEP2_172"] = new AFBtauLEP2(myEW,172.);
    //thobs["AFBtauLEP2_183"] = new AFBtauLEP2(myEW,183.);
    //thobs["AFBtauLEP2_189"] = new AFBtauLEP2(myEW,189.);
    //thobs["AFBtauLEP2_192"] = new AFBtauLEP2(myEW,192.);
    //thobs["AFBtauLEP2_196"] = new AFBtauLEP2(myEW,196.);
    //thobs["AFBtauLEP2_200"] = new AFBtauLEP2(myEW,200.);
    //thobs["AFBtauLEP2_202"] = new AFBtauLEP2(myEW,202.);
    //thobs["AFBtauLEP2_205"] = new AFBtauLEP2(myEW,205.);
    //thobs["AFBtauLEP2_207"] = new AFBtauLEP2(myEW,207.);
    
    //thobs["AFBbottomLEP2_133"] = new AFBbottomLEP2(myEW,133.);
    //thobs["AFBbottomLEP2_167"] = new AFBbottomLEP2(myEW,167.);
    //thobs["AFBbottomLEP2_183"] = new AFBbottomLEP2(myEW,183.);
    //thobs["AFBbottomLEP2_189"] = new AFBbottomLEP2(myEW,189.);
    //thobs["AFBbottomLEP2_192"] = new AFBbottomLEP2(myEW,192.);
    //thobs["AFBbottomLEP2_196"] = new AFBbottomLEP2(myEW,196.);
    //thobs["AFBbottomLEP2_200"] = new AFBbottomLEP2(myEW,200.);
    //thobs["AFBbottomLEP2_202"] = new AFBbottomLEP2(myEW,202.);
    //thobs["AFBbottomLEP2_205"] = new AFBbottomLEP2(myEW,205.);
    //thobs["AFBbottomLEP2_207"] = new AFBbottomLEP2(myEW,207.);
    
    //thobs["AFBcharmLEP2_133"] = new AFBcharmLEP2(myEW,133.);    
    //thobs["AFBcharmLEP2_167"] = new AFBcharmLEP2(myEW,167.);    
    //thobs["AFBcharmLEP2_183"] = new AFBcharmLEP2(myEW,183.);    
    //thobs["AFBcharmLEP2_189"] = new AFBcharmLEP2(myEW,189.);    
    
    //thobs["RbottomLEP2_133"] = new RbottomLEP2(myEW,133.);    
    //thobs["RbottomLEP2_167"] = new RbottomLEP2(myEW,167.);    
    //thobs["RbottomLEP2_183"] = new RbottomLEP2(myEW,183.);    
    //thobs["RbottomLEP2_189"] = new RbottomLEP2(myEW,189.);    
    //thobs["RbottomLEP2_192"] = new RbottomLEP2(myEW,192.);    
    //thobs["RbottomLEP2_196"] = new RbottomLEP2(myEW,196.);    
    //thobs["RbottomLEP2_200"] = new RbottomLEP2(myEW,200.);    
    //thobs["RbottomLEP2_202"] = new RbottomLEP2(myEW,202.);    
    //thobs["RbottomLEP2_205"] = new RbottomLEP2(myEW,205.);    
    //thobs["RbottomLEP2_207"] = new RbottomLEP2(myEW,207.);    
    
    //thobs["RcharmLEP2_183"] = new RcharmLEP2(myEW,183.); 
    //thobs["RcharmLEP2_189"] = new RcharmLEP2(myEW,189.); 
    //thobs["RcharmLEP2_196"] = new RcharmLEP2(myEW,196.); 
    //thobs["RcharmLEP2_200"] = new RcharmLEP2(myEW,200.); 
    //thobs["RcharmLEP2_205"] = new RcharmLEP2(myEW,205.); 
    //thobs["RcharmLEP2_207"] = new RcharmLEP2(myEW,207.); 

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
