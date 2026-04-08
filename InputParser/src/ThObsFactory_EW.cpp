/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "EWObservables.h"
#include "masses.h"
#include "MtMSbar.h"
#include "alpha_s.h"
using namespace ThObsConst;

void ThObsFactory::registerEWObservables()
{
    //-----  StandardModel observables  -----
    obsThFactory["MtMSbar"] = [](const StandardModel& SM) { return new MtMSbar(SM); };
    obsThFactory["alpha_s_LO"] = [=](const StandardModel& SM) { return new alpha_s(SM, LO); };
    obsThFactory["alpha_s_NLO"] = [=](const StandardModel& SM) { return new alpha_s(SM, NLO); };
    obsThFactory["alpha_s_FULLNLO"] = [=](const StandardModel& SM) { return new alpha_s(SM, FULLNLO); };
    obsThFactory["alpha_s_NNLO"] = [=](const StandardModel& SM) { return new alpha_s(SM, NNLO); };
    obsThFactory["alpha_s_FULLNNLO"] = [=](const StandardModel& SM) { return new alpha_s(SM, FULLNNLO); };
    obsThFactory["up_mass"] = [](const StandardModel& SM) { return new up_mass(SM); };
    obsThFactory["down_mass"] = [](const StandardModel& SM) { return new down_mass(SM); };
    obsThFactory["strange_mass"] = [](const StandardModel& SM) { return new strange_mass(SM); };
    obsThFactory["charm_mass"] = [](const StandardModel& SM) { return new charm_mass(SM); };
    obsThFactory["bottom_mass"] = [](const StandardModel& SM) { return new bottom_mass(SM); };
    obsThFactory["top_mass"] = [](const StandardModel& SM) { return new top_mass(SM); };
    obsThFactory["mtpole"] = [](const StandardModel& SM) { return new mtpole(SM); };
    obsThFactory["electron_mass"] = [](const StandardModel& SM) { return new electron_mass(SM); };
    obsThFactory["muon_mass"] = [](const StandardModel& SM) { return new muon_mass(SM); };
    obsThFactory["tau_mass"] = [](const StandardModel& SM) { return new tau_mass(SM); };
    //-----  Electroweak precision observables  -----
    obsThFactory["alphaMz"] = [](const StandardModel& SM) { return new AlphaEmMz(SM); };
    obsThFactory["Dalpha_5h_Mz"] = [](const StandardModel& SM) { return new DAlpha5hadMz(SM); };
    obsThFactory["Mw"] = [](const StandardModel& SM) { return new Mw(SM); };
    obsThFactory["GammaW"] = [](const StandardModel& SM) { return new GammaW(SM); };
    obsThFactory["BrWlepton"] = [](const StandardModel& SM) { return new BrWlepton(SM); };
    obsThFactory["BrWelectron"] = [](const StandardModel& SM) { return new BrWelectron(SM); };
    obsThFactory["BrWmuon"] = [](const StandardModel& SM) { return new BrWmuon(SM); };
    obsThFactory["BrWtau"] = [](const StandardModel& SM) { return new BrWtau(SM); };
    obsThFactory["BrWhadrons"] = [](const StandardModel& SM) { return new BrWhadrons(SM); };
    obsThFactory["RWc"] = [](const StandardModel& SM) { return new RWc(SM); };
    obsThFactory["RW_mu_e"] = [](const StandardModel& SM) { return new RWmue(SM); };
    obsThFactory["RW_tau_e"] = [](const StandardModel& SM) { return new RWtaue(SM); };
    obsThFactory["RW_tau_mu"] = [](const StandardModel& SM) { return new RWtaumu(SM); };
    obsThFactory["RZ_mu_e"] = [](const StandardModel& SM) { return new RZmue(SM); };
    obsThFactory["RZ_tau_e"] = [](const StandardModel& SM) { return new RZtaue(SM); };
    obsThFactory["RZ_tau_mu"] = [](const StandardModel& SM) { return new RZtaumu(SM); };
    obsThFactory["GammaZ"] = [](const StandardModel& SM) { return new GammaZ(SM); };
    obsThFactory["GammaZhad"] = [](const StandardModel& SM) { return new GammaZhad(SM); };
    obsThFactory["sigmaHadron"] = [](const StandardModel& SM) { return new sigmaHadron(SM); };
    obsThFactory["sin2thetaEff"] = [](const StandardModel& SM) { return new sin2thetaEff(SM); };
    obsThFactory["sin2thetaEffelectron"] = [](const StandardModel& SM) { return new sin2thetaEffel(SM); };
    obsThFactory["sin2thetaEffmuon"] = [](const StandardModel& SM) { return new sin2thetaEffmu(SM); };
    obsThFactory["PtauPol"] = [](const StandardModel& SM) { return new PtauPol(SM); };
    obsThFactory["Alepton"] = [](const StandardModel& SM) { return new Alepton(SM); };
    obsThFactory["Aelectron"] = [](const StandardModel& SM) { return new Aelectron(SM); };
    obsThFactory["Amuon"] = [](const StandardModel& SM) { return new Amuon(SM); };
    obsThFactory["Atau"] = [](const StandardModel& SM) { return new Atau(SM); };
    obsThFactory["Astrange"] = [](const StandardModel& SM) { return new Astrange(SM); };
    obsThFactory["Acharm"] = [](const StandardModel& SM) { return new Acharm(SM); };
    obsThFactory["Abottom"] = [](const StandardModel& SM) { return new Abottom(SM); };
    obsThFactory["AFBlepton"] = [](const StandardModel& SM) { return new AFBlepton(SM); };
    obsThFactory["AFBelectron"] = [](const StandardModel& SM) { return new AFBelectron(SM); };
    obsThFactory["AFBmuon"] = [](const StandardModel& SM) { return new AFBmuon(SM); };
    obsThFactory["AFBtau"] = [](const StandardModel& SM) { return new AFBtau(SM); };
    obsThFactory["AFBstrange"] = [](const StandardModel& SM) { return new AFBstrange(SM); };
    obsThFactory["AFBcharm"] = [](const StandardModel& SM) { return new AFBcharm(SM); };
    obsThFactory["AFBbottom"] = [](const StandardModel& SM) { return new AFBbottom(SM); };
    obsThFactory["Nneutrinos"] = [](const StandardModel& SM) { return new Nneutrinos(SM); };
    obsThFactory["Rlepton"] = [](const StandardModel& SM) { return new Rlepton(SM); };
    obsThFactory["Relectron"] = [](const StandardModel& SM) { return new Relectron(SM); };
    obsThFactory["Rmuon"] = [](const StandardModel& SM) { return new Rmuon(SM); };
    obsThFactory["Rtau"] = [](const StandardModel& SM) { return new Rtau(SM); };
    obsThFactory["Rneutrinos"] = [](const StandardModel& SM) { return new Rneutrinos(SM); };
    obsThFactory["Rinv"] = [](const StandardModel& SM) { return new Rinv(SM); };
    obsThFactory["Ruc"] = [](const StandardModel& SM) { return new Ruc(SM); };
    obsThFactory["Rstrange"] = [](const StandardModel& SM) { return new Rstrange(SM); };
    obsThFactory["Rcharm"] = [](const StandardModel& SM) { return new Rcharm(SM); };
    obsThFactory["Rbottom"] = [](const StandardModel& SM) { return new Rbottom(SM); };
    //
    obsThFactory["GammaZee"] = [](const StandardModel& SM) { return new GammaZee(SM); };
    obsThFactory["GammaZmumu"] = [](const StandardModel& SM) { return new GammaZmumu(SM); };
    obsThFactory["GammaZtautau"] = [](const StandardModel& SM) { return new GammaZtautau(SM); };
    obsThFactory["GammaZinv"] = [](const StandardModel& SM) { return new GammaZinv(SM); };
    obsThFactory["GammaZuu"] = [](const StandardModel& SM) { return new GammaZuu(SM); };
    obsThFactory["GammaZcc"] = [](const StandardModel& SM) { return new GammaZcc(SM); };
    obsThFactory["GammaZdd"] = [](const StandardModel& SM) { return new GammaZdd(SM); };
    obsThFactory["GammaZss"] = [](const StandardModel& SM) { return new GammaZss(SM); };
    obsThFactory["GammaZbb"] = [](const StandardModel& SM) { return new GammaZbb(SM); };
    //----- Higgs observables: Decay widths
    obsThFactory["GammaHtobb"] = [](const StandardModel& SM) { return new Htobb(SM); };
    obsThFactory["GammaHtocc"] = [](const StandardModel& SM) { return new Htocc(SM); };
    obsThFactory["GammaHtoss"] = [](const StandardModel& SM) { return new Htoss(SM); };
    obsThFactory["GammaHtotautau"] = [](const StandardModel& SM) { return new Htotautau(SM); };
    obsThFactory["GammaHtomumu"] = [](const StandardModel& SM) { return new Htomumu(SM); };
    obsThFactory["GammaHtoWW"] = [](const StandardModel& SM) { return new HtoWW(SM); };
    obsThFactory["GammaHtoZZ"] = [](const StandardModel& SM) { return new HtoZZ(SM); };
    obsThFactory["GammaHtogaga"] = [](const StandardModel& SM) { return new Htogaga(SM); };
    obsThFactory["GammaHtoZga"] = [](const StandardModel& SM) { return new HtoZga(SM); };
    obsThFactory["GammaHtogg"] = [](const StandardModel& SM) { return new Htogg(SM); };
    obsThFactory["GammaH"] = [](const StandardModel& SM) { return new Hwidth(SM); };
}
