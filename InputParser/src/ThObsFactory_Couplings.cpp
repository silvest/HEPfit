/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "NP_couplings.h"
#include "EWObservables.h"
#include "OtherThObservables.h"
#include <boost/lexical_cast.hpp>
using namespace ThObsConst;

void ThObsFactory::registerCouplingObservables()
{
    //-----  Observables for particle couplings -----
    //-----  EM constant correction  ----------
    obsThFactory["deltae"] = [=](const StandardModel& SM) { return new deltae(SM, muEW); };
    //-----  Triple gauge coupling observables (scale independent definition -> muEW)  -----
    obsThFactory["deltag1Z"] = [=](const StandardModel& SM) { return new deltag1Z(SM, muEW); };
    obsThFactory["deltag1gamma"] = [=](const StandardModel& SM) { return new deltag1gamma(SM, muEW); };
    obsThFactory["deltaKgamma"] = [=](const StandardModel& SM) { return new deltaKgamma(SM, muEW); };
    obsThFactory["lambdaZ"] = [=](const StandardModel& SM) { return new lambdaZ(SM, muEW); };
    //-----  Zff couplings observables: relative corrections (scale independent definition -> muEW)  ----------
    //-----  Z couplings with neutrinos ---------
    obsThFactory["deltagZveveL"] = [=](const StandardModel& SM) { return new deltagZveveL(SM, muEW); };
    obsThFactory["deltagZvmuvmuL"] = [=](const StandardModel& SM) { return new deltagZvmuvmuL(SM, muEW); };
    obsThFactory["deltagZvtavtaL"] = [=](const StandardModel& SM) { return new deltagZvtavtaL(SM, muEW); };
    //-----  Z couplings with leptons ---------
    obsThFactory["deltagZeeL"] = [=](const StandardModel& SM) { return new deltagZeeL(SM, muEW); };
    obsThFactory["deltagZeeR"] = [=](const StandardModel& SM) { return new deltagZeeR(SM, muEW); };
    obsThFactory["deltagZmumuL"] = [=](const StandardModel& SM) { return new deltagZmumuL(SM, muEW); };
    obsThFactory["deltagZmumuR"] = [=](const StandardModel& SM) { return new deltagZmumuR(SM, muEW); };
    obsThFactory["deltagZtataL"] = [=](const StandardModel& SM) { return new deltagZtataL(SM, muEW); };
    obsThFactory["deltagZtataR"] = [=](const StandardModel& SM) { return new deltagZtataR(SM, muEW); };
    //-----  Z couplings with up sector quarks ---------
    obsThFactory["deltagZuuL"] = [=](const StandardModel& SM) { return new deltagZuuL(SM, muEW); };
    obsThFactory["deltagZuuR"] = [=](const StandardModel& SM) { return new deltagZuuR(SM, muEW); };
    obsThFactory["deltagZuuV"] = [=](const StandardModel& SM) { return new deltagZuuV(SM, muEW); };
    obsThFactory["deltagZuuA"] = [=](const StandardModel& SM) { return new deltagZuuA(SM, muEW); };
    obsThFactory["deltagZccL"] = [=](const StandardModel& SM) { return new deltagZccL(SM, muEW); };
    obsThFactory["deltagZccR"] = [=](const StandardModel& SM) { return new deltagZccR(SM, muEW); };
    obsThFactory["deltagZttL"] = [=](const StandardModel& SM) { return new deltagZttL(SM, muEW); };
    obsThFactory["deltagZttR"] = [=](const StandardModel& SM) { return new deltagZttR(SM, muEW); };
    obsThFactory["deltagZttV"] = [=](const StandardModel& SM) { return new deltagZttV(SM, muEW); };
    obsThFactory["deltagZttA"] = [=](const StandardModel& SM) { return new deltagZttA(SM, muEW); };
    //-----  Z couplings with down sector quarks ---------
    obsThFactory["deltagZddL"] = [=](const StandardModel& SM) { return new deltagZddL(SM, muEW); };
    obsThFactory["deltagZddR"] = [=](const StandardModel& SM) { return new deltagZddR(SM, muEW); };
    obsThFactory["deltagZddV"] = [=](const StandardModel& SM) { return new deltagZddV(SM, muEW); };
    obsThFactory["deltagZddA"] = [=](const StandardModel& SM) { return new deltagZddA(SM, muEW); };
    obsThFactory["deltagZssL"] = [=](const StandardModel& SM) { return new deltagZssL(SM, muEW); };
    obsThFactory["deltagZssR"] = [=](const StandardModel& SM) { return new deltagZssR(SM, muEW); };
    obsThFactory["deltagZbbL"] = [=](const StandardModel& SM) { return new deltagZbbL(SM, muEW); };
    obsThFactory["deltagZbbR"] = [=](const StandardModel& SM) { return new deltagZbbR(SM, muEW); };
    //-----  Zff couplings observables: absolute corrections (scale independent definition -> muEW)  ----------
    //-----  Z couplings with leptons ---------
    obsThFactory["delgZeL"] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::ELECTRON, muEW); };
    obsThFactory["delgZeR"] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::ELECTRON, muEW); };
    obsThFactory["delgZmuL"] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::MU, muEW); };
    obsThFactory["delgZmuR"] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::MU, muEW); };
    obsThFactory["delgZtaL"] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::TAU, muEW); };
    obsThFactory["delgZtaR"] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::TAU, muEW); };
    //-----  Z couplings with up sector quarks ---------
    obsThFactory["delgZuL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::UP, muEW); };
    obsThFactory["delgZuR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::UP, muEW); };
    obsThFactory["delgZcL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::CHARM, muEW); };
    obsThFactory["delgZcR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::CHARM, muEW); };
    obsThFactory["delgZtL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::TOP, muEW); };
    obsThFactory["delgZtR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::TOP, muEW); };
    //-----  Z couplings with down sector quarks ---------
    obsThFactory["delgZdL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::DOWN, muEW); };
    obsThFactory["delgZdR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::DOWN, muEW); };
    obsThFactory["delgZsL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::STRANGE, muEW); };
    obsThFactory["delgZsR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::STRANGE, muEW); };
    obsThFactory["delgZbL"] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::BOTTOM, muEW); };
    obsThFactory["delgZbR"] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::BOTTOM, muEW); };
    //-----  Wff couplings observables (scale independent definition -> muEW) ----------
    obsThFactory["deltaUWeve"] = [=](const StandardModel& SM) { return new deltaUWeve(SM, muEW); };
    obsThFactory["deltaUWmuvmu"] = [=](const StandardModel& SM) { return new deltaUWmuvmu(SM, muEW); };
    obsThFactory["deltaUWtavta"] = [=](const StandardModel& SM) { return new deltaUWtavta(SM, muEW); };
    obsThFactory["deltaVudL"] = [=](const StandardModel& SM) { return new deltaVudL(SM, muEW); };
    obsThFactory["deltaVudR"] = [=](const StandardModel& SM) { return new deltaVudR(SM, muEW); };
    obsThFactory["deltaVcsL"] = [=](const StandardModel& SM) { return new deltaVcsL(SM, muEW); };
    obsThFactory["deltaVcsR"] = [=](const StandardModel& SM) { return new deltaVcsR(SM, muEW); };
    obsThFactory["deltaVtbL"] = [=](const StandardModel& SM) { return new deltaVtbL(SM, muEW); };
    obsThFactory["deltaVtbR"] = [=](const StandardModel& SM) { return new deltaVtbR(SM, muEW); };
    //
    // Energy dependent definitions of the above
    for (int i = 0; i < 16; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_see[i]);
        
    //-----  EM constant correction  ----------
        obsThFactory["deltae_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltae(SM, sqrt_s_eeff[i]); };
    //-----  Triple gauge coupling observables (scale dependent definition)  -----
        obsThFactory["deltag1Z_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltag1Z(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltag1gamma_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltag1gamma(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaKgamma_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaKgamma(SM, sqrt_s_eeff[i]); };
        obsThFactory["lambdaZ_" + sqrt_s_str] = [=](const StandardModel& SM) { return new lambdaZ(SM, sqrt_s_eeff[i]); };
    //-----  Zff couplings observables: relative corrections (scale dependent definition)  ----------
    //-----  Z couplings with neutrinos ---------
        obsThFactory["deltagZveveL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZveveL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZvmuvmuL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZvmuvmuL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZvtavtaL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZvtavtaL(SM, sqrt_s_eeff[i]); };
    //-----  Z couplings with leptons ---------
        obsThFactory["deltagZeeL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZeeL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZeeR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZeeR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZmumuL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZmumuL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZmumuR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZmumuR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZtataL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZtataL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZtataR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZtataR(SM, sqrt_s_eeff[i]); };
    //-----  Z couplings with up sector quarks ---------
        obsThFactory["deltagZuuL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZuuL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZuuR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZuuR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZuuV_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZuuV(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZuuA_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZuuA(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZccL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZccL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZccR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZccR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZttL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZttL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZttR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZttR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZttV_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZttV(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZttA_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZttA(SM, sqrt_s_eeff[i]); };
    //-----  Z couplings with down sector quarks ---------
        obsThFactory["deltagZddL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZddL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZddR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZddR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZddV_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZddV(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZddA_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZddA(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZssL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZssL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZssR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZssR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZbbL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZbbL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltagZbbR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltagZbbR(SM, sqrt_s_eeff[i]); };
    //-----  Zff couplings observables: absolute corrections (scale dependent definition)  ----------
    //-----  Z couplings with leptons ---------
        obsThFactory["delgZeL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::ELECTRON, sqrt_s_eeff[i]); };
        obsThFactory["delgZeR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::ELECTRON, sqrt_s_eeff[i]); };
        obsThFactory["delgZmuL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::MU, sqrt_s_eeff[i]); };
        obsThFactory["delgZmuR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::MU, sqrt_s_eeff[i]); };
        obsThFactory["delgZtaL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlL(SM, StandardModel::TAU, sqrt_s_eeff[i]); };
        obsThFactory["delgZtaR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZlR(SM, StandardModel::TAU, sqrt_s_eeff[i]); };
    //-----  Z couplings with up sector quarks ---------
        obsThFactory["delgZuL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::UP, sqrt_s_eeff[i]); };
        obsThFactory["delgZuR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::UP, sqrt_s_eeff[i]); };
        obsThFactory["delgZcL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::CHARM, sqrt_s_eeff[i]); };
        obsThFactory["delgZcR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::CHARM, sqrt_s_eeff[i]); };
        obsThFactory["delgZtL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::TOP, sqrt_s_eeff[i]); };
        obsThFactory["delgZtR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::TOP, sqrt_s_eeff[i]); };
    //-----  Z couplings with down sector quarks ---------
        obsThFactory["delgZdL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::DOWN, sqrt_s_eeff[i]); };
        obsThFactory["delgZdR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::DOWN, sqrt_s_eeff[i]); };
        obsThFactory["delgZsL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::STRANGE, sqrt_s_eeff[i]); };
        obsThFactory["delgZsR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::STRANGE, sqrt_s_eeff[i]); };
        obsThFactory["delgZbL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqL(SM, StandardModel::BOTTOM, sqrt_s_eeff[i]); };
        obsThFactory["delgZbR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new delgZqR(SM, StandardModel::BOTTOM, sqrt_s_eeff[i]); };
    //-----  Wff couplings observables (scale dependent) ----------
        obsThFactory["deltaUWeve_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaUWeve(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaUWmuvmu_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaUWmuvmu(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaUWtavta_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaUWtavta(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVudL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVudL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVudR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVudR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVcsL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVcsL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVcsR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVcsR(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVtbL_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVtbL(SM, sqrt_s_eeff[i]); };
        obsThFactory["deltaVtbR_" + sqrt_s_str] = [=](const StandardModel& SM) { return new deltaVtbR(SM, sqrt_s_eeff[i]); };
    }
    //
    //-----  Zff EFFECTIVE couplings observables: relative corrections (derived from Af and Gamma(Z->ff)  ----------
    //-----  Z couplings with neutrinos ---------
    obsThFactory["deltagZveveLeff"] = [](const StandardModel& SM) { return new deltagEffZveveL(SM); };
    obsThFactory["deltagZvmuvmuLeff"] = [](const StandardModel& SM) { return new deltagEffZvmuvmuL(SM); };
    obsThFactory["deltagZvtavtaLeff"] = [](const StandardModel& SM) { return new deltagEffZvtavtaL(SM); };
    //-----  Z couplings with leptons ---------
    obsThFactory["deltagZeeLeff"] = [](const StandardModel& SM) { return new deltagEffZeeL(SM); };
    obsThFactory["deltagZeeReff"] = [](const StandardModel& SM) { return new deltagEffZeeR(SM); };
    obsThFactory["deltagZmumuLeff"] = [](const StandardModel& SM) { return new deltagEffZmumuL(SM); };
    obsThFactory["deltagZmumuReff"] = [](const StandardModel& SM) { return new deltagEffZmumuR(SM); };
    obsThFactory["deltagZtataLeff"] = [](const StandardModel& SM) { return new deltagEffZtataL(SM); };
    obsThFactory["deltagZtataReff"] = [](const StandardModel& SM) { return new deltagEffZtataR(SM); };
    //-----  Z couplings with up sector quarks ---------
    obsThFactory["deltagZccLeff"] = [](const StandardModel& SM) { return new deltagEffZccL(SM); };
    obsThFactory["deltagZccReff"] = [](const StandardModel& SM) { return new deltagEffZccR(SM); };
    //-----  Z couplings with down sector quarks ---------
    obsThFactory["deltagZssLeff"] = [](const StandardModel& SM) { return new deltagEffZssL(SM); };
    obsThFactory["deltagZssReff"] = [](const StandardModel& SM) { return new deltagEffZssR(SM); };
    obsThFactory["deltagZbbLeff"] = [](const StandardModel& SM) { return new deltagEffZbbL(SM); };
    obsThFactory["deltagZbbReff"] = [](const StandardModel& SM) { return new deltagEffZbbR(SM); };
    //
    //-----  W mass correction  ----------
    obsThFactory["deltaMW"] = [](const StandardModel& SM) { return new deltaMW(SM); };
    //-----  Hff couplings observables  ----------
    obsThFactory["gHmumueff"] = [](const StandardModel& SM) { return new gHmumueff(SM); };
    obsThFactory["gHtataeff"] = [](const StandardModel& SM) { return new gHtataeff(SM); };
    obsThFactory["gHcceff"] = [](const StandardModel& SM) { return new gHcceff(SM); };
    obsThFactory["gHbbeff"] = [](const StandardModel& SM) { return new gHbbeff(SM); };
    obsThFactory["deltagHee"] = [=](const StandardModel& SM) { return new deltagHee(SM, muMH); };//boost::factory<deltagHee*>();
    obsThFactory["deltagHmumu"] = [=](const StandardModel& SM) { return new deltagHmumu(SM, muMH); };//boost::factory<deltagHmumu*>();
    obsThFactory["deltagHtata"] = [=](const StandardModel& SM) { return new deltagHtata(SM, muMH); };//boost::factory<deltagHtata*>();
    obsThFactory["deltagHuu"] = [=](const StandardModel& SM) { return new deltagHuu(SM, muMH); };//boost::factory<deltagHuu*>();
    obsThFactory["deltagHcc"] = [=](const StandardModel& SM) { return new deltagHcc(SM, muMH); };//boost::factory<deltagHcc*>();
    obsThFactory["deltagHtt"] = [=](const StandardModel& SM) { return new deltagHtt(SM, muMH); };//boost::factory<deltagHtt*>();
    obsThFactory["deltagHdd"] = [=](const StandardModel& SM) { return new deltagHdd(SM, muMH); };//boost::factory<deltagHdd*>();
    obsThFactory["deltagHss"] = [=](const StandardModel& SM) { return new deltagHss(SM, muMH); };//boost::factory<deltagHss*>();
    obsThFactory["deltagHbb"] = [=](const StandardModel& SM) { return new deltagHbb(SM, muMH); };//boost::factory<deltagHbb*>();
    //-----  HGG couplings observables  ----------
    obsThFactory["gHGGeff"] = [](const StandardModel& SM) { return new gHGGeff(SM); };
    obsThFactory["deltagHGG"] = [=](const StandardModel& SM) { return new deltagHGG(SM, muMH); };//boost::factory<deltagHGG*>();
    //-----  HZZ couplings observables  ----------
    obsThFactory["gHZZeff"] = [](const StandardModel& SM) { return new gHZZeff(SM); };
    obsThFactory["gHZZ4feff"] = [](const StandardModel& SM) { return new gHZZ4feff(SM); };
    obsThFactory["deltagHZZ"] = [=](const StandardModel& SM) { return new deltagHZZ(SM, muMH); };//boost::factory<deltagHZZ*>();
    obsThFactory["gHZZ1"] = [=](const StandardModel& SM) { return new gHZZ1(SM, muMH); };//boost::factory<gHZZ1*>();
    obsThFactory["gHZZ2"] = [=](const StandardModel& SM) { return new gHZZ2(SM, muMH); };//boost::factory<gHZZ2*>();
    //-----  HAA couplings observables  ----------
    obsThFactory["gHAAeff"] = [](const StandardModel& SM) { return new gHAAeff(SM); };
    obsThFactory["deltagHAA"] = [=](const StandardModel& SM) { return new deltagHAA(SM, muMH); };//boost::factory<deltagHAA*>();
    //-----  HZA couplings observables  ----------
    obsThFactory["gHZAeff"] = [](const StandardModel& SM) { return new gHZAeff(SM); };
    obsThFactory["deltagHZA"] = [=](const StandardModel& SM) { return new deltagHZA(SM, muMH); };//boost::factory<deltagHZA*>();
    obsThFactory["gHZA2"] = [=](const StandardModel& SM) { return new gHZA2(SM, muMH); };//boost::factory<gHZA2*>();
    //-----  HWW couplings observables  ----------
    obsThFactory["gHWWeff"] = [](const StandardModel& SM) { return new gHWWeff(SM); };
    obsThFactory["gHWW4feff"] = [](const StandardModel& SM) { return new gHWW4feff(SM); };
    obsThFactory["deltagHWW"] = [=](const StandardModel& SM) { return new deltagHWW(SM, muMH); };//boost::factory<deltagHWW*>();
    obsThFactory["gHWW1"] = [=](const StandardModel& SM) { return new gHWW1(SM, muMH); };//boost::factory<gHWW1*>();
    obsThFactory["gHWW2"] = [=](const StandardModel& SM) { return new gHWW2(SM, muMH); };//boost::factory<gHWW2*>();
    //-----  HHH couplings observables  ----------
    obsThFactory["deltalHHH"] = [=](const StandardModel& SM) { return new deltalHHH(SM, muMH); };//boost::factory<deltalHHH*>();
    obsThFactory["deltalHHH246"] = [=](const StandardModel& SM) { return new deltalHHH(SM, 240.0); };//boost::factory<deltalHHH*>();
    obsThFactory["deltalHHH1000"] = [=](const StandardModel& SM) { return new deltalHHH(SM, 1000.0); };//boost::factory<deltalHHH*>();
    //-----  Other Higgs couplings observables  ----------
    obsThFactory["gHWZeff_Ratio"] = [](const StandardModel& SM) { return new gHWZeff(SM); };
    obsThFactory["gHbWeff_Ratio"] = [](const StandardModel& SM) { return new gHbWeff(SM); };
    obsThFactory["gHtaWeff_Ratio"] = [](const StandardModel& SM) { return new gHtaWeff(SM); };
    //-----  VVV couplings observables  ----------
    obsThFactory["deltag1ZEff"] = [](const StandardModel& SM) { return new deltag1ZEff(SM); };
    obsThFactory["deltaKgammaEff"] = [](const StandardModel& SM) { return new deltaKgammaEff(SM); };
    //-----  Basic interactions of the so-called Higgs basis  ----------
    obsThFactory["deltayt_HB"] = [=](const StandardModel& SM) { return new deltaytHB(SM, muMH); };////boost::factory<deltaytHB*>();
    obsThFactory["deltayb_HB"] = [=](const StandardModel& SM) { return new deltaybHB(SM, muMH); };////boost::factory<deltaybHB*>();
    obsThFactory["deltaytau_HB"] = [=](const StandardModel& SM) { return new deltaytauHB(SM, muMH); };////boost::factory<deltaytauHB*>();
    obsThFactory["deltayc_HB"] = [=](const StandardModel& SM) { return new deltaycHB(SM, muMH); };////boost::factory<deltaycHB*>();
    obsThFactory["deltaymu_HB"] = [=](const StandardModel& SM) { return new deltaymuHB(SM, muMH); };////boost::factory<deltaymuHB*>();
    obsThFactory["deltacZ_HB"] = [=](const StandardModel& SM) { return new deltacZHB(SM, muMH); };////boost::factory<deltacZHB*>();
    obsThFactory["cZBox_HB"] = [=](const StandardModel& SM) { return new cZBoxHB(SM, muMH); };////boost::factory<cZBoxHB*>();
    obsThFactory["cZZ_HB"] = [=](const StandardModel& SM) { return new cZZHB(SM, muMH); };////boost::factory<cZZHB*>();
    obsThFactory["cZga_HB"] = [=](const StandardModel& SM) { return new cZgaHB(SM, muMH); };////boost::factory<cZgaHB*>();
    obsThFactory["cgaga_HB"] = [=](const StandardModel& SM) { return new cgagaHB(SM, muMH); };////boost::factory<cgagaHB*>();
    obsThFactory["cgg_HB"] = [=](const StandardModel& SM) { return new cggHB(SM, muMH); };////boost::factory<cggHB*>();
    obsThFactory["cggEff_HB"] = [=](const StandardModel& SM) { return new cggEffHB(SM, muMH); };////boost::factory<cggEffHB*>();
    obsThFactory["lambz_HB"] = [=](const StandardModel& SM) { return new lambzHB(SM, muMH); };////boost::factory<lambzHB*>();
    //-----  Other useful observables to work with new physics  ----------
    //-----  Oblique Parameters ---------
    obsThFactory["oblSpar"] = [](const StandardModel& SM) { return new oblS(SM); };
    obsThFactory["oblTpar"] = [](const StandardModel& SM) { return new oblT(SM); };
    obsThFactory["oblWpar"] = [](const StandardModel& SM) { return new oblW(SM); };
    obsThFactory["oblYpar"] = [](const StandardModel& SM) { return new oblY(SM); };

    //-----  (Relative) Deviations of SM inputs with respect to reference value  ---------
    obsThFactory["deltaalphaMz"] = [](const StandardModel& SM) { return new dalphaMzRef(SM); };
    obsThFactory["deltaalphaSMz"] = [](const StandardModel& SM) { return new dalphaSMzRef(SM); };
    obsThFactory["deltaMz"] = [](const StandardModel& SM) { return new dMzRef(SM); };
    obsThFactory["deltaMh"] = [](const StandardModel& SM) { return new dMHRef(SM); };
    obsThFactory["deltamt"] = [](const StandardModel& SM) { return new dmtRef(SM); };

    //-----  Combinations of Warsaw basis coefficients constrained by EWPO  ----------
    obsThFactory["CEWHL1_11"] = [=](const StandardModel& SM) { return new CEWHL111(SM, muEW); };
    obsThFactory["CEWHL1_22"] = [=](const StandardModel& SM) { return new CEWHL122(SM, muEW); };
    obsThFactory["CEWHL1_33"] = [=](const StandardModel& SM) { return new CEWHL133(SM, muEW); };
    obsThFactory["CEWHL3_11"] = [=](const StandardModel& SM) { return new CEWHL311(SM, muEW); };
    obsThFactory["CEWHL3_22"] = [=](const StandardModel& SM) { return new CEWHL322(SM, muEW); };
    obsThFactory["CEWHL3_33"] = [=](const StandardModel& SM) { return new CEWHL333(SM, muEW); };
    obsThFactory["CEWHQ1_11"] = [=](const StandardModel& SM) { return new CEWHQ111(SM, muEW); };
    obsThFactory["CEWHQ1_22"] = [=](const StandardModel& SM) { return new CEWHQ122(SM, muEW); };
    obsThFactory["CEWHQ1_33"] = [=](const StandardModel& SM) { return new CEWHQ133(SM, muEW); };
    obsThFactory["CEWHQ3_11"] = [=](const StandardModel& SM) { return new CEWHQ311(SM, muEW); };
    obsThFactory["CEWHQ3_22"] = [=](const StandardModel& SM) { return new CEWHQ322(SM, muEW); };
    obsThFactory["CEWHQ3_33"] = [=](const StandardModel& SM) { return new CEWHQ333(SM, muEW); };
    obsThFactory["CEWHQd_33"] = [=](const StandardModel& SM) { return new CEWHQd33(SM, muEW); };
    obsThFactory["CEWHQu_33"] = [=](const StandardModel& SM) { return new CEWHQu33(SM, muEW); };
    obsThFactory["CEWHe_11"] = [=](const StandardModel& SM) { return new CEWHe11(SM, muEW); };
    obsThFactory["CEWHe_22"] = [=](const StandardModel& SM) { return new CEWHe22(SM, muEW); };
    obsThFactory["CEWHe_33"] = [=](const StandardModel& SM) { return new CEWHe33(SM, muEW); };
    obsThFactory["CEWHu_11"] = [=](const StandardModel& SM) { return new CEWHu11(SM, muEW); };
    obsThFactory["CEWHu_22"] = [=](const StandardModel& SM) { return new CEWHu22(SM, muEW); };
    obsThFactory["CEWHu_33"] = [=](const StandardModel& SM) { return new CEWHu33(SM, muEW); };
    obsThFactory["CEWHd_11"] = [=](const StandardModel& SM) { return new CEWHd11(SM, muEW); };
    obsThFactory["CEWHd_22"] = [=](const StandardModel& SM) { return new CEWHd22(SM, muEW); };
    obsThFactory["CEWHd_33"] = [=](const StandardModel& SM) { return new CEWHd33(SM, muEW); };
    obsThFactory["CEWll"] = [=](const StandardModel& SM) { return new CEWll(SM, muEW); };
    
    //-----  The same, in the quark mass basis  ----------
    obsThFactory["CEWHQ1_uu"] = [=](const StandardModel& SM) { return new CEWHQ1uu(SM, muEW); };
    obsThFactory["CEWHQ1_cc"] = [=](const StandardModel& SM) { return new CEWHQ1cc(SM, muEW); };
    obsThFactory["CEWHQ1_tt"] = [=](const StandardModel& SM) { return new CEWHQ1tt(SM, muEW); };
    //
    obsThFactory["CEWHQ1_dd"] = [=](const StandardModel& SM) { return new CEWHQ1dd(SM, muEW); };
    obsThFactory["CEWHQ1_ss"] = [=](const StandardModel& SM) { return new CEWHQ1ss(SM, muEW); };
    obsThFactory["CEWHQ1_bb"] = [=](const StandardModel& SM) { return new CEWHQ1bb(SM, muEW); };
    //
    obsThFactory["CEWHQ3_uu"] = [=](const StandardModel& SM) { return new CEWHQ3uu(SM, muEW); };
    obsThFactory["CEWHQ3_cc"] = [=](const StandardModel& SM) { return new CEWHQ3cc(SM, muEW); };
    obsThFactory["CEWHQ3_tt"] = [=](const StandardModel& SM) { return new CEWHQ3tt(SM, muEW); };
    //
    obsThFactory["CEWHQ3_dd"] = [=](const StandardModel& SM) { return new CEWHQ3dd(SM, muEW); };
    obsThFactory["CEWHQ3_ss"] = [=](const StandardModel& SM) { return new CEWHQ3ss(SM, muEW); };
    obsThFactory["CEWHQ3_bb"] = [=](const StandardModel& SM) { return new CEWHQ3bb(SM, muEW); };
    //    
    obsThFactory["CEWHu_uu"] = [=](const StandardModel& SM) { return new CEWHuuu(SM, muEW); };
    obsThFactory["CEWHu_cc"] = [=](const StandardModel& SM) { return new CEWHucc(SM, muEW); };
    obsThFactory["CEWHu_tt"] = [=](const StandardModel& SM) { return new CEWHutt(SM, muEW); };
    //
    obsThFactory["CEWHd_dd"] = [=](const StandardModel& SM) { return new CEWHddd(SM, muEW); };
    obsThFactory["CEWHd_ss"] = [=](const StandardModel& SM) { return new CEWHdss(SM, muEW); };
    obsThFactory["CEWHd_bb"] = [=](const StandardModel& SM) { return new CEWHdbb(SM, muEW); };
    
    
    //-----  Combinations of Warsaw basis coefficients constrained by EWPO  ----------
    obsThFactory["CEWHL1_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHL111(SM, 1000.); };
    obsThFactory["CEWHL1_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHL122(SM, 1000.); };
    obsThFactory["CEWHL1_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHL133(SM, 1000.); };
    obsThFactory["CEWHL3_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHL311(SM, 1000.); };
    obsThFactory["CEWHL3_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHL322(SM, 1000.); };
    obsThFactory["CEWHL3_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHL333(SM, 1000.); };
    obsThFactory["CEWHQ1_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ111(SM, 1000.); };
    obsThFactory["CEWHQ1_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ122(SM, 1000.); };
    obsThFactory["CEWHQ1_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ133(SM, 1000.); };
    obsThFactory["CEWHQ3_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ311(SM, 1000.); };
    obsThFactory["CEWHQ3_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ322(SM, 1000.); };
    obsThFactory["CEWHQ3_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ333(SM, 1000.); };
    obsThFactory["CEWHQd_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHQd33(SM, 1000.); };
    obsThFactory["CEWHQu_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHQu33(SM, 1000.); };
    obsThFactory["CEWHe_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHe11(SM, 1000.); };
    obsThFactory["CEWHe_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHe22(SM, 1000.); };
    obsThFactory["CEWHe_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHe33(SM, 1000.); };
    obsThFactory["CEWHu_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHu11(SM, 1000.); };
    obsThFactory["CEWHu_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHu22(SM, 1000.); };
    obsThFactory["CEWHu_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHu33(SM, 1000.); };
    obsThFactory["CEWHd_11_1TeV"] = [=](const StandardModel& SM) { return new CEWHd11(SM, 1000.); };
    obsThFactory["CEWHd_22_1TeV"] = [=](const StandardModel& SM) { return new CEWHd22(SM, 1000.); };
    obsThFactory["CEWHd_33_1TeV"] = [=](const StandardModel& SM) { return new CEWHd33(SM, 1000.); };
    obsThFactory["CEWll_1TeV"] = [=](const StandardModel& SM) { return new CEWll(SM, 1000.); };
    
    //-----  The same, in the quark mass basis  ----------
    obsThFactory["CEWHQ1_uu_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1uu(SM, 1000.); };
    obsThFactory["CEWHQ1_cc_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1cc(SM, 1000.); };
    obsThFactory["CEWHQ1_tt_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1tt(SM, 1000.); };
    //
    obsThFactory["CEWHQ1_dd_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1dd(SM, 1000.); };
    obsThFactory["CEWHQ1_ss_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1ss(SM, 1000.); };
    obsThFactory["CEWHQ1_bb_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ1bb(SM, 1000.); };
    //
    obsThFactory["CEWHQ3_uu_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3uu(SM, 1000.); };
    obsThFactory["CEWHQ3_cc_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3cc(SM, 1000.); };
    obsThFactory["CEWHQ3_tt_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3tt(SM, 1000.); };
    //
    obsThFactory["CEWHQ3_dd_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3dd(SM, 1000.); };
    obsThFactory["CEWHQ3_ss_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3ss(SM, 1000.); };
    obsThFactory["CEWHQ3_bb_1TeV"] = [=](const StandardModel& SM) { return new CEWHQ3bb(SM, 1000.); };
    //    
    obsThFactory["CEWHu_uu_1TeV"] = [=](const StandardModel& SM) { return new CEWHuuu(SM, 1000.); };
    obsThFactory["CEWHu_cc_1TeV"] = [=](const StandardModel& SM) { return new CEWHucc(SM, 1000.); };
    obsThFactory["CEWHu_tt_1TeV"] = [=](const StandardModel& SM) { return new CEWHutt(SM, 1000.); };
    //
    obsThFactory["CEWHd_dd_1TeV"] = [=](const StandardModel& SM) { return new CEWHddd(SM, 1000.); };
    obsThFactory["CEWHd_ss_1TeV"] = [=](const StandardModel& SM) { return new CEWHdss(SM, 1000.); };
    obsThFactory["CEWHd_bb_1TeV"] = [=](const StandardModel& SM) { return new CEWHdbb(SM, 1000.); };


    //-----  Auxiliary observables to work with new physics  ----------
    obsThFactory["AuxObsNP1"] = [](const StandardModel& SM) { return new AuxObsNP1(SM); };
    obsThFactory["AuxObsNP2"] = [](const StandardModel& SM) { return new AuxObsNP2(SM); };
    obsThFactory["AuxObsNP3"] = [](const StandardModel& SM) { return new AuxObsNP3(SM); };
    obsThFactory["AuxObsNP4"] = [](const StandardModel& SM) { return new AuxObsNP4(SM); };
    obsThFactory["AuxObsNP5"] = [](const StandardModel& SM) { return new AuxObsNP5(SM); };
    obsThFactory["AuxObsNP6"] = [](const StandardModel& SM) { return new AuxObsNP6(SM); };
    obsThFactory["AuxObsNP7"] = [](const StandardModel& SM) { return new AuxObsNP7(SM); };
    obsThFactory["AuxObsNP8"] = [](const StandardModel& SM) { return new AuxObsNP8(SM); };
    obsThFactory["AuxObsNP9"] = [](const StandardModel& SM) { return new AuxObsNP9(SM); };
    obsThFactory["AuxObsNP10"] = [](const StandardModel& SM) { return new AuxObsNP10(SM); };
    obsThFactory["AuxObsNP11"] = [](const StandardModel& SM) { return new AuxObsNP11(SM); };
    obsThFactory["AuxObsNP12"] = [](const StandardModel& SM) { return new AuxObsNP12(SM); };
    obsThFactory["AuxObsNP13"] = [](const StandardModel& SM) { return new AuxObsNP13(SM); };
    obsThFactory["AuxObsNP14"] = [](const StandardModel& SM) { return new AuxObsNP14(SM); };
    obsThFactory["AuxObsNP15"] = [](const StandardModel& SM) { return new AuxObsNP15(SM); };
    obsThFactory["AuxObsNP16"] = [](const StandardModel& SM) { return new AuxObsNP16(SM); };
    obsThFactory["AuxObsNP17"] = [](const StandardModel& SM) { return new AuxObsNP17(SM); };
    obsThFactory["AuxObsNP18"] = [](const StandardModel& SM) { return new AuxObsNP18(SM); };
    obsThFactory["AuxObsNP19"] = [](const StandardModel& SM) { return new AuxObsNP19(SM); };
    obsThFactory["AuxObsNP20"] = [](const StandardModel& SM) { return new AuxObsNP20(SM); };
    obsThFactory["AuxObsNP21"] = [](const StandardModel& SM) { return new AuxObsNP21(SM); };
    obsThFactory["AuxObsNP22"] = [](const StandardModel& SM) { return new AuxObsNP22(SM); };
    obsThFactory["AuxObsNP23"] = [](const StandardModel& SM) { return new AuxObsNP23(SM); };
    obsThFactory["AuxObsNP24"] = [](const StandardModel& SM) { return new AuxObsNP24(SM); };
    obsThFactory["AuxObsNP25"] = [](const StandardModel& SM) { return new AuxObsNP25(SM); };
    obsThFactory["AuxObsNP26"] = [](const StandardModel& SM) { return new AuxObsNP26(SM); };
    obsThFactory["AuxObsNP27"] = [](const StandardModel& SM) { return new AuxObsNP27(SM); };
    obsThFactory["AuxObsNP28"] = [](const StandardModel& SM) { return new AuxObsNP28(SM); };
    obsThFactory["AuxObsNP29"] = [](const StandardModel& SM) { return new AuxObsNP29(SM); };
    obsThFactory["AuxObsNP30"] = [](const StandardModel& SM) { return new AuxObsNP30(SM); };

}
