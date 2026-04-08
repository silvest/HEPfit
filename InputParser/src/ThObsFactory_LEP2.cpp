/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObsFactoryConstants.h"
#include "EWObservables.h"
#include "NP_couplings.h"
#include <boost/lexical_cast.hpp>
using namespace ThObsConst;

void ThObsFactory::registerLEP2Observables()
{
    //-----  LEP-II two-fermion processes  -----

    for (int i = 0; i < 12; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s[i]);
        obsThFactory["sigmaqLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2sigmaHadron(SM, sqrt_s_LEP2eeff[i]); };
        obsThFactory["sigmamuLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2sigmaMu(SM, sqrt_s_LEP2eeff[i]); };
        obsThFactory["sigmatauLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2sigmaTau(SM, sqrt_s_LEP2eeff[i]); };
        obsThFactory["AFBmuLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2AFBmu(SM, sqrt_s_LEP2eeff[i]); };
        obsThFactory["AFBtauLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2AFBtau(SM, sqrt_s_LEP2eeff[i]); };
    }

    for (int i = 0; i < 8; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_sDiffll[i]);

        for (int j = 0; j < 5; j++) {
        std::string cos_str = boost::lexical_cast<std::string, double>(fabs(10.*cos_Diffll[j]));
        std::string cosee_str = boost::lexical_cast<std::string, double>(fabs(10.*cos_DiffeeInp[j]));

        obsThFactory["dsigmadcosmuLEP2_" + sqrt_s_str + "_m0" + cos_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosMu(SM, sqrt_sDiffll[i], cos_Diffll[j], cosmin_Diffll[j], cosmax_Diffll[j]); };
        obsThFactory["dsigmadcostauLEP2_" + sqrt_s_str + "_m0" + cos_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosTau(SM, sqrt_sDiffll[i], cos_Diffll[j], cosmin_Diffll[j], cosmax_Diffll[j]); };

        if (i>0) {obsThFactory["dsigmadcoseLEP2_" + sqrt_s_str + "_m0" + cosee_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosElectron(SM, sqrt_sDiffll[i], cos_Diffee[j], cosmin_Diffee[j], cosmax_Diffee[j]); };}
        }

        for (int j = 5; j < 10; j++) {
        std::string cos_str = boost::lexical_cast<std::string, double>(10.*cos_Diffll[j]);
        std::string cosee_str = boost::lexical_cast<std::string, double>(fabs(10.*cos_DiffeeInp[j]));

        obsThFactory["dsigmadcosmuLEP2_" + sqrt_s_str + "_0" + cos_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosMu(SM, sqrt_sDiffll[i], cos_Diffll[j], cosmin_Diffll[j], cosmax_Diffll[j]); };
        obsThFactory["dsigmadcostauLEP2_" + sqrt_s_str + "_0" + cos_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosTau(SM, sqrt_sDiffll[i], cos_Diffll[j], cosmin_Diffll[j], cosmax_Diffll[j]); };

        if (i>0) {obsThFactory["dsigmadcoseLEP2_" + sqrt_s_str + "_0" + cosee_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosElectron(SM, sqrt_sDiffll[i], cos_Diffee[j], cosmin_Diffee[j], cosmax_Diffee[j]); };}
        }

        if (i>0) {
        for (int j = 10; j < 15; j++) {
        std::string cosee_str = boost::lexical_cast<std::string, double>(fabs(10.*cos_DiffeeInp[j]));
        obsThFactory["dsigmadcoseLEP2_" + sqrt_s_str + "_0" + cosee_str] = [=](const StandardModel& SM) { return new LEP2dsigmadcosElectron(SM, sqrt_sDiffll[i], cos_Diffee[j], cosmin_Diffee[j], cosmax_Diffee[j]); };
        }
        }

    }

    //-----  e+ e- two-fermion processes  -----

    for (int i = 0; i < 16; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_see[i]);

        // Unpolarized
        // Total cross sections
        obsThFactory["sigmaeeee_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0., 0., sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0., 0., sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["Reeee_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0., 0., sqrt_s_eeff[i]); };

        obsThFactory["Reess_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["Reecc_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["Reebb_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0., 0., sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0., 0., sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0., 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0., 0., sqrt_s_eeff[i]); };
                
        // Polarized: Pe-: -80% Pe+: 0%
        // Total cross sections
        obsThFactory["sigmaeeee_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_m80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, -0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, -0.8, 0., sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reeee_m80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, -0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["Reess_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reecc_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reebb_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, -0.8, 0., sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_m80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, -0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, -0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_m80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, -0.8, 0., sqrt_s_eeff[i]); };

        // Polarized: Pe-: 80% Pe+: 0%
        // Total cross sections
        obsThFactory["sigmaeeee_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_p80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0.8, 0., sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reeee_p80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["Reess_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reecc_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["Reebb_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0.8, 0., sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_p80_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0.8, 0., sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0.8, 0., sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_p80_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0.8, 0., sqrt_s_eeff[i]); };

        // Polarized: Pe-: -80% Pe+: +20%
        // Total cross sections
        obsThFactory["sigmaeeee_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_m80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_m80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["Reess_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_m80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, -0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_m80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, -0.8, 0.2, sqrt_s_eeff[i]); };

        // Polarized: Pe-: 80% Pe+: -20%
        // Total cross sections
        obsThFactory["sigmaeeee_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_p80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0.8, -0.2, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_p80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["Reess_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0.8, -0.2, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_p80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_p80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0.8, -0.2, sqrt_s_eeff[i]); };
               
        // Polarized: Pe-: -80% Pe+: -20%
        // Total cross sections
        obsThFactory["sigmaeeee_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_m80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_m80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["Reess_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_m80_m20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, -0.8, -0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_m80_m20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, -0.8, -0.2, sqrt_s_eeff[i]); };

        // Polarized: Pe-: 80% Pe+: 20%
        // Total cross sections
        obsThFactory["sigmaeeee_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_p80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_p80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["Reess_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_p80_p20_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0.8, 0.2, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_p80_p20_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0.8, 0.2, sqrt_s_eeff[i]); };

        // Polarized: Pe-: -80% Pe+: +30%
        // Total cross sections
        obsThFactory["sigmaeeee_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_m80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_m80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["Reess_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_m80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, -0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_m80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, -0.8, 0.3, sqrt_s_eeff[i]); };

        // Polarized: Pe-: 80% Pe+: -30%
        // Total cross sections
        obsThFactory["sigmaeeee_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_p80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0.8, -0.3, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_p80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["Reess_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0.8, -0.3, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_p80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_p80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0.8, -0.3, sqrt_s_eeff[i]); };
                
        // Polarized: Pe-: -80% Pe+: -30%
        // Total cross sections
        obsThFactory["sigmaeeee_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_m80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_m80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["Reess_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_m80_m30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, -0.8, -0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_m80_m30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, -0.8, -0.3, sqrt_s_eeff[i]); };

        // Polarized: Pe-: 80% Pe+: 30%
        // Total cross sections
        obsThFactory["sigmaeeee_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaElectron(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeeee_p80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmatsubElectron(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeemumu_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaMu(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeetautau_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaTau(SM, 0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["sigmaeeqq_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaHadron(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeess_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaStrange(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeecc_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaCharm(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["sigmaeebb_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffsigmaBottom(SM, 0.8, 0.3, sqrt_s_eeff[i]); };

        //Ratios
        obsThFactory["Reeee_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectron(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reeee_p80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRelectrontsub(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reemumu_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRmuon(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reetautau_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRtau(SM, 0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["Reess_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRstrange(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reecc_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRcharm(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["Reebb_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffRbottom(SM, 0.8, 0.3, sqrt_s_eeff[i]); };

        //Asymmetries
        obsThFactory["AFBeeee_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectron(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeeee_p80_p30_tsub_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBelectrontsub(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeemumu_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBmu(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeetautau_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBtau(SM, 0.8, 0.3, sqrt_s_eeff[i]); };

        obsThFactory["AFBeess_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBstrange(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeecc_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBcharm(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
        obsThFactory["AFBeebb_p80_p30_" + sqrt_s_str] = [=](const StandardModel& SM) { return new eeffAFBbottom(SM, 0.8, 0.3, sqrt_s_eeff[i]); };
                
    }

    //-----  Wilson coefficients of Top operators in the conventions of the LHC Top WG ----------
    //-----  (Energy-independent definition) ----------
    //
    obsThFactory["cQQ1_TWG"] = [=](const StandardModel& SM) { return new TWGcQQ1(SM, muEW); };
    obsThFactory["cQQ8_TWG"] = [=](const StandardModel& SM) { return new TWGcQQ8(SM, muEW); };
    obsThFactory["ctt1_TWG"] = [=](const StandardModel& SM) { return new TWGctt1(SM, muEW); };
    obsThFactory["cQt1_TWG"] = [=](const StandardModel& SM) { return new TWGcQt1(SM, muEW); };
    obsThFactory["cQt8_TWG"] = [=](const StandardModel& SM) { return new TWGcQt8(SM, muEW); };
    //
    obsThFactory["cQq31_TWG"] = [=](const StandardModel& SM) { return new TWGcQq31(SM, muEW); };
    obsThFactory["cQq38_TWG"] = [=](const StandardModel& SM) { return new TWGcQq38(SM, muEW); };
    obsThFactory["cQq11_TWG"] = [=](const StandardModel& SM) { return new TWGcQq11(SM, muEW); };
    obsThFactory["cQq18_TWG"] = [=](const StandardModel& SM) { return new TWGcQq18(SM, muEW); };
    obsThFactory["cQu1_TWG"] = [=](const StandardModel& SM) { return new TWGcQu1(SM, muEW); };
    obsThFactory["cQu8_TWG"] = [=](const StandardModel& SM) { return new TWGcQu8(SM, muEW); };
    obsThFactory["cQd1_TWG"] = [=](const StandardModel& SM) { return new TWGcQd1(SM, muEW); };
    obsThFactory["cQd8_TWG"] = [=](const StandardModel& SM) { return new TWGcQd8(SM, muEW); };
    obsThFactory["ctq1_TWG"] = [=](const StandardModel& SM) { return new TWGctq1(SM, muEW); };
    obsThFactory["ctq8_TWG"] = [=](const StandardModel& SM) { return new TWGctq8(SM, muEW); };
    obsThFactory["ctu1_TWG"] = [=](const StandardModel& SM) { return new TWGctu1(SM, muEW); };
    obsThFactory["ctu8_TWG"] = [=](const StandardModel& SM) { return new TWGctu8(SM, muEW); };
    obsThFactory["ctd1_TWG"] = [=](const StandardModel& SM) { return new TWGctd1(SM, muEW); };
    obsThFactory["ctd8_TWG"] = [=](const StandardModel& SM) { return new TWGctd8(SM, muEW); };
    //
    obsThFactory["ctH_TWG"] = [=](const StandardModel& SM) { return new TWGctH(SM, muEW); };
    obsThFactory["cHQm_TWG"] = [=](const StandardModel& SM) { return new TWGcHQm(SM, muEW); };
    obsThFactory["cHQp_TWG"] = [=](const StandardModel& SM) { return new TWGcHQp(SM, muEW); };
    obsThFactory["cHQ3_TWG"] = [=](const StandardModel& SM) { return new TWGcHQ3(SM, muEW); };
    obsThFactory["cHt_TWG"] = [=](const StandardModel& SM) { return new TWGcHt(SM, muEW); };
    obsThFactory["cHb_TWG"] = [=](const StandardModel& SM) { return new TWGcHb(SM, muEW); };
    obsThFactory["cHtb_TWG"] = [=](const StandardModel& SM) { return new TWGcHtb(SM, muEW); };
    //
    obsThFactory["ctW_TWG"] = [=](const StandardModel& SM) { return new TWGctW(SM, muEW); };
    obsThFactory["ImctW_TWG"] = [=](const StandardModel& SM) { return new TWGImctW(SM, muEW); };
    obsThFactory["ctZ_TWG"] = [=](const StandardModel& SM) { return new TWGctZ(SM, muEW); };
    obsThFactory["ImctZ_TWG"] = [=](const StandardModel& SM) { return new TWGImctZ(SM, muEW); };
    obsThFactory["ctG_TWG"] = [=](const StandardModel& SM) { return new TWGctG(SM, muEW); };
    //
    obsThFactory["cbW_TWG"] = [=](const StandardModel& SM) { return new TWGcbW(SM, muEW); };
    //
    obsThFactory["cQlM_TWG"] = [=](const StandardModel& SM) { return new TWGcQlM(SM, muEW); };
    obsThFactory["cQlP_TWG"] = [=](const StandardModel& SM) { return new TWGcQlP(SM, muEW); };
    obsThFactory["cQl3_TWG"] = [=](const StandardModel& SM) { return new TWGcQl3(SM, muEW); };
    obsThFactory["cQe_TWG"] = [=](const StandardModel& SM) { return new TWGcQe(SM, muEW); };
    obsThFactory["ctl_TWG"] = [=](const StandardModel& SM) { return new TWGctl(SM, muEW); };
    obsThFactory["cte_TWG"] = [=](const StandardModel& SM) { return new TWGcte(SM, muEW); };
    obsThFactory["ctlS_TWG"] = [=](const StandardModel& SM) { return new TWGctlS(SM, muEW); };
    obsThFactory["ctlT_TWG"] = [=](const StandardModel& SM) { return new TWGctlT(SM, muEW); };

    //-----  Wilson coefficients of Top operators in the conventions of the LHC Top WG ----------
    //-----  (Energy-dependent definition: Only for energies above the ttbar threshold ) --------

    for (int i = 6; i < 16; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_see[i]);

        obsThFactory["cQQ1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQQ1(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQQ8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQQ8(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctt1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctt1(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQt1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQt1(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQt8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQt8(SM, sqrt_s_eeff[i]); };
        //
        obsThFactory["cQq31_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQq31(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQq38_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQq38(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQq11_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQq11(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQq18_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQq18(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQu1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQu1(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQu8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQu8(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQd1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQd1(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQd8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQd8(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctq1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctq1(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctq8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctq8(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctu1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctu1(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctu8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctu8(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctd1_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctd1(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctd8_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctd8(SM, sqrt_s_eeff[i]); };
        //
        obsThFactory["ctH_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctH(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHQm_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHQm(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHQp_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHQp(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHQ3_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHQ3(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHt_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHt(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHb_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHb(SM, sqrt_s_eeff[i]); };
        obsThFactory["cHtb_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcHtb(SM, sqrt_s_eeff[i]); };
        //
        obsThFactory["ctW_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctW(SM, sqrt_s_eeff[i]); };
        obsThFactory["ImctW_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGImctW(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctZ_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctZ(SM, sqrt_s_eeff[i]); };
        obsThFactory["ImctZ_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGImctZ(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctG_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctG(SM, sqrt_s_eeff[i]); };
        //
        obsThFactory["cbW_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcbW(SM, sqrt_s_eeff[i]); };
        //
        obsThFactory["cQlM_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQlM(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQlP_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQlP(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQl3_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQl3(SM, sqrt_s_eeff[i]); };
        obsThFactory["cQe_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcQe(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctl_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctl(SM, sqrt_s_eeff[i]); };
        obsThFactory["cte_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGcte(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctlS_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctlS(SM, sqrt_s_eeff[i]); };
        obsThFactory["ctlT_TWG_" + sqrt_s_str] = [=](const StandardModel& SM) { return new TWGctlT(SM, sqrt_s_eeff[i]); };

    }
    
    //-----  (Extra temporary definition, evaluated at 10 TeV) ----------
    //
    obsThFactory["cQQ1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQQ1(SM, 10000.); };
    obsThFactory["cQQ8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQQ8(SM, 10000.); };
    obsThFactory["ctt1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctt1(SM, 10000.); };
    obsThFactory["cQt1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQt1(SM, 10000.); };
    obsThFactory["cQt8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQt8(SM, 10000.); };
    //
    obsThFactory["cQq31_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQq31(SM, 10000.); };
    obsThFactory["cQq38_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQq38(SM, 10000.); };
    obsThFactory["cQq11_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQq11(SM, 10000.); };
    obsThFactory["cQq18_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQq18(SM, 10000.); };
    obsThFactory["cQu1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQu1(SM, 10000.); };
    obsThFactory["cQu8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQu8(SM, 10000.); };
    obsThFactory["cQd1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQd1(SM, 10000.); };
    obsThFactory["cQd8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQd8(SM, 10000.); };
    obsThFactory["ctq1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctq1(SM, 10000.); };
    obsThFactory["ctq8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctq8(SM, 10000.); };
    obsThFactory["ctu1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctu1(SM, 10000.); };
    obsThFactory["ctu8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctu8(SM, 10000.); };
    obsThFactory["ctd1_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctd1(SM, 10000.); };
    obsThFactory["ctd8_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctd8(SM, 10000.); };
    //
    obsThFactory["ctH_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctH(SM, 10000.); };
    obsThFactory["cHQm_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHQm(SM, 10000.); };
    obsThFactory["cHQp_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHQp(SM, 10000.); };
    obsThFactory["cHQ3_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHQ3(SM, 10000.); };
    obsThFactory["cHt_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHt(SM, 10000.); };
    obsThFactory["cHb_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHb(SM, 10000.); };
    obsThFactory["cHtb_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcHtb(SM, 10000.); };
    //
    obsThFactory["ctW_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctW(SM, 10000.); };
    obsThFactory["ImctW_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGImctW(SM, 10000.); };
    obsThFactory["ctZ_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctZ(SM, 10000.); };
    obsThFactory["ImctZ_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGImctZ(SM, 10000.); };
    obsThFactory["ctG_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctG(SM, 10000.); };
    //
    obsThFactory["cbW_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcbW(SM, 10000.); };
    //
    obsThFactory["cQlM_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQlM(SM, 10000.); };
    obsThFactory["cQlP_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQlP(SM, 10000.); };
    obsThFactory["cQl3_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQl3(SM, 10000.); };
    obsThFactory["cQe_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcQe(SM, 10000.); };
    obsThFactory["ctl_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctl(SM, 10000.); };
    obsThFactory["cte_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGcte(SM, 10000.); };
    obsThFactory["ctlS_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctlS(SM, 10000.); };
    obsThFactory["ctlT_TWG_10TeV"] = [=](const StandardModel& SM) { return new TWGctlT(SM, 10000.); };

    /* BEGIN: REMOVE FROM THE PACKAGE */

    for (int i = 0; i < 10; i++) {
        std::string sqrt_s_str = boost::lexical_cast<std::string, double>(sqrt_s_HF[i]);
        obsThFactory["AFBbottomLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2AFBbottom(SM, sqrt_s_LEP2_HF[i]); };
        obsThFactory["AFBcharmLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2AFBcharm(SM, sqrt_s_LEP2_HF[i]); };
        obsThFactory["RbottomLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2Rbottom(SM, sqrt_s_LEP2_HF[i]); };
        obsThFactory["RcharmLEP2_" + sqrt_s_str] = [=](const StandardModel& SM) { return new LEP2Rcharm(SM, sqrt_s_LEP2_HF[i]); };
    }
    /* END: REMOVE FROM THE PACKAGE */


    //-----  Low Energy EW observables  -----
    // Parity violation
    obsThFactory["QWeMoller"] = [](const StandardModel& SM) { return new QWe(SM); };
    obsThFactory["QWproton"] = [](const StandardModel& SM) { return new QWp(SM); };
    obsThFactory["QWCs133_55"] = [=](const StandardModel& SM) { return new QWAPV(SM, 55, 78); };
    obsThFactory["QWTl205_81"] = [=](const StandardModel& SM) { return new QWAPV(SM, 81, 124); };
    // Neutrino scattering
    obsThFactory["gL2_nuN"] = [](const StandardModel& SM) { return new gLnuN2(SM); };
    obsThFactory["gR2_nuN"] = [](const StandardModel& SM) { return new gRnuN2(SM); };
    obsThFactory["gV_nue"] = [](const StandardModel& SM) { return new gVnue(SM); };
    obsThFactory["gA_nue"] = [](const StandardModel& SM) { return new gAnue(SM); };
    // Temporary observable: for testing
    obsThFactory["amuongminus2"] = [](const StandardModel& SM) { return new agminus2muon(SM); };

    //-----  Lepton decays  -----
    // LFU tests in Tau decays
    obsThFactory["gmuge_TauLFU"] = [](const StandardModel& SM) { return new gmugeTauLFU(SM); };
    obsThFactory["gtaugmu_TauLFU"] = [](const StandardModel& SM) { return new gtaugmuTauLFU(SM); };
    obsThFactory["gtauge_TauLFU"] = [](const StandardModel& SM) { return new gtaugeTauLFU(SM); };

    obsThFactory["gtaugmuPi_TauLFU"] = [](const StandardModel& SM) { return new gtaugmuPiTauLFU(SM); };
    obsThFactory["gtaugmuK_TauLFU"] = [](const StandardModel& SM) { return new gtaugmuKTauLFU(SM); };

}
