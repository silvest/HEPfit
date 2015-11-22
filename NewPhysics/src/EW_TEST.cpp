/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_TEST.h"
#include "NPEpsilons.h"
#include "NPSTU.h"
#include "EWSMApproximateFormulae.h"
#include <stdexcept>

EW_TEST::EW_TEST(const std::string mode_i, const std::string obsname_i, const StandardModel& SM_i)
: ThObservable(SM_i), mode(mode_i), obsname(obsname_i), SM(SM_i)
{
    if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
        myEW_ABC = new EW_ABC(dynamic_cast<const NPEpsilons &> (SM));
    else if (mode.compare("BURGESS") == 0)
        myEW_BURGESS = new EW_BURGESS(dynamic_cast<const NPSTU &> (SM));
    else if (mode.compare("CHMN") == 0)
        myEW_CHMN = new EW_CHMN(SM);
    else
        throw std::runtime_error("EW_TEST::EW_TEST(): Wrong mode " + mode);
}

double EW_TEST::computeThValue()
{
    if (obsname.compare("Mw") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->Mw(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->Mw(true);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->Mw(SM.getTrueSM().Mw());
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->Mw();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("GammaW") == 0) {
        if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->GammaW(SM.getTrueSM().GammaW());
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->GammaW();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("GammaZ") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->GammaZ(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->GammaZ(true);
        else if (mode.compare("BURGESS") == 0) {
            double Gamma_Z_SM;
            if (SM.IsFlagNoApproximateGammaZ())
                Gamma_Z_SM = SM.getTrueSM().Gamma_Z();
            else
                Gamma_Z_SM = SM.getTrueSM().getMyApproximateFormulae()->X_extended("GammaZ");
            Gamma_Z_SM += SM.getDelGammaZ(); /* Theoretical uncertainty */
            return myEW_BURGESS->GammaZ(Gamma_Z_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->GammaZ();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("sigmaHadron") == 0) {
        if (mode.compare("ABC") == 0)
            return ( myEW_ABC->sigma0_had(false) * SM.GeVminus2_to_nb);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->sigma0_had(true);
        else if (mode.compare("BURGESS") == 0) {
            double sigma_had_SM;
            if (SM.IsFlagNoApproximateGammaZ())
                sigma_had_SM = SM.getTrueSM().sigma0_had();
            else
                sigma_had_SM = SM.getTrueSM().getMyApproximateFormulae()->X_extended("sigmaHadron");
            return ( myEW_BURGESS->sigmaHadron(sigma_had_SM,
                    SM.Gamma_Z(), SM.Gamma_had(), SM.GammaZ(SM.getLeptons(StandardModel::ELECTRON)))
                    * SM.GeVminus2_to_nb);
        } else if (mode.compare("CHMN") == 0)
            return ( myEW_CHMN->sigma0_had() * SM.GeVminus2_to_nb);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("sin2thetaEff") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->sin2thetaEff(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->sin2thetaEff(true);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->sin2thetaEff(SM.sin2thetaEff(SM.getLeptons(SM.ELECTRON)));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->sin2thetaEff();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("PtauPol") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->A_l(SM.TAU, false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->A_l(SM.TAU, true);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->PtauPol(SM.A_f(SM.getLeptons(SM.TAU)));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_l(SM.TAU);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Alepton") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->A_l(SM.ELECTRON, false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->A_l(SM.ELECTRON, true);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->Alepton(SM.A_f(SM.getLeptons(SM.ELECTRON)));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Acharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->A_q(SM.CHARM);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->Acharm(SM.A_f(SM.getQuarks(SM.CHARM)), SM.A_f(SM.getLeptons(SM.ELECTRON)));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_q(SM.CHARM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Abottom") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->A_b();
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->Abottom(SM.A_f(SM.getQuarks(SM.BOTTOM)), SM.A_f(SM.getLeptons(SM.ELECTRON)));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_q(SM.BOTTOM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("AFBlepton") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->AFB_l(SM.ELECTRON, false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->AFB_l(SM.ELECTRON, true);
        else if (mode.compare("BURGESS") == 0) {
            double AFB_l_SM = 3.0 / 4.0 * SM.A_f(SM.getLeptons(SM.ELECTRON)) * SM.A_f(SM.getLeptons(SM.ELECTRON));
            return myEW_BURGESS->AFBlepton(AFB_l_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->AFB_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("AFBcharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->AFB_c();
        else if (mode.compare("BURGESS") == 0) {
            double AFB_c_SM = 3.0 / 4.0 * SM.A_f(SM.getLeptons(SM.ELECTRON)) * SM.A_f(SM.getQuarks(SM.CHARM));
            return myEW_BURGESS->AFBcharm(AFB_c_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->AFB_q(SM.CHARM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("AFBbottom") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->AFB_b();
        else if (mode.compare("BURGESS") == 0) {
            double AFB_b_SM = 3.0 / 4.0 * SM.A_f(SM.getLeptons(SM.ELECTRON)) * SM.A_f(SM.getQuarks(SM.BOTTOM));
            return myEW_BURGESS->AFBbottom(AFB_b_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->AFB_q(SM.BOTTOM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Rlepton") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->R_l(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->R_l(true);
        else if (mode.compare("BURGESS") == 0) {
            double R0_l_SM;
            if (SM.IsFlagNoApproximateGammaZ())
                R0_l_SM = SM.getMyApproximateFormulae()->X_extended("R0_lepton");
            else
                R0_l_SM = SM.Gamma_had() / SM.GammaZ(SM.getLeptons(SM.ELECTRON));
            return myEW_BURGESS->Rlepton(R0_l_SM, SM.Gamma_had(), SM.GammaZ(SM.getLeptons(SM.ELECTRON)));
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->R_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Rcharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->R_c();
        else if (mode.compare("BURGESS") == 0) {
            double R0_c_SM;
            if (SM.IsFlagNoApproximateGammaZ())
                R0_c_SM = SM.getMyApproximateFormulae()->X_extended("R0_charm");
            else
                R0_c_SM = SM.GammaZ(SM.getQuarks(SM.CHARM)) / SM.Gamma_had();
            return myEW_BURGESS->Rcharm(R0_c_SM, SM.Gamma_had());
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->R_c();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Rbottom") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->R_b(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->R_b(true);
        else if (mode.compare("BURGESS") == 0) {
            double R0_b_SM;
            if (SM.IsFlagNoApproximateGammaZ())
                R0_b_SM = SM.getMyApproximateFormulae()->X_extended("R0_bottom");
            else
                R0_b_SM = SM.GammaZ(SM.getQuarks(SM.BOTTOM)) / SM.Gamma_had();
            return myEW_BURGESS->Rbottom(R0_b_SM, SM.Gamma_had(), SM.GammaZ(SM.getQuarks(SM.BOTTOM)));
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->R_b();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else
        throw std::runtime_error("EW_TEST::computeThValue(): Wrong obsname " + obsname);
}

