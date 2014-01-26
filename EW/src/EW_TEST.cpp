/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <EWSM.h>
#include "EW_TEST.h"


EW_TEST::EW_TEST(const std::string mode_i, const std::string obsname_i, const EW& EW_i)
: ThObservable(EW_i), mode(mode_i), obsname(obsname_i), myEW(EW_i), SM(EW_i.getSM())
{
    if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
       myEW_ABC = new EW_ABC(myEW.getSM());
    else if (mode.compare("BURGESS") == 0)
       myEW_BURGESS = new EW_BURGESS(myEW.getSM());
    else if (mode.compare("CHMN") == 0 )
       myEW_CHMN = new EW_CHMN(myEW.getSM());
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
            return myEW_BURGESS->Mw(SM.getEWSM()->Mw_SM());
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->Mw();
        else 
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("GammaW") == 0) {
        if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->GammaW(SM.getEWSM()->GammaW_SM());
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
            if (SM.IsFlagNoApproximateGammaZ() && SM.ModelName() != "NPEpsilons")
                Gamma_Z_SM = myEW.Gamma_Z();
            else
                Gamma_Z_SM = SM.getEWSM()->getMyApproximateFormulae()->X_extended("GammaZ", SM.DeltaAlphaL5q());
            Gamma_Z_SM += SM.getDelGammaZ();/* Theoretical uncertainty */
            return myEW_BURGESS->GammaZ(Gamma_Z_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->GammaZ();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("sigmaHadron") == 0) {
        if (mode.compare("ABC") == 0)
            return ( myEW_ABC->sigma0_had(false)*GeVminus2_to_nb );
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->sigma0_had(true);
        else if (mode.compare("BURGESS") == 0) {
            double sigma_had_SM;
            if (SM.IsFlagNoApproximateSigmaH() && SM.ModelName() != "NPEpsilons")
                sigma_had_SM = myEW.sigma0_had();
            else
                sigma_had_SM = SM.getEWSM()->getMyApproximateFormulae()->X_extended("sigmaHadron", SM.DeltaAlphaL5q());
            return ( myEW_BURGESS->sigmaHadron(sigma_had_SM,
                       myEW.Gamma_Z(), myEW.Gamma_had(), myEW.Gamma_l(SM.ELECTRON))
                    *GeVminus2_to_nb );
        } else if (mode.compare("CHMN") == 0)
            return ( myEW_CHMN->sigma0_had()*GeVminus2_to_nb );
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("sin2thetaEff") == 0) {
        if (mode.compare("ABC") == 0)
            return myEW_ABC->sin2thetaEff(false);
        else if (mode.compare("ABC2") == 0)
            return myEW_ABC->sin2thetaEff(true);
        else if (mode.compare("BURGESS") == 0)
            return myEW_BURGESS->sin2thetaEff(myEW.sin2thetaEff(SM.ELECTRON));
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
            return myEW_BURGESS->PtauPol(myEW.A_l(SM.TAU));
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
            return myEW_BURGESS->Alepton(myEW.A_l(SM.ELECTRON));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Acharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0) 
            return myEW_ABC->A_q(SM.CHARM);
        else if (mode.compare("BURGESS") == 0) 
            return myEW_BURGESS->Acharm(myEW.A_q(SM.CHARM), myEW.A_l(SM.ELECTRON));
        else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->A_q(SM.CHARM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Abottom") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->A_b();
        else if (mode.compare("BURGESS") == 0) 
            return myEW_BURGESS->Abottom(myEW.A_q(SM.BOTTOM), myEW.A_l(SM.ELECTRON));
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
            double AFB_l_SM = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_l(SM.ELECTRON);
            return myEW_BURGESS->AFBlepton(AFB_l_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->AFB_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("AFBcharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->AFB_c();
        else if (mode.compare("BURGESS") == 0)  {
            double AFB_c_SM = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.CHARM);
            return myEW_BURGESS->AFBcharm(AFB_c_SM);
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->AFB_q(SM.CHARM);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("AFBbottom") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->AFB_b();
        else if (mode.compare("BURGESS") == 0) {
            double AFB_b_SM = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
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
            if (!SM.IsFlagNoApproximateRl() && SM.ModelName() != "NPEpsilons")
                R0_l_SM = SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_lepton", SM.DeltaAlphaL5q());
            else
                R0_l_SM = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
            return myEW_BURGESS->Rlepton(R0_l_SM, myEW.Gamma_had(), myEW.Gamma_l(SM.ELECTRON));
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->R_l(SM.ELECTRON);
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else if (obsname.compare("Rcharm") == 0) {
        if (mode.compare("ABC") == 0 || mode.compare("ABC2") == 0)
            return myEW_ABC->R_c();
        else if (mode.compare("BURGESS") == 0) {
            double R0_c_SM;
            if (!SM.IsFlagNoApproximateRc() && SM.ModelName() != "NPEpsilons")
                R0_c_SM = SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_charm", SM.DeltaAlphaL5q());
            else
                R0_c_SM = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
            return myEW_BURGESS->Rcharm(R0_c_SM, myEW.Gamma_had());
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
            if (!SM.IsFlagNoApproximateRb() && SM.ModelName() != "NPEpsilons")
                R0_b_SM = SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_bottom", SM.DeltaAlphaL5q());
            else
                R0_b_SM = myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();
            return myEW_BURGESS->Rbottom(R0_b_SM, myEW.Gamma_had(), myEW.Gamma_q(SM.BOTTOM));
        } else if (mode.compare("CHMN") == 0)
            return myEW_CHMN->R_b();
        else
            throw std::runtime_error("EW_TEST::computeThValue(): Undefined");
    } else
        throw std::runtime_error("EW_TEST::computeThValue(): Wrong obsname " + obsname);
}

