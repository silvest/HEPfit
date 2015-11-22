/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <stdexcept>
#include "EWSMTwoFermionsLEP2_Hollik.h"

EWSMTwoFermionsLEP2_Hollik::EWSMTwoFermionsLEP2_Hollik(const StandardModel& SM_i)
: SM(SM_i), myOneLoopEW_HV(SM_i), PV(true)
{
    bUseHollik = false;
    //bUseHollik = true; // for test (use the self-energies in Hollik's paper)
}


//////////////////////////////////////////////////////////////////////// 

double EWSMTwoFermionsLEP2_Hollik::sigma_l(const StandardModel::lepton l,
        const double s, const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK, const bool bQED) const
{
    double Ncf = 1.0;
    double mf = myOneLoopEW_HV.ml(l);
    double betaf = sqrt(1.0 - 4.0 * mf * mf / s);

    return ( 4.0 * M_PI * SM.getAle() * SM.getAle() / (3.0 * s) * Ncf * betaf
            * (G1_l(l, s, Mw, GammaZ, bDP, bWEAK, bQED)
            + 2.0 * mf * mf / s * G2_l(l, s, Mw, bDP)));
}

double EWSMTwoFermionsLEP2_Hollik::sigma_q(const QCD::quark q,
        const double s, const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK, const bool bQED) const
{
    double Ncf = 3.0;
    double mf = myOneLoopEW_HV.mq(q, sqrt(s));
    double betaf = sqrt(1.0 - 4.0 * mf * mf / s);

    return ( 4.0 * M_PI * SM.getAle() * SM.getAle() / (3.0 * s) * Ncf * betaf
            * (G1_q(q, s, Mw, GammaZ, bDP, bWEAK, bQED)
            + 2.0 * mf * mf / s * G2_q(q, s, Mw, bDP)));
}

double EWSMTwoFermionsLEP2_Hollik::AFB_l(const StandardModel::lepton l,
        const double s, const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK, const bool bQED) const
{
    double mf = myOneLoopEW_HV.ml(l);
    double betaf = sqrt(1.0 - 4.0 * mf * mf / s);

    return ( 3.0 / 4.0 * betaf * G3_l(l, s, Mw, GammaZ, bDP, bWEAK, bQED)
            / (G1_l(l, s, Mw, GammaZ, bDP, bWEAK, bQED)
            + 2.0 * mf * mf / s * G2_l(l, s, Mw, bDP)));
}

double EWSMTwoFermionsLEP2_Hollik::AFB_q(const QCD::quark q,
        const double s, const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK, const bool bQED) const
{
    double mf = myOneLoopEW_HV.mq(q, sqrt(s));
    double betaf = sqrt(1.0 - 4.0 * mf * mf / s);

    return ( 3.0 / 4.0 * betaf * G3_q(q, s, Mw, GammaZ, bDP, bWEAK, bQED)
            / (G1_q(q, s, Mw, GammaZ, bDP, bWEAK, bQED)
            + 2.0 * mf * mf / s * G2_q(q, s, Mw, bDP)));
}


////////////////////////////////////////////////////////////////////////

double EWSMTwoFermionsLEP2_Hollik::sigma_l_old(const StandardModel::lepton l, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bQED) const
{
    double mf = SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double I3f = SM.getLeptons(l).getIsospin();

    return sigma_f_old(s, Mw, GammaZ, mf, Qf, I3f, 1.0, bDP, bQED);
}

double EWSMTwoFermionsLEP2_Hollik::sigma_q_old(const QCD::quark q, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bQED) const
{
    double mf = myOneLoopEW_HV.mq(q, sqrt(s));
    double Qf = SM.getQuarks(q).getCharge();
    double I3f = SM.getQuarks(q).getIsospin();

    return sigma_f_old(s, Mw, GammaZ, mf, Qf, I3f, 3.0, bDP, bQED);
}

double EWSMTwoFermionsLEP2_Hollik::sigma_f_old(const double s, const double Mw, const double GammaZ,
        const double mf, const double Qf, const double I3f,
        const double Ncf,
        const bool bDP, const bool bQED) const
{
    double betaf = sqrt(1.0 - 4.0 * mf * mf / s);
    double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
    double I3e = SM.getLeptons(SM.ELECTRON).getIsospin();
    double Mz = SM.getMz();
    double cW2 = Mw * Mw / Mz / Mz, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);
    double ve = -(I3e - 2.0 * Qe * sW2) / (2.0 * sW * cW);
    double ae = -I3e / (2.0 * sW * cW);
    double vf = -(I3f - 2.0 * Qf * sW2) / (2.0 * sW * cW);
    double af = -I3f / (2.0 * sW * cW);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf; //, betaf2 = betaf*betaf;
    double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;

    double mu = sqrt(s); // The renormalized self-energies are independent of the scale.  
    gslpp::complex chiG = chi_gamma(mu, s, Mw, bDP);
    gslpp::complex chiZ = chi_Z(mu, s, Mw, bDP);
    gslpp::complex chiGZ = chi_gammaZ(mu, s, Mw, bDP);

    double G1 = Qe2 * Qf2 * chiG.abs2()
            + 2.0 * ve * vf * Qe * Qf * (chiZ * chiG.conjugate()).real()
            //+ (ve2 + ae2)*(vf2 + betaf2*af2)*chiZ.abs2() // with betaf2
            + (ve2 + ae2)*(vf2 + af2) * chiZ.abs2() // without betaf2
            + (Qe2 * (vf2 + af2) + 2.0 * ve * vf * Qe * Qf + (ve2 + ae2) * Qf2) * chiGZ.abs2()
            + 2.0 * (vf * Qe2 * Qf + ve * Qe * Qf2)*(chiG * chiGZ.conjugate()).real()
            + 2.0 * (ve * (vf2 + af2) * Qe + (ve2 + ae2) * vf * Qf)*(chiZ * chiGZ.conjugate()).real();
    double G2 = Qe2 * Qf2 * chiG.abs2()
            + 2.0 * ve * vf * Qe * Qf * (chiZ * chiG.conjugate()).real()
            + (ve2 + ae2) * vf2 * chiZ.abs2();

    // QED corrections
    if (bQED) {
        G1 += Qf * Qf * C11V(s, mf, Qf) * chiG.abs2()
                + 2.0 * Qe * Qf * ((ve * vf * C12V(s, GammaZ, mf, Qf) + ae * af * C12A(s, mf, Qf))
                * chiG * chiGZ.conjugate()).real()
                + ((ve2 + ae2)*(vf2 + af2) * C22V(s, GammaZ, mf, Qf)
                + 4.0 * ve * ae * vf * af * C22A(s, mf, Qf)) * chiZ.abs2();
    }

    return ( 4.0 * M_PI * SM.getAle() * SM.getAle() / (3.0 * s) * Ncf * betaf * (G1 + 2.0 * mf * mf / s * G2));
}


//////////////////////////////////////////////////////////////////////// 

void EWSMTwoFermionsLEP2_Hollik::TEST(const double s, const double Mw) const
{
    //-----------------------------------
    // Test for renormalization scale dependence of the renormalized self-energies
    std::cout << "TEST1 (mu=Mz)   : "
            << Sigma_hat_ZZ(SM.getMz(), s, Mw) << " "
            << Sigma_hat_gZ(SM.getMz(), s, Mw) << " "
            << Sigma_hat_gg(SM.getMz(), s, Mw) << std::endl;
    std::cout << "TEST2 (mu=2*Mz) : "
            << Sigma_hat_ZZ(2.0 * SM.getMz(), s, Mw) << " "
            << Sigma_hat_gZ(2.0 * SM.getMz(), s, Mw) << " "
            << Sigma_hat_gg(2.0 * SM.getMz(), s, Mw) << std::endl;
    //-----------------------------------
}


//////////////////////////////////////////////////////////////////////// 
// Renormalized self-energies

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Sigma_hat_ZZ(const double mu, const double s,
        const double Mw) const
{
    double Mw2 = Mw*Mw;
    double Mz = SM.getMz(), Mz2 = Mz*Mz;
    double cW2 = Mw2 / Mz2, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);

    // Bosonic contributions to self-energies
    gslpp::complex Sigma_WW_Mw2, Sigma_ZZ_s, Sigma_ZZ_Mz2, Sigma_Zg_0, Pi_gg_0;
    if (!bUseHollik) {
        Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos(mu, Mw2, Mw);
        Sigma_ZZ_s = myOneLoopEW_HV.SigmaZZ_bos(mu, s, Mw);
        Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos(mu, Mz2, Mw);
        Sigma_Zg_0 = myOneLoopEW_HV.SigmaZgamma_bos(mu, 0.0, Mw);
        Pi_gg_0 = myOneLoopEW_HV.PiGammaGamma_bos(mu, 0.0, Mw);
    } else {
        //-- TEST (use the self-energies in Hollik's paper) --
        Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos_Hollik(mu, Mw2, Mw);
        Sigma_ZZ_s = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, s, Mw);
        Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, Mz2, Mw);
        Sigma_Zg_0 = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, 0.0, Mw);
        Pi_gg_0 = myOneLoopEW_HV.PiGammaGamma_bos_Hollik(mu, 0.0, Mw);
    }

    // Fermionic contributions to self-energies
    double muForMq = sqrt(s); // renormalization scale for the running quark mass
    Sigma_WW_Mw2 += myOneLoopEW_HV.SigmaWW_fer(mu, muForMq, Mw2);
    Sigma_ZZ_s += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, s, Mw);
    Sigma_ZZ_Mz2 += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, Mz2, Mw);
    Sigma_Zg_0 += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, 0.0, Mw);
    Pi_gg_0 += myOneLoopEW_HV.PiGammaGamma_fer(mu, muForMq, 0.0);

    // Refactoring
    Sigma_WW_Mw2 *= SM.getAle() / 4.0 / M_PI / sW2;
    Sigma_ZZ_s *= SM.getAle() / 4.0 / M_PI / sW2 / cW2;
    Sigma_ZZ_Mz2 *= SM.getAle() / 4.0 / M_PI / sW2 / cW2;
    Sigma_Zg_0 *= SM.getAle() / 4.0 / M_PI / sW / cW;
    Pi_gg_0 *= SM.getAle() / 4.0 / M_PI;

    // Counter terms for the mass renormalization
    double deltaMw2 = Sigma_WW_Mw2.real();
    double deltaMz2 = Sigma_ZZ_Mz2.real();

    // Counter terms for the wave-function renormalization
    gslpp::complex deltaZz = -Pi_gg_0 + (cW2 - sW2) / sW2 * (deltaMz2 / Mz2 - deltaMw2 / Mw2
            + 2.0 * sW / cW * Sigma_Zg_0 / Mz2);

    return ( Sigma_ZZ_s + (s - Mz2) * deltaZz - deltaMz2);
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Sigma_hat_gZ(const double mu, const double s,
        const double Mw) const
{
    double Mw2 = Mw*Mw;
    double Mz = SM.getMz(), Mz2 = Mz*Mz;
    double cW2 = Mw2 / Mz2, cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);

    // Bosonic contributions to self-energies
    gslpp::complex Sigma_WW_Mw2, Sigma_ZZ_Mz2, Sigma_Zg_s, Sigma_Zg_0;
    if (!bUseHollik) {
        Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos(mu, Mw2, Mw);
        Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos(mu, Mz2, Mw);
        Sigma_Zg_s = myOneLoopEW_HV.SigmaZgamma_bos(mu, s, Mw);
        Sigma_Zg_0 = myOneLoopEW_HV.SigmaZgamma_bos(mu, 0.0, Mw);
    } else {
        //-- TEST (use the self-energies in Hollik's paper) --
        Sigma_WW_Mw2 = myOneLoopEW_HV.SigmaWW_bos_Hollik(mu, Mw2, Mw);
        Sigma_ZZ_Mz2 = myOneLoopEW_HV.SigmaZZ_bos_Hollik(mu, Mz2, Mw);
        Sigma_Zg_s = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, s, Mw);
        Sigma_Zg_0 = myOneLoopEW_HV.SigmaZgamma_bos_Hollik(mu, 0.0, Mw);
    }

    // Fermionic contributions to self-energies
    double muForMq = sqrt(s); // renormalization scale for the running quark mass
    Sigma_WW_Mw2 += myOneLoopEW_HV.SigmaWW_fer(mu, muForMq, Mw2);
    Sigma_ZZ_Mz2 += myOneLoopEW_HV.SigmaZZ_fer(mu, muForMq, Mz2, Mw);
    Sigma_Zg_s += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, s, Mw);
    Sigma_Zg_0 += myOneLoopEW_HV.SigmaZgamma_fer(mu, muForMq, 0.0, Mw);

    // Refactoring
    Sigma_WW_Mw2 *= SM.getAle() / 4.0 / M_PI / sW2;
    Sigma_ZZ_Mz2 *= SM.getAle() / 4.0 / M_PI / sW2 / cW2;
    Sigma_Zg_s *= SM.getAle() / 4.0 / M_PI / sW / cW;
    Sigma_Zg_0 *= SM.getAle() / 4.0 / M_PI / sW / cW;

    // Counter terms for the mass renormalization
    double deltaMw2 = Sigma_WW_Mw2.real();
    double deltaMz2 = Sigma_ZZ_Mz2.real();

    return ( Sigma_Zg_s - Sigma_Zg_0
            + s * (cW / sW * (deltaMz2 / Mz2 - deltaMw2 / Mw2) + 2.0 * Sigma_Zg_0 / Mz2));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Sigma_hat_gg(const double mu, const double s,
        const double Mw) const
{
    // Bosonic contributions to self-energies
    gslpp::complex Sigma_gg_s, Pi_gg_0;
    if (!bUseHollik) {
        Sigma_gg_s = myOneLoopEW_HV.SigmaGammaGamma_bos(mu, s, Mw);
        Pi_gg_0 = myOneLoopEW_HV.PiGammaGamma_bos(mu, 0.0, Mw);
    } else {
        //-- TEST (use the self-energies in Hollik's paper) --
        Sigma_gg_s = myOneLoopEW_HV.SigmaGammaGamma_bos_Hollik(mu, s, Mw);
        Pi_gg_0 = myOneLoopEW_HV.PiGammaGamma_bos_Hollik(mu, 0.0, Mw);
    }

    // Fermionic contributions to self-energies
    double muForMq = sqrt(s); // renormalization scale for the running quark mass
    Sigma_gg_s += myOneLoopEW_HV.SigmaGammaGamma_fer(mu, muForMq, s);
    Pi_gg_0 += myOneLoopEW_HV.PiGammaGamma_fer(mu, muForMq, 0.0);

    // Refactoring
    Sigma_gg_s *= SM.getAle() / 4.0 / M_PI;
    Pi_gg_0 *= SM.getAle() / 4.0 / M_PI;

    return ( Sigma_gg_s - s * Pi_gg_0);
}


//////////////////////////////////////////////////////////////////////// 
// Tree-level couplings for neutral-current interactions

double EWSMTwoFermionsLEP2_Hollik::vl(const StandardModel::lepton l, const double Mw) const
{
    double cW = Mw / SM.getMz(), sW2 = 1.0 - cW*cW, sW = sqrt(sW2);
    return ( -(SM.getLeptons(l).getIsospin()
            - 2.0 * SM.getLeptons(l).getCharge() * sW2) / (2.0 * sW * cW));
}

double EWSMTwoFermionsLEP2_Hollik::vq(const QCD::quark q, const double Mw) const
{
    double cW = Mw / SM.getMz(), sW2 = 1.0 - cW*cW, sW = sqrt(sW2);
    return ( -(SM.getQuarks(q).getIsospin()
            - 2.0 * SM.getQuarks(q).getCharge() * sW2) / (2.0 * sW * cW));
}

double EWSMTwoFermionsLEP2_Hollik::al(const StandardModel::lepton l, const double Mw) const
{
    double cW = Mw / SM.getMz(), sW = sqrt(1.0 - cW * cW);
    return ( -SM.getLeptons(l).getIsospin() / (2.0 * sW * cW));
}

double EWSMTwoFermionsLEP2_Hollik::aq(const QCD::quark q, const double Mw) const
{
    double cW = Mw / SM.getMz(), sW = sqrt(1.0 - cW * cW);
    return ( -SM.getQuarks(q).getIsospin() / (2.0 * sW * cW));
}


//////////////////////////////////////////////////////////////////////// 
// Dressed gauge-boson propagators

gslpp::complex EWSMTwoFermionsLEP2_Hollik::chi_Z(const double mu, const double s,
        const double Mw, const bool bDP) const
{
    gslpp::complex chi;
    if (bDP) {
        double Mz = SM.getMz();
        chi = s / (s - Mz * Mz + Sigma_hat_ZZ(mu, s, Mw));
    } else {
        double Mw2 = Mw*Mw;
        double Mz = SM.getMz(), Mz2 = Mz*Mz;
        double cW2 = Mw2 / Mz2;
        double sW2 = 1.0 - cW2;
        double GammaZ_tree = 7.0 * SM.getAle() * Mz / 16.0 / sW2 / cW2;
        gslpp::complex denom = gslpp::complex(s - Mz*Mz, Mz*GammaZ_tree, false);
        chi = s / denom;
    }

    return chi;
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::chi_gamma(const double mu, const double s,
        const double Mw, const bool bDP) const
{
    gslpp::complex chi;
    if (bDP)
        chi = s / (s + Sigma_hat_gg(mu, s, Mw));
    else
        chi = 1.0;

    return chi;
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::chi_gammaZ(const double mu, const double s,
        const double Mw, const bool bDP) const
{
    gslpp::complex chi;
    if (bDP) {
        // O(alpha) approximation
        chi = Sigma_hat_gZ(mu, s, Mw) / s * chi_Z(mu, s, Mw, bDP);
    } else
        chi = 0.0;

    return chi;
}


//////////////////////////////////////////////////////////////////////// 
// Renormalized vertex form factors for the Z-f-f vertex (non-QED part)  

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FVZ_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    double cW = Mw / SM.getMz(), cW2 = cW*cW;
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);
    double v_l = vl(l, Mw), a_l = al(l, Mw);

    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return ( -SM.getAle() / (16.0 * M_PI * sW * cW)
                    * (Lambda2(s, SM.getMz()) / (4.0 * cW2 * sW2)
                    + Lambda2(s, Mw)*(2.0 * sW2 - 1.0) / (2.0 * sW2)
                    + Lambda3(s, Mw)*3.0 * cW2 / sW2));
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( SM.getAle() / (4.0 * M_PI)
                    * (v_l * (v_l * v_l + 3.0 * a_l * a_l) * Lambda2(s, SM.getMz()) + FL_l(l, s, Mw)));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::FL_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FAZ_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    double cW = Mw / SM.getMz(), cW2 = cW*cW;
    double sW2 = 1.0 - cW2, sW = sqrt(sW2);
    double v_l = vl(l, Mw), a_l = al(l, Mw);

    switch (l) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            return ( -SM.getAle() / (16.0 * M_PI * sW * cW)
                    * (Lambda2(s, SM.getMz()) / (4.0 * cW2 * sW2)
                    + Lambda2(s, Mw)*(2.0 * sW2 - 1.0) / (2.0 * sW2)
                    + Lambda3(s, Mw)*3.0 * cW2 / sW2));
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            return ( SM.getAle() / (4.0 * M_PI)
                    * (a_l * (3.0 * v_l * v_l + a_l * a_l) * Lambda2(s, SM.getMz()) + FL_l(l, s, Mw)));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::FL_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FVZ_q(const QCD::quark q,
        const double s, const double Mw) const
{
    double v_q = vq(q, Mw), a_q = aq(q, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (v_q * (v_q * v_q + 3.0 * a_q * a_q) * Lambda2(s, SM.getMz()) + FL_q(q, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FAZ_q(const QCD::quark q,
        const double s, const double Mw) const
{
    double v_q = vq(q, Mw), a_q = aq(q, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (a_q * (3.0 * v_q * v_q + a_q * a_q) * Lambda2(s, SM.getMz()) + FL_q(q, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FL_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    if ((l == StandardModel::NEUTRINO_1) || (l == StandardModel::NEUTRINO_2)
            || (l == StandardModel::NEUTRINO_3))
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::FL_l()");

    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW, sW = sqrt(sW2), sW3 = sW2*sW;

    return ( -1.0 / (8.0 * sW3 * cW) * Lambda2(s, Mw)
            + 3.0 * cW / (4.0 * sW3) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FL_q(const QCD::quark q,
        const double s, const double Mw) const
{
    switch (q) {
        case QCD::UP:
        case QCD::CHARM:
            return FL_u(s, Mw);
        case QCD::DOWN:
        case QCD::STRANGE:
            return FL_d(s, Mw);
        case QCD::BOTTOM:
        {
            double cW = Mw / SM.getMz();
            double sW2 = 1.0 - cW*cW, sW = sqrt(sW2);
            return ( Fb(s, Mw) + Fc(s, Mw) + Fd(s, Mw)
                    + Fe(s, Mw) + Ff(s, Mw) + Fg(s, Mw)
                    - (2.0 / 3.0 * sW2 - 1.0) / (4.0 * sW * cW) * deltaZL_fin(s, Mw));
        }
        case QCD::TOP:
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::FL_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FL_u(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW, sW = sqrt(sW2), sW3 = sW2*sW;

    return ( (1.0 - 2.0 / 3.0 * sW2) / (8.0 * sW3 * cW) * Lambda2(s, Mw)
            - 3.0 * cW / (4.0 * sW3) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FL_d(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW, sW = sqrt(sW2), sW3 = sW2*sW;

    return ( -(1.0 - 4.0 / 3.0 * sW2) / (8.0 * sW3 * cW) * Lambda2(s, Mw)
            + 3.0 * cW / (4.0 * sW3) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Lambda2(const double s, const double M) const
{
    if (s <= 0.0)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::Lambda2()");

    double w = M * M / s;
    double real = -7.0 / 2.0 - 2.0 * w - (2.0 * w + 3.0) * log(w)
            + 2.0 * (1.0 + w)*(1.0 + w)
            *(log(w) * log((1.0 + w) / w) - Polylog.Li2(-1.0 / w).real());
    double imag = -M_PI * (3.0 + 2.0 * w - 2.0 * (1.0 + w)*(1.0 + w) * log((1.0 + w) / w));
    return gslpp::complex(real, imag, false);
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Lambda3(const double s, const double M) const
{
    if (s <= 0.0)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::Lambda3()");

    double w = M * M / s;
    if (s < 4.0 * M * M) {
        double sqrt_tmp = sqrt(4.0 * w - 1.0);
        double atan_tmp = atan(1.0 / sqrt_tmp);
        return ( 5.0 / 6.0 - 2.0 * w / 3.0 + 2.0 / 3.0 * (2.0 * w + 1.0) * sqrt_tmp * atan_tmp
                - 8.0 / 3.0 * w * (w + 2.0) * atan_tmp * atan_tmp);
    } else {
        double sqrt_tmp = sqrt(1.0 - 4.0 * w);
        gslpp::complex log_tmp;
        if (sqrt_tmp > 1.0)
            log_tmp = log((sqrt_tmp - 1.0) / (sqrt_tmp + 1.0));
        else
            log_tmp = gslpp::complex(log(-(sqrt_tmp - 1.0) / (sqrt_tmp + 1.0)), M_PI, false);
        return ( 5.0 / 6.0 - 2.0 * w / 3.0 - (2.0 * w + 1.0) / 3.0 * sqrt_tmp * log_tmp
                + 2.0 / 3.0 * w * (w + 2.0) * log_tmp * log_tmp);
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::deltaZL_fin(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mb = myOneLoopEW_HV.mq(QCD::BOTTOM, sqrt(s));
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( (2.0 + mt * mt / Mw / Mw) / (2.0 * sW2)
            * (B1bar_Hollik(mb*mb, mt, Mw)
            + mb * mb * B1barPrime_Hollik(mb*mb, mt, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Fb(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));
    double vt = vq(QCD::TOP, Mw);
    double at = aq(QCD::TOP, Mw);

    return ( (vt + at) / (4.0 * sW2)
            * (-3.0 / 2.0 + 2.0 * log(Mw / mt) + 4.0 * C2zero_Hollik(s, mt, Mw)
            - 2.0 * s * (C2plus_Hollik(s, mt, Mw) - C2minus_Hollik(s, mt, Mw))
            + 4.0 * s * C1plus_Hollik(s, mt, Mw) - 2.0 * s * C0_Hollik(s, mt, Mw))
            - (vt - at) / (4.0 * sW2)*2.0 * mt * mt * C0_Hollik(s, mt, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Fc(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW, sW = sqrt(sW2), sW3 = sW2*sW;
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( cW / (4.0 * sW3)
            * (-3.0 / 2.0 + 12.0 * C2zero_Hollik(s, Mw, mt)
            - 2.0 * s * (C2plus_Hollik(s, Mw, mt) - C2minus_Hollik(s, Mw, mt))
            + 4.0 * s * C1plus_Hollik(s, Mw, mt)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Fd(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));
    double vt = vq(QCD::TOP, Mw);
    double at = aq(QCD::TOP, Mw);

    return ( (vt - at) / (4.0 * sW2) * mt * mt / Mw / Mw
            * (-3.0 / 4.0 + log(Mw / mt) + 2.0 * C2zero_Hollik(s, mt, Mw)
            - s * (C2plus_Hollik(s, mt, Mw) - C2minus_Hollik(s, mt, Mw)))
            - (vt + at) / (4.0 * sW2) * mt * mt / Mw / Mw * mt * mt * C0_Hollik(s, mt, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Fe(const double s, const double Mw) const
{
    double cW2 = Mw * Mw / SM.getMz() / SM.getMz(), cW = sqrt(cW2);
    double sW2 = 1.0 - cW2, sW = sqrt(sW2), sW3 = sW2*sW;
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( -(sW2 - sW2) / (8.0 * sW3 * cW) * mt * mt / Mw / Mw
            * (-1.0 / 4.0 - 2.0 * C2zero_Hollik(s, Mw, mt)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Ff(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW, sW = sqrt(sW2);
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( mt * mt / (4.0 * sW * cW) * C0_Hollik(s, Mw, mt));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Fg(const double s, const double Mw) const
{
    return ( Ff(s, Mw));
}


//////////////////////////////////////////////////////////////////////// 
// Renormalized vertex form factors for the gamma-f-f vertex (non-QED part)

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FVgamma_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    double Q_l = SM.getLeptons(l).getCharge();
    double v_l = vl(l, Mw), a_l = al(l, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (Q_l * (v_l * v_l + a_l * a_l) * Lambda2(s, SM.getMz()) + GL_l(l, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FAgamma_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    double Q_l = SM.getLeptons(l).getCharge();
    double v_l = vl(l, Mw), a_l = al(l, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (2.0 * Q_l * v_l * a_l * Lambda2(s, SM.getMz()) + GL_l(l, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FVgamma_q(const QCD::quark q,
        const double s, const double Mw) const
{
    double Q_q = SM.getQuarks(q).getCharge();
    double v_q = vq(q, Mw), a_q = aq(q, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (Q_q * (v_q * v_q + a_q * a_q) * Lambda2(s, SM.getMz()) + GL_q(q, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::FAgamma_q(const QCD::quark q,
        const double s, const double Mw) const
{
    double Q_q = SM.getQuarks(q).getCharge();
    double v_q = vq(q, Mw), a_q = aq(q, Mw);

    return ( SM.getAle() / (4.0 * M_PI)
            * (2.0 * Q_q * v_q * a_q * Lambda2(s, SM.getMz()) + GL_q(q, s, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::GL_l(const StandardModel::lepton l,
        const double s, const double Mw) const
{
    if ((l == StandardModel::NEUTRINO_1) || (l == StandardModel::NEUTRINO_2)
            || (l == StandardModel::NEUTRINO_3))
        return gslpp::complex(0.0, 0.0, false);

    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW;

    return ( -3.0 / (4.0 * sW2) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::GL_q(const QCD::quark q,
        const double s, const double Mw) const
{
    switch (q) {
        case QCD::UP:
        case QCD::CHARM:
            return GL_u(s, Mw);
        case QCD::DOWN:
        case QCD::STRANGE:
            return GL_d(s, Mw);
        case QCD::BOTTOM:
            return ( Gb(s, Mw) + Gc(s, Mw) + Gd(s, Mw)
                    + Ge(s, Mw) + Gf(s, Mw) + Gg(s, Mw)
                    - 1.0 / 6.0 * deltaZL_fin(s, Mw));
        case QCD::TOP:
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::GL_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::GL_u(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW;

    return ( -Lambda2(s, Mw) / (12.0 * sW2) + 3.0 / (4.0 * sW2) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::GL_d(const double s, const double Mw) const
{
    double cW = Mw / SM.getMz();
    double sW2 = 1.0 - cW*cW;

    return ( Lambda2(s, Mw) / (6.0 * sW2) - 3.0 / (4.0 * sW2) * Lambda3(s, Mw));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Gb(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( 1 / (6.0 * sW2)
            * (-3.0 / 2.0 + 2.0 * log(Mw / mt) + 4.0 * C2zero_Hollik(s, mt, Mw)
            - 2.0 * s * (C2plus_Hollik(s, mt, Mw) - C2minus_Hollik(s, mt, Mw))
            + 4.0 * s * C1plus_Hollik(s, mt, Mw) - 2.0 * s * C0_Hollik(s, mt, Mw)
            - 2.0 * mt * mt * C0_Hollik(s, mt, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Gc(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( -1.0 / (4.0 * sW2)
            * (-3.0 / 2.0 + 12.0 * C2zero_Hollik(s, Mw, mt)
            - 2.0 * s * (C2plus_Hollik(s, Mw, mt) - C2minus_Hollik(s, Mw, mt))
            + 4.0 * s * C1plus_Hollik(s, Mw, mt)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Gd(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( 1.0 / (6.0 * sW2) * mt * mt / Mw / Mw
            * (-3.0 / 4.0 + log(Mw / mt) + 2.0 * C2zero_Hollik(s, mt, Mw)
            - s * (C2plus_Hollik(s, mt, Mw) - C2minus_Hollik(s, mt, Mw))
            - mt * mt * C0_Hollik(s, mt, Mw)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Ge(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( -1.0 / (4.0 * sW2) * mt * mt / Mw / Mw
            * (-1.0 / 4.0 + 2.0 * C2zero_Hollik(s, Mw, mt)));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Gf(const double s, const double Mw) const
{
    double sW2 = 1.0 - Mw * Mw / SM.getMz() / SM.getMz();
    double mt = myOneLoopEW_HV.mq(QCD::TOP, sqrt(s));

    return ( mt * mt / (4.0 * sW2) * C0_Hollik(s, Mw, mt));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::Gg(const double s, const double Mw) const
{
    return ( Gf(s, Mw));
}


//////////////////////////////////////////////////////////////////////// 
// Born + dressed propagators + non-QED vertex corrections 

gslpp::complex EWSMTwoFermionsLEP2_Hollik::V_e(const int j, const double s,
        const double Mw, const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 3:
        case 6:
            return SM.getLeptons(SM.ELECTRON).getCharge();
        case 2:
        case 4:
        case 8:
            return vl(SM.ELECTRON, Mw);
        case 5:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVgamma_l(SM.ELECTRON, s, Mw);
        case 7:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVZ_l(SM.ELECTRON, s, Mw);
        case 9:
            return ( vl(SM.ELECTRON, Mw) * vl(SM.ELECTRON, Mw)
                    + al(SM.ELECTRON, Mw) * al(SM.ELECTRON, Mw));
        case 10:
            return ( 2.0 * vl(SM.ELECTRON, Mw) * al(SM.ELECTRON, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::V_e()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::A_e(const int j, const double s,
        const double Mw, const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 3:
        case 6:
            return 0.0;
        case 2:
        case 4:
        case 8:
            return al(SM.ELECTRON, Mw);
        case 5:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAgamma_l(SM.ELECTRON, s, Mw);
        case 7:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAZ_l(SM.ELECTRON, s, Mw);
        case 9:
            return ( 2.0 * vl(SM.ELECTRON, Mw) * al(SM.ELECTRON, Mw));
        case 10:
            return ( vl(SM.ELECTRON, Mw) * vl(SM.ELECTRON, Mw)
                    + al(SM.ELECTRON, Mw) * al(SM.ELECTRON, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::A_e()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::V_l(const int j, const StandardModel::lepton l,
        const double s, const double Mw,
        const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 4:
        case 5:
            return SM.getLeptons(l).getCharge();
        case 2:
        case 3:
        case 7:
            return vl(l, Mw);
        case 6:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVgamma_l(l, s, Mw);
        case 8:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVZ_l(l, s, Mw);
        case 9:
            return ( vl(l, Mw) * vl(l, Mw) + al(l, Mw) * al(l, Mw));
        case 10:
            return ( 2.0 * vl(l, Mw) * al(l, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::V_l()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::V_q(const int j, const QCD::quark q,
        const double s, const double Mw,
        const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 4:
        case 5:
            return SM.getQuarks(q).getCharge();
        case 2:
        case 3:
        case 7:
            return vq(q, Mw);
        case 6:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVgamma_q(q, s, Mw);
        case 8:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FVZ_q(q, s, Mw);
        case 9:
            return ( vq(q, Mw) * vq(q, Mw) + aq(q, Mw) * aq(q, Mw));
        case 10:
            return ( 2.0 * vq(q, Mw) * aq(q, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::V_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::A_l(const int j, const StandardModel::lepton l,
        const double s, const double Mw,
        const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 4:
        case 5:
            return 0.0;
        case 2:
        case 3:
        case 7:
            return al(l, Mw);
        case 6:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAgamma_l(l, s, Mw);
        case 8:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAZ_l(l, s, Mw);
        case 9:
            return ( 2.0 * vl(l, Mw) * al(l, Mw));
        case 10:
            return ( vl(l, Mw) * vl(l, Mw) + al(l, Mw) * al(l, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::A_l()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::A_q(const int j, const QCD::quark q,
        const double s, const double Mw,
        const bool bWEAK) const
{
    switch (j) {
        case 1:
        case 4:
        case 5:
            return 0.0;
        case 2:
        case 3:
        case 7:
            return aq(q, Mw);
        case 6:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAgamma_q(q, s, Mw);
        case 8:
            if (!bWEAK) return gslpp::complex(0.0, 0.0, false);
            return FAZ_q(q, s, Mw);
        case 9:
            return ( 2.0 * vq(q, Mw) * aq(q, Mw));
        case 10:
            return ( vq(q, Mw) * vq(q, Mw) + aq(q, Mw) * aq(q, Mw));
        case 11:
            return ( 1.0 / 4.0 / (1.0 - Mw * Mw / SM.getMz() / SM.getMz()));
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::A_q()");
    }
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::chi(const int j, const double s,
        const double Mw, const bool bDP) const
{
    double mu = Mw; // The result is independent of the renormalization scale.  

    switch (j) {
        case 1:
            return chi_gamma(mu, s, Mw, bDP);
        case 2:
            return chi_Z(mu, s, Mw, bDP);
        case 3:
            return chi_gammaZ(mu, s, Mw, bDP);
        case 4:
            return chi_gammaZ(mu, s, Mw, bDP);
        case 5:
            return chi_gamma(mu, s, Mw, bDP);
        case 6:
            return chi_gamma(mu, s, Mw, bDP);
        case 7:
            return chi_Z(mu, s, Mw, bDP);
        case 8:
            return chi_Z(mu, s, Mw, bDP);
        case 9:
            // box contribution. add codes!!!
        case 10:
            // box contribution. add codes!!!
        case 11:
            // box contribution. add codes!!!
        default:
            throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::chi()");
    }
}

double EWSMTwoFermionsLEP2_Hollik::G1_l(const StandardModel::lepton l, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK,
        const bool bQED) const
{
    int j, k;
    double G1 = 0.0;
    for (j = 1; j <= 8; j++) {
        for (k = 1; k <= 8; k++) {
            G1 += ((V_e(j, s, Mw, bWEAK) * V_e(k, s, Mw, bWEAK).conjugate()
                    + A_e(j, s, Mw, bWEAK) * A_e(k, s, Mw, bWEAK).conjugate())
                    * (V_l(j, l, s, Mw, bWEAK) * V_l(k, l, s, Mw, bWEAK).conjugate()
                    + A_l(j, l, s, Mw, bWEAK) * A_l(k, l, s, Mw, bWEAK).conjugate())
                    * chi(j, s, Mw, bDP) * chi(k, s, Mw, bDP).conjugate()).real();
        }
    }

    // QED corrections
    if (bQED) {
        double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = myOneLoopEW_HV.ml(l);
        ;
        double ve = vl(SM.ELECTRON, Mw), ae = al(SM.ELECTRON, Mw);
        double vf = vl(l, Mw), af = al(l, Mw);
        double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;
        gslpp::complex chi1 = chi(1, s, Mw, bDP);
        gslpp::complex chi2 = chi(2, s, Mw, bDP);

        G1 += Qf * Qf * C11V(s, mf, Qf) * chi1.abs2()
                + 2.0 * Qe * Qf * ((ve * vf * C12V(s, GammaZ, mf, Qf) + ae * af * C12A(s, mf, Qf))
                * chi1 * chi2.conjugate()).real()
                + ((ve2 + ae2)*(vf2 + af2) * C22V(s, GammaZ, mf, Qf)
                + 4.0 * ve * ae * vf * af * C22A(s, mf, Qf)) * chi2.abs2();
    }

    return G1;
}

double EWSMTwoFermionsLEP2_Hollik::G1_q(const QCD::quark q, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK,
        const bool bQED) const
{
    int j, k;
    double G1 = 0.0;
    for (j = 1; j <= 8; j++) {
        for (k = 1; k <= 8; k++) {
            G1 += ((V_e(j, s, Mw, bWEAK) * V_e(k, s, Mw, bWEAK).conjugate()
                    + A_e(j, s, Mw, bWEAK) * A_e(k, s, Mw, bWEAK).conjugate())
                    * (V_q(j, q, s, Mw, bWEAK) * V_q(k, q, s, Mw, bWEAK).conjugate()
                    + A_q(j, q, s, Mw, bWEAK) * A_q(k, q, s, Mw, bWEAK).conjugate())
                    * chi(j, s, Mw, bDP) * chi(k, s, Mw, bDP).conjugate()).real();
        }
    }

    // QED corrections
    if (bQED) {
        double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = myOneLoopEW_HV.mq(q, sqrt(s));
        ;
        double ve = vl(SM.ELECTRON, Mw), ae = al(SM.ELECTRON, Mw);
        double vf = vq(q, Mw), af = aq(q, Mw);
        double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;
        gslpp::complex chi1 = chi(1, s, Mw, bDP);
        gslpp::complex chi2 = chi(2, s, Mw, bDP);

        G1 += Qf * Qf * C11V(s, mf, Qf) * chi1.abs2()
                + 2.0 * Qe * Qf * ((ve * vf * C12V(s, GammaZ, mf, Qf) + ae * af * C12A(s, mf, Qf))
                * chi1 * chi2.conjugate()).real()
                + ((ve2 + ae2)*(vf2 + af2) * C22V(s, GammaZ, mf, Qf)
                + 4.0 * ve * ae * vf * af * C22A(s, mf, Qf)) * chi2.abs2();
    }

    return G1;
}

double EWSMTwoFermionsLEP2_Hollik::G2_l(const StandardModel::lepton l, const double s,
        const double Mw, const bool bDP) const

{
    double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
    double Qf = SM.getLeptons(l).getCharge();
    double ve = vl(SM.ELECTRON, Mw);
    double ae = al(SM.ELECTRON, Mw);
    double vf = vl(l, Mw);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf;
    gslpp::complex chi1 = chi(1, s, Mw, bDP);
    gslpp::complex chi2 = chi(2, s, Mw, bDP);

    return ( Qe2 * Qf2 * chi1.abs2()
            + 2.0 * ve * vf * Qe * Qf * (chi2 * chi1.conjugate()).real()
            + (ve2 + ae2) * vf2 * chi2.abs2());
}

double EWSMTwoFermionsLEP2_Hollik::G2_q(const QCD::quark q, const double s,
        const double Mw, const bool bDP) const
{
    double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
    double Qf = SM.getQuarks(q).getCharge();
    double ve = vl(SM.ELECTRON, Mw);
    double ae = al(SM.ELECTRON, Mw);
    double vf = vq(q, Mw);
    double Qe2 = Qe*Qe, Qf2 = Qf*Qf;
    double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf;
    gslpp::complex chi1 = chi(1, s, Mw, bDP);
    gslpp::complex chi2 = chi(2, s, Mw, bDP);

    return ( Qe2 * Qf2 * chi1.abs2()
            + 2.0 * ve * vf * Qe * Qf * (chi2 * chi1.conjugate()).real()
            + (ve2 + ae2) * vf2 * chi2.abs2());
}

double EWSMTwoFermionsLEP2_Hollik::G3_l(const StandardModel::lepton l, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK,
        const bool bQED) const
{
    int j, k;
    double G3 = 0.0;
    for (j = 1; j <= 8; j++) {
        for (k = 1; k <= 8; k++) {
            G3 += ((V_e(j, s, Mw, bWEAK) * A_e(k, s, Mw, bWEAK).conjugate()
                    + A_e(j, s, Mw, bWEAK) * V_e(k, s, Mw, bWEAK).conjugate())
                    * (V_l(j, l, s, Mw, bWEAK) * A_l(k, l, s, Mw, bWEAK).conjugate()
                    + A_l(j, l, s, Mw, bWEAK) * V_l(k, l, s, Mw, bWEAK).conjugate())
                    * chi(j, s, Mw, bDP) * chi(k, s, Mw, bDP).conjugate()).real();
        }
    }

    // QED corrections
    if (bQED) {
        double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = myOneLoopEW_HV.ml(l);
        ;
        double ve = vl(SM.ELECTRON, Mw), ae = al(SM.ELECTRON, Mw);
        double vf = vl(l, Mw), af = al(l, Mw);
        double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;
        gslpp::complex chi1 = chi(1, s, Mw, bDP);
        gslpp::complex chi2 = chi(2, s, Mw, bDP);

        G3 += Qf * Qf * C11A(s, mf, Qf) * chi1.abs2()
                + 2.0 * Qe * Qf * ((ae * af * C12V(s, GammaZ, mf, Qf) + ve * vf * C12A(s, mf, Qf))
                * chi1 * chi2.conjugate()).real()
                + (4.0 * ve * ae * vf * af * C22V(s, GammaZ, mf, Qf)
                + (ve2 + ae2)*(vf2 + af2) * C22A(s, mf, Qf)) * chi2.abs2();
    }

    return G3;
}

double EWSMTwoFermionsLEP2_Hollik::G3_q(const QCD::quark q, const double s,
        const double Mw, const double GammaZ,
        const bool bDP, const bool bWEAK,
        const bool bQED) const
{
    int j, k;
    double G3 = 0.0;
    for (j = 1; j <= 8; j++) {
        for (k = 1; k <= 8; k++) {
            G3 += ((V_e(j, s, Mw, bWEAK) * A_e(k, s, Mw, bWEAK).conjugate()
                    + A_e(j, s, Mw, bWEAK) * V_e(k, s, Mw, bWEAK).conjugate())
                    * (V_q(j, q, s, Mw, bWEAK) * A_q(k, q, s, Mw, bWEAK).conjugate()
                    + A_q(j, q, s, Mw, bWEAK) * V_q(k, q, s, Mw, bWEAK).conjugate())
                    * chi(j, s, Mw, bDP) * chi(k, s, Mw, bDP).conjugate()).real();
        }
    }

    // QED corrections
    if (bQED) {
        double Qe = SM.getLeptons(SM.ELECTRON).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = myOneLoopEW_HV.mq(q, sqrt(s));
        ;
        double ve = vl(SM.ELECTRON, Mw), ae = al(SM.ELECTRON, Mw);
        double vf = vq(q, Mw), af = aq(q, Mw);
        double ve2 = ve*ve, ae2 = ae*ae, vf2 = vf*vf, af2 = af*af;
        gslpp::complex chi1 = chi(1, s, Mw, bDP);
        gslpp::complex chi2 = chi(2, s, Mw, bDP);

        G3 += Qf * Qf * C11A(s, mf, Qf) * chi1.abs2()
                + 2.0 * Qe * Qf * ((ae * af * C12V(s, GammaZ, mf, Qf) + ve * vf * C12A(s, mf, Qf))
                * chi1 * chi2.conjugate()).real()
                + (4.0 * ve * ae * vf * af * C22V(s, GammaZ, mf, Qf)
                + (ve2 + ae2)*(vf2 + af2) * C22A(s, mf, Qf)) * chi2.abs2();
    }

    return G3;
}


//////////////////////////////////////////////////////////////////////// 
// QED corrections    

double EWSMTwoFermionsLEP2_Hollik::delta() const
{
    return ( 1.0 - 0.85 * 0.85); // sqrt{s'} > 0.85*sqrt{s}

}

double EWSMTwoFermionsLEP2_Hollik::Bf(const double s, const double mf) const
{
    return ( log(s / mf / mf) - 1.0);

}

double EWSMTwoFermionsLEP2_Hollik::gamma_delta(const double s, const double mf, const double Qf) const
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();

    return ( 2.0 * SM.getAle() / M_PI
            * (Bf(s, me) + Qf * Qf * Bf(s, mf)) * log(delta()));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::gamma_delta_int(const double s, const double GammaZ,
        const double mf, const double Qf) const
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    gslpp::complex M2 = gslpp::complex(Mz*Mz, -Mz*GammaZ, false);
    double d = delta();

    return ( 2.0 * SM.getAle() / M_PI
            * (Bf(s, me) * log(d * (s - M2) / (s - s * d - M2)) + Qf * Qf * Bf(s, mf) * log(d)));
}

double EWSMTwoFermionsLEP2_Hollik::gamma_delta_res(const double s, const double GammaZ,
        const double mf, const double Qf) const
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    gslpp::complex M2 = gslpp::complex(Mz*Mz, -Mz*GammaZ, false);
    double d = delta();

    return ( 2.0 * SM.getAle() / M_PI
            * (Bf(s, me) * log((d * (s - M2) / (s - s * d - M2)).abs())
            + Qf * Qf * Bf(s, mf) * log(d)));
}

double EWSMTwoFermionsLEP2_Hollik::gamma_tail(const double s, const double GammaZ) const
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();
    double Mz = SM.getMz();
    double d = delta();

    return ( 2.0 * SM.getAle() / M_PI
            * Bf(s, me)*(s - Mz * Mz) / Mz / GammaZ
            * (atan((Mz * Mz - s + s * d) / Mz / GammaZ) - atan((Mz * Mz - s) / Mz / GammaZ)));
}

double EWSMTwoFermionsLEP2_Hollik::gamma_fin(const double s, const double mf, const double Qf) const
{
    double me = SM.getLeptons(SM.ELECTRON).getMass();

    return ( 3.0 * SM.getAle() / 2.0 / M_PI * (Bf(s, me) + Qf * Qf * Bf(s, mf))
            + SM.getAle() / M_PI * (1 + Qf * Qf)*(M_PI * M_PI / 3.0 - 1.0 / 2.0));
}

double EWSMTwoFermionsLEP2_Hollik::C11V(const double s, const double mf, const double Qf) const
{
    return ( gamma_delta(s, mf, Qf) + gamma_fin(s, mf, Qf));
}

double EWSMTwoFermionsLEP2_Hollik::C11A(const double s, const double mf, const double Qf) const
{
    return 0.0;
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C12V(const double s, const double GammaZ,
        const double mf, const double Qf) const
{
    return ( gamma_delta_int(s, GammaZ, mf, Qf).conjugate() + gamma_fin(s, mf, Qf));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C12A(const double s, const double mf, const double Qf) const
{
    return gslpp::complex(0.0, 0.0, false);
}

double EWSMTwoFermionsLEP2_Hollik::C22V(const double s, const double GammaZ,
        const double mf, const double Qf) const
{
    return ( gamma_delta_res(s, GammaZ, mf, Qf) + gamma_tail(s, GammaZ) + gamma_fin(s, mf, Qf));
}

double EWSMTwoFermionsLEP2_Hollik::C22A(const double s, const double mf, const double Qf) const
{
    return 0.0;
}


//////////////////////////////////////////////////////////////////////// 
// Loop functions

gslpp::complex EWSMTwoFermionsLEP2_Hollik::B0bar_Hollik(const double s, const double m1,
        const double m2) const
{
    double mu = sqrt(s); // The result is independent of the renormalization scale. 
    return ( PV.B0(mu*mu, s, m1*m1, m2 * m2) + log(m1 * m2 / mu / mu));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::B1bar_Hollik(const double s, const double m1,
        const double m2) const
{
    double mu = sqrt(s); // The result is independent of the renormalization scale. 
    return ( PV.B1(mu*mu, s, m1*m1, m2 * m2) - log(m1 * m2 / mu / mu) / 2.0);
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::B1barPrime_Hollik(const double s, const double m1,
        const double m2) const
{
    double mu = sqrt(s); // The result is independent of the renormalization scale. 
    return ( PV.B1p(mu*mu, s, m1*m1, m2 * m2));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C0_Hollik(const double s, const double M,
        const double Mprime) const
{
    return ( -PV.C0(s, M*M, Mprime*Mprime, M * M));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C1plus_Hollik(const double s, const double M,
        const double Mprime) const
{
    double mb = myOneLoopEW_HV.mq(QCD::BOTTOM, sqrt(s));
    if (s == 4.0 * mb * mb)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::C1plus()");

    return ( (log(Mprime / M) + B0bar_Hollik(s, M, M)
            - B0bar_Hollik(mb*mb, M, Mprime)
            + (Mprime * Mprime - M * M + mb * mb) * C0_Hollik(s, M, Mprime))
            / (4.0 * mb * mb - s));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C2zero_Hollik(const double s, const double M,
        const double Mprime) const
{
    double mb = myOneLoopEW_HV.mq(QCD::BOTTOM, sqrt(s));
    return ( (B0bar_Hollik(s, M, M) + 1.0) / 4.0
            + (M * M - Mprime * Mprime - mb * mb) / 2.0 * C1plus_Hollik(s, M, Mprime)
            + Mprime * Mprime / 2.0 * C0_Hollik(s, M, Mprime));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C2plus_Hollik(const double s, const double M,
        const double Mprime) const
{
    double mb = myOneLoopEW_HV.mq(QCD::BOTTOM, sqrt(s));
    if (s == 4.0 * mb * mb)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::C2plus()");

    return ( (B0bar_Hollik(s, M, M) / 2.0
            + (B1bar_Hollik(mb*mb, Mprime, M) - 1.0 / 4.0) / 2.0
            + (Mprime * Mprime - M * M + mb * mb) * C1plus_Hollik(s, M, Mprime)
            - C2zero_Hollik(s, M, Mprime)) / (4.0 * mb * mb - s));
}

gslpp::complex EWSMTwoFermionsLEP2_Hollik::C2minus_Hollik(const double s, const double M,
        const double Mprime) const
{
    if (s == 0.0)
        throw std::runtime_error("Error in EWSMTwoFermionsLEP2_Hollik::C1plus()");

    double mb = myOneLoopEW_HV.mq(QCD::BOTTOM, sqrt(s));
    return ( (-(B1bar_Hollik(mb*mb, Mprime, M) - 1.0 / 4.0) / 2.0
            - C2zero_Hollik(s, M, Mprime)) / s);
}


