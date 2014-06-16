/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPbase.h"

NPbase::NPbase()
: StandardModel()
{
}

bool NPbase::PostUpdate()
{
    bool SMup = trueSM.PostUpdate();
    trueSM = *this;
    return (SMup);
}

double NPbase::Mw() const
{
    double myMw = trueSM.Mw();

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    myMw *= 1.0 - alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - s2 / 2.0 / (c2 - s2) * DeltaGF();

    //std::cout << "Mw: c_S=" << - alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Mw: c_T=" << - alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Mw: c_U=" << - alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return myMw;
}

double NPbase::GammaW() const
{
    double Gamma_W = trueSM.GammaW();

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    Gamma_W *= 1.0 - 3.0 * alpha / 4.0 / (c2 - s2)
            *(obliqueS() - 2.0 * c2 * obliqueT() - (c2 - s2) * obliqueU() / 2.0 / s2)
            - (1.0 + c2) / 2.0 / (c2 - s2) * DeltaGF();

    //std::cout << "Gw: c_S=" << - 3.0*alpha/4.0/(c2-s2) << std::endl;
    //std::cout << "Gw: c_T=" << - 3.0*alpha/4.0/(c2-s2)*(- 2.0*c2) << std::endl;
    //std::cout << "Gw: c_U=" << - 3.0*alpha/4.0/(c2-s2)*(- (c2-s2)/2.0/s2) << std::endl;

    return Gamma_W;
}

double NPbase::deltaGV_f(const Particle p) const
{
    if (p.is("TOP")) return 0.;

    /* SM values */
    double alpha = trueSM.alphaMz();
    double sW2SM = trueSM.sW2();
    double cW2SM = trueSM.cW2();
    double gVSM = trueSM.gV_f(p).real();
    double gASM = trueSM.gA_f(p).real();

    return ( gVSM * (alpha * obliqueT() - DeltaGF()) / 2.0
            + (gVSM - gASM) / 4.0 / sW2SM / (cW2SM - sW2SM)
            *(alpha * (obliqueS() - 4.0 * cW2SM * sW2SM * obliqueT())
            + 4.0 * cW2SM * sW2SM * DeltaGF()));
}

complex NPbase::gV_f(const Particle p) const
{
    return ( trueSM.gV_f(p) + deltaGV_f(p));
}

double NPbase::deltaGA_f(const Particle p) const
{
    if (p.is("TOP")) return 0.;
    /* SM values */
    double alpha = trueSM.alphaMz();
    double gASM = trueSM.gA_f(p).real();

    return ( gASM * (alpha * obliqueT() - DeltaGF()) / 2.0);
}

complex NPbase::gA_f(const Particle p) const
{
    return ( trueSM.gA_f(p) + deltaGA_f(p));
}

complex NPbase::rhoZ_f(const Particle p) const
{
    return ( gA_f(p) * gA_f(p) / p.getIsospin() / p.getIsospin());

}

complex NPbase::kappaZ_f(const Particle p) const
{
    return ( (1.0 - gV_f(p) / gA_f(p)) / (4.0 * fabs(p.getCharge()) * sW2()));
}

////////////////////////////////////////////////////////////////////////

double NPbase::deltaGamma_Z() const
{
    double deltaGamma_Z = 0.;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p] != 0.0 || delGAl[p] != 0.0
                || delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double deltaGl[6], deltaGq[6];
        double delGammaZ = 0.0;
        for (int p = 0; p < 6; ++p) {
            gVf = trueSM.gV_f(leptons[p]).real();
            gAf = trueSM.gA_f(leptons[p]).real();
            deltaGl[p] = 2.0 * (gVf * delGVl[p] + gAf * delGAl[p]);

            gVf = trueSM.gV_f(quarks[p]).real();
            gAf = trueSM.gA_f(quarks[p]).real();
            deltaGq[p] = 2.0 * (gVf * delGVq[p] + gAf * delGAq[p]);

            delGammaZ += deltaGl[p] + 3.0 * deltaGq[p];
        }

        double sW2_SM = trueSM.sW2();
        double cW2_SM = trueSM.cW2();
        deltaGamma_Z = alphaMz() * Mz / 12.0 / sW2_SM / cW2_SM
                * delGammaZ;
    }

    return deltaGamma_Z;
}

double NPbase::Gamma_Z() const
{
    return (trueSM.Gamma_Z() + deltaGamma_Z());
}

double NPbase::deltaSigmaHadron() const
{
    double sigma_had = 0.;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p = 0; p < 6; ++p) {
        delGVl[p] = deltaGV_f(leptons[p]);
        delGAl[p] = deltaGA_f(leptons[p]);
        delGVq[p] = deltaGV_f(quarks[p]);
        delGAq[p] = deltaGA_f(quarks[p]);
        if (delGVl[p] != 0.0 || delGAl[p] != 0.0
                || delGVq[p] != 0.0 || delGAq[p] != 0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double Gl[6], deltaGl[6], Gq[6], deltaGq[6];
        double Gq_sum = 0.0, delGq_sum = 0.0;
        double Gf_sum = 0.0, delGf_sum = 0.0;
        for (int p = 0; p < 6; ++p) {
            gVf = trueSM.gV_f(leptons[p]).real();
            gAf = trueSM.gA_f(leptons[p]).real();
            Gl[p] = gVf * gVf + gAf*gAf;
            deltaGl[p] = 2.0 * (gVf * delGVl[p] + gAf * delGAl[p]);

            gVf = trueSM.gV_f(quarks[p]).real();
            gAf = trueSM.gA_f(quarks[p]).real();
            Gq[p] = gVf * gVf + gAf*gAf;
            deltaGq[p] = 2.0 * (gVf * delGVq[p] + gAf * delGAq[p]);

            Gq_sum += 3.0 * Gq[p];
            Gf_sum += Gl[p] + 3.0 * Gq[p];
            delGq_sum += 3.0 * deltaGq[p];
            delGf_sum += deltaGl[p] + 3.0 * deltaGq[p];
        }

        sigma_had = 12.0 * M_PI / Mz / Mz
                * Gl[ELECTRON] * Gq_sum / Gf_sum / Gf_sum
                * (deltaGl[ELECTRON] / Gl[ELECTRON]
                + delGq_sum / Gq_sum - 2.0 * delGf_sum / Gf_sum);
    }

    return sigma_had;
}

double NPbase::sigma0_had() const
{
    return (trueSM.sigma0_had() + deltaSigmaHadron());
}

double NPbase::deltaSin2thetaEff_e() const
{
    double sin2_theta_eff = 0.;
    double delGVf = deltaGV_f(leptons[ELECTRON]);
    double delGAf = deltaGA_f(leptons[ELECTRON]);
    if (delGVf != 0.0 || delGAf != 0.0) {
        double gVf = trueSM.gV_f(leptons[ELECTRON]).real();
        double gAf = trueSM.gA_f(leptons[ELECTRON]).real();
        double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

        sin2_theta_eff = -delGVfOverGAf / 4.0;
    }
    return sin2_theta_eff;
}

double NPbase::sin2thetaEff(const Particle p) const
{
    if (p.is("ELECTRON"))
        return (trueSM.sin2thetaEff(p) + deltaSin2thetaEff_e());
    else
        return (trueSM.sin2thetaEff(p));
}

double NPbase::deltaA_f(const Particle p) const
{
    double dAf = 0.;
    double delGVf = deltaGV_f(p);
    double delGAf = deltaGA_f(p);
    if (delGVf != 0.0 || delGAf != 0.0) {
        double gVf = trueSM.gV_f(p).real();
        double gAf = trueSM.gA_f(p).real();
        double Gf = gVf * gVf + gAf*gAf;
        double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

        dAf = -2.0 * (gVf * gVf - gAf * gAf) * gAf * gAf / Gf / Gf*delGVfOverGAf;
    }

    return dAf;
}

double NPbase::A_f(const Particle p) const
{
    return (trueSM.A_f(p) + deltaA_f(p));
}

double NPbase::deltaAFB(const Particle p) const
{
    double dAFB = 0.;
    double delGVf = deltaGV_f(p);
    double delGAf = deltaGA_f(p);
    if (p.is("LEPTON")) {
        if (delGVf != 0.0 || delGAf != 0.0) {
            double gVe = trueSM.gV_f(p).real();
            double gAe = trueSM.gA_f(p).real();
            double Ge = gVe * gVe + gAe*gAe;
            double delGVeOverGAe = (gAe * delGVf - gVe * delGAf) / gAe / gAe;
            dAFB = -6.0 * gVe * gAe * (gVe * gVe - gAe * gAe) * gAe * gAe / Ge / Ge / Ge*delGVeOverGAe;
        }
    } else {
        double delGVe = deltaGV_f(leptons[ELECTRON]);
        double delGAe = deltaGA_f(leptons[ELECTRON]);
        if (delGVe != 0.0 || delGAe != 0.0 || delGVf != 0.0 || delGAf != 0.0) {
            double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
            double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
            double Ge = gVe * gVe + gAe*gAe;
            double delGVeOverGAe = (gAe * delGVe - gVe * delGAe) / gAe / gAe;
            //
            double gVf = trueSM.gV_f(p).real();
            double gAf = trueSM.gA_f(p).real();
            double Gf = gVf * gVf + gAf*gAf;
            double delGVfOverGAf = (gAf * delGVf - gVf * delGAf) / gAf / gAf;

            dAFB = -(3.0 * gVf * gAf * (gVe * gVe - gAe * gAe) * gAe * gAe / Gf / Ge / Ge * delGVeOverGAe
                    + 3.0 * gVe * gAe * (gVf * gVf - gAf * gAf) * gAf * gAf / Ge / Gf / Gf * delGVfOverGAf);
        }
    }

    return dAFB;
}

double NPbase::AFB(const Particle p) const
{
    return (trueSM.AFB(p) + deltaAFB(p));
}

double NPbase::deltaR0_f(const Particle p) const
{
    double dR0_f = 0., delGVe = 0., delGAe = 0., deltaGe = 0., Ge = 0.;
    bool nonZeroNP = false;
    if (p.is("LEPTON")) {
        delGVe = deltaGV_f(leptons[ELECTRON]);
        delGAe = deltaGA_f(leptons[ELECTRON]);
        if (delGVe != 0.0 || delGAe != 0.0) nonZeroNP = true;
    }

    double delGVq[6], delGAq[6];
    for (int q = 0; q < 6; ++q) {
        delGVq[q] = deltaGV_f(quarks[q]);
        delGAq[q] = deltaGA_f(quarks[q]);
        if (delGVq[q] != 0.0 || delGAq[q] != 0.0) nonZeroNP = true;
    }

    if (nonZeroNP) {
        double CF = 1.;
        if (p.is("LEPTON")) {
            double gVe = trueSM.gV_f(leptons[ELECTRON]).real();
            double gAe = trueSM.gA_f(leptons[ELECTRON]).real();
            Ge = gVe * gVe + gAe*gAe;
            deltaGe = 2.0 * (gVe * delGVe + gAe * delGAe);
            CF = 3.;
        }
        double Gq[6], deltaGq[6];
        double gVq, gAq;
        double Gq_sum = 0.0, delGq_sum = 0.0;
        for (int q = 0; q < 6; ++q) {
            gVq = trueSM.gV_f(quarks[q]).real();
            gAq = trueSM.gA_f(quarks[q]).real();
            Gq[q] = gVq * gVq + gAq*gAq;
            deltaGq[q] = 2.0 * (gVq * delGVq[q] + gAq * delGAq[q]);

            Gq_sum += CF * Gq[q];
            delGq_sum += CF * deltaGq[q];
        }
        if (p.is("LEPTON"))
            dR0_f = delGq_sum / Ge - Gq_sum * deltaGe / Ge / Ge;
        else
            dR0_f = deltaGq[p.getIndex() - 6] / Gq_sum
                - Gq[p.getIndex() - 6] * delGq_sum / Gq_sum / Gq_sum;
    }
    return dR0_f;
}

double NPbase::R0_f(const Particle p) const
{
    return (trueSM.R0_f(p) + deltaR0_f(p));
}
