/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LEP2oblique.h"
#include <cmath>
#include <complex>

LEP2oblique::LEP2oblique(const StandardModel& SM_i) : SM(SM_i)
{
}

double LEP2oblique::Rates::deltaAFB() const
{
    return (deltaForwardBackward - (forwardBackward / sigma) * deltaSigma) / sigma;
}

LEP2oblique::Rates LEP2oblique::rates(const Particle& finalState, double s, double mass,
                                      double colors, const double oblique[]) const
{
    // The SM supplied by NPSTUVWXY is its pure SM reference. Match the photon
    // and neutral-current Born couplings to the LEP2 running-alpha and G_mu inputs.
    const double mz = SM.getMz(), mw = SM.Mw();
    const double sw2 = SM.sW2(mw), cw2 = SM.cW2(mw);
    const double sw = std::sqrt(sw2), cw = std::sqrt(cw2);
    const double e = std::sqrt(4. * M_PI * SM.ale_OS(std::sqrt(s), FULLNLO));
    const double gz = std::sqrt(4. * std::sqrt(2.) * SM.getGF()) * mz;
    typedef std::complex<double> Complex;
    const Complex chi = s / Complex(s - mz*mz, s * SM.Gamma_Z() / mz);
    const double contact = s / (mw * mw);

    // Barbieri et al., hep-ph/0405040, Eqs. (6), (10), (11).
    // Propagators below include a factor s; Z couples to T3 - sw2 Q.
    const double d1 = oblique[That] - oblique[W] + 2.*sw/cw*oblique[X] - sw2/cw2*oblique[Y];
    const double d3 = oblique[Shat] - oblique[W] + oblique[X]/(sw*cw) - oblique[Y];
    const double eZZ = cw2*oblique[W] - 2.*sw*cw*oblique[X] + sw2*oblique[Y];
    const double eAA = sw2*oblique[W] + 2.*sw*cw*oblique[X] + cw2*oblique[Y];
    const double eZA = (cw2-sw2)*oblique[X] + sw*cw*(oblique[W]-oblique[Y]);
    const Complex dZZ = chi*d1 - contact*eZZ;
    // At fixed alpha, G_F and M_Z, delta(s_eff^2) includes the shift of the
    // on-shell angle induced by delta M_W. The delta epsilon_2 terms cancel.
    double deltaSin2 = sw2 * (d3 - cw2*d1) / (cw2 - sw2);
    double deltaAlpha = 0.;
    if (SM.IsFlagMWinput()) {
        // M_W replaces alpha(M_Z) as an input. Hold the on-shell angle fixed
        // and include the predicted relative shift of the hard photon coupling squared.
        const double d2 = oblique[Uhat] - oblique[V] - oblique[W] + 2.*sw/cw*oblique[X];
        deltaSin2 = cw2 * (d1 - d2) - sw2 * d3;
        deltaAlpha = (cw2*d1 - 2.*sw2*d3 - (cw2-sw2)*d2) / sw2;
    }
    const Complex dZA = -chi * (gz/e) * deltaSin2 - contact*eZA;
    const double dAA = deltaAlpha - contact*eAA;

    const Particle& electron = SM.getLeptons(QCD::ELECTRON);
    const double qe = e * electron.getCharge(), qf = e * finalState.getCharge();
    const double af = gz * finalState.getIsospin() / 2.;
    const double vf = af - gz * sw2 * finalState.getCharge();
    const double r = mass*mass/s, beta2 = 1. - 4.*r;
    double total = 0., deltaTotal = 0., odd = 0., deltaOdd = 0.;
    for (int chirality : {1, -1}) { // left, right incoming electron
        const double t3 = chirality == 1 ? electron.getIsospin() : 0.;
        const double ze = gz * (t3 - sw2 * electron.getCharge());
        const Complex vec = qe*qf + ze*vf*chi;
        const Complex ax = ze*af*chi;
        const Complex dVec = qe*qf*dAA + (ze*qf + qe*vf)*dZA + ze*vf*dZZ;
        const Complex dAx = qe*af*dZA + ze*af*dZZ;

        total += (1.+2.*r)*std::norm(vec) + beta2*std::norm(ax);
        deltaTotal += 2.*std::real((1.+2.*r)*std::conj(vec)*dVec + beta2*std::conj(ax)*dAx);
        odd += chirality * std::real(vec*std::conj(ax));
        deltaOdd += chirality * std::real(dVec*std::conj(ax) + vec*std::conj(dAx));
    }
    const double beta = std::sqrt(beta2);
    const double normalization = colors*beta/(24.*M_PI*s);
    const double fbNormalization = normalization*1.5*beta;
    return {normalization*total, normalization*deltaTotal, fbNormalization*odd, fbNormalization*deltaOdd};
}

double LEP2oblique::sigma_l_LEP2_NP(QCD::lepton l, double s, double ml, const double oblique[]) const
{
    return rates(SM.getLeptons(l), s, ml, 1., oblique).deltaSigma;
}

double LEP2oblique::sigma_q_LEP2_NP(QCD::quark q, double s, double mq, const double oblique[]) const
{
    return rates(SM.getQuarks(q), s, mq, 3., oblique).deltaSigma;
}

double LEP2oblique::AFB_l_LEP2_NP(QCD::lepton l, double s, double ml, const double oblique[]) const
{
    return rates(SM.getLeptons(l), s, ml, 1., oblique).deltaAFB();
}

double LEP2oblique::AFB_q_LEP2_NP(QCD::quark q, double s, double mq, const double oblique[]) const
{
    return rates(SM.getQuarks(q), s, mq, 3., oblique).deltaAFB();
}

double LEP2oblique::R_q_LEP2_NP(QCD::quark q, double s, double mq, const double oblique[]) const
{
    const Rates selected = rates(SM.getQuarks(q), s, mq, 3., oblique);
    double hadronic = 0., deltaHadronic = 0.;
    for (const auto flavor : {QCD::UP, QCD::DOWN, QCD::STRANGE, QCD::CHARM, QCD::BOTTOM}) {
        const Rates channel = flavor == q ? selected
                : rates(SM.getQuarks(flavor), s, SM.getmq(flavor, std::sqrt(s)), 3., oblique);
        hadronic += channel.sigma;
        deltaHadronic += channel.deltaSigma;
    }
    return (selected.deltaSigma - (selected.sigma / hadronic)*deltaHadronic) / hadronic;
}
