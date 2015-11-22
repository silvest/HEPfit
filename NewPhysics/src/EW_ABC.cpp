/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_ABC.h"

EW_ABC::EW_ABC(const NPEpsilons& NPE_i)
: NPE(NPE_i)
{
}

double EW_ABC::Mw(const bool bAlternative) const
{
    if (!bAlternative)
        return NPE.Mw_eps(eps1(), eps2(), eps3());
    else {
        double delta_alpha = (NPE.alphaMz() - 1.0 / 128.90) / NPE.getAle();
        double cW2_Born = 0.768905 * (1.0 - 0.40 * delta_alpha);
        double cW2 = cW2_Born * (1.0 + 1.43 * eps1() - 1.00 * eps2() - 0.86 * eps3());
        return ( sqrt(cW2) * NPE.getMz());
    }
}

double EW_ABC::Gamma_l(StandardModel::lepton l) const
{
    double Qf = NPE.getLeptons(l).getCharge();
    double RQED = 1.0 + 3.0 * NPE.alphaMz() / 4.0 / M_PI * Qf*Qf;
    double mf = NPE.getLeptons(l).getMass();
    double beta = sqrt(1.0 - 4.0 * mf * mf / NPE.getMz() / NPE.getMz());
    double factor = NPE.getGF() * NPE.getMz() * NPE.getMz() * NPE.getMz() / 6.0 / M_PI / sqrt(2.0);
    return ( factor * beta * ((3.0 - beta * beta) / 2.0 * gVl(l).abs2()
            + beta * beta * gAl(l).abs2()) * RQED);
}

double EW_ABC::Gamma_q(QCD::quark q) const
{
    if (q == QCD::BOTTOM || q == QCD::TOP)
        throw std::runtime_error("Error in EW_ABC::Gamma_q()");
    double Qf = NPE.getQuarks(q).getCharge();
    double RQED = 1.0 + 3.0 * NPE.alphaMz() / 4.0 / M_PI * Qf*Qf;
    double a = NPE.Als(NPE.getMz(), FULLNNLO) / M_PI;
    double RQCD = 1.0 + 1.2 * a - 1.1 * a * a - 13.0 * a * a*a;
    double mf = 0.0;
    if (q == QCD::CHARM)
        mf = 1.67; // pole mass (PDG2012)
    else
        mf = 0.0;
    //mf = NPE.Mrun(NPE.getMz(),NPE.getQuarks(q).getMass_scale(),NPE.getQuarks(q).getMass(),FULLNNLO);
    double beta = sqrt(1.0 - 4.0 * mf * mf / NPE.getMz() / NPE.getMz());
    double Nc = 3.0;
    double factor = NPE.getGF() * NPE.getMz() * NPE.getMz() * NPE.getMz() / 6.0 / M_PI / sqrt(2.0);

    return ( factor * beta
            * ((3.0 - beta * beta) / 2.0 * gVq(q).abs2()
            + beta * beta * gAq(q).abs2()) * Nc * RQCD * RQED);
}

double EW_ABC::Gamma_b() const
{
    double RQED = 1.0 + NPE.alphaMz() / 12.0 / M_PI;
    double a = NPE.Als(NPE.getMz(), FULLNNLO) / M_PI;
    double RQCD = 1.0 + 1.2 * a - 1.1 * a * a - 13.0 * a * a*a;
    //double mb = 4.78; // pole mass (PDG2012)
    double mb = 4.7; // used by Altarelli et al (1998)
    double beta = sqrt(1.0 - 4.0 * mb * mb / NPE.getMz() / NPE.getMz());
    double Nc = 3.0;
    double factor = NPE.getGF() * NPE.getMz() * NPE.getMz() * NPE.getMz() / 6.0 / M_PI / sqrt(2.0);

    return ( factor * beta
            * ((3.0 - beta * beta) / 2.0 * gVb().abs2()
            + beta * beta * gAb().abs2()) * Nc * RQCD * RQED);
}

double EW_ABC::Gamma_had() const
{
    return ( Gamma_q(NPE.UP) + Gamma_q(NPE.DOWN)
            + Gamma_q(NPE.CHARM) + Gamma_q(NPE.STRANGE)
            + Gamma_b());
}

double EW_ABC::GammaZ(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_l(NPE.NEUTRINO_1) + Gamma_l(NPE.NEUTRINO_2)
            + Gamma_l(NPE.NEUTRINO_3) + Gamma_l(NPE.ELECTRON)
            + Gamma_l(NPE.MU) + Gamma_l(NPE.TAU)
            + Gamma_had());
    else {
        double delta_als = (NPE.Als(NPE.getMz(), FULLNNLO) - 0.119) / M_PI;
        double delta_alpha = (NPE.alphaMz() - 1.0 / 128.90) / NPE.getAle();
        double Gamma_T0 = 2.48946 * (1.0 + 0.73 * delta_als - 0.35 * delta_alpha);
        return ( Gamma_T0 * (1.0 + 1.35 * eps1() - 0.46 * eps3() + 0.35 * epsb()));
    }
}

double EW_ABC::R_l(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_had() / Gamma_l(NPE.ELECTRON));
    else {
        double delta_als = (NPE.Als(NPE.getMz(), FULLNNLO) - 0.119) / M_PI;
        double delta_alpha = (NPE.alphaMz() - 1.0 / 128.90) / NPE.getAle();
        double R_0 = 20.8228 * (1.0 + 1.05 * delta_als - 0.28 * delta_alpha);
        return ( R_0 * (1.0 + 0.28 * eps1() - 0.36 * eps3() + 0.50 * epsb()));
    }
}

double EW_ABC::R_c() const
{
    return ( Gamma_q(NPE.CHARM) / Gamma_had());
}

double EW_ABC::R_b(const bool bAlternative) const
{
    if (!bAlternative)
        return ( Gamma_b() / Gamma_had());
    else {
        double R_b0 = 0.2182355;
        return ( R_b0 * (1.0 - 0.06 * eps1() + 0.07 * eps3() + 1.79 * epsb()));
    }
}

double EW_ABC::sigma0_had(const bool bAlternative) const
{
    if (!bAlternative)
        return ( 12.0 * M_PI / NPE.getMz() / NPE.getMz()
            * Gamma_l(NPE.ELECTRON) * Gamma_had()
            / GammaZ(false) / GammaZ(false));
    else {
        double delta_als = (NPE.Als(NPE.getMz(), FULLNNLO) - 0.119) / M_PI;
        double delta_alpha = (NPE.alphaMz() - 1.0 / 128.90) / NPE.getAle();
        double sigma_h0 = 41.420 * (1.0 - 0.41 * delta_als + 0.03 * delta_alpha);
        return ( sigma_h0 * (1.0 - 0.03 * eps1() + 0.04 * eps3() - 0.20 * epsb()));
    }
}

double EW_ABC::A_l(StandardModel::lepton l, const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(l).real();
        return ( 2.0 * x / (1.0 + x * x));
    } else {
        double delta_alpha = (NPE.alphaMz() - 1.0 / 128.90) / NPE.getAle();
        double x0 = 0.075619 - 1.32 * delta_alpha;
        double x = x0 * (1.0 + 17.6 * eps1() - 22.9 * eps3());
        return ( 2.0 * x / (1.0 + x * x));
    }
}

double EW_ABC::A_q(QCD::quark q) const
{
    if (q == QCD::BOTTOM || q == QCD::TOP)
        throw std::runtime_error("Error in EW_ABC::A_q()");
    double x = gVq_over_gAq(q).real();
    return ( 2.0 * x / (1.0 + x * x));
}

double EW_ABC::A_b() const
{
    double x = gVb_over_gAb().real();
    return ( 2.0 * x / (1.0 + x * x));
}

double EW_ABC::AFB_l(StandardModel::lepton l, const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(l).real();
        return ( 3.0 * x * x / (1.0 + x * x) / (1.0 + x * x));
    } else {
        double delta_als = (NPE.Als(NPE.getMz(), FULLNNLO) - 0.119) / M_PI;
        double AFB_l_Born = 0.01696 * (1.0 - 34.0 * delta_als);
        return ( AFB_l_Born * (1.0 + 34.72 * eps1() - 45.15 * eps3()));
    }
}

double EW_ABC::AFB_c() const
{
    double x = gVl_over_gAl(NPE.ELECTRON).real();
    double xc = gVq_over_gAq(NPE.CHARM).real();
    return ( 3.0 * x * xc / (1.0 + x * x) / (1.0 + xc * xc));
}

double EW_ABC::AFB_b() const
{
    double x = gVl_over_gAl(NPE.ELECTRON).real();
    double xb = gVb_over_gAb().real();
    return ( 3.0 * x * xb / (1.0 + x * x) / (1.0 + xb * xb));
}

double EW_ABC::sin2thetaEff(const bool bAlternative) const
{
    if (!bAlternative) {
        double x = gVl_over_gAl(NPE.ELECTRON).real();
        return ( (1.0 - x) / 4.0);
    } else {
        double delta_als = (NPE.Als(NPE.getMz(), FULLNNLO) - 0.119) / M_PI;
        double x0 = 0.075619 - 1.32 * delta_als;
        double x = x0 * (1.0 + 17.6 * eps1() - 22.9 * eps3());
        return ( (1.0 - x) / 4.0);
    }
}


////////////////////////////////////////////////////////////////////////

double EW_ABC::eps1() const
{
    return NPE.epsilon1();
}

double EW_ABC::eps2() const
{
    return NPE.epsilon2();
}

double EW_ABC::eps3() const
{
    return NPE.epsilon3();
}

double EW_ABC::epsb() const
{
    return NPE.epsilonb();
}


////////////////////////////////////////////////////////////////////////

gslpp::complex EW_ABC::gVl(StandardModel::lepton l) const
{
    return NPE.gV_f_eps(NPE.getLeptons(l), eps1(), eps3());
}

gslpp::complex EW_ABC::gAl(StandardModel::lepton l) const
{
    return NPE.gA_f_eps(NPE.getLeptons(l), eps1());
}

gslpp::complex EW_ABC::gVl_over_gAl(StandardModel::lepton l) const
{
    return ( gVl(l) / gAl(l));
}

gslpp::complex EW_ABC::gVq(QCD::quark q) const
{
    if (q == QCD::BOTTOM || q == QCD::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq()");
    return NPE.gV_f_eps(NPE.getQuarks(q), eps1(), eps3());
}

gslpp::complex EW_ABC::gAq(QCD::quark q) const
{
    if (q == QCD::BOTTOM || q == QCD::TOP)
        throw std::runtime_error("Error in EW_ABC::gAq()");
    return NPE.gA_f_eps(NPE.getQuarks(q), eps1());
}

gslpp::complex EW_ABC::gVq_over_gAq(QCD::quark q) const
{
    if (q == QCD::BOTTOM || q == QCD::TOP)
        throw std::runtime_error("Error in EW_ABC::gVq_over_gAq()");
    return ( gVq(q) / gAq(q));
}

gslpp::complex EW_ABC::gVb() const
{
    return NPE.gV_f_eps(NPE.getQuarks(NPE.BOTTOM), eps1(), eps3(), epsb());
}

gslpp::complex EW_ABC::gAb() const
{
    return NPE.gA_f_eps(NPE.getQuarks(NPE.BOTTOM), eps1(), epsb());
}

gslpp::complex EW_ABC::gVb_over_gAb() const
{
    return ( gVb() / gAb());
}

