/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <stdexcept>
#include "EWSMOneLoopEW_HV.h"
#include "EWSMOneLoopEW.h"

EWSMOneLoopEW_HV::EWSMOneLoopEW_HV(const StandardModel& SM_i)
: SM(SM_i), PV(true)
{
}


//////////////////////////////////////////////////////////////////////// 

gslpp::complex EWSMOneLoopEW_HV::SigmaWW_bos(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2 / Mz2, cW4 = cW2*cW2;
    double sW2 = 1.0 - cW2, sW4 = sW2*sW2;
    double w = -s / Mw2, h = -s / mh2;
    double wh = mh2 / Mw2, zh = mh2 / Mz2;

    gslpp::complex Sigma(0.0, 0.0, false);
    if (s == 0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaWW_bos(s=0)");
    } else {
        /* Loop functions */
        double A0_Mw2 = PV.A0(mu*mu, Mw2);
        double A0_mh2 = PV.A0(mu*mu, mh2);
        double A0_Mz2 = PV.A0(mu*mu, Mz2);
        gslpp::complex B0_s_Mw2_Mz2 = PV.B0(mu*mu, s, Mw2, Mz2);
        gslpp::complex B0_s_mh2_Mw2 = PV.B0(mu*mu, s, mh2, Mw2);
        gslpp::complex B0_s_0_Mw2 = PV.B0(mu*mu, s, 0.0, Mw2);

        Sigma = Mw2 / 12.0
                * (-(sW4 / cW4 * (1.0 + 8.0 * cW2) / w - 10.0 / cW2 + 54.0
                + 16.0 * cW2 + (1.0 - 40.0 * cW2) * w) * B0_s_Mw2_Mz2
                - ((1.0 - wh)*(1.0 - wh) / w + 2.0 * wh - 10.0 + w) * B0_s_mh2_Mw2
                - 8.0 * sW2 * (1.0 / w + 2.0 - 5.0 * w) * B0_s_0_Mw2
                + ((wh + 1.0 / cW2 - 2.0) / w + 36.0 / wh - 14.0) * A0_Mw2 / Mw2
                - (1.0 / h - 1.0 / w - 7.0) * A0_mh2 / Mw2
                - (sW2 / cW2 * (1.0 + 8.0 * cW2) / w - 18.0 / zh - 1.0 + 16.0 * cW2)
                * A0_Mz2 / Mw2
                + 12.0 / wh * (1.0 / cW4 + 2.0)
                - 2.0 * (1.0 / cW2 + 18.0 + wh - 2.0 * w / 3.0));
    }
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaWW_fer(const double mu, const double muForMq,
        const double s) const
{
    double ml[6], mq[6], ml2[6], mq2[6];
    for (int i = 0; i < 6; i++) {
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((QCD::quark) i, muForMq);
        ml2[i] = ml[i] * ml[i];
        mq2[i] = mq[i] * mq[i];
    }

    gslpp::complex Sigma(0.0, 0.0, false);
    if (s == 0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaWW_fer(s=0)");
    } else {
        /* Loop functions */
        gslpp::complex B1_s_ml2_mlprime2[3], B1_s_mq2_mqprime2[3];
        gslpp::complex B1_s_mlprime2_ml2[3], B1_s_mqprime2_mq2[3];
        gslpp::complex Bf_s_mlprime2_ml2[3], Bf_s_mqprime2_mq2[3];
        for (int gen = 0; gen < 3; gen++) {
            B1_s_ml2_mlprime2[gen] = PV.B1(mu*mu, s, ml2[2 * gen], ml2[2 * gen + 1]);
            B1_s_mq2_mqprime2[gen] = PV.B1(mu*mu, s, mq2[2 * gen], mq2[2 * gen + 1]);
            B1_s_mlprime2_ml2[gen] = PV.B1(mu*mu, s, ml2[2 * gen + 1], ml2[2 * gen]);
            B1_s_mqprime2_mq2[gen] = PV.B1(mu*mu, s, mq2[2 * gen + 1], mq2[2 * gen]);
            Bf_s_mlprime2_ml2[gen] = PV.Bf(mu*mu, s, ml2[2 * gen + 1], ml2[2 * gen]);
            Bf_s_mqprime2_mq2[gen] = PV.Bf(mu*mu, s, mq2[2 * gen + 1], mq2[2 * gen]);
        }

        double mf2, mfprime2;
        for (int gen = 0; gen < 3; gen++) {
            mf2 = ml2[2 * gen];
            mfprime2 = ml2[2 * gen + 1];
            if (s != 0.0) Sigma += -s * Bf_s_mlprime2_ml2[gen];
            Sigma += mfprime2 * B1_s_ml2_mlprime2[gen] + mf2 * B1_s_mlprime2_ml2[gen];
            //
            mf2 = mq2[2 * gen];
            mfprime2 = mq2[2 * gen + 1];
            if (s != 0.0) Sigma += 3.0 * (-s * Bf_s_mqprime2_mq2[gen]);
            Sigma += 3.0 * (mfprime2 * B1_s_mq2_mqprime2[gen] + mf2 * B1_s_mqprime2_mq2[gen]);
        }
    }
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZZ_bos(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2 / Mz2, cW4 = cW2*cW2, cW6 = cW4*cW2;
    double z = -s / Mz2, h = -s / mh2;
    double wh = mh2 / Mw2, zh = mh2 / Mz2;

    gslpp::complex Sigma(0.0, 0.0, false);
    if (s == 0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaZZ_bos(s=0)");
    } else {
        /* Loop functions */
        double A0_Mw2 = PV.A0(mu*mu, Mw2);
        double A0_mh2 = PV.A0(mu*mu, mh2);
        double A0_Mz2 = PV.A0(mu*mu, Mz2);
        gslpp::complex B0_s_Mw2_Mw2 = PV.B0(mu*mu, s, Mw2, Mw2);
        gslpp::complex B0_s_mh2_Mz2 = PV.B0(mu*mu, s, mh2, Mz2);

        Sigma = Mz2 / 12.0
                * ((4.0 * cW2 * (5.0 - 8.0 * cW2 - 12.0 * cW4)
                - (1.0 - 4.0 * cW2 - 36.0 * cW4) * z) * B0_s_Mw2_Mw2
                - ((1.0 - zh)*(1.0 - zh) / z + 2.0 * zh - 10.0 + z) * B0_s_mh2_Mz2
                - (1.0 / h - 1.0 / z - 7.0) * A0_mh2 / Mz2
                + 2.0 * (18.0 / wh + 1.0 + 8.0 * cW2 - 24.0 * cW4) * A0_Mw2 / Mz2
                + (1.0 / h - 1.0 / z + 18.0 / zh + 1.0) * A0_Mz2 / Mz2
                + 2.0 * (6.0 * (1.0 + 2.0 * cW4) / zh - zh - 1.0 - 2.0 * cW2 + 8.0 * cW4
                - 24.0 * cW6 - 2.0 / 3.0 * (1.0 - 2.0 * cW2) * z));
    }
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZZ_fer(const double mu, const double muForMq,
        const double s, const double Mw) const
{
    double ml[6], mq[6], ml2[6], mq2[6];
    for (int i = 0; i < 6; i++) {
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((QCD::quark) i, muForMq);
        ml2[i] = ml[i] * ml[i];
        mq2[i] = mq[i] * mq[i];
    }

    gslpp::complex Sigma(0.0, 0.0, false);
    if (s == 0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaZZ_fer(s=0)");
    } else {
        /* Loop functions */
        gslpp::complex Bf_s_ml2_ml2[6], Bf_s_mq2_mq2[6];
        gslpp::complex B0_s_ml2_ml2[6], B0_s_mq2_mq2[6];
        for (int i = 0; i < 6; i++) {
            Bf_s_ml2_ml2[i] = PV.Bf(mu*mu, s, ml2[i], ml2[i]);
            Bf_s_mq2_mq2[i] = PV.Bf(mu*mu, s, mq2[i], mq2[i]);
            B0_s_ml2_ml2[i] = PV.B0(mu*mu, s, ml2[i], ml2[i]);
            B0_s_mq2_mq2[i] = PV.B0(mu*mu, s, mq2[i], mq2[i]);
        }

        double mf2, vf2, af2;
        for (int i = 0; i < 6; i++) {
            mf2 = ml2[i];
            vf2 = pow(vl((StandardModel::lepton) i, Mw), 2.0);
            af2 = pow(al((StandardModel::lepton) i), 2.0);
            Sigma += -(vf2 + af2) * s * Bf_s_ml2_ml2[i];
            Sigma += -2.0 * af2 * mf2 * B0_s_ml2_ml2[i];
            //
            mf2 = mq2[i];
            vf2 = pow(vq((QCD::quark) i, Mw), 2.0);
            af2 = pow(aq((QCD::quark) i), 2.0);
            Sigma += -3.0 * (vf2 + af2) * s * Bf_s_mq2_mq2[i];
            Sigma += -3.0 * 2.0 * af2 * mf2 * B0_s_mq2_mq2[i];
        }
    }
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaGammaGamma_bos(const double mu, const double s,
        const double Mw) const
{
    return ( s * PiGammaGamma_bos(mu, s, Mw));
}

gslpp::complex EWSMOneLoopEW_HV::SigmaGammaGamma_fer(const double mu, const double muForMq,
        const double s) const
{
    return ( s * PiGammaGamma_fer(mu, muForMq, s));
}

gslpp::complex EWSMOneLoopEW_HV::PiGammaGamma_bos(const double mu, const double s,
        const double Mw) const
{
    double Mw2 = Mw*Mw;
    double w = -s / Mw2;

    gslpp::complex Pi(0.0, 0.0, false);
    if (s == 0.0) {
        Pi = 3.0 * log(Mw2 / mu / mu) - 2.0 / 3.0;
    } else {
        /* Loop functions */
        double A0_Mw2 = PV.A0(mu*mu, Mw2);
        gslpp::complex B0_s_Mw2_Mw2 = PV.B0(mu*mu, s, Mw2, Mw2);

        Pi = -Mw2 / s * ((4.0 - 3.0 * w) * B0_s_Mw2_Mw2 + 4.0 * (A0_Mw2 / Mw2 + 1.0));
    }
    return Pi;
}

gslpp::complex EWSMOneLoopEW_HV::PiGammaGamma_fer_l(const double mu, const double s,
        const StandardModel::lepton l) const
{
    // Neutrinos do not contribute, since Qf=0.
    if ((l == StandardModel::NEUTRINO_1) || (l == StandardModel::NEUTRINO_2)
            || (l == StandardModel::NEUTRINO_3))
        return 0.0;

    double mf = this->ml((StandardModel::lepton) l);
    double Qf = SM.getLeptons(l).getCharge();

    /* Loop functions */
    gslpp::complex Bf_s_mf2_mf2;
    if (mf == 0.0)
        Bf_s_mf2_mf2 = 0.0;
    else
        Bf_s_mf2_mf2 = PV.Bf(mu*mu, s, mf*mf, mf * mf);

    return ( -4.0 * Qf * Qf * Bf_s_mf2_mf2);
}

gslpp::complex EWSMOneLoopEW_HV::PiGammaGamma_fer_q(const double mu, const double muForMq,
        const double s, const QCD::quark q) const
{
    double mf = this->mq((QCD::quark) q, muForMq);
    double Qf = SM.getQuarks(q).getCharge();

    /* Loop functions */
    gslpp::complex Bf_s_mf2_mf2;
    if (mf == 0.0)
        Bf_s_mf2_mf2 = 0.0;
    else
        Bf_s_mf2_mf2 = PV.Bf(mu*mu, s, mf*mf, mf * mf);

    return ( -4.0 * 3.0 * Qf * Qf * Bf_s_mf2_mf2);
}

gslpp::complex EWSMOneLoopEW_HV::PiGammaGamma_fer(const double mu, const double muForMq,
        const double s) const
{
    gslpp::complex Pi(0.0, 0.0, false);
    for (int i = 0; i < 6; i++) {
        Pi += PiGammaGamma_fer_l(mu, s, (StandardModel::lepton) i);
        Pi += PiGammaGamma_fer_q(mu, muForMq, s, (QCD::quark) i);
    }
    return Pi;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZgamma_bos(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double cW2 = Mw2 / Mz2;
    double w = -s / Mw2;

    gslpp::complex Sigma(0.0, 0.0, false);
    if (s == 0.0) {
        Sigma = 2.0 * Mw2 * log(Mw2 / mu / mu);
    } else {
        /* Loop functions */
        double A0_Mw2 = PV.A0(mu*mu, Mw2);
        gslpp::complex B0_s_Mw2_Mw2 = PV.B0(mu*mu, s, Mw2, Mw2);

        Sigma = (4.0 * (1.0 / 3.0 + cW2) / w - 1.0 / 6.0 - 3.0 * cW2) * s * B0_s_Mw2_Mw2
                + (2.0 / 3.0 - 4.0 * cW2) * Mw2 * (A0_Mw2 / Mw2 + 1.0) - s / 9.0;
    }
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZgamma_fer(const double mu, const double muForMq,
        const double s, const double Mw) const
{
    double ml[6], mq[6], ml2[6], mq2[6];
    for (int i = 0; i < 6; i++) {
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((QCD::quark) i, muForMq);
        ml2[i] = ml[i] * ml[i];
        mq2[i] = mq[i] * mq[i];
    }
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double sW2 = 1.0 - Mw2 / Mz2;

    /* Loop functions */
    gslpp::complex Bf_s_ml2_ml2[6], Bf_s_mq2_mq2[6];
    for (int i = 0; i < 6; i++) {
        if (i == 0 || i == 2 || i == 4)
            Bf_s_ml2_ml2[i] = 0.0; // Neutrinos do not contribute, since Ql=0.
        else
            Bf_s_ml2_ml2[i] = PV.Bf(mu*mu, s, ml2[i], ml2[i]);
        Bf_s_mq2_mq2[i] = PV.Bf(mu*mu, s, mq2[i], mq2[i]);
    }

    gslpp::complex Pi(0.0, 0.0, false);
    double Ql, Qq;
    for (int i = 0; i < 6; i++) {
        Ql = SM.getLeptons((StandardModel::lepton) i).getCharge();
        Pi += -(fabs(Ql) - 4.0 * sW2 * Ql * Ql) * Bf_s_ml2_ml2[i];
        //
        Qq = SM.getQuarks((QCD::quark) i).getCharge();
        Pi += -3.0 * (fabs(Qq) - 4.0 * sW2 * Qq * Qq) * Bf_s_mq2_mq2[i];
    }
    return ( s * Pi);
}


//////////////////////////////////////////////////////////////////////// 

gslpp::complex EWSMOneLoopEW_HV::F_Hollik(const double s, const double m1,
        const double m2) const
{
    double m12 = m1*m1, m22 = m2*m2;
    double mu = SM.getMz(); // The result is independent of mu. 

    if (m1 != 0.0 && m2 != 0.0) {
        if (m1 == m2)
            return ( PV.B0(mu * mu, s, m12, m12) + log(m1 * m1 / mu / mu));
        else
            return ( PV.B0(mu*mu, s, m12, m22) + log(m1 * m2 / mu / mu) - 1.0
                + (m12 + m22) / (m12 - m22) * log(m1 / m2));
    } else if (m1 == 0.0 && m2 != 0.0)
        return ( PV.B0(mu * mu, s, 0.0, m22) + log(m2 * m2 / mu / mu) - 1.0);
    else if (m1 != 0.0 && m2 == 0.0)
        return ( PV.B0(mu * mu, s, m12, 0.0) + log(m1 * m1 / mu / mu) - 1.0);

    else
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::F_Hollik()");
}

gslpp::complex EWSMOneLoopEW_HV::Fprime_Hollik(const double muIR, const double s,
        const double m1, const double m2) const
{
    return ( PV.B0p(muIR*muIR, s, m1*m1, m2 * m2));
}

gslpp::complex EWSMOneLoopEW_HV::SigmaWW_bos_Hollik(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2 / Mz2;
    double sW2 = 1.0 - cW2;
    double w = Mw2, z = Mz2, h = mh2;

    if (mu <= 0.0 || s == 0.0 || w == 0.0 || z == 0.0 || h == 0.0 || h == w || z == w)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaWW_bos_Hollik()");

    double DeltaW = -log(Mw2 / mu / mu);

    gslpp::complex Sigma = -(19.0 / 2.0 * s + 3.0 * w * (1.0 - sW2 / cW2)) * DeltaW / 3.0 // correct, Hollik (90)
            //- (19.0/2.0*s + 3.0*w*(1.0 - sW2/cW2))*DeltaW // incorrect, Consoli, Hollik, Jegerlehner (89)
            + (sW2 * sW2 * z - cW2 / 3.0 * (7.0 * (z + w) + 10.0 * s - 2.0 * (z - w)*(z - w) / s)
            - (w + z - s / 2.0 - (z - w)*(z - w) / 2.0 / s) / 6.0)
            * F_Hollik(s, Mz, Mw)
            + sW2 / 3.0 * (-4.0 * w - 10.0 * s + 2.0 * w * w / s) * F_Hollik(s, 0.0, Mw)
            + (5.0 * w - h + s / 2.0 + (h - w)*(h - w) / 2.0 / s) / 6.0
            * F_Hollik(s, mh, Mw)
            + (cW2 / 3.0 * (7.0 * z + 7.0 * w + 10.0 * s - 4.0 * (z - w))
            - sW2 * sW2 * z + (2.0 * w - s / 2.0) / 6.0)*3.0 * z / (z - w) * log(z / w) // Hollik (90)
            //   - sW2*sW2*z +(2.0*w - s/2.0)/6.0)*z/(z - w)*log(z/w) // Consoli, Hollik, Jegerlehner (89)
            - (2.0 / 3.0 * w + s / 12.0) * h / (h - w) * log(h / w)
            - cW2 / 3.0 * (7.0 * z + 7.0 * w + 32.0 / 3.0 * s) + sW2 * sW2 * z
            + (5.0 / 3.0 * s + 4.0 * w - z - h) / 6.0
            - sW2 / 3.0 * (4.0 * w + 32.0 / 3.0 * s);
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZZ_bos_Hollik(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2 / Mz2;
    double sW2 = 1.0 - cW2;
    double w = Mw2, z = Mz2, h = mh2;

    if (mu <= 0.0 || s == 0.0 || w == 0.0 || z == 0.0 || h == 0.0 || h == z)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaZZ_bos_Hollik()");

    double DeltaW = -log(Mw2 / mu / mu);

    gslpp::complex Sigma = ((3.0 - 19.0 / 6.0 / sW2 + 1.0 / 6.0 / cW2) * s
            + (4.0 + 1.0 / cW2 - 1.0 / sW2) * Mz2) * DeltaW * sW2 * cW2
            + ((-cW2 * cW2 * (40.0 * s + 80.0 * w)
            + (cW2 - sW2)*(cW2 - sW2)*(8.0 * w + s) + 12.0 * w)
            * F_Hollik(s, Mw, Mw)
            + (10.0 * z - 2.0 * h + s + (h - z)*(h - z) / s)
            * F_Hollik(s, mh, Mz)
            - 2.0 * h * log(h / w) - 2.0 * z * log(z / w)
            + (10.0 * z - 2.0 * h + s)
            *(1.0 - (h + z) / (h - z) * log(mh / Mz) - log(mh * Mz / w))
            + 2.0 / 3.0 * s * (1.0 + (cW2 - sW2)*(cW2 - sW2) - 4.0 * cW2)
            ) / 12.0;
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaGammaGamma_bos_Hollik(const double mu, const double s,
        const double Mw) const
{
    double Mw2 = Mw*Mw;
    double w = Mw2;

    if (mu <= 0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaGammaGamma_bos_Hollik()");

    double DeltaW = -log(Mw2 / mu / mu);

    gslpp::complex Sigma = -3.0 * s * DeltaW - (3.0 * s + 4.0 * w) * F_Hollik(s, Mw, Mw);
    return Sigma;
}

gslpp::complex EWSMOneLoopEW_HV::PiGammaGamma_bos_Hollik(const double mu, const double s,
        const double Mw) const
{
    double Mw2 = Mw*Mw;
    double w = Mw2;

    if (mu <= 0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::PiGammaGamma_bos_Hollik()");

    double DeltaW = -log(Mw2 / mu / mu);
    double muIR = Mw; // relevant only for Fprime_Hollik(muIR, Mw*Mw, 0, Mw)

    gslpp::complex Pi;
    if (s == 0.0) {
        //-- Pi(0) = dSigma(s)/ds|_{s=0} --
        Pi = -3.0 * DeltaW - 3.0 * F_Hollik(s, Mw, Mw)
                - (3.0 * s + 4.0 * w) * Fprime_Hollik(muIR, s, Mw, Mw);
    } else {
        //-- Pi(s) = Sigma(s)/s --
        Pi = SigmaGammaGamma_bos_Hollik(mu, s, Mw) / s;
    }
    return Pi;
}

gslpp::complex EWSMOneLoopEW_HV::SigmaZgamma_bos_Hollik(const double mu, const double s,
        const double Mw) const
{
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double cW2 = Mw2 / Mz2;
    double w = Mw2;

    if (mu <= 0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaZgamma_bos_Hollik()");

    double DeltaW = -log(Mw2 / mu / mu);

    gslpp::complex Sigma = ((3.0 * cW2 + 1.0 / 6.0) * s + 2.0 * w) * DeltaW
            + ((3.0 * cW2 + 1.0 / 6.0) * s + (4.0 * cW2 + 4.0 / 3.0) * w)
            * F_Hollik(s, Mw, Mw)
            + s / 9.0;
    return ( -Sigma); // The minus sign is attributed to the different definition of s_W.
}




