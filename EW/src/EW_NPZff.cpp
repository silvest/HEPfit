/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "EW_NPZff.h"
#include "Rbottom.h"

EW_NPZff::EW_NPZff(const StandardModel& SM_i)
: SM(SM_i)
{
}

double EW_NPZff::GammaZ(const double GammaZ_SM) const
{
    double Gamma_Z = GammaZ_SM;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p=0; p<6; ++p) {
        delGVl[p] = SM.deltaGVl((StandardModel::lepton)p);
        delGAl[p] = SM.deltaGAl((StandardModel::lepton)p);
        delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
        delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
        if (delGVl[p]!=0.0 || delGAl[p]!=0.0
                || delGVq[p]!=0.0 || delGAq[p]!=0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double deltaGl[6], deltaGq[6];
        double delGammaZ = 0.0;
        for (int p=0; p<6; ++p) {
            gVf = SM.StandardModel::gVl((StandardModel::lepton)p).real();
            gAf = SM.StandardModel::gAl((StandardModel::lepton)p).real();
            deltaGl[p] = 2.0*(gVf*delGVl[p] + gAf*delGAl[p]);

            gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
            gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
            deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

            delGammaZ += deltaGl[p] + 3.0*deltaGq[p];
        }

        double sW2_SM = SM.StandardModel::sW2();
        double cW2_SM = SM.StandardModel::cW2();
        Gamma_Z += SM.alphaMz()*SM.getMz()/12.0/sW2_SM/cW2_SM
                   * delGammaZ;
    }
    return Gamma_Z;
}

double EW_NPZff::sigmaHadron(const double sigmaHadron_SM) const
{
    double sigma_had = sigmaHadron_SM;
    bool nonZeroNP = false;

    double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
    for (int p=0; p<6; ++p) {
        delGVl[p] = SM.deltaGVl((StandardModel::lepton)p);
        delGAl[p] = SM.deltaGAl((StandardModel::lepton)p);
        delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
        delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
        if (delGVl[p]!=0.0 || delGAl[p]!=0.0
                || delGVq[p]!=0.0 || delGAq[p]!=0.0)
            nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double Gl[6], deltaGl[6], Gq[6], deltaGq[6];
        double Gq_sum = 0.0, delGq_sum = 0.0;
        double Gf_sum = 0.0, delGf_sum = 0.0;
        for (int p=0; p<6; ++p) {
            gVf = SM.StandardModel::gVl((StandardModel::lepton)p).real();
            gAf = SM.StandardModel::gAl((StandardModel::lepton)p).real();
            Gl[p] = gVf*gVf + gAf*gAf;
            deltaGl[p] = 2.0*(gVf*delGVl[p] + gAf*delGAl[p]);

            gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
            gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
            Gq[p] = gVf*gVf + gAf*gAf;
            deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

            Gq_sum += 3.0*Gq[p];
            Gf_sum += Gl[p] + 3.0*Gq[p];
            delGq_sum += 3.0*deltaGq[p];
            delGf_sum += deltaGl[p] + 3.0*deltaGq[p];
        }

        sigma_had += 12.0*M_PI/SM.getMz()/SM.getMz()
                     *Gl[(int)SM.ELECTRON]*Gq_sum/Gf_sum/Gf_sum
                     *( deltaGl[(int)SM.ELECTRON]/Gl[(int)SM.ELECTRON]
                        + delGq_sum/Gq_sum - 2.0*delGf_sum/Gf_sum );
    }
    return sigma_had;
}

double EW_NPZff::sin2thetaEff(const double sin2thetaEff_SM) const
{
    double sin2_theta_eff = sin2thetaEff_SM;
    double delGVf = SM.deltaGVl(SM.ELECTRON);
    double delGAf = SM.deltaGAl(SM.ELECTRON);
    if (delGVf!=0.0 || delGAf!=0.0) {
        double gVf = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAf = SM.StandardModel::gAl(SM.ELECTRON).real();
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        sin2_theta_eff -= delGVfOverGAf/4.0;
    }
    return sin2_theta_eff;
}

double EW_NPZff::PtauPol(const double PtauPol_SM) const
{
    double P_tau_pol = PtauPol_SM;
    double delGVf = SM.deltaGVl(SM.TAU);
    double delGAf = SM.deltaGAl(SM.TAU);
    if (delGVf!=0.0 || delGAf!=0.0) {
        double gVf = SM.StandardModel::gVl(SM.TAU).real();
        double gAf = SM.StandardModel::gAl(SM.TAU).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        P_tau_pol -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
    }
    return P_tau_pol;
}

double EW_NPZff::Alepton(const double Alepton_SM) const
{
    double A_l = Alepton_SM;
    double delGVf = SM.deltaGVl(SM.ELECTRON);
    double delGAf = SM.deltaGAl(SM.ELECTRON);
    if (delGVf!=0.0 || delGAf!=0.0) {
        double gVf = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAf = SM.StandardModel::gAl(SM.ELECTRON).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        A_l -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
    }
    return A_l;
}

double EW_NPZff::Acharm(const double Acharm_SM) const
{
    double A_c = Acharm_SM;
    double delGVf = SM.deltaGVq(SM.CHARM);
    double delGAf = SM.deltaGAq(SM.CHARM);
    if (delGVf!=0.0 || delGAf!=0.0) {
        double gVf = SM.StandardModel::gVq(SM.CHARM).real();
        double gAf = SM.StandardModel::gAq(SM.CHARM).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        A_c -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
    }
    return A_c;
}

double EW_NPZff::Abottom(const double Abottom_SM) const
{
    double A_b = Abottom_SM;
    double delGVf = SM.deltaGVq(SM.BOTTOM);
    double delGAf = SM.deltaGAq(SM.BOTTOM);
    if (delGVf!=0.0 || delGAf!=0.0) {
        double gVf = SM.StandardModel::gVq(SM.BOTTOM).real();
        double gAf = SM.StandardModel::gAq(SM.BOTTOM).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        A_b -= 2.0*(gVf*gVf - gAf*gAf)*gAf*gAf/Gf/Gf*delGVfOverGAf;
    }
    return A_b;
}

double EW_NPZff::AFBlepton(const double AFBlepton_SM) const
{
    double AFB_l = AFBlepton_SM;
    double delGVe = SM.deltaGVl(SM.ELECTRON);
    double delGAe = SM.deltaGAl(SM.ELECTRON);
    if (delGVe!=0.0 || delGAe!=0.0) {
        double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
        double Ge = gVe*gVe + gAe*gAe;
        double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;

        AFB_l -= 6.0*gVe*gAe*(gVe*gVe - gAe*gAe)*gAe*gAe/Ge/Ge/Ge*delGVeOverGAe;
    }
    return AFB_l;
}

double EW_NPZff::AFBcharm(const double AFBcharm_SM) const
{
    double AFB_c = AFBcharm_SM;
    double delGVe = SM.deltaGVl(SM.ELECTRON);
    double delGAe = SM.deltaGAl(SM.ELECTRON);
    double delGVf = SM.deltaGVq(SM.CHARM);
    double delGAf = SM.deltaGAq(SM.CHARM);
    if (delGVe!=0.0 || delGAe!=0.0 || delGVf!=0.0 || delGAf!=0.0) {
        double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
        double Ge = gVe*gVe + gAe*gAe;
        double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;
        //
        double gVf = SM.StandardModel::gVq(SM.CHARM).real();
        double gAf = SM.StandardModel::gAq(SM.CHARM).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        AFB_c -= 3.0*gVf*gAf*(gVe*gVe - gAe*gAe)*gAe*gAe/Gf/Ge/Ge*delGVeOverGAe
                 + 3.0*gVe*gAe*(gVf*gVf - gAf*gAf)*gAf*gAf/Ge/Gf/Gf*delGVfOverGAf;
    } 
    return AFB_c;
}

double EW_NPZff::AFBbottom(const double AFBbottom_SM) const
{
    double AFB_b = AFBbottom_SM;
    double delGVe = SM.deltaGVl(SM.ELECTRON);
    double delGAe = SM.deltaGAl(SM.ELECTRON);
    double delGVf = SM.deltaGVq(SM.BOTTOM);
    double delGAf = SM.deltaGAq(SM.BOTTOM);
    if (delGVe!=0.0 || delGAe!=0.0 || delGVf!=0.0 || delGAf!=0.0) {
        double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
        double Ge = gVe*gVe + gAe*gAe;
        double delGVeOverGAe = (gAe*delGVe - gVe*delGAe)/gAe/gAe;
        //
        double gVf = SM.StandardModel::gVq(SM.BOTTOM).real();
        double gAf = SM.StandardModel::gAq(SM.BOTTOM).real();
        double Gf = gVf*gVf + gAf*gAf;
        double delGVfOverGAf = (gAf*delGVf - gVf*delGAf)/gAf/gAf;

        AFB_b -= 3.0*gVf*gAf*(gVe*gVe - gAe*gAe)*gAe*gAe/Gf/Ge/Ge*delGVeOverGAe
                 + 3.0*gVe*gAe*(gVf*gVf - gAf*gAf)*gAf*gAf/Ge/Gf/Gf*delGVfOverGAf;
    }
    return AFB_b;
}

double EW_NPZff::Rlepton(const double Rlepton_SM) const
{
    double R0_l = Rlepton_SM;
    bool nonZeroNP = false;

    double delGVe = SM.deltaGVl(SM.ELECTRON);
    double delGAe = SM.deltaGAl(SM.ELECTRON);
    if (delGVe!=0.0 || delGAe!=0.0) nonZeroNP = true;

    double delGVq[6], delGAq[6];
    for (int p=0; p<6; ++p) {
        delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
        delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
        if (delGVq[p]!=0.0 || delGAq[p]!=0.0) nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
        double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
        double Ge = gVe*gVe + gAe*gAe;
        double deltaGe = 2.0*(gVe*delGVe + gAe*delGAe);

        double Gq[6], deltaGq[6];
        double gVq, gAq;
        double Gq_sum = 0.0, delGq_sum = 0.0;
        for (int p=0; p<6; ++p) {
            gVq = SM.StandardModel::gVq((StandardModel::quark)p).real();
            gAq = SM.StandardModel::gAq((StandardModel::quark)p).real();
            Gq[p] = gVq*gVq + gAq*gAq;
            deltaGq[p] = 2.0*(gVq*delGVq[p] + gAq*delGAq[p]);

            Gq_sum += 3.0*Gq[p];
            delGq_sum += 3.0*deltaGq[p];
        }

        R0_l += delGq_sum/Ge - Gq_sum*deltaGe/Ge/Ge;
    }
    return R0_l;
}

double EW_NPZff::Rcharm(const double Rcharm_SM) const
{
    double R0_c = Rcharm_SM;
    bool nonZeroNP = false;
    double delGVq[6], delGAq[6];
    for (int p=0; p<6; ++p) {
        delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
        delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
        if (delGVq[p]!=0.0 || delGAq[p]!=0.0) nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double Gq[6], deltaGq[6];
        double Gq_sum = 0.0, delGq_sum = 0.0;
        for (int p=0; p<6; ++p) {
            gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
            gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
            Gq[p] = gVf*gVf + gAf*gAf;
            deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

            Gq_sum += Gq[p]; /* without the color factor */
            delGq_sum += deltaGq[p]; /* without the color factor */
        }

        R0_c += deltaGq[(int)SM.CHARM]/Gq_sum
                - Gq[(int)SM.CHARM]*delGq_sum/Gq_sum/Gq_sum;
    }
    return R0_c;
}

double EW_NPZff::Rbottom(const double Rbottom_SM) const
{
    double R0_b = Rbottom_SM;
    bool nonZeroNP = false;
    double delGVq[6], delGAq[6];
    for (int p=0; p<6; ++p) {
        delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
        delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
        if (delGVq[p]!=0.0 || delGAq[p]!=0.0) nonZeroNP = true;
    }

    if (nonZeroNP) {
        double gVf, gAf;
        double Gq[6], deltaGq[6];
        double Gq_sum = 0.0, delGq_sum = 0.0;
        for (int p=0; p<6; ++p) {
            gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
            gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
            Gq[p] = gVf*gVf + gAf*gAf;
            deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

            Gq_sum += Gq[p]; /* without the color factor */
            delGq_sum += deltaGq[p]; /* without the color factor */
        }

        R0_b += deltaGq[(int)SM.BOTTOM]/Gq_sum
                - Gq[(int)SM.BOTTOM]*delGq_sum/Gq_sum/Gq_sum;
    }
    return R0_b;
}