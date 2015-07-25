/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_BURGESS.h"

EW_BURGESS::EW_BURGESS(const NPSTU& NP_i)
: NP(NP_i)
{
}

double EW_BURGESS::Mw(const double Mw_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();
    double obliqueU = NP.obliqueU();

    return ( Mw_SM * (1.0 - 0.00723 / 2.0 * obliqueS + 0.0111 / 2.0 * obliqueT
            + 0.00849 / 2.0 * obliqueU));
}

double EW_BURGESS::GammaW(const double GammaW_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();
    double obliqueU = NP.obliqueU();

    return ( GammaW_SM * (1.0 - 0.00723 * obliqueS + 0.0111 * obliqueT
            + 0.00849 * obliqueU));
}

double EW_BURGESS::GammaZ(const double GammaZ_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( GammaZ_SM - 0.00961 * obliqueS + 0.0263 * obliqueT);
}

double EW_BURGESS::sigmaHadron(const double sigmaHadron_SM, const double GammaZ_SM,
        const double GammaHad_SM, const double Gamma_l_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double delta_l = -0.000192 * obliqueS + 0.000790 * obliqueT;
    double delta_had = -0.00901 * obliqueS + 0.0200 * obliqueT;
    double delta_Z = -0.00961 * obliqueS + 0.0263 * obliqueT;
    return ( sigmaHadron_SM
            * (1.0 + delta_l / Gamma_l_SM + delta_had / GammaHad_SM
            - 2.0 * delta_Z / GammaZ_SM));
}

double EW_BURGESS::sin2thetaEff(const double sin2thetaEff_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( sin2thetaEff_SM + 0.00362 * obliqueS - 0.00256 * obliqueT);
}

double EW_BURGESS::PtauPol(const double PtauPol_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( PtauPol_SM - 0.0284 * obliqueS + 0.0201 * obliqueT);
}

double EW_BURGESS::Alepton(const double Alepton_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( Alepton_SM - 0.0284 * obliqueS + 0.0201 * obliqueT);
}

double EW_BURGESS::Acharm(const double Acharm_SM, const double Alepton_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double AFB_c = 3.0 / 4.0 * Alepton_SM*Acharm_SM;
    double delta_AFB_c = -0.0147 * obliqueS + 0.0104 * obliqueT;
    double delta_A_l = -0.0284 * obliqueS + 0.0201 * obliqueT;
    return ( Acharm_SM * (1.0 + delta_AFB_c / AFB_c - delta_A_l / Alepton_SM));
}

double EW_BURGESS::Abottom(const double Abottom_SM, const double Alepton_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double AFB_b = 3.0 / 4.0 * Alepton_SM*Abottom_SM;
    double delta_AFB_b = -0.0188 * obliqueS + 0.0131 * obliqueT;
    double delta_A_l = -0.0284 * obliqueS + 0.0201 * obliqueT;
    return ( Abottom_SM * (1.0 + delta_AFB_b / AFB_b - delta_A_l / Alepton_SM));
}

double EW_BURGESS::AFBlepton(const double AFBlepton_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( AFBlepton_SM - 0.00677 * obliqueS + 0.00479 * obliqueT);
}

double EW_BURGESS::AFBcharm(const double AFBcharm_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( AFBcharm_SM - 0.0147 * obliqueS + 0.0104 * obliqueT);
}

double EW_BURGESS::AFBbottom(const double AFBbottom_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    return ( AFBbottom_SM - 0.0188 * obliqueS + 0.0131 * obliqueT);
}

double EW_BURGESS::Rlepton(const double Rlepton_SM, const double GammaHad_SM,
        const double Gamma_l_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double delta_had = -0.00901 * obliqueS + 0.0200 * obliqueT;
    double delta_l = -0.000192 * obliqueS + 0.000790 * obliqueT;
    return ( Rlepton_SM * (1.0 + delta_had / GammaHad_SM - delta_l / Gamma_l_SM));
}

double EW_BURGESS::Rcharm(const double Rcharm_SM, const double GammaHad_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double delta_c_over_Gamma_c = -0.00649 * obliqueS + 0.0124 * obliqueT;
    double delta_had = -0.00901 * obliqueS + 0.0200 * obliqueT;
    return ( Rcharm_SM * (1.0 + delta_c_over_Gamma_c - delta_had / GammaHad_SM));
}

double EW_BURGESS::Rbottom(const double Rbottom_SM, const double GammaHad_SM,
        const double Gamma_b_SM) const
{
    double obliqueS = NP.obliqueS();
    double obliqueT = NP.obliqueT();

    double delta_b = -0.00171 * obliqueS + 0.00416 * obliqueT;
    double delta_had = -0.00901 * obliqueS + 0.0200 * obliqueT;
    return ( Rbottom_SM * (1.0 + delta_b / Gamma_b_SM - delta_had / GammaHad_SM));
}

