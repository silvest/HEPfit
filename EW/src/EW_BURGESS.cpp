/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EW_BURGESS.h"

double EW_BURGESS::Mw(const double Mw_SM) const
{
    return ( Mw_SM * (1.0 - 0.00723/2.0*SM.obliqueS() + 0.0111/2.0*SM.obliqueT()
                      + 0.00849/2.0*SM.obliqueU()) );
}

double EW_BURGESS::GammaW(const double GammaW_SM) const
{
    return ( GammaW_SM * (1.0 - 0.00723*SM.obliqueS() + 0.0111*SM.obliqueT()
                          + 0.00849*SM.obliqueU()) );
}

double EW_BURGESS::GammaZ(const double GammaZ_SM) const
{
    return ( GammaZ_SM - 0.00961*SM.obliqueS() + 0.0263*SM.obliqueT() );
}

double EW_BURGESS::sigmaHadron(const double sigmaHadron_SM, const double GammaZ_SM,
                               const double GammaHad_SM, const double Gamma_l_SM) const
{
    double delta_l = - 0.000192*SM.obliqueS() + 0.000790*SM.obliqueT();
    double delta_had = - 0.00901*SM.obliqueS() + 0.0200*SM.obliqueT();
    double delta_Z = - 0.00961*SM.obliqueS() + 0.0263*SM.obliqueT();
    return ( sigmaHadron_SM
             * (1.0 + delta_l/Gamma_l_SM + delta_had/GammaHad_SM
                - 2.0*delta_Z/GammaZ_SM) );
}

double EW_BURGESS::sin2thetaEff(const double sin2thetaEff_SM) const
{
    return ( sin2thetaEff_SM + 0.00362*SM.obliqueS() - 0.00256*SM.obliqueT() );
}

double EW_BURGESS::PtauPol(const double PtauPol_SM) const
{
    return ( PtauPol_SM - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT() );
}

double EW_BURGESS::Alepton(const double Alepton_SM) const
{
    return ( Alepton_SM - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT() );
}

double EW_BURGESS::Acharm(const double Acharm_SM, const double Alepton_SM) const
{
    double AFB_c = 3.0/4.0*Alepton_SM*Acharm_SM;
    double delta_AFB_c = - 0.0147*SM.obliqueS() + 0.0104*SM.obliqueT();
    double delta_A_l = - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT();
    return ( Acharm_SM * (1.0 + delta_AFB_c/AFB_c - delta_A_l/Alepton_SM) );
}

double EW_BURGESS::Abottom(const double Abottom_SM, const double Alepton_SM) const
{
    double AFB_b = 3.0/4.0*Alepton_SM*Abottom_SM;
    double delta_AFB_b = - 0.0188*SM.obliqueS() + 0.0131*SM.obliqueT();
    double delta_A_l = - 0.0284*SM.obliqueS() + 0.0201*SM.obliqueT();
    return ( Abottom_SM * (1.0 + delta_AFB_b/AFB_b - delta_A_l/Alepton_SM) );
}

double EW_BURGESS::AFBlepton(const double AFBlepton_SM) const
{
    return ( AFBlepton_SM - 0.00677*SM.obliqueS() + 0.00479*SM.obliqueT() );
}

double EW_BURGESS::AFBcharm(const double AFBcharm_SM) const
{
    return ( AFBcharm_SM  - 0.0147*SM.obliqueS() + 0.0104*SM.obliqueT() );
}

double EW_BURGESS::AFBbottom(const double AFBbottom_SM) const
{
    return ( AFBbottom_SM - 0.0188*SM.obliqueS() + 0.0131*SM.obliqueT() );
}

double EW_BURGESS::Rlepton(const double Rlepton_SM, const double GammaHad_SM,
                           const double Gamma_l_SM) const
{
    double delta_had = - 0.00901*SM.obliqueS() + 0.0200*SM.obliqueT();
    double delta_l = - 0.000192*SM.obliqueS() + 0.000790*SM.obliqueT();
    return ( Rlepton_SM * (1.0 + delta_had/GammaHad_SM - delta_l/Gamma_l_SM) );
}

double EW_BURGESS::Rcharm(const double Rcharm_SM, const double GammaHad_SM) const
{
    double delta_c_over_Gamma_c = - 0.00649*SM.obliqueS() + 0.0124*SM.obliqueT();
    double delta_had = - 0.00901*SM.obliqueS() + 0.0200*SM.obliqueT();
    return ( Rcharm_SM * (1.0 + delta_c_over_Gamma_c - delta_had/GammaHad_SM) );
}

double EW_BURGESS::Rbottom(const double Rbottom_SM, const double GammaHad_SM,
                           const double Gamma_b_SM) const
{
    double delta_b = - 0.00171*SM.obliqueS() + 0.00416*SM.obliqueT();
    double delta_had = - 0.00901*SM.obliqueS() + 0.0200*SM.obliqueT();
    return ( Rbottom_SM * (1.0 + delta_b/Gamma_b_SM - delta_had/GammaHad_SM) );
}

