/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include "FeynHiggsWrapper.h"
#include "EWSUSY.h"

const double EWSUSY::Mw_unphysical = 2.0;
const double EWSUSY::RenormalizationScaleFactor = 1.0;
//const double EWSUSY::RenormalizationScaleFactor = 2.0; // for debug

EWSUSY::EWSUSY(const SUSY& SUSY_in)
: EWSM(SUSY_in), PV(true), mySUSY(SUSY_in),
        Yu(3,3,0.0), Yd(3,3,0.0), Yl(3,3,0.0),
        Au(3,3,0.0), Ad(3,3,0.0), Al(3,3,0.0),
        Zm(2,2,0.), Zp(2,2,0.), ZN(4,4,0.),
        ZU(6,6,0.), ZD(6,6,0.), ZL(6,6,0.), Zne(6,6,0.),
        ZR(2,2,0.), ZH(2,2,0.)
{
}

void EWSUSY::SetRosiekParameters()
{
    Yu = mySUSY.getYu();
    Yd = - mySUSY.getYd();
    Yl = - mySUSY.getYe();

    Au = - mySUSY.getTUhat().transpose();
    Ad = mySUSY.getTDhat().transpose();
    Al = mySUSY.getTEhat().transpose();

    Zm = mySUSY.getU().hconjugate();
    Zp = mySUSY.getV().hconjugate();
    ZN = mySUSY.getN().hconjugate();

    ZU = mySUSY.getRu().hconjugate();
    ZD = mySUSY.getRd().transpose();
    ZL = mySUSY.getRl().transpose();
    Zne = mySUSY.getRn().hconjugate();

    double sinAlpha = mySUSY.getSaeff().real(); /* Correct? */
    double cosAlpha = sqrt(1.0 - sinAlpha*sinAlpha); /* -Pi/2 < alpha < 0 */
    ZR.assign(0,0, cosAlpha);
    ZR.assign(0,1, - sinAlpha);
    ZR.assign(1,0, sinAlpha);
    ZR.assign(1,1, cosAlpha);

    ZH.assign(0,0, mySUSY.getSinb());
    ZH.assign(0,1, - mySUSY.getCosb());
    ZH.assign(1,0, mySUSY.getCosb());
    ZH.assign(1,1, mySUSY.getSinb());

    /* particle massses */
    for (int I=0; I<3; ++I) {
        /* up-type quarks */
        m_u[I] = getMyCache()->mq((StandardModel::quark)(2*I), mySUSY.getMz(), FULLNNLO);
        /* down-type quarks */
        m_d[I] = getMyCache()->mq((StandardModel::quark)(2*I + 1), mySUSY.getMz(), FULLNNLO);
        /* charged leptons */
        m_l[I] = mySUSY.getLeptons((StandardModel::lepton)(2*I + 1)).getMass();
    }
    /* H^0_i = (H^0, h^0, A^0, G^0) */
    mH02[0] = mySUSY.getMHh()*mySUSY.getMHh();
    mH02[1] = mySUSY.getMHl()*mySUSY.getMHl();
    mH02[2] = mySUSY.getMHa()*mySUSY.getMHa();
    mH02[3] = mySUSY.getMz()*mySUSY.getMz(); /* mass squared of the neutral Goldstone boson */
    for (int k=0; k<6; ++k) {
        Msu2[k] = mySUSY.getMsu2()(k);
        Msd2[k] = mySUSY.getMsd2()(k);
        Mse2[k] = mySUSY.getMse2()(k);
        Msn2[k] = mySUSY.getMsn2()(k);
    }
    for (int i=0; i<2; ++i)
        mC[i] = mySUSY.getMch()(i);
    for (int j=0; j<4; ++j)
        mN[j] = mySUSY.getMneu()(j);
}

complex EWSUSY::FA(const double mu, const double p2,
                   const double mi, const double mj,
                   const complex cV_aij, const complex cV_bji,
                   const complex cA_aij, const complex cA_bji) const
{
    double mu2 = mu*mu, mi2 = mi*mi, mj2 = mj*mj;

    /* PV functions */
    double A0i = PV.A0(mu2, mi2);
    double A0j = PV.A0(mu2, mj2);
    complex B0 = PV.B0(mu2, p2, mi2, mj2);
    complex B00 = PV.B00(mu2, p2, mi2, mj2);

    return ( -2.0*(cV_aij*cV_bji + cA_aij*cA_bji)
               *(4.0*B00 + A0i + A0j + (p2 - mi*mi - mj*mj)*B0)
             -4.0*(cV_aij*cV_bji - cA_aij*cA_bji)*mi*mj*B0 );
}

complex EWSUSY::dFA(const double mu, const double p2,
                    const double mi, const double mj,
                    const complex cV_aij, const complex cV_bji,
                    const complex cA_aij, const complex cA_bji) const
{
    double mu2 = mu*mu, mi2 = mi*mi, mj2 = mj*mj;

    /* PV functions */
    complex B0 = PV.B0(mu2, p2, mi2, mj2);
    complex B0p = PV.B0p(mu2, p2, mi2, mj2);
    complex B00p = PV.B00p(mu2, p2, mi2, mj2);

    if (mi == mj && cA_aij == 0.0 && cA_bji == 0.0)
        return ( -2.0*cV_aij*cV_bji*(4.0*B00p + p2*B0p + B0) );
    else
        return ( -2.0*(cV_aij*cV_bji + cA_aij*cA_bji)
                  *(4.0*B00p + (p2 - mi*mi - mj*mj)*B0p + B0)
                 -4.0*(cV_aij*cV_bji - cA_aij*cA_bji)*mi*mj*B0p );
}

complex EWSUSY::PiT_Z(const double mu, const double p2, const double Mw_i) const
{
    double mu2 = mu*mu;
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double Mw2 = Mw_i*Mw_i;
    double mHp2[2] = {mySUSY.getMHp()*mySUSY.getMHp(), Mw2}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2sq = e2/sW2; /* g2 squared */
    double e_4sc = e/4.0/sW/cW;
    double e_2sc = 2.0*e_4sc;

    complex PiT_f = complex(0.0, 0.0, false);
    complex PiT_sf = complex(0.0, 0.0, false);
    complex PiT_ch = complex(0.0, 0.0, false);
    complex PiT_WZH = complex(0.0, 0.0, false);
    double a0;
    complex b0, b00;
    complex cV_Zij, cV_Zji, cA_Zij, cA_Zji;
    matrix<double> Id6 = matrix<double>::Id(6);
    matrix<double> Id2 = matrix<double>::Id(2);

    /* neutrino loops */
    b0 = PV.B0(mu2, p2, 0.0, 0.0);
    b00 = PV.B00(mu2, p2, 0.0, 0.0);
    PiT_f += - 3.0/4.0*g2sq/cW2*(4.0*b00 + p2*b0);

    /* other SM fermion loops */
    complex cV_Zee = - e_4sc*(1.0 - 4.0*sW2);
    complex cA_Zee = - e_4sc;
    complex cV_Zdd = - e_4sc*(1.0 - 4.0/3.0*sW2);
    complex cA_Zdd = - e_4sc;
    complex cV_Zuu = e_4sc*(1.0 - 8.0/3.0*sW2);
    complex cA_Zuu = e_4sc;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        PiT_f += FA(mu, p2, m_l[I], m_l[I], cV_Zee, cV_Zee, cA_Zee, cA_Zee);

        /* down-type quarks */
        PiT_f += Nc*FA(mu, p2, m_d[I], m_d[I], cV_Zdd, cV_Zdd, cA_Zdd, cA_Zdd);

        /* up-type quarks */
        PiT_f += Nc*FA(mu, p2, m_u[I], m_u[I], cV_Zuu, cV_Zuu, cA_Zuu, cA_Zuu);
    }

    /* sneutrino loops */
    complex VZsnsn_II = e_2sc;
    complex VZZsnsn_II = e2/2.0/sW2/cW2;
    for (int I=0; I<3; ++I) {  /* I=0-2 for left-handed sneutrinos */
        b00 = PV.B00(mu2, p2, Msn2[I], Msn2[I]);
        PiT_sf += 4.0*VZsnsn_II.abs2()*b00;
        a0 = PV.A0(mu2, Msn2[I]);
        PiT_sf += VZZsnsn_II*a0;
    }

    /* charged-slepton loops */
    complex VZLL_mn, VZZLL_nn;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VZLL_mn = complex(0.0, 0.0, false);
            for (int I=0; I<3; ++I) /* sum over left-handed sleptons */
                VZLL_mn += - e_2sc*ZL(I,n)*ZL(I,m).conjugate();
            VZLL_mn += - e_2sc*(- 2.0*sW2*Id6(m,n));
            b00 = PV.B00(mu2, p2, Mse2[m], Mse2[n]);
            PiT_sf += 4.0*VZLL_mn.abs2()*b00;
        }
        VZZLL_nn = complex(0.0, 0.0, false);
        VZZLL_nn += 2.0*e2/cW2*sW2;
        for (int I=0; I<3; ++I) /* sum over left-handed sleptons */
            VZZLL_nn += 2.0*e2/cW2*(1.0 - 4.0*sW2)/4.0/sW2*ZL(I,n)*ZL(I,n).conjugate();
        a0 = PV.A0(mu2, Mse2[n]);
        PiT_sf += VZZLL_nn*a0;
    }

    /* down-type squark loops */
    complex VZDD_mn, VZZDD_nn;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VZDD_mn = complex(0.0, 0.0, false);
            for (int I=0; I<3; ++I) /* sum over left-handed squarks */
                VZDD_mn += - e_2sc*ZD(I,n)*ZD(I,m).conjugate();
            VZDD_mn += - e_2sc*(- 2.0/3.0*sW2*Id6(m,n));
            b00 = PV.B00(mu2, p2, Msd2[m], Msd2[n]);
            PiT_sf += 4.0*Nc*VZDD_mn.abs2()*b00;
        }
        VZZDD_nn = complex(0.0, 0.0, false);
        VZZDD_nn += 2.0*e2/3.0/cW2*sW2/3.0;
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VZZDD_nn += 2.0*e2/3.0/cW2*(3.0 - 4.0*sW2)/4.0/sW2*ZD(I,n)*ZD(I,n).conjugate();
        a0 = PV.A0(mu2, Msd2[n]);
        PiT_sf += Nc*VZZDD_nn*a0;
    }

    /* up-type squark loops */
    complex VZUU_mn, VZZUU_nn;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VZUU_mn = complex(0.0, 0.0, false);
            for (int I=0; I<3; ++I) /* sum over left-handed squarks */
                VZUU_mn += e_2sc*ZU(I,m).conjugate()*ZU(I,n);
            VZUU_mn += e_2sc*(- 4.0/3.0*sW2*Id6(m,n));
            b00 = PV.B00(mu2, p2, Msu2[m], Msu2[n]);
            PiT_sf += 4.0*Nc*VZUU_mn.abs2()*b00;
        }
        VZZUU_nn = complex(0.0, 0.0, false);
        VZZUU_nn += 2.0*e2/3.0/cW2*4.0*sW2/3.0;
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VZZUU_nn += 2.0*e2/3.0/cW2*(3.0 - 8.0*sW2)/4.0/sW2*ZU(I,n).conjugate()*ZU(I,n);
        a0 = PV.A0(mu2, Msu2[n]);
        PiT_sf += Nc*VZZUU_nn*a0;
    }

    /* chargino loops */
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j) {
            cV_Zij = e_4sc*(  Zp(0,j).conjugate()*Zp(0,i)
                            + Zm(0,j)*Zm(0,i).conjugate()
                            + 2.0*(cW2 - sW2)*Id2(j,i) );
            cV_Zji = cV_Zij.conjugate();
            cA_Zij = e_4sc*(  Zp(0,j).conjugate()*Zp(0,i)
                            - Zm(0,j)*Zm(0,i).conjugate() );
            cA_Zji = cA_Zij.conjugate();
            PiT_ch += FA(mu, p2, mC[i], mC[j], cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }

    /* neutralino loops */
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j) {
            cV_Zij = - e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                              - ZN(2,j).conjugate()*ZN(2,i)
                              - ZN(3,j)*ZN(3,i).conjugate()
                              + ZN(2,j)*ZN(2,i).conjugate() );
            cV_Zji = cV_Zij.conjugate();
            cA_Zij = - e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                              - ZN(2,j).conjugate()*ZN(2,i)
                              + ZN(3,j)*ZN(3,i).conjugate()
                              - ZN(2,j)*ZN(2,i).conjugate() );
            cA_Zji = cA_Zij.conjugate();
            PiT_ch += 0.5*FA(mu, p2, mN[i], mN[j], cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }

    /* charged-Higgs loops */
    double cot_2thW = (cW2 - sW2)/(2.0*sW*cW);
    for (int i=0; i<2; ++i) {
        b00 = PV.B00(mu2, p2, mHp2[i], mHp2[i]);
        a0 = PV.A0(mu2, mHp2[i]);
        PiT_WZH += 2.0*e2*cot_2thW*cot_2thW*(2.0*b00 + a0);
    }

    /* neutral-Higgs loops */
    double AM_ij;
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j) {
            AM_ij = ZR(0,i)*ZH(0,j) - ZR(1,i)*ZH(1,j);
            b00 = PV.B00(mu2, p2, mH02[i], mH02[j+2]);
            PiT_WZH += g2sq/cW2*AM_ij*AM_ij*b00;
        }
    for (int j=0; j<4; ++j) {
        a0 = PV.A0(mu2, mH02[j]);
        PiT_WZH += g2sq/4.0/cW2*a0;
    }

    /* W-boson - charged-Goldstone-boson loop*/
    b0 = PV.B0(mu2, p2, Mw2, Mw2);
    PiT_WZH += - 2.0*g2sq*sW2*sW2*Mz*Mz*b0;

    /* Z-boson - Higgs loops */
    double CR_i;
    for (int i=0; i<2; ++i) {
        CR_i = mySUSY.v1()*ZR(0,i) + mySUSY.v2()*ZR(1,i);
        b0 = PV.B0(mu2, p2, Mz*Mz, mH02[i]);
        /* Mw^2/v^2 is substituted for g2^2/4 compared to the expression in the
         * paper, in order to ensure the cancellation of the UV divergences in
         * the case where Mw is not the tree-level value. */
        PiT_WZH += - g2sq*Mw2/mySUSY.v()/mySUSY.v()/cW2/cW2*CR_i*CR_i*b0;
    }

    /* W-boson loops */
    a0 = PV.A0(mu2, Mw2);
    b0 = PV.B0(mu2, p2, Mw2, Mw2);
    b00 = PV.B00(mu2, p2, Mw2, Mw2);
    /* typo in the paper: a0^2 --> a0 in the first term */
    PiT_WZH += 2.0*g2sq*cW2*(2.0*a0 + (2.0*p2 + Mw2)*b0 + 4.0*b00);

    /* Sum of all contributions */
    complex PiT = PiT_f + PiT_sf + PiT_ch + PiT_WZH;

    return ( PiT/16.0/M_PI/M_PI );
}


complex EWSUSY::PiT_W(const double mu, const double p2, const double Mw_i) const
{
    double mu2 = mu*mu;
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double Mw2 = Mw_i*Mw_i;
    double mHp2[2] = {mySUSY.getMHp()*mySUSY.getMHp(), Mw2}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2sq = e2/sW2; /* g2 squared */
    double e_2s = e/2.0/sW;
    double e_sq2s = e/sqrt(2.0)/sW;
    double e_2sq2s = e_sq2s/2.0;
    double e2_2s2 = e2/2.0/sW2;

    complex PiT_f = complex(0.0, 0.0, false);
    complex PiT_sf = complex(0.0, 0.0, false);
    complex PiT_ch = complex(0.0, 0.0, false);
    complex PiT_WZH = complex(0.0, 0.0, false);
    double a0;
    complex b0, b00;

    /* SM fermion loops */
    complex cV_Wen = e_2sq2s;
    complex cA_Wen = e_2sq2s;
    complex cV_Wne = cV_Wen.conjugate();
    complex cA_Wne = cA_Wen.conjugate();
    complex cV_Wdu = e_2sq2s; /* no CKM */
    complex cA_Wdu = e_2sq2s; /* no CKM */
    complex cV_Wud = cV_Wdu.conjugate();
    complex cA_Wud = cA_Wdu.conjugate();
    for (int I=0; I<3; ++I) {
        /* leptons */
        PiT_f += FA(mu, p2, m_l[I], 0.0, cV_Wen, cV_Wne, cA_Wen, cA_Wne);

        /* quarks (no CKM) */
        PiT_f += Nc*FA(mu, p2, m_d[I], m_u[I], cV_Wdu, cV_Wud, cA_Wdu, cA_Wud);
    }

    /* slepton loops */
    complex VWsnL_In, VWWsnsn_II, VWWLL_nn;
    for (int I=0; I<3; ++I) {  /* I=0-2 for left-handed sneutrinos */
        for (int n=0; n<6; ++n) {
            VWsnL_In = complex(0.0, 0.0, false);
            for (int J=0; J<3; ++J) /* sum over left-handed sleptons */
                VWsnL_In += e_sq2s*Zne(J,I)*ZL(J,n);
            b00 = PV.B00(mu2, p2, Msn2[I], Mse2[n]);
            PiT_sf += 4.0*VWsnL_In.abs2()*b00;
        }
        VWWsnsn_II = e2_2s2;
        a0 = PV.A0(mu2, Msn2[I]);
        PiT_sf += VWWsnsn_II*a0;
    }
    for (int n=0; n<6; ++n) {
        VWWLL_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed sleptons */
            VWWLL_nn += e2_2s2*ZL(I,n)*ZL(I,n).conjugate();
        a0 = PV.A0(mu2, Mse2[n]);
        PiT_sf += VWWLL_nn*a0;
    }

    /* squark loops (no CKM) */
    complex VWDU_nm, VWWDD_nn, VWWUU_nn;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VWDU_nm = complex(0.0, 0.0, false);
            for (int I=0; I<3; ++I) /* sum over left-handed squarks */
                VWDU_nm += e_sq2s*ZD(I,n).conjugate()*ZU(I,m).conjugate();
            b00 = PV.B00(mu2, p2, Msd2[n], Msu2[m]);
            PiT_sf += 4.0*Nc*VWDU_nm.abs2()*b00;
        }
        VWWDD_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VWWDD_nn += e2_2s2*ZD(I,n)*ZD(I,n).conjugate();
        a0 = PV.A0(mu2, Msd2[n]);
        PiT_sf += Nc*VWWDD_nn*a0;
        VWWUU_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VWWUU_nn += e2_2s2*ZU(I,n).conjugate()*ZU(I,n);
        a0 = PV.A0(mu2, Msu2[n]);
        PiT_sf += Nc*VWWUU_nn*a0;
    }

    /* chargino - neutralino loops */
    complex cV_Wij, cV_Wji, cA_Wij, cA_Wji;
    for (int i=0; i<2; ++i)
        for (int j=0; j<4; ++j) {
            /* W^+ + neutralino(j) -> chi^+(i) */
            /* W^+ + chi^-(i) -> neutralino(j) */
            cV_Wji = - e_2s*(  ZN(1,j)*Zp(0,i).conjugate()
                             - ZN(3,j)*Zp(1,i).conjugate()/sqrt(2.0)
                             + ZN(1,j).conjugate()*Zm(0,i)
                             + ZN(2,j).conjugate()*Zm(1,i)/sqrt(2.0) );
            cV_Wij = cV_Wji.conjugate();
            cA_Wji = - e_2s*(  ZN(1,j)*Zp(0,i).conjugate()
                             - ZN(3,j)*Zp(1,i).conjugate()/sqrt(2.0)
                             - ZN(1,j).conjugate()*Zm(0,i)
                             - ZN(2,j).conjugate()*Zm(1,i)/sqrt(2.0) );
            cA_Wij = cA_Wji.conjugate();
            PiT_ch += FA(mu, p2, mC[i], mN[j], cV_Wij, cV_Wji, cA_Wij, cA_Wji);
        }

    /* Higgs loops */
    double AM_ij;
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            AM_ij = ZR(0,i)*ZH(0,j) - ZR(1,i)*ZH(1,j);
            b00 = PV.B00(mu2, p2, mH02[i], mHp2[j]);
            PiT_WZH += g2sq*AM_ij*AM_ij*b00;
        }
        b00 = PV.B00(mu2, p2, mH02[i+2], mHp2[i]);
        PiT_WZH += g2sq*b00;
    }
    for (int i=0; i<4; ++i) {
        a0 = PV.A0(mu2, mH02[i]);
        PiT_WZH += g2sq/4.0*a0;
    }
    for (int i=0; i<2; ++i) {
        a0 = PV.A0(mu2, mHp2[i]);
        PiT_WZH += g2sq/2.0*a0;
    }

    /* photon - charged-Goldstone-boson loops */
    b0 = PV.B0(mu2, p2, Mw2, 0.0);
    PiT_WZH += - e2*Mw2*b0;

    /* W-boson - Higgs loops */
    double CR_i;
    for (int i=0; i<2; ++i) {
        CR_i = mySUSY.v1()*ZR(0,i) + mySUSY.v2()*ZR(1,i);
        b0 = PV.B0(mu2, p2, Mw2, mH02[i]);
        /* Mw^2/v^2 is substituted for g2^2/4 compared to the expression in the
         * paper, in order to ensure the cancellation of the UV divergences in
         * the case where Mw is not the tree-level value. */
        PiT_WZH += - g2sq*Mw2/mySUSY.v()/mySUSY.v()*CR_i*CR_i*b0;
    }

    /* Z-boson - charged-Goldstone-boson loops */
    b0 = PV.B0(mu2, p2, Mw2, Mz*Mz);
    PiT_WZH += - e2*sW2*Mz*Mz*b0;

    /* gauge-boson loops */
    a0 = PV.A0(mu2, Mw2);
    b0 = PV.B0(mu2, p2, Mw2, Mz*Mz);
    b00 = PV.B00(mu2, p2, Mz*Mz, Mw2);
    PiT_WZH += g2sq*cW2*(2.0*PV.A0(mu2, Mz*Mz) - a0
                         + (4.0*p2 + Mz*Mz + Mw2)*b0 + 8.0*b00);
    //
    PiT_WZH += 3.0*g2sq*a0;
    //
    b0 = PV.B0(mu2, p2, Mw2, 0.0);
    b00 = PV.B00(mu2, p2, Mw2, 0.0);
    PiT_WZH += e2*((4.0*p2 + Mw2)*b0 - a0 + 8.0*b00);

    /* Sum of all contributions */
    complex PiT = PiT_f + PiT_sf + PiT_ch + PiT_WZH;

    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiT_AZ(const double mu, const double p2, const double Mw_i) const
{
    double mu2 = mu*mu;
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double Mw2 = Mw_i*Mw_i;
    double mHp2[2] = {mySUSY.getMHp()*mySUSY.getMHp(), Mw2}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2 = e/sW;
    double e_4sc = e/4.0/sW/cW;
    double e_2sc = 2.0*e_4sc;
    double e2_sc = e2/sW/cW;

    complex PiT_f = complex(0.0, 0.0, false);
    complex PiT_sf = complex(0.0, 0.0, false);
    complex PiT_ch = complex(0.0, 0.0, false);
    complex PiT_WZH = complex(0.0, 0.0, false);
    double a0;
    complex b0, b00;

    /* SM fermion loops */
    complex cV_Aee = - e;
    complex cA_Aee = 0.0;
    complex cV_Add = - e/3.0;
    complex cA_Add = 0.0;
    complex cV_Auu = 2.0/3.0*e;
    complex cA_Auu = 0.0;
    complex cV_Zee = - e_4sc*(1.0 - 4.0*sW2);
    complex cA_Zee = - e_4sc;
    complex cV_Zdd = - e_4sc*(1.0 - 4.0/3.0*sW2);
    complex cA_Zdd = - e_4sc;
    complex cV_Zuu = e_4sc*(1.0 - 8.0/3.0*sW2);
    complex cA_Zuu = e_4sc;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        PiT_f += FA(mu, p2, m_l[I], m_l[I], cV_Aee, cV_Zee, cA_Aee, cA_Zee);

        /* down-type quarks */
        PiT_f += Nc*FA(mu, p2, m_d[I], m_d[I], cV_Add, cV_Zdd, cA_Add, cA_Zdd);

        /* up-type quarks */
        PiT_f += Nc*FA(mu, p2, m_u[I], m_u[I], cV_Auu, cV_Zuu, cA_Auu, cA_Zuu);
    }

    /* charged-slepton loops */
    complex VZLL_nn, VAZLL_nn;
    for (int n=0; n<6; ++n) {
        VZLL_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed sleptons */
            VZLL_nn += - e_2sc*ZL(I,n)*ZL(I,n).conjugate();
        VZLL_nn += - e_2sc*(- 2.0*sW2);
        b00 = PV.B00(mu2, p2, Mse2[n], Mse2[n]);
        /* typo in the paper: e^2 --> e */
        PiT_sf += - 4.0*e*VZLL_nn*b00;

        VAZLL_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed sleptons */
            VAZLL_nn += e2_sc*ZL(I,n)*ZL(I,n).conjugate();
        VAZLL_nn += e2_sc*(- 2.0*sW2);
        a0 = PV.A0(mu2, Mse2[n]);
        PiT_sf += VAZLL_nn*a0;
    }

    /* down-type squark loops */
    complex VZDD_nn, VAZDD_nn;
    for (int n=0; n<6; ++n) {
        VZDD_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VZDD_nn += - e_2sc*ZD(I,n)*ZD(I,n).conjugate();
        VZDD_nn += - e_2sc*(- 2.0/3.0*sW2);
        b00 = PV.B00(mu2, p2, Msd2[n], Msd2[n]);
        /* typo in the paper: e^2 --> e */
        PiT_sf += - 4.0*e*Nc/3.0*VZDD_nn*b00;

        VAZDD_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VAZDD_nn += e2_sc/3.0*ZD(I,n)*ZD(I,n).conjugate();
        VAZDD_nn += e2_sc/3.0*(- 2.0/3.0*sW2);
        a0 = PV.A0(mu2, Msd2[n]);
        PiT_sf += Nc*VAZDD_nn*a0;
    }

    /* up-type squark loops */
    complex VZUU_nn, VAZUU_nn;
    for (int n=0; n<6; ++n) {
        VZUU_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VZUU_nn += e_2sc*ZU(I,n).conjugate()*ZU(I,n);
        VZUU_nn += e_2sc*(- 4.0/3.0*sW2);
        b00 = PV.B00(mu2, p2, Msu2[n], Msu2[n]);
        /* typo in the paper: e^2 --> e */
        PiT_sf += 4.0*e*Nc*2.0/3.0*VZUU_nn*b00;

        VAZUU_nn = complex(0.0, 0.0, false);
        for (int I=0; I<3; ++I) /* sum over left-handed squarks */
            VAZUU_nn += 2.0*e2_sc/3.0*ZU(I,n).conjugate()*ZU(I,n);
        VAZUU_nn += 2.0*e2_sc/3.0*(- 4.0/3.0*sW2);
        a0 = PV.A0(mu2, Msu2[n]);
        PiT_sf += Nc*VAZUU_nn*a0;
    }

    /* chargino loops */
    complex cV_Aii = e;
    complex cA_Aii = 0.0;
    complex cV_Zii, cA_Zii;
    for (int i=0; i<2; ++i) {
        cV_Zii = e_4sc*(  Zp(0,i).conjugate()*Zp(0,i)
                        + Zm(0,i)*Zm(0,i).conjugate() + 2.0*(cW2 - sW2) );
        cA_Zii = e_4sc*(  Zp(0,i).conjugate()*Zp(0,i)
                        - Zm(0,i)*Zm(0,i).conjugate() );
        PiT_ch += FA(mu, p2, mC[i], mC[i], cV_Aii, cV_Zii, cA_Aii, cA_Zii);
    }

    /* W-boson - charged-Goldstone-boson loops */
    b0 = PV.B0(mu2, p2, Mw2, Mw2);
    PiT_WZH += 2.0*e2*cW*sW*Mz*Mz*b0;

    /* W-boson loops */
    a0 = PV.A0(mu2, Mw2);
    //b0 = PV.B0(mu2, p2, Mw2, Mw2); /* same as the above */
    b00 = PV.B00(mu2, p2, Mw2, Mw2);
    PiT_WZH += 2.0*e*g2*cW*(2.0*a0 + (2.0*p2 + Mw2)*b0 + 4.0*b00);

    /* charged-Higgs loops */
    double cot_2thW = (cW2 - sW2)/(2.0*sW*cW);
    for (int i=0; i<2; ++i) {
        a0 = PV.A0(mu2, mHp2[i]);
        b00 = PV.B00(mu2, p2, mHp2[i], mHp2[i]);
        PiT_WZH += 2.0*e2*cot_2thW*(2.0*b00 + a0);
    }

    /* Sum of all contributions */
    complex PiT = PiT_f + PiT_sf + PiT_ch + PiT_WZH;

    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiTp_A(const double mu, const double p2, const double Mw_i) const
{
    double mu2 = mu*mu;
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double Mw2 = Mw_i*Mw_i;
    double mHp2[2] = {mySUSY.getMHp()*mySUSY.getMHp(), Mw2}; /* H^+_i = (H^+, G^+) */

    complex PiTp_f = complex(0.0, 0.0, false);
    complex PiTp_sf = complex(0.0, 0.0, false);
    complex PiTp_ch = complex(0.0, 0.0, false);
    complex PiTp_WZH = complex(0.0, 0.0, false);
    complex b0, b0p, b00p;

    /* SM fermion loops */
    complex cV_Aee = - e;
    complex cA_Aee = 0.0;
    complex cV_Add = - e/3.0;
    complex cA_Add = 0.0;
    complex cV_Auu = 2.0/3.0*e;
    complex cA_Auu = 0.0;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        PiTp_f += dFA(mu, p2, m_l[I], m_l[I], cV_Aee, cV_Aee, cA_Aee, cA_Aee);

        /* down-type quarks */
        PiTp_f += Nc*dFA(mu, p2, m_d[I], m_d[I], cV_Add, cV_Add, cA_Add, cA_Add);

        /* up-type quarks */
        PiTp_f += Nc*dFA(mu, p2, m_u[I], m_u[I], cV_Auu, cV_Auu, cA_Auu, cA_Auu);
    }

    /* charged-slepton loops */
    for (int n=0; n<6; ++n) {
        b00p = PV.B00p(mu2, p2, Mse2[n], Mse2[n]);
        PiTp_sf += 4.0*e2*b00p;
    }

    /* down-type squark loops */
    for (int n=0; n<6; ++n) {
        b00p = PV.B00p(mu2, p2, Msd2[n], Msd2[n]);
        PiTp_sf += 4.0*e2*Nc/3.0/3.0*b00p;
    }

    /* up-type squark loops */
    for (int n=0; n<6; ++n) {
        b00p = PV.B00p(mu2, p2, Msu2[n], Msu2[n]);
        PiTp_sf += 4.0*e2*Nc*2.0/3.0*2.0/3.0*b00p;
    }

    /* chargino loops */
    complex cV_Aii = e;
    complex cA_Aii = 0.0;
    for (int i=0; i<2; ++i)
        PiTp_ch += dFA(mu, p2, mC[i], mC[i], cV_Aii, cV_Aii, cA_Aii, cA_Aii);

    /* charged-Higgs loops */
    for (int i=0; i<2; ++i) {
        b00p = PV.B00p(mu2, p2, mHp2[i], mHp2[i]);
        PiTp_WZH += 4.0*e2*b00p;
    }

    /* W-boson loops */
    /* The Mw_i*Mw_i*b0p term, adding the corresponding contribution from
     * the W-G loop below, differs from the one in the paper. */
    b0 = PV.B0(mu2, p2, Mw2, Mw2);
    b0p = PV.B0p(mu2, p2, Mw2, Mw2);
    b00p = PV.B00p(mu2, p2, Mw2, Mw2);
    PiTp_WZH += 2.0*e2*( (2.0*p2 + Mw2)*b0p + 2.0*b0 + 4.0*b00p);

    /* W-boson - charged-Goldstone-boson loop */
    //b0p = PV.B0p(mu2, p2, Mw2, Mw2); /* Same as the above */
    PiTp_WZH += - 2.0*e2*Mw2*b0p;

    /* Sum of all contributions */
    complex PiTp = PiTp_f + PiTp_sf + PiTp_ch + PiTp_WZH;

    return ( PiTp/16.0/M_PI/M_PI );
}

double EWSUSY::PiThat_W_0(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i * RenormalizationScaleFactor;

    double Mz = mySUSY.getMz();
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);

    double PiThat = 0.0;

    /* W self-energy */
    PiThat += PiT_W(mu, 0.0, Mw_i).real();

    /* W-mass counter term */
    double delMw2 = PiT_W(mu, Mw_i*Mw_i, Mw_i).real();
    PiThat -= delMw2;

    /* counter term for e: (del e)/e */
    double dele_over_e = PiTp_A(mu, 0.0, Mw_i).real()/2.0
                         + sW/cW*PiT_AZ(mu, 0.0, Mw_i).real()/Mz/Mz;
    PiThat += 2.0*Mw_i*Mw_i*dele_over_e;

    /* counter term for sW: (del sW)/sW */
    double delSw_overSw = - cW2/2.0/sW2
                            *( PiT_W(mu, Mw_i*Mw_i, Mw_i).real()/Mw_i/Mw_i
                              - PiT_Z(mu, Mz*Mz, Mw_i).real()/Mz/Mz );
    PiThat -= 2.0*Mw_i*Mw_i*delSw_overSw;

    /* remaining counter terms,
     * usually denoted by 2/(sW*cW)*PiT_AZ(0)/Mz/Mz. */
    PiThat += - 2.0*Mw_i*Mw_i/(sW*cW)*PiT_AZ(mu, 0.0, Mw_i).real()/Mz/Mz;

    return PiThat;
}

double EWSUSY::DeltaR_rem_SM(const double Mw_i) const
{
    double cW2 = Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    double sW2 = 1.0 - cW2;

    /* renormalized vertex corrections + box */
    return ( mySUSY.getAle()/4.0/M_PI/sW2
             *(6.0 + (7.0 - 4.0*sW2)/2.0/sW2*log(cW2)) );
}

double EWSUSY::DeltaR_boxLL_SUSY(const double Mw_i) const
{
    int M = 1; // MU
    int N = 0; // ELECTRON
    int J = 1; // NEUTRINO_2
    int I = 0; // NEUTRINO_1

    complex a11 = complex(0.0, 0.0, false);
    complex a12 = complex(0.0, 0.0, false);

    /* charged-lepton - sneutrino - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    complex FF = F(sqrt(Mse2[k]), sqrt(Msn2[k]), mC[i], mN[j]);
                    a11 += 0.5
                           *L_esnC(M, K, i, Mw_i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *mC[i]*mN[j]*FF;
                    a11 += 0.5
                           *L_eLN(M, k, j, Mw_i)
                           *L_nsnN(I, K, j, Mw_i)
                           *L_nLC(J, k, i, Mw_i).conjugate()
                           *L_esnC(N, K, i, Mw_i).conjugate()
                           *mC[i]*mN[j]*FF;
                }

    /* charged-lepton - charged-lepton - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int l=0; l<6; ++l)
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a11 +=  L_eLN(M, k, j, Mw_i)
                           *L_nLC(J, k, i, Mw_i).conjugate()
                           *L_nLC(I, l, i, Mw_i)
                           *L_eLN(N, l, j, Mw_i).conjugate()
                           *H(sqrt(Mse2[k]), sqrt(Mse2[l]), mC[i], mN[j]);
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a11 +=  L_esnC(M, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_nsnN(I, L, j, Mw_i)
                           *L_esnC(N, L, i, Mw_i).conjugate()
                           *H(sqrt(Msn2[K]), sqrt(Msn2[L]), mC[i], mN[j]);
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {
                    a12 += 0.5
                           *L_esnC(M, K, i, Mw_i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *L_esnC(N, K, j, Mw_i).conjugate()
                           *mC[i]*mC[j]*F(sqrt(Mse2[k]), sqrt(Msn2[K]), mC[i], mC[j]);
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    a12 += 0.5
                           *L_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F(sqrt(Mse2[k]), sqrt(Msn2[K]), mN[i], mN[j]);
                    a12 +=  L_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *H(sqrt(Mse2[k]), sqrt(Msn2[K]), mN[i], mN[j]);
               }

    complex a1 = (a11 + a12)/16.0/M_PI/M_PI;

    double sW2 = 1.0 - Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    return ( - sW2*Mw_i*Mw_i/2.0/M_PI/mySUSY.getAle()*a1.real() );
}

double EWSUSY::DeltaR_boxLR_SUSY(const double Mw_i) const
{
    int M = 1; // MU
    int N = 0; // ELECTRON
    int J = 1; // NEUTRINO_2
    int I = 0; // NEUTRINO_1

    complex a21 = complex(0.0, 0.0, false);
    complex a22 = complex(0.0, 0.0, false);

    /* charged-lepton - sneutrino - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    complex HH = H(sqrt(Mse2[k]), sqrt(Msn2[K]), mC[i], mN[j]);
                    a21 += - 2.0
                             *R_esnC(M, K, i)
                             *L_nLC(I, k, i, Mw_i)
                             *L_nsnN(J, K, j, Mw_i).conjugate()
                             *R_eLN(N, k, j, Mw_i).conjugate()
                             *HH;
                    a21 += - 2.0
                             *R_eLN(M, k, j, Mw_i)
                             *L_nsnN(I, K, j, Mw_i)
                             *L_nLC(J, k, i, Mw_i).conjugate()
                             *R_esnC(N, K, i).conjugate()
                             *HH;
                }

    /* charged-lepton - charged-lepton - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int l=0; l<6; ++l)
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a21 += - 2.0
                             *R_eLN(M, k, j, Mw_i)
                             *L_nLC(J, k, i, Mw_i).conjugate()
                             *L_nLC(I, l, i, Mw_i)
                             *R_eLN(N, l, j, Mw_i).conjugate()
                             *H(sqrt(Mse2[k]), sqrt(Mse2[l]), mC[i], mN[j]);
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a21 += - 2.0
                             *R_esnC(M, K, i)
                             *L_nsnN(J, K, j, Mw_i).conjugate()
                             *L_nsnN(I, L, j, Mw_i)
                             *R_esnC(N, L, i).conjugate()
                             *H(sqrt(Msn2[K]), sqrt(Msn2[L]), mC[i], mN[j]);
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    a22 += R_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F(sqrt(Mse2[k]), sqrt(Msn2[K]), mN[i], mN[j]);
                    a22 += 2.0
                           *R_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *H(sqrt(Mse2[k]), sqrt(Msn2[K]), mN[i], mN[j]);
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {
                    a22 += 2.0
                           *R_esnC(M, K, i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *R_esnC(N, K, j).conjugate()
                           *H(sqrt(Mse2[k]), sqrt(Msn2[K]), mC[i], mC[j]);
                }

    complex a2 = (a21 + a22)/16.0/M_PI/M_PI;

    double sW2 = 1.0 - Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    return ( sW2*Mw_i*Mw_i/4.0/M_PI/mySUSY.getAle()*a2.real() );
}

complex EWSUSY::v(const double mu, const StandardModel::lepton M,
                  const StandardModel::lepton J, const double Mw_i) const
{
    int intM, intJ;
    switch (M) {
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            intM = ((int)M - StandardModel::ELECTRON)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::v(): Wrong argument!");
    }
    switch (J) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            intJ = ((int)J - StandardModel::NEUTRINO_1)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::v(): Wrong argument!");
    }

    complex v = complex(0.0, 0.0, false);
    complex b0, ff;
    complex CL_ji, CR_ji; /* chargino-neutralino-W couplings */

    /* charged-slepton - chargino - neutralino loops */
    for (int k=0; k<6; ++k)
        for (int j=0; j<4; ++j)
            for (int i=0; i<2; ++i) {
                CL_ji = ZN(1,j)*Zp(0,i).conjugate()
                        - ZN(3,j)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,j).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu*mu, 0.0, mC[i]*mC[i], mN[j]*mN[j]);
                ff = f(sqrt(Mse2[k]), mC[i], mN[j]);
                v += L_nLC(intJ, k, i, Mw_i).conjugate()*L_eLN(intM, k, j, Mw_i)
                     *( sqrt(2.0)*CL_ji*mC[i]*mN[j]*ff
                        - CR_ji/sqrt(2.0)*(b0 - 0.5 + Mse2[k]*ff) );
            }

    /* sneutrino - neutralino - chargino loops */
    for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
        for (int j=0; j<4; ++j)
            for (int i=0; i<2; ++i) {
                CL_ji = ZN(1,j)*Zp(0,i).conjugate()
                        - ZN(3,j)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,j).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu*mu, 0.0, mC[i]*mC[i], mN[j]*mN[j]);
                ff = f(sqrt(Msn2[K]), mC[i], mN[j]);
                v += L_nsnN(intJ, K, j, Mw_i).conjugate()*L_esnC(intM, K, i, Mw_i)
                     *( - sqrt(2.0)*CR_ji*mC[i]*mN[j]*ff
                        + CL_ji/sqrt(2.0)*(b0 - 0.5 + Msn2[K]*ff) );
            }

    /* sneutrino - charged-slepton - neutralino loops */
    matrix<complex> ZneT_ZL = Zne.transpose()*ZL;
    for (int i=0; i<6; ++i)
        for (int j=0; j<4; ++j)
            for (int K=0; K<3; ++K) {  /* K=0-2 for left-handed sneutrinos */
                b0 = PV.B0(mu*mu, 0.0, Mse2[i], Msn2[K]);
                ff = f(mN[j], sqrt(Mse2[i]), sqrt(Msn2[K]));
                v += 0.5*L_nsnN(intJ, K, j, Mw_i).conjugate()*L_eLN(intM, i, j, Mw_i)
                     *ZneT_ZL(K, i).conjugate()*(b0 + 0.5 + mN[j]*mN[j]*ff);
            }

    return ( v/16.0/M_PI/M_PI );
}

complex EWSUSY::delta_v(const double mu, const StandardModel::lepton M,
                        const StandardModel::lepton J, const double Mw_i) const
{
    int intM, intJ;
    switch (M) {
        case StandardModel::ELECTRON:
        case StandardModel::MU:
        case StandardModel::TAU:
            intM = ((int)M - StandardModel::ELECTRON)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::delta_v(): Wrong argument!");
    }
    switch (J) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            intJ = ((int)J - StandardModel::NEUTRINO_1)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::delta_v(): Wrong argument!");
    }

    complex delv = complex(0.0, 0.0, false);
    double muIR = mu; /* fictional scale, since B0p(0,m1^2,m2^2) is IR finite */
    complex b0p, b0;

    /* charged-slepton - neutralino loops */
    for (int k=0; k<6; ++k)
        for (int j=0; j<4; ++j) {
            b0p = PV.B0p(muIR*muIR, 0.0, Mse2[k], mN[j]*mN[j]);
            b0 = PV.B0(mu*mu, 0.0, Mse2[k], mN[j]*mN[j]);
            delv += 0.5*L_eLN(intM, k, j, Mw_i)*L_eLN(intJ, k, j, Mw_i).conjugate()
                    *( (Mse2[k] - mN[j]*mN[j])*b0p - b0 );
        }

    /* sneutrino - chargino loops */
    for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
        for (int i=0; i<2; ++i) {
            b0p = PV.B0p(muIR*muIR, 0.0, Msn2[K], mC[i]*mC[i]);
            b0 = PV.B0(mu*mu, 0.0, Msn2[K], mC[i]*mC[i]);
            delv += 0.5*L_esnC(intM, K, i, Mw_i)*L_esnC(intJ, K, i, Mw_i).conjugate()
                    *( (Msn2[K] - mC[i]*mC[i])*b0p - b0 );
        }

    return ( delv/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_vertex_SUSY(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i * RenormalizationScaleFactor;

    return ( v(mu, mySUSY.ELECTRON, mySUSY.NEUTRINO_1, Mw_i).real()
            + delta_v(mu, mySUSY.ELECTRON, mySUSY.NEUTRINO_1, Mw_i).real()
            + v(mu, mySUSY.MU, mySUSY.NEUTRINO_2, Mw_i).real()
            + delta_v(mu, mySUSY.MU, mySUSY.NEUTRINO_2, Mw_i).real() );
}

complex EWSUSY::Sigma_nu_0(const double mu, const StandardModel::lepton I,
                           const StandardModel::lepton J, const double Mw_i) const
{
    int intI, intJ;
    switch (I) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            intI = ((int)I - StandardModel::NEUTRINO_1)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::Sigma_nu(): Wrong argument!");
    }
    switch (J) {
        case StandardModel::NEUTRINO_1:
        case StandardModel::NEUTRINO_2:
        case StandardModel::NEUTRINO_3:
            intJ = ((int)J - StandardModel::NEUTRINO_1)/2;
            break;
        default:
            throw std::runtime_error("EWSUSY::Sigma_nu(): Wrong argument!");
    }

    complex Sigma = complex(0.0, 0.0, false);
    double muIR = mu; /* fictional scale, since B0p(0,m1,m2) is IR finite */
    complex b0p, b0;

    /* charged-slepton - chargino loops */
    for (int k=0; k<6; ++k)
        for (int i=0; i<2; ++i) {
            b0p = PV.B0p(muIR*muIR, 0.0, Mse2[k], mC[i]*mC[i]);
            b0 = PV.B0(mu*mu, 0.0, Mse2[k], mC[i]*mC[i]);
            Sigma += 0.5*L_nLC(intI, k, i, Mw_i)*L_nLC(intJ, k, i, Mw_i).conjugate()
                     *( (Mse2[k] - mC[i]*mC[i])*b0p - b0 );
        }

    /* sneutrino - neutralino loops */
    for (int K=0; K<3; ++K)  /* K=0-2 for left-handed sneutrinos */
        for (int j=0; j<4; ++j) {
            b0p = PV.B0p(muIR*muIR, 0.0, Msn2[K], mN[j]*mN[j]);
            b0 = PV.B0(mu*mu, 0.0, Msn2[K], mN[j]*mN[j]);
            Sigma += 0.5*L_nsnN(intI, K, j, Mw_i)*L_nsnN(intJ, K, j, Mw_i).conjugate()
                     *( (Msn2[K] - mN[j]*mN[j])*b0p - b0 );
        }

    return ( Sigma/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_neutrino_SUSY(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i * RenormalizationScaleFactor;

    return ( ( Sigma_nu_0(mu, mySUSY.NEUTRINO_1, mySUSY.NEUTRINO_1, Mw_i).real()
              - delta_v(mu, mySUSY.ELECTRON, mySUSY.NEUTRINO_1, Mw_i).real()
              + Sigma_nu_0(mu, mySUSY.NEUTRINO_2, mySUSY.NEUTRINO_2, Mw_i).real()
              - delta_v(mu, mySUSY.MU, mySUSY.NEUTRINO_2, Mw_i).real() )/2.0 );
}

double EWSUSY::DeltaR_TOTAL_EW1(const double Mw_i) const
{
    double DeltaR = 0.0;

    /* SM+SUSY renormalized W self energy */
    DeltaR += - PiThat_W_0(Mw_i)/Mw_i/Mw_i;

    /* SM renormalized vertex + box */
    DeltaR += DeltaR_rem_SM(Mw_i);

    /* SUSY box corrections */
    DeltaR += DeltaR_boxLL_SUSY(Mw_i);
    DeltaR += DeltaR_boxLR_SUSY(Mw_i);

    /* SUSY renormalized vertex corrections */
    DeltaR += DeltaR_vertex_SUSY(Mw_i);

    /* SUSY renormalized neutrino wave function */
    DeltaR += DeltaR_neutrino_SUSY(Mw_i);

    /* Debug */
    //std::cout << "MSSM WSE = " << - PiThat_W_0(Mw_i)/Mw_i/Mw_i << std::endl;
    //std::cout << "SM VC+Box = " << DeltaR_rem_SM(Mw_i) << std::endl;
    //std::cout << "SUSY BoxLL = " << DeltaR_boxLL_SUSY(Mw_i) << std::endl;
    //std::cout << "SUSY BoxLR = " << DeltaR_boxLR_SUSY(Mw_i) << std::endl;
    //std::cout << "SUSY VC = " << DeltaR_vertex_SUSY(Mw_i) << std::endl;
    //std::cout << "SUSY nuSE = " << DeltaR_neutrino_SUSY(Mw_i) << std::endl;

    return DeltaR;
}

double EWSUSY::DeltaAlphaL5q_SM_EW1() const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = mySUSY.getMz() * RenormalizationScaleFactor;

    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double Nc = mySUSY.getNc();

    double DelA_l = 0.0, DelA_d = 0.0, DelA_u = 0.0;

    /* SM fermion loops */
    complex cV_Aee = - e;
    complex cA_Aee = 0.0;
    complex cV_Add = - e/3.0;
    complex cA_Add = 0.0;
    complex cV_Auu = 2.0/3.0*e;
    complex cA_Auu = 0.0;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        DelA_l += FA(mu, Mz2, m_l[I], m_l[I], cV_Aee, cV_Aee, cA_Aee, cA_Aee).real()/Mz2;
        DelA_l -= dFA(mu, 0.0, m_l[I], m_l[I], cV_Aee, cV_Aee, cA_Aee, cA_Aee).real();

        /* down-type quarks */
        DelA_d += Nc*FA(mu, Mz2, m_d[I], m_d[I], cV_Add, cV_Add, cA_Add, cA_Add).real()/Mz2;
        DelA_d -= Nc*dFA(mu, 0.0, m_d[I], m_d[I], cV_Add, cV_Add, cA_Add, cA_Add).real();

        /* up-type quarks, not including top quark */
        if (I!=3) {
            DelA_u += Nc*FA(mu, Mz2, m_u[I], m_u[I], cV_Auu, cV_Auu, cA_Auu, cA_Auu).real()/Mz2;
            DelA_u -= Nc*dFA(mu, 0.0, m_u[I], m_u[I], cV_Auu, cV_Auu, cA_Auu, cA_Auu).real();
        }
    }

    /* Debug */
    //std::cout << "EWSUSY(l) " << DelA_l/16.0/M_PI/M_PI << std::endl;
    //std::cout << "EWSM(l)   " << getMyOneLoopEW()->DeltaAlpha_l(Mz2) << std::endl;
    //std::cout << "EWSUSY(q) " << (DelA_d + DelA_u)/16.0/M_PI/M_PI << std::endl;
    //std::cout << "EWSM(had) " << getMyOneLoopEW()->DeltaAlpha_5q(Mz2) << std::endl;

    return ( (DelA_l + DelA_d + DelA_u)/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_SUSY_EW1(const double Mw_i) const
{
    double cW2 = Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    double sW2 = 1.0 - cW2;

    /* SM one-loop contributions */
    double DeltaAlphaL5q_EW1 = DeltaAlphaL5q_SM_EW1();
    double DeltaRho_EW1 = getMyOneLoopEW()->DeltaRho(Mw_i);
    double DeltaR_rem_EW1 = getMyOneLoopEW()->DeltaR_rem(Mw_i);
    double DeltaR_SM_EW1 = DeltaAlphaL5q_EW1 - cW2/sW2*DeltaRho_EW1 + DeltaR_rem_EW1;

    /* Debug */
    //std::cout << std::endl;
    //std::cout << "DeltaAlphaL5q_EW1 = " << DeltaAlphaL5q_EW1 << std::endl;
    //std::cout << "-cW2/sW2*DeltaRho_EW1 = " << - cW2/sW2*DeltaRho_EW1 << std::endl;
    //std::cout << "DeltaR_rem_EW1 = " << DeltaR_rem_EW1 << std::endl;
    //std::cout << "DeltaR_SM_EW1 = " << DeltaR_SM_EW1 << std::endl;

    return ( DeltaR_TOTAL_EW1(Mw_i) - DeltaR_SM_EW1 );
}

double EWSUSY::Mw_MSSM_TMP(const double Mw_i) const
{
    if (getSchemeMw()!=EWSM::NORESUM)
        throw std::runtime_error("EWSUSY::Mw_SUSY(): Scheme for Mw is not applicable");

    double cW2 = Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    double sW2 = 1.0 - cW2;
    if (sW2 < 0.0)
        throw std::runtime_error("EWSUSY::Mw_SUSY(): negative sW2");

    /* SM contributions to Delta r */
    double dAleL5q = DeltaAlphaL5q();
    double DeltaRho[EWSM::orders_EW_size], DeltaRho_sum = 0.0;
    double DeltaR_rem[EWSM::orders_EW_size], DeltaR_rem_sum = 0.0;
    ComputeDeltaRho(Mw_i, DeltaRho);
    ComputeDeltaR_rem(Mw_i, DeltaR_rem);
    for (int j=0; j<EWSM::orders_EW_size; ++j) {
        /* excluding EW two-loop contributions, which will be added below */
        if (j!=(int)EWSM::EW2) {
            DeltaRho_sum += DeltaRho[(EWSM::orders_EW)j];
            DeltaR_rem_sum += DeltaR_rem[(EWSM::orders_EW)j];
        }
    }

    /* Full EW one-loop contribution (without the full DeltaAlphaL5q) */
    double DeltaR_EW1 = - cW2/sW2*DeltaRho[EWSM::EW1] + DeltaR_rem[EWSM::EW1];

    /* Full EW two-loop contribution with reducible corrections */
    double DeltaR_EW2 = dAleL5q*dAleL5q + 2.0*dAleL5q*DeltaR_EW1
    + getMyApproximateFormulae()->DeltaR_TwoLoopEW_rem(Mw_i);

    /* R = 1 + Delta r */
    double R = 1.0 + dAleL5q - cW2/sW2*DeltaRho_sum + DeltaR_rem_sum + DeltaR_EW2;

    /* SUSY contribution */
    R += DeltaR_SUSY_EW1(Mw_i);

    /* the W-boson mass in the complex pole scheme */
    double tmp = 4.0*M_PI*mySUSY.getAle()/sqrt(2.0)/mySUSY.getGF()/Mzbar()/Mzbar();
    if (tmp*R > 1.0) throw std::runtime_error("EWSUSY::Mw(): Negative (1-tmp*R)");
    double Mwbar = Mzbar()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R));

    /* complex-pole/fixed-width scheme --> experimental/running-width scheme */
    double Mw_exp = MwFromMwbar(Mwbar);

    if (Mw_exp >= mySUSY.getMz()) {
        std::cout << "WARNING: Mw > Mz in EWSUSY::Mw_MSSM_TMP" << std::endl;
        //std::cout << Mw_i << std::endl;
        //std::cout << Mw_exp << std::endl;
        //std::cout << - PiThat_W_0(Mw_i)/Mw_i/Mw_i << std::endl;
        //double mu = Mw_i, Mz = mySUSY.getMz();
        //std::cout << "  " << PiT_W(mu, 0.0, Mw_i).real() << std::endl;
        //std::cout << "  " << PiT_W(mu, Mw_i*Mw_i, Mw_i).real() << std::endl;
        //std::cout << "  " << PiT_W(mu, Mw_i*Mw_i, Mw_i).real() << std::endl;
        //std::cout << "  " << PiT_Z(mu, Mz*Mz, Mw_i).real() << std::endl;
        //std::cout << "  " << PiT_AZ(mu, 0.0, Mw_i).real() << std::endl;
        //std::cout << "  " << PiT_AZ(mu, 0.0, Mw_i).real() << std::endl;
        //std::cout << "  " << PiTp_A(mu, 0.0, Mw_i).real() << std::endl;
        //std::cout << DeltaR_rem_SM(Mw_i) << std::endl;
        //std::cout << DeltaR_boxLL_SUSY(Mw_i) << std::endl;
        //std::cout << DeltaR_boxLR_SUSY(Mw_i) << std::endl;
        //std::cout << DeltaR_vertex_SUSY(Mw_i) << std::endl;
        //std::cout << DeltaR_neutrino_SUSY(Mw_i) << std::endl;

        return Mw_unphysical;
    } else
        return Mw_exp;
}

double EWSUSY::Mw_MSSM() const
{
    /* initial value for Mw */
    double Mw_org = mySUSY.getMyFH()->getMw_FHinput();

    double Mw = Mw_MSSM_TMP(Mw_org);
    //std::cout << std::endl << std::setprecision(12)
    //          << "EWSUSY::Mw_MSSM(): Mw_org = " << Mw_org
    //          << "  Mw_new = " << Mw << std::endl;

    if (Mw == Mw_unphysical) return Mw_unphysical;

    /* iterations */
    while (fabs(Mw - Mw_org) > EWSM::Mw_error) {
        Mw_org = Mw;
        Mw = Mw_MSSM_TMP(Mw);
        //std::cout << std::setprecision(12)
        //          << "EWSUSY::Mw_MSSM(): Mw_org = " << Mw_org
        //          << "  Mw_new = " << Mw << std::endl;

        if (Mw == Mw_unphysical) return Mw_unphysical;
    }

    return Mw;
}

complex EWSUSY::L_esnC(const int N, const int K, const int i, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();
    double sW = sqrt(1.0 - cW*cW);

    return ( e/sW*Zp(0,i)*Zne(N,K).conjugate() );
}

complex EWSUSY::R_esnC(const int N, const int K, const int i) const
{
    return ( Yl(N,N)*Zne(N,K).conjugate()*Zm(1,i).conjugate() );
}

complex EWSUSY::L_nLC(const int I, const int k, const int i, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();
    double sW = sqrt(1.0 - cW*cW);

    return ( e/sW*ZL(I,k)*Zm(0,i) + Yl(I,I)*ZL(I+3,k)*Zm(1,i) );
}

complex EWSUSY::L_nsnN(const int J, const int K, const int j, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();
    double sW = sqrt(1.0 - cW*cW);

    return ( - e/sqrt(2.0)/sW/cW*Zne(J,K).conjugate()*(ZN(0,j)*sW - ZN(1,j)*cW) );
}

complex EWSUSY::L_eLN(const int N, const int k, const int j, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();
    double sW = sqrt(1.0 - cW*cW);

    return ( - e/sqrt(2.0)/sW/cW*ZL(N,k)*(ZN(0,j)*sW + ZN(1,j)*cW)
             - Yl(N,N)*ZL(N+3,k)*ZN(2,j) );
}

complex EWSUSY::R_eLN(const int N, const int k, const int j, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();

    return ( e*sqrt(2.0)/cW*ZL(N+3,k)*ZN(0,j).conjugate()
             - Yl(N,N)*ZL(N,k)*ZN(2,j).conjugate() );
}

complex EWSUSY::F(const double m1, const double m2, const double m3,
                  const double m4) const
{
    return PV.D0(0.0, 0.0, m1*m1, m2*m2, m3*m3, m4*m4);
}

complex EWSUSY::H(const double m1, const double m2, const double m3,
                  const double m4) const
{
    return PV.D00(0.0, 0.0, m1*m1, m2*m2, m3*m3, m4*m4);
}

complex EWSUSY::f(const double m1, const double m2, const double m3) const
{
    return ( - PV.C0(0.0, m1*m1, m2*m2, m3*m3) );
}


