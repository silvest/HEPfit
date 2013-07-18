/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <stdexcept>
#include "EWSUSY.h"

EWSUSY::EWSUSY(const SUSY& SUSY_in)
: mySUSY(SUSY_in), myEWSMCache(SUSY_in), myEWSMOneLoopEW(myEWSMCache),
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

    Au = - mySUSY.getTU();
    Ad = mySUSY.getTD();
    Al = mySUSY.getTE();

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
        m_u[I] = myEWSMCache.mq((StandardModel::quark)(2*I), mySUSY.getMz(), FULLNNLO);
        /* down-type quarks */
        m_d[I] = myEWSMCache.mq((StandardModel::quark)(2*I + 1), mySUSY.getMz(), FULLNNLO);
        /* charged leptons */
        m_l[I] = mySUSY.getLeptons((StandardModel::lepton)(2*I + 1)).getMass();
    }
    /* H^0_i = (H^0, h^0, A^0, G^0) */
    mH0[0] = mySUSY.getMHh();
    mH0[1] = mySUSY.getMHl();
    mH0[2] = mySUSY.getMHa();
    mH0[3] = mySUSY.getMz(); /* mass of the neutral Goldstone boson */
    for (int k=0; k<6; ++k) {
        Msu[k] = sqrt(mySUSY.getMsu2()(k));
        Msd[k] = sqrt(mySUSY.getMsd2()(k));
        Mse[k] = sqrt(mySUSY.getMse2()(k));
    }
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        Msn[K] = sqrt(mySUSY.getMsn2()(K));
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
    /* PV functions */
    double A0i = PV.A0(mu, mi);
    double A0j = PV.A0(mu, mj);
    complex B0 = PV.B0(mu, p2, mi, mj);
    complex B22 = PV.B22(mu, p2, mi, mj);

    return ( -2.0*(cV_aij*cV_bji + cA_aij*cA_bji) 
              *(4.0*B22 + A0i + A0j + (p2 - mi*mi - mj*mj)*B0)
             -4.0*(cV_aij*cV_bji - cA_aij*cA_bji)*mi*mj*B0 );
}

complex EWSUSY::dFA(const double mu, const double p2,
                   const double mi, const double mj,
                   const complex cV_aij, const complex cV_bji,
                   const complex cA_aij, const complex cA_bji) const
{
    /* PV functions */
    complex B0 = PV.B0(mu, p2, mi, mj);
    complex B0p = PV.B0p(mu, p2, mi, mj);
    complex B22p = PV.B22p(mu, p2, mi, mj);

    if (mi == mj && cA_aij == 0.0 && cA_bji == 0.0)
        return ( -2.0*cV_aij*cV_bji*(4.0*B22p + p2*B0p + B0) );
    else
        return ( -2.0*(cV_aij*cV_bji + cA_aij*cA_bji)
                  *(4.0*B22p + (p2 - mi*mi - mj*mj)*B0p + B0)
                 -4.0*(cV_aij*cV_bji - cA_aij*cA_bji)*mi*mj*B0p );
}

complex EWSUSY::PiT_Z(const double mu, const double p2, const double Mw_i) const
{
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double mHp[2] = {mySUSY.getMHp(), Mw_i}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2sq = e2/sW2; /* g2 squared */
    double e_4sc = e/4.0/sW/cW;
    double e_2sc = 2.0*e_4sc;

    complex PiT = complex(0.0, 0.0, false);
    double a0;
    complex b0, b22;
    complex cV_Zij, cV_Zji, cA_Zij, cA_Zji;
    matrix<double> Id6 = matrix<double>::Id(6);
    matrix<double> Id2 = matrix<double>::Id(2);
    
    /* neutrino loops */
    b0 = PV.B0(mu, p2, 0.0, 0.0);
    b22 = PV.B22(mu, p2, 0.0, 0.0);
    PiT += - 3.0/4.0*g2sq/cW2*(4.0*b22 + p2*b0);

    /* other SM fermion loops */
    complex cV_Zee = - e_4sc*(1.0 - 4.0*sW2);
    complex cA_Zee = - e_4sc;
    complex cV_Zdd = - e_4sc*(1.0 - 4.0/3.0*sW2);
    complex cA_Zdd = - e_4sc;
    complex cV_Zuu = e_4sc*(1.0 - 8.0/3.0*sW2);
    complex cA_Zuu = e_4sc;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        PiT += FA(mu, p2, m_l[I], m_l[I], cV_Zee, cV_Zee, cA_Zee, cA_Zee);

        /* down-type quarks */
        PiT += Nc*FA(mu, p2, m_d[I], m_d[I], cV_Zdd, cV_Zdd, cA_Zdd, cA_Zdd);

        /* up-type quarks */
        PiT += Nc*FA(mu, p2, m_u[I], m_u[I], cV_Zuu, cV_Zuu, cA_Zuu, cA_Zuu);
    }

    /* sneutrino loops */
    complex VZsnsn_II = e_2sc;
    complex VZZsnsn_II = e2/2.0/sW2/cW2;
    for (int I=0; I<3; ++I) {  /* I=0-3 for left-handed sneutrinos */
        b22 = PV.B22(mu, p2, Msn[I], Msn[I]);
        PiT += 4.0*VZsnsn_II.abs2()*b22;

        a0 = PV.A0(mu, Msn[I]);
        PiT += VZZsnsn_II*a0;
    }

    /* charged-slepton loops */
    complex VZLL_mn, VZZLL_nn;
    matrix<complex> ZLhc_ZL = ZL.hconjugate()*ZL; 
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VZLL_mn = - e_2sc*(ZLhc_ZL(m,n) - 2.0*sW2*Id6(m,n));
            b22 = PV.B22(mu, p2, Mse[m], Mse[n]);
            PiT += 4.0*VZLL_mn.abs2()*b22;
        }
        VZZLL_nn = 2.0*e2/cW2*(sW2 + (1.0 - 4.0*sW2)/4.0/sW2*ZLhc_ZL(n,n));
        a0 = PV.A0(mu, Mse[n]);
        PiT += VZZLL_nn*a0;
    }

    /* down-type squark loops */
    complex VZDD_mn, VZZDD_nn;
    matrix<complex> ZDhc_ZD = ZD.hconjugate()*ZD; 
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            VZDD_mn = - e_2sc*(ZDhc_ZD(m,n) - 2.0/3.0*sW2*Id6(m,n));
            b22 = PV.B22(mu, p2, Msd[m], Msd[n]);
            PiT += 4.0*Nc*VZDD_mn.abs2()*b22;
        }
        VZZDD_nn = 2.0*e2/3.0/cW2*(sW2/3.0 + (3.0 - 4.0*sW2)/4.0/sW2*ZDhc_ZD(n,n));
        a0 = PV.A0(mu, Msd[n]);
        PiT += Nc*VZZDD_nn*a0;
    }

    /* up-type squark loops */
    complex VZUU_mn, VZZUU_nn;
    matrix<complex> ZUhc_ZU = ZU.hconjugate()*ZU;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            /* (n,m) instead of (m,n) */
            VZUU_mn = e_2sc*(ZUhc_ZU(n,m) - 4.0/3.0*sW2*Id6(m,n));
            b22 = PV.B22(mu, p2, Msu[m], Msu[n]);
            PiT += 4.0*Nc*VZUU_mn.abs2()*b22;
        }
        VZZUU_nn = 2.0*e2/3.0/cW2*(4.0*sW2/3.0 + (3.0 - 8.0*sW2)/4.0/sW2*ZUhc_ZU(n,n));
        a0 = PV.A0(mu, Msu[n]);
        PiT += Nc*VZZUU_nn*a0;
    }

    /* chargino loops */
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j) {
            cV_Zij = e_4sc*(Zp(0,j).conjugate()*Zp(0,i) 
                            + Zm(0,j)*Zm(0,i).conjugate()
                            + 2.0*(cW2 - sW2)*Id2(j,i));
            cV_Zji = cV_Zij.conjugate();
            cA_Zij = e_4sc*(Zp(0,j).conjugate()*Zp(0,i)
                            - Zm(0,j)*Zm(0,i).conjugate());
            cA_Zji = cA_Zij.conjugate();
            PiT += FA(mu, p2 ,mC[i], mC[j], cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }

    /* neutralino loops */
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j) {
            cV_Zij = e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                            - ZN(2,j).conjugate()*ZN(2,i)
                            - ZN(3,j)*ZN(3,i).conjugate()
                            + ZN(2,j)*ZN(2,i).conjugate());
            cV_Zji = cV_Zij.conjugate();
            cA_Zij = e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                            - ZN(2,j).conjugate()*ZN(2,i)
                            + ZN(3,j)*ZN(3,i).conjugate()
                            - ZN(2,j)*ZN(2,i).conjugate());
            cA_Zji = cA_Zij.conjugate();
            PiT += 0.5*FA(mu, p2 ,mN[i], mN[j], cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }
    
    /* charged-Higgs loops */
    double cot_2thW = (cW2 - sW2)/(2.0*sW*cW);
    for (int i=0; i<2; ++i) {
        b22 = PV.B22(mu, p2, mHp[i], mHp[i]);
        a0 = PV.A0(mu, mHp[i]);
        PiT += 2.0*e2*cot_2thW*cot_2thW*(2.0*b22 + a0);
    }
    
    /* neutral-Higgs loops */
    double AM_ij;
    for (int i=0; i<2; ++i)
        for (int j=0; j<2; ++j) {
            AM_ij = ZR(0,i)*ZH(0,j) - ZR(1,i)*ZH(1,j);
            b22 = PV.B22(mu, p2, mH0[i], mH0[j+2]);
            PiT += g2sq/cW2*AM_ij*AM_ij*b22;
        }
    for (int j=0; j<4; ++j) {
        a0 = PV.A0(mu, mH0[j]);
        PiT += g2sq/4.0/cW2*a0;
    }

    /* W-boson - charged-Goldstone-boson loop*/
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    PiT += - 2.0*g2sq*sW2*sW2*Mz*Mz*b0;

    /* Z-boson - Higgs loops */
    double CR_i;
    for (int i=0; i<2; ++i) {
        CR_i = mySUSY.v1()*ZR(0,i) + mySUSY.v2()*ZR(1,i);
        b0 = PV.B0(mu, p2, Mz, mH0[i]);
        PiT += - g2sq*g2sq/4.0/cW2/cW2*CR_i*CR_i*b0;
    }

    /* W-boson loops */
    a0 = PV.A0(mu, Mw_i);
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    b22 = PV.B22(mu, p2, Mw_i, Mw_i);
    /* a0^2 --> a0 in the first term */
    PiT += 2.0*g2sq*cW2*(2.0*a0 + (2.0*p2 + Mw_i*Mw_i)*b0 + 4.0*b22);

    return ( PiT/16.0/M_PI/M_PI );
}


complex EWSUSY::PiT_W(const double mu, const double p2, const double Mw_i) const
{
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double mHp[2] = {mySUSY.getMHp(), Mw_i}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2sq = e2/sW2; /* g2 squared */
    double e_2s = e/2.0/sW;
    double e_sq2s = e/sqrt(2.0)/sW;
    double e_2sq2s = e_sq2s/2.0;
    double e2_2s2 = e2/2.0/sW2;

    complex PiT = complex(0.0, 0.0, false);
    double a0;
    complex b0, b22;
    
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
        PiT += FA(mu, p2, m_l[I], 0.0, cV_Wen, cV_Wne, cA_Wen, cA_Wne);

        /* quarks (no CKM) */
        PiT += Nc*FA(mu, p2, m_d[I], m_u[I], cV_Wdu, cV_Wud, cA_Wdu, cA_Wud);
    }

    /* slepton loops */
    complex VWsnL_In, VWWsnsn_II, VWWLL_nn;
    matrix<complex> ZLT_Zne = ZL.transpose()*Zne;
    for (int I=0; I<3; ++I) {  /* I=0-3 for left-handed sneutrinos */
        for (int n=0; n<6; ++n) {
            VWsnL_In = e_sq2s*ZLT_Zne(n,I);
            b22 = PV.B22(mu, p2, Msn[I], Mse[n]);
            PiT += 4.0*VWsnL_In.abs2()*b22;
        }
        VWWsnsn_II = e2_2s2;
        a0 = PV.A0(mu, Msn[I]);
        PiT += VWWsnsn_II*a0;
    }
    for (int n=0; n<6; ++n) {
        VWWLL_nn = complex(0.0, 0.0, false);
        for (int I=0; I<6; ++I)
            VWWLL_nn += e2_2s2*ZL(I,n)*ZL(I,n).conjugate();
        a0 = PV.A0(mu, Mse[n]);
        PiT += VWWLL_nn*a0;
    }

    /* squark loops (no CKM) */
    complex VWUD_mn, VWWDD_nn, VWWUU_nn;
    matrix<complex> ZDT_ZU = ZD.transpose()*ZU;
    for (int n=0; n<6; ++n) {
        for (int m=0; m<6; ++m) {
            /* VWDU^dagger = VWUD */
            VWUD_mn = e_sq2s*ZDT_ZU(n,m);
            b22 = PV.B22(mu, p2, Msu[m], Msd[n]);
            PiT += 4.0*Nc*VWUD_mn.abs2()*b22;
        }
        VWWDD_nn = complex(0.0, 0.0, false);
        for (int I=0; I<6; ++I)
            VWWDD_nn += e2_2s2*ZD(I,n)*ZD(I,n).conjugate();
        a0 = PV.A0(mu, Msd[n]);
        PiT += Nc*VWWDD_nn*a0;
        VWWUU_nn = complex(0.0, 0.0, false);
        for (int I=0; I<6; ++I)
            VWWUU_nn += e2_2s2*ZU(I,n).conjugate()*ZU(I,n);
        a0 = PV.A0(mu, Msu[n]);
        PiT += Nc*VWWUU_nn*a0;
    }

    /* chargino - neutralino loops */
    complex cV_Wij, cV_Wji, cA_Wij, cA_Wji;
    for (int i=0; i<2; ++i)
        for (int j=0; j<4; ++j) {
            /* indices: chargino=i, neutralino=j */
            cV_Wij = - e_2s*(  ZN(1,j).conjugate()*Zp(0,i)
                             - ZN(3,j).conjugate()*Zp(1,i)/sqrt(2.0)
                             + ZN(1,j)*Zp(0,i).conjugate()
                             + ZN(2,j)*Zp(1,i).conjugate()/sqrt(2.0) );
            cV_Wji = cV_Wij.conjugate();
            cA_Wij = - e_2s*(  ZN(1,j).conjugate()*Zp(0,i)
                             - ZN(3,j).conjugate()*Zp(1,i)/sqrt(2.0)
                             - ZN(1,j)*Zp(0,i).conjugate()
                             - ZN(2,j)*Zp(1,i).conjugate()/sqrt(2.0) );
            cA_Wji = cA_Wij.conjugate();
            PiT += FA(mu, p2 ,mC[i], mN[j], cV_Wij, cV_Wji, cA_Wij, cA_Wji);
        }

    /* Higgs loops */
    double AM_ij;
    for (int i=0; i<2; ++i) {
        for (int j=0; j<2; ++j) {
            AM_ij = ZR(0,i)*ZH(0,j) - ZR(1,i)*ZH(1,j);
            b22 = PV.B22(mu, p2, mH0[i], mHp[j]);
            PiT += g2sq*AM_ij*AM_ij*b22;
        }
        b22 = PV.B22(mu, p2, mH0[i+2], mHp[i]);
        PiT += g2sq*b22;
    }
    for (int i=0; i<4; ++i) {
        a0 = PV.A0(mu, mH0[i]);
        PiT += g2sq/4.0*a0;
    }
    for (int i=0; i<2; ++i) {
        a0 = PV.A0(mu, mHp[i]);
        PiT += g2sq/2.0*a0;
    }

    /* photon - charged-Goldstone-boson loops */
    b0 = PV.B0(mu, p2, Mw_i, 0.0);
    PiT += - e2*Mw_i*Mw_i*b0;

    /* W-boson - Higgs loops */
    double CR_i;
    for (int i=0; i<2; ++i) {
        CR_i = mySUSY.v1()*ZR(0,i) + mySUSY.v2()*ZR(1,i);
        b0 = PV.B0(mu, p2, Mw_i, mH0[i]);
        PiT += - g2sq*g2sq/4.0*CR_i*CR_i*b0;
    }

    /* Z-boson - charged-Goldstone-boson loops */
    b0 = PV.B0(mu, p2, Mw_i, Mz);
    PiT += - e2*sW2*Mz*Mz*b0;

    /* gauge-boson loops */
    a0 = PV.A0(mu, Mw_i);
    b0 = PV.B0(mu, p2, Mw_i, Mz);
    b22 = PV.B22(mu, p2, Mz, Mw_i);
    PiT += g2sq*cW2*(2.0*PV.A0(mu, Mz) - a0
                     + (4.0*p2 + Mz*Mz + Mw_i*Mw_i)*b0 + 8.0*b22);
    //
    PiT += 3.0*g2sq*a0;
    //
    b0 = PV.B0(mu, p2, Mw_i, 0.0);
    b22 = PV.B22(mu, p2, Mw_i, 0.0);
    PiT += e2*((4.0*p2 + Mw_i*Mw_i)*b0 - a0 + 8.0*b22);
    
    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiT_AZ(const double mu, const double p2, const double Mw_i) const
{
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double mHp[2] = {mySUSY.getMHp(), Mw_i}; /* H^+_i = (H^+, G^+) */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2 = e/sW;
    double e_4sc = e/4.0/sW/cW;
    double e_2sc = 2.0*e_4sc;
    double e2_sc = e2/sW/cW;
    
    complex PiT = complex(0.0, 0.0, false);
    double a0;
    complex b0, b22;

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
        PiT += FA(mu, p2, m_l[I], m_l[I], cV_Aee, cV_Zee, cA_Aee, cA_Zee);

        /* down-type quarks */
        PiT += Nc*FA(mu, p2, m_d[I], m_d[I], cV_Add, cV_Zdd, cA_Add, cA_Zdd);

        /* up-type quarks */
        PiT += Nc*FA(mu, p2, m_u[I], m_u[I], cV_Auu, cV_Zuu, cA_Auu, cA_Zuu);
    }

    /* charged-slepton loops */
    complex VZLL_nn, VAZLL_nn;
    matrix<complex> ZLhc_ZL = ZL.hconjugate()*ZL;
    for (int n=0; n<6; ++n) {
        VZLL_nn = - e_2sc*(ZLhc_ZL(n,n) - 2.0*sW2);
        b22 = PV.B22(mu, p2, Mse[n], Mse[n]);
        /* e^2 --> e */
        PiT += - 4.0*e*VZLL_nn*b22;

        VAZLL_nn = e2_sc*(ZLhc_ZL(n,n) - 2.0*sW2);
        a0 = PV.A0(mu, Mse[n]);
        PiT += VAZLL_nn*a0;
    }

    /* down-type squark loops */
    complex VZDD_nn, VAZDD_nn;
    matrix<complex> ZDhc_ZD = ZD.hconjugate()*ZD;
    for (int n=0; n<6; ++n) {
        VZDD_nn = - e_2sc*(ZDhc_ZD(n,n) - 2.0/3.0*sW2);
        b22 = PV.B22(mu, p2, Msd[n], Msd[n]);
        /* e^2 --> e */
        PiT += - 4.0*e*Nc/3.0*VZDD_nn*b22;

        VAZDD_nn = e2_sc/3.0*(ZDhc_ZD(n,n) - 2.0/3.0*sW2);
        a0 = PV.A0(mu, Msd[n]);
        PiT += Nc*VAZDD_nn*a0;
    }

    /* up-type squark loops */
    complex VZUU_nn, VAZUU_nn;
    matrix<complex> ZUhc_ZU = ZU.hconjugate()*ZU;
    for (int n=0; n<6; ++n) {
        VZUU_nn = e_2sc*(ZUhc_ZU(n,n) - 4.0/3.0*sW2);
        b22 = PV.B22(mu, p2, Msu[n], Msu[n]);
        /* e^2 --> e */
        PiT += 4.0*e*Nc*2.0/3.0*VZUU_nn*b22;

        VAZUU_nn = 2.0*e2_sc/3.0*(ZUhc_ZU(n,n) - 4.0/3.0*sW2);
        a0 = PV.A0(mu, Msu[n]);
        PiT += Nc*VAZUU_nn*a0;
    }

    /* chargino loops */
    complex cV_Aii = e;
    complex cA_Aii = 0.0;
    complex cV_Zii, cA_Zii;
    for (int i=0; i<2; ++i) {
        cV_Zii = e_4sc*( Zp(0,i).conjugate()*Zp(0,i)
                         + Zm(0,i)*Zm(0,i).conjugate() + 2.0*(cW2 - sW2) );
        cA_Zii = e_4sc*( Zp(0,i).conjugate()*Zp(0,i)
                         - Zm(0,i)*Zm(0,i).conjugate() );
        PiT += FA(mu, p2 ,mC[i], mC[i], cV_Aii, cV_Zii, cA_Aii, cA_Zii);
    }

    /* W-boson - charged-Goldstone-boson loops */
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    PiT += 2.0*e2*cW*sW*Mz*Mz*b0;

    /* W-boson loops */
    a0 = PV.A0(mu, Mw_i);
    //b0 = PV.B0(mu, p2, Mw_i, Mw_i); /* same as the above */
    b22 = PV.B22(mu, p2, Mw_i, Mw_i);
    PiT += 2.0*e*g2*cW*(2.0*a0 + (2.0*p2 + Mw_i*Mw_i)*b0 + 4.0*b22);

    /* charged-Higgs loops */
    double cot_2thW = (cW2 - sW2)/(2.0*sW*cW);
    for (int i=0; i<2; ++i) {
        a0 = PV.A0(mu, mHp[i]);
        b22 = PV.B22(mu, p2, mHp[i], mHp[i]);
        PiT += 2.0*e2*cot_2thW*(2.0*b22 + a0);
    }

    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiTp_A(const double mu, const double p2, const double Mw_i) const
{
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double mHp[2] = {mySUSY.getMHp(), Mw_i}; /* H^+_i = (H^+, G^+) */
    
    complex PiTp = complex(0.0, 0.0, false);
    complex b0, b0p, b22p;

    /* SM fermion loops */
    complex cV_Aee = - e;
    complex cA_Aee = 0.0;
    complex cV_Add = - e/3.0;
    complex cA_Add = 0.0;
    complex cV_Auu = 2.0/3.0*e;
    complex cA_Auu = 0.0;
    for (int I=0; I<3; ++I) {
        /* charged leptons */
        PiTp += dFA(mu, p2, m_l[I], m_l[I], cV_Aee, cV_Aee, cA_Aee, cA_Aee);

        /* down-type quarks */
        PiTp += Nc*dFA(mu, p2, m_d[I], m_d[I], cV_Add, cV_Add, cA_Add, cA_Add);

        /* up-type quarks */
        PiTp += Nc*dFA(mu, p2, m_u[I], m_u[I], cV_Auu, cV_Auu, cA_Auu, cA_Auu);
    }

    /* charged-slepton loops */
    for (int n=0; n<6; ++n) {
        b22p = PV.B22p(mu, p2, Mse[n], Mse[n]);
        PiTp += 4.0*e2*b22p;
    }

    /* down-type squark loops */
    for (int n=0; n<6; ++n) {
        b22p = PV.B22p(mu, p2, Msd[n], Msd[n]);
        PiTp += 4.0*e2*Nc/3.0/3.0*b22p;
    }

    /* up-type squark loops */
    for (int n=0; n<6; ++n) {
        b22p = PV.B22p(mu, p2, Msu[n], Msu[n]);
        PiTp += 4.0*e2*Nc*2.0/3.0*2.0/3.0*b22p;
    }

    /* chargino loops */
    complex cV_Aii = e;
    complex cA_Aii = 0.0;
    for (int i=0; i<2; ++i)
        PiTp += FA(mu, p2 ,mC[i], mC[i], cV_Aii, cV_Aii, cA_Aii, cA_Aii);

    /* charged-Higgs loops */
    for (int i=0; i<2; ++i) {
        b22p = PV.B22p(mu, p2, mHp[i], mHp[i]);
        PiTp += 4.0*e2*b22p;
    }

    /* W-boson loops */
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    b0p = PV.B0p(mu, p2, Mw_i, Mw_i);
    b22p = PV.B22(mu, p2, Mw_i, Mw_i);
    PiTp += 4.0*e2*( (p2 + 2.0*Mw_i*Mw_i)*b0p + b0 + 2.0*b22p);

    return ( PiTp/16.0/M_PI/M_PI );
}

double EWSUSY::PiThat_W_0(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i;

    double Mz = mySUSY.getMz();
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);

    double PiThat = 0.0;

    /* see e.g. Eqs. (3.16), (3.17), (3.23) in [Hollik, Fortsch.Phys. 38 (1990) 165] */

    /* W self-energy */
    PiThat += PiT_W(mu, 0.0, Mw_i).real();

    /* W-mass counter term */
    double delMw2 = PiT_W(mu, Mw_i*Mw_i, Mw_i).real();
    PiThat -= delMw2;

    /* W wave-function counter term */
    double delZ2w = - PiTp_A(mu, 0.0, Mw_i).real()
                    + 2.0*cW/sW*PiT_AZ(mu, 0.0, Mw_i).real()/Mz/Mz
                    - cW2/sW2*(PiT_W(mu, Mw_i*Mw_i, Mw_i).real()/Mw_i/Mw_i
                               - PiT_Z(mu, Mz*Mz, Mw_i).real()/Mz/Mz);
    PiThat -= delZ2w*Mw_i*Mw_i;

    return PiThat;
}

double EWSUSY::DeltaR_rem_SM(const double Mw_i) const
{
    double cW2 = Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    double sW2 = 1.0 - cW2;

    /* renormalized vertex corrections + box*/
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
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    complex FF = F(Mse[k], Msn[K], mC[i], mN[j]);
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
                           *H(Mse[k], Mse[l], mC[i], mN[j]);
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a11 +=  L_esnC(M, K, j, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, L, i, Mw_i)
                           *L_esnC(N, L, j, Mw_i).conjugate()
                           *H(Msn[K], Msn[L], mC[i], mN[j]);
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {
                    a12 += 0.5
                           *L_esnC(M, K, i, Mw_i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *L_esnC(N, K, j, Mw_i).conjugate()
                           *mC[i]*mC[j]*F(Mse[k], Msn[K], mC[i], mC[j]);
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    a12 += 0.5
                           *L_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F(Mse[k], Msn[K], mN[i], mN[j]);
                    a12 +=  L_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *H(Mse[k], Msn[K], mN[i], mN[j]);
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
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    complex HH = H(Mse[k], Msn[K], mC[i], mN[j]);
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
                             *H(Mse[k], Mse[l], mC[i], mN[j]);
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    a21 += - 2.0
                             *R_esnC(M, K, i)
                             *L_nsnN(J, K, j, Mw_i).conjugate()
                             *L_nsnN(I, L, j, Mw_i)
                             *R_esnC(N, L, i).conjugate()
                             *H(Msn[K], Msn[L], mC[i], mN[j]);
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    a22 += R_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F(Mse[k], Msn[K], mN[i], mN[j]);
                    a22 += 2.0
                           *R_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *H(Mse[k], Msn[K], mN[i], mN[j]);
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {
                    a22 += 2.0
                           *R_esnC(M, K, i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *R_esnC(N, K, j).conjugate()
                           *H(Mse[k], Msn[K], mC[i], mC[j]);
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
                        - ZN(3,i)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,i).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu, 0.0, mC[i], mN[j]);
                ff = f(Mse[k], mC[i], mN[j]);
                v += L_nLC(intJ, k, i, Mw_i).conjugate()*L_eLN(intM, k, j, Mw_i)
                     *( sqrt(2.0)*CL_ji*mC[i]*mN[j]*ff
                        - CR_ji/sqrt(2.0)*(b0 - 0.5 + Mse[k]*Mse[k]*ff) );
            }

    /* sneutrino - neutralino - chargino loops */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int j=0; j<4; ++j)
            for (int i=0; i<2; ++i) {
                CL_ji = ZN(1,j)*Zp(0,i).conjugate()
                        - ZN(3,i)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,i).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu, 0.0, mC[i], mN[j]);
                ff = f(Msn[K], mC[i], mN[j]);
                v += L_nsnN(intJ, K, j, Mw_i).conjugate()*L_esnC(intM, K, i, Mw_i)
                     *( - sqrt(2.0)*CR_ji*mC[i]*mN[j]*ff
                        + CL_ji/sqrt(2.0)*(b0 - 0.5 + Msn[K]*Msn[K]*ff) );
            }
    
    /* sneutrino - charged-slepton - neutralino loops */
    matrix<complex> ZneT_ZL = Zne.transpose()*ZL;
    for (int i=0; i<6; ++i)
        for (int j=0; j<4; ++j)
            for (int K=0; K<3; ++K) {  /* K=0-3 for left-handed sneutrinos */
                b0 = PV.B0(mu, 0.0, Mse[i], Msn[K]);
                ff = f(mN[j], Mse[i], Msn[K]);
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
    double muIR = mu; /* fictional scale, since B0p(0,m1,m2) is IR finite */
    complex b0p, b0;

    /* charged-slepton - neutralino loops */
    for (int k=0; k<6; ++k)
        for (int j=0; j<4; ++j) {
            b0p = PV.B0p(muIR, 0.0, Mse[k], mN[j]);
            b0 = PV.B0(mu, 0.0, Mse[k], mN[j]);
            delv += 0.5*L_eLN(intM, k, j, Mw_i)*L_eLN(intJ, k, j, Mw_i).conjugate()
                    *( (Mse[k]*Mse[k] - mN[j]*mN[j])*b0p - b0 );
        }

    /* sneutrino - chargino loops */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int i=0; i<2; ++i) {
            b0p = PV.B0p(muIR, 0.0, Msn[K], mC[i]);
            b0 = PV.B0(mu, 0.0, Msn[K], mC[i]);
            delv += 0.5*L_esnC(intM, K, i, Mw_i)*L_esnC(intJ, K, i, Mw_i).conjugate()
                    *( (Msn[K]*Msn[K] - mC[i]*mC[i])*b0p - b0 );
        }

    return ( delv/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_vertex_SUSY(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i;

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
            b0p = PV.B0p(muIR, 0.0, Mse[k], mC[i]);
            b0 = PV.B0(mu, 0.0, Mse[k], mC[i]);
            Sigma += 0.5*L_nLC(intI, k, i, Mw_i)*L_nLC(intJ, k, i, Mw_i).conjugate()
                     *( (Mse[k]*Mse[k] - mC[i]*mC[i])*b0p - b0 );
        }

    /* sneutrino - neutralino loops */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int j=0; j<4; ++j) {
            b0p = PV.B0p(muIR, 0.0, Msn[K], mN[j]);
            b0 = PV.B0(mu, 0.0, Msn[K], mN[j]);
            Sigma += 0.5*L_nsnN(intI, K, j, Mw_i)*L_nsnN(intJ, K, j, Mw_i).conjugate()
                     *( (Msn[K]*Msn[K] - mN[j]*mN[j])*b0p - b0 );
        }

    return ( Sigma/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_neutrino_SUSY(const double Mw_i) const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = Mw_i;
    
    return ( ( Sigma_nu_0(mu, mySUSY.NEUTRINO_1, mySUSY.NEUTRINO_1, Mw_i).real()
              - delta_v(mu, mySUSY.ELECTRON, mySUSY.NEUTRINO_1, Mw_i).real()
              + Sigma_nu_0(mu, mySUSY.NEUTRINO_2, mySUSY.NEUTRINO_2, Mw_i).real()
              - delta_v(mu, mySUSY.MU, mySUSY.NEUTRINO_2, Mw_i).real() )/2.0 );
}

double EWSUSY::DeltaR_MSSM_EW1(const double Mw_i) const
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

    return DeltaR;
}

double EWSUSY::DeltaAlphaL5q_SM_EW1() const
{
    /* Renormalization scale (varied for checking the cancellation of UV divergences */
    double mu = mySUSY.getMz();

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
    //std::cout << "EWSM(l)   " << myEWSMOneLoopEW.DeltaAlpha_l(Mz2) << std::endl;
    //std::cout << "EWSUSY(q) " << (DelA_d + DelA_u)/16.0/M_PI/M_PI << std::endl;
    //std::cout << "EWSM(had) " << myEWSMOneLoopEW.DeltaAlpha_5q(Mz2) << std::endl;

    return ( (DelA_l + DelA_d + DelA_u)/16.0/M_PI/M_PI );
}

double EWSUSY::DeltaR_SUSY_EW1(const double Mw_i) const
{
    double cW2 = Mw_i*Mw_i/mySUSY.getMz()/mySUSY.getMz();
    double sW2 = 1.0 - cW2;

    /* SM one-loop contributions */
    double DeltaAlphaL5q_EW1 = DeltaAlphaL5q_SM_EW1();
    double DeltaRho_EW1 = - cW2/sW2*myEWSMOneLoopEW.DeltaRho(Mw_i);
    double DeltaR_rem_EW1 = myEWSMOneLoopEW.DeltaR_rem(Mw_i);
    double DeltaR_SM_EW1 = DeltaAlphaL5q_EW1 - cW2/sW2*DeltaRho_EW1 + DeltaR_rem_EW1;

    return ( DeltaR_MSSM_EW1(Mw_i) - DeltaR_SM_EW1 );
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
    return PV.D0(0.0, 0.0, m1, m2, m3, m4);
}

complex EWSUSY::H(const double m1, const double m2, const double m3,
                  const double m4) const
{
    if ( fabs(m1 - m2) > 1.0e-7 ) {
        double m12 = m1*m1;
        double m22 = m2*m2;
        return ( 0.25/(m12 - m22)*(m12*f(m1, m3, m4) - m22*f(m2, m3, m4)) );
    } else
        return  ( 0.5*f(m1, m3, m4) );
}

complex EWSUSY::f(const double m1, const double m2, const double m3) const
{
    return ( - PV.C0(0.0, m1, m2, m3) );
}


