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
: mySUSY(SUSY_in),
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

complex EWSUSY::PiT_Z(const double mu, const double p2, const double Mw_i) const
{
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double e = sqrt(e2);
    double Mz = mySUSY.getMz();
    double Nc = mySUSY.getNc();

    /* variables depending on Mw_i */
    double cW = Mw_i/Mz;
    double cW2 = cW*cW;
    double sW2 = 1.0 - cW2;
    double sW = sqrt(sW2);
    double g2sq = e2/sW2; /* g2 squared */
    double e_4sc = e/4.0/sW/cW;
    double e_2sc = 2.0*e_4sc;

    complex PiT = complex(0.0, 0.0, false);
    double mI, mm, mn, mi, mj;
    double a0;
    complex b0, b22;
    complex VZss, VZZss, cV_Zij, cV_Zji, cA_Zij, cA_Zji;
    matrix<double> Id6 = matrix<double>::Id(6);

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
        mI = mySUSY.Ml_Q((StandardModel::lepton)(2*I + 1));
        PiT += FA(mu, p2 ,mI, mI, cV_Zee, cV_Zee, cA_Zee, cA_Zee);

        /* down-type quarks */
        mI = mySUSY.Mq_Q((StandardModel::quark)(2*I + 1));
        PiT += Nc*FA(mu, p2 ,mI, mI, cV_Zdd, cV_Zdd, cA_Zdd, cA_Zdd);

        /* up-type quarks */
        mI = mySUSY.Mq_Q((StandardModel::quark)(2*I));
        PiT += Nc*FA(mu, p2 ,mI, mI, cV_Zuu, cV_Zuu, cA_Zuu, cA_Zuu);
    }

    /* sneutrino loops */
    VZss = e_2sc;
    VZZss = e2/2.0/sW2/cW2;
    for (int I=0; I<3; ++I) {  /* I=0-3 for left-handed sneutrinos */
        mI = sqrt(mySUSY.getMsn2()(I));

        b22 = PV.B22(mu, p2, mI, mI);
        PiT += 4.0*VZss.abs2()*b22;

        a0 = PV.A0(mu, mI);
        PiT += VZZss*a0;
    }

    /* charged slepton loops */
    matrix<complex> ZLhc_ZL = ZL.hconjugate()*ZL; 
    for (int m=0; m<6; ++m) {
        mm = sqrt(mySUSY.getMse2()(m));
        for (int n=0; n<6; ++n) {
            mn = sqrt(mySUSY.getMse2()(n));
            VZss = - e_2sc*(ZLhc_ZL(m,n) - 2.0*sW2*Id6(m,n));
            b22 = PV.B22(mu, p2, mm, mn);
            PiT += 4.0*VZss.abs2()*b22;
        }
        VZZss = 2.0*e2/cW2*(sW2 + (1.0 - 4.0*sW2)/4.0/sW2*ZLhc_ZL(m,m));
        a0 = PV.A0(mu, mm);
        PiT += VZZss*a0;
    }

    /* down-type squark loops */
    matrix<complex> ZDhc_ZD = ZD.hconjugate()*ZD; 
    for (int m=0; m<6; ++m) {
        mm = sqrt(mySUSY.getMsd2()(m));
        for (int n=0; n<6; ++n) {
            VZss = - e_2sc*(ZDhc_ZD(m,n) - 2.0/3.0*sW2*Id6(m,n));
            mn = sqrt(mySUSY.getMsd2()(n));
            b22 = PV.B22(mu, p2, mm, mn);
            PiT += 4.0*Nc*VZss.abs2()*b22;
        }
        VZZss = 2.0*e2/3.0/cW2*(sW2/3.0 + (3.0 - 4.0*sW2)/4.0/sW2*ZDhc_ZD(m,m));
        a0 = PV.A0(mu, mm);
        PiT += Nc*VZZss*a0;
    }

    /* up-type squark loops */
    matrix<complex> ZUhc_ZU = ZU.hconjugate()*ZU;
    for (int m=0; m<6; ++m) {
        mm = sqrt(mySUSY.getMsu2()(m));
        for (int n=0; n<6; ++n) {
            /* (n,m) instead of (m,n) */
            VZss = e_2sc*(ZUhc_ZU(n,m) - 4.0/3.0*sW2*Id6(m,n));
            mn = sqrt(mySUSY.getMsu2()(n));
            b22 = PV.B22(mu, p2, mm, mn);
            PiT += 4.0*Nc*VZss.abs2()*b22;
        }
        VZZss = 2.0*e2/3.0/cW2*(4.0*sW2/3.0 + (3.0 - 8.0*sW2)/4.0/sW2*ZUhc_ZU(m,m));
        a0 = PV.A0(mu, mm);
        PiT += Nc*VZZss*a0;
    }

    /* chargino loops */
    matrix<double> Id2 = matrix<double>::Id(2);
    for (int i=0; i<2; ++i) {
        mi = mySUSY.getMch()(i);
        for (int j=0; j<2; ++j) {
            mj = mySUSY.getMch()(j);
            cV_Zij = e_4sc*(Zp(0,j).conjugate()*Zp(0,i) 
                            + Zm(0,j)*Zm(0,i).conjugate() + 2.0*(cW2 - sW2)*Id2(j,i));
            cV_Zji = e_4sc*(Zp(0,i).conjugate()*Zp(0,j)
                            + Zm(0,i)*Zm(0,j).conjugate() + 2.0*(cW2 - sW2)*Id2(i,j));
            cA_Zij = e_4sc*(Zp(0,j).conjugate()*Zp(0,i)
                            - Zm(0,j)*Zm(0,i).conjugate());
            cA_Zji = e_4sc*(Zp(0,i).conjugate()*Zp(0,j)
                            - Zm(0,i)*Zm(0,j).conjugate());
            PiT += FA(mu, p2 ,mi, mj, cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }
    }

    /* neutralino loops */
    for (int i=0; i<4; ++i) {
        mi = mySUSY.getMneu()(i);
        for (int j=0; j<4; ++j) {
            mj = mySUSY.getMneu()(j);
            cV_Zij = e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                            - ZN(2,j).conjugate()*ZN(2,i)
                            - ZN(3,j)*ZN(3,i).conjugate()
                            + ZN(2,j)*ZN(2,i).conjugate());
            cV_Zji = e_4sc*(  ZN(3,i).conjugate()*ZN(3,j)
                            - ZN(2,i).conjugate()*ZN(2,j)
                            - ZN(3,i)*ZN(3,j).conjugate()
                            + ZN(2,i)*ZN(2,j).conjugate());
            cA_Zij = e_4sc*(  ZN(3,j).conjugate()*ZN(3,i)
                            - ZN(2,j).conjugate()*ZN(2,i)
                            + ZN(3,j)*ZN(3,i).conjugate()
                            - ZN(2,j)*ZN(2,i).conjugate());
            cA_Zji = e_4sc*(  ZN(3,i).conjugate()*ZN(3,j)
                            - ZN(2,i).conjugate()*ZN(2,j)
                            + ZN(3,i)*ZN(3,j).conjugate()
                            - ZN(2,i)*ZN(2,j).conjugate());
            PiT += 0.5*FA(mu, p2 ,mi, mj, cV_Zij, cV_Zji, cA_Zij, cA_Zji);
        }
    }
    
    /* charged-Higgs loops */
    double cot_2thW = (cW2 - sW2)/(2.0*sW*cW);
    mi = mySUSY.getMHp();
    b22 = PV.B22(mu, p2, mi, mi);
    a0 = PV.A0(mu, mi);
    PiT += 4.0*e2*cot_2thW*cot_2thW*(2.0*b22 + a0);

    /* neutral-Higgs loops */
    double AM_ij;
    for (int i=0; i<2; ++i) {
        if (i==0) mi = mySUSY.getMHh();
        else mi = mySUSY.getMHl();
        for (int j=0; j<2; ++j) {
            if (j==0) mj = mySUSY.getMHa();
            else mj = Mz; /* mass of the neutral Goldstone boson */
            AM_ij = ZR(0,i)*ZH(0,j) - ZR(1,i)*ZH(1,j);
            b22 = PV.B22(mu, p2, mi, mj);
            PiT += g2sq/cW2*AM_ij*AM_ij*b22;
        }
    }
    for (int j=0; j<4; ++j) {
        mj = mySUSY.getMH0()(j);
        a0 = PV.A0(mu, mj);
        PiT += g2sq/4.0/cW2*a0;
    }

    /* W-boson loops */
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    PiT += - 2.0*g2sq*sW2*sW2*Mz*Mz*b0;
    a0 = PV.A0(mu, Mw_i);
    b0 = PV.B0(mu, p2, Mw_i, Mw_i);
    b22 = PV.B22(mu, p2, Mw_i, Mw_i);
    PiT += 2.0*g2sq*cW2*(2.0*a0*a0 + (2.0*p2 + Mw_i*Mw_i)*b0 + 4.0*b22);

    /* Z-boson and Higgs loops */
    double CR_i;
    for (int i=0; i<2; ++i) {
        if (i==0) mi = mySUSY.getMHh();
        else mi = mySUSY.getMHl();
        CR_i = mySUSY.v1()*ZR(0,i) + mySUSY.v2()*ZR(1,i);
        b0 = PV.B0(mu, p2, mySUSY.getMz(), mi);
        PiT += - g2sq/4.0/cW2/cW2*CR_i*CR_i*b0;
    }
    
    return ( PiT/16.0/M_PI/M_PI );
}


complex EWSUSY::PiT_W(const double mu, const double p2, const double Mw_i) const
{
    complex PiT = complex(0.0, 0.0, false);


    /* Write codes! */

    
    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiT_AZ(const double mu, const double p2, const double Mw_i) const
{
    complex PiT = complex(0.0, 0.0, false);


    /* Write codes! */


    return ( PiT/16.0/M_PI/M_PI );
}

complex EWSUSY::PiTp_A(const double mu, const double p2, const double Mw_i) const
{
    complex PiTp = complex(0.0, 0.0, false);


    /* Write codes! */


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
    complex F, H;
    double ML[6] = {sqrt(mySUSY.getMse2()(0)), sqrt(mySUSY.getMse2()(1)),
                    sqrt(mySUSY.getMse2()(2)), sqrt(mySUSY.getMse2()(3)),
                    sqrt(mySUSY.getMse2()(4)), sqrt(mySUSY.getMse2()(5))};
    double Msn[3] = {sqrt(mySUSY.getMsn2()(0)), sqrt(mySUSY.getMsn2()(1))};
    double mC[2] = {mySUSY.getMch()(0), mySUSY.getMch()(1)};
    double mN[4] = {mySUSY.getMneu()(0), mySUSY.getMneu()(1),
                    mySUSY.getMneu()(2), mySUSY.getMneu()(3)};

    /* charged-lepton - sneutrino - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {
                    F = PV.D0(0.0, 0.0, ML[k], Msn[K], mC[i], mN[j]);
                    a11 += 0.5
                           *L_esnC(M, K, i, Mw_i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *mC[i]*mN[j]*F;
                    a11 += 0.5
                           *L_eLN(M, k, j, Mw_i)
                           *L_nsnN(I, K, j, Mw_i)
                           *L_nLC(J, k, i, Mw_i).conjugate()
                           *L_esnC(N, K, i, Mw_i).conjugate()
                           *mC[i]*mN[j]*F;
                }
    
    /* charged-lepton - charged-lepton - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int l=0; l<6; ++l)
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {


                    H = 0.0; /* Write codes! */


                    a11 +=  L_eLN(M, k, j, Mw_i)
                           *L_nLC(J, k, i, Mw_i).conjugate()
                           *L_nLC(I, l, i, Mw_i)
                           *L_eLN(N, l, j, Mw_i).conjugate()
                           *H;
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {


                    H = 0.0; /* Write codes! */

                    
                    a11 +=  L_esnC(M, K, j, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, L, i, Mw_i)
                           *L_esnC(N, L, j, Mw_i).conjugate()
                           *H;
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {
                    F = PV.D0(0.0, 0.0, ML[k], Msn[K], mC[i], mC[j]);

                    a12 += 0.5
                           *L_esnC(M, K, i, Mw_i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *L_esnC(N, K, j, Mw_i).conjugate()
                           *mC[i]*mC[j]*F;
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    F = PV.D0(0.0, 0.0, ML[k], Msn[K], mN[i], mN[j]);
                    a12 += 0.5
                           *L_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F;


                    H = 0.0; /* Write codes! */


                    a12 +=  L_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *L_eLN(N, k, j, Mw_i).conjugate()
                           *H;
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
    complex F, H;
    double ML[6] = {sqrt(mySUSY.getMse2()(0)), sqrt(mySUSY.getMse2()(1)),
                    sqrt(mySUSY.getMse2()(2)), sqrt(mySUSY.getMse2()(3)),
                    sqrt(mySUSY.getMse2()(4)), sqrt(mySUSY.getMse2()(5))};
    double Msn[3] = {sqrt(mySUSY.getMsn2()(0)), sqrt(mySUSY.getMsn2()(1))};
    double mC[2] = {mySUSY.getMch()(0), mySUSY.getMch()(1)};
    double mN[4] = {mySUSY.getMneu()(0), mySUSY.getMneu()(1),
                    mySUSY.getMneu()(2), mySUSY.getMneu()(3)};

    /* charged-lepton - sneutrino - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {


                    H = 0.0; /* Write codes! */


                    a21 += - 2.0
                             *R_esnC(M, K, i)
                             *L_nLC(I, k, i, Mw_i)
                             *L_nsnN(J, K, j, Mw_i).conjugate()
                             *R_eLN(N, k, j, Mw_i).conjugate()
                             *H;
                    a21 += - 2.0
                             *R_eLN(M, k, j, Mw_i)
                             *L_nsnN(I, K, j, Mw_i)
                             *L_nLC(J, k, i, Mw_i).conjugate()
                             *R_esnC(N, K, i).conjugate()
                             *H;
                }

    /* charged-lepton - charged-lepton - chargino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int l=0; l<6; ++l)
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {


                    H = 0.0; /* Write codes! */


                    a21 += - 2.0
                             *R_eLN(M, k, j, Mw_i)
                             *L_nLC(J, k, i, Mw_i).conjugate()
                             *L_nLC(I, l, i, Mw_i)
                             *R_eLN(N, l, j, Mw_i).conjugate()
                             *H;
                }

    /* sneutrino - sneutrino - chargino - neutralino loop */
    for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
        for (int L=0; L<3; ++L)  /* L=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<4; ++j) {


                    H = 0.0; /* Write codes! */


                    a21 += - 2.0
                             *R_esnC(M, K, i)
                             *L_nsnN(J, K, j, Mw_i).conjugate()
                             *L_nsnN(I, L, j, Mw_i)
                             *R_esnC(N, L, i).conjugate()
                             *H;
                }

    /* charged-lepton - sneutrino - neutralino - neutralino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j) {
                    F = PV.D0(0.0, 0.0, ML[k], Msn[K], mN[i], mN[j]);
                    a22 += R_eLN(M, k, i, Mw_i)
                           *L_nsnN(J, K, i, Mw_i).conjugate()
                           *L_nsnN(I, K, j, Mw_i)
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *mN[i]*mN[j]*F;

                    
                    H = 0.0; /* Write codes! */


                    a22 += 2.0
                           *R_eLN(M, k, i, Mw_i)
                           *L_nsnN(I, K, i, Mw_i)
                           *L_nsnN(J, K, j, Mw_i).conjugate()
                           *R_eLN(N, k, j, Mw_i).conjugate()
                           *H;
                }

    /* charged-lepton - sneutrino - chargino - chargino loop */
    for (int k=0; k<6; ++k)
        for (int K=0; K<3; ++K)  /* K=0-3 for left-handed sneutrinos */
            for (int i=0; i<2; ++i)
                for (int j=0; j<2; ++j) {


                    H = 0.0; /* Write codes! */

                    
                    a22 += 2.0
                           *R_esnC(M, K, i)
                           *L_nLC(I, k, i, Mw_i)
                           *L_nLC(J, k, j, Mw_i).conjugate()
                           *R_esnC(N, K, j).conjugate()
                           *H;
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
    complex b0, f;
    complex CL_ji, CR_ji; /* chargino-neutralino-W couplings */
    double MLi, MLk, MsnK, mNj, mCi;
    
    /* charged-slepton - chargino - neutralino loops */
    for (int k=0; k<6; ++k) {
        MLk = sqrt(mySUSY.getMse2()(k));
        for (int j=0; j<4; ++j) {
            mNj = mySUSY.getMneu()(j);
            for (int i=0; i<2; ++i) {
                mCi = mySUSY.getMch()(i);
                CL_ji = ZN(1,j)*Zp(0,i).conjugate()
                        - ZN(3,i)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,i).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu, 0.0, mCi, mNj);
                f = - PV.C0(0.0, MLk, mCi, mNj);
                v += L_nLC(intJ, k, i, Mw_i).conjugate()*L_eLN(intM, k, j, Mw_i)
                     *( sqrt(2.0)*CL_ji*mCi*mNj*f
                        - CR_ji/sqrt(2.0)*(b0 - 0.5 + MLk*MLk*f) );
            }
        }
    }

    /* sneutrino - neutralino - chargino loops */
    for (int K=0; K<3; ++K) {  /* K=0-3 for left-handed sneutrinos */
        MsnK = sqrt(mySUSY.getMsn2()(K));
        for (int j=0; j<4; ++j) {
            mNj = mySUSY.getMneu()(j);
            for (int i=0; i<2; ++i) {
                mCi = mySUSY.getMch()(i);
                CL_ji = ZN(1,j)*Zp(0,i).conjugate()
                        - ZN(3,i)*Zp(1,i).conjugate()/sqrt(2.0);
                CR_ji = ZN(1,j).conjugate()*Zm(0,i)
                        + ZN(2,i).conjugate()*Zm(1,i)/sqrt(2.0);
                b0 = PV.B0(mu, 0.0, mCi, mNj);
                f = - PV.C0(0.0, MsnK, mCi, mNj);
                v += L_nsnN(intJ, K, j, Mw_i).conjugate()*L_esnC(intM, K, i, Mw_i)
                     *( - sqrt(2.0)*CR_ji*mCi*mNj*f
                        + CL_ji/sqrt(2.0)*(b0 - 0.5 + MsnK*MsnK*f) );
            }
        }
    }
    
    /* sneutrino - charged-slepton - neutralino loops */
    matrix<complex> ZneT_ZL = Zne.transpose()*ZL;
    for (int i=0; i<6; ++i) {
        MLi = sqrt(mySUSY.getMse2()(i));
        for (int j=0; j<4; ++j) {
            mNj = mySUSY.getMneu()(j);
            for (int K=0; K<3; ++K) {  /* K=0-3 for left-handed sneutrinos */
                MsnK = sqrt(mySUSY.getMsn2()(K));
                b0 = PV.B0(mu, 0.0, MLi, MsnK);
                f = - PV.C0(0.0, mNj, MLi, MsnK);
                v += 0.5*L_nsnN(intJ, K, j, Mw_i).conjugate()*L_eLN(intM, i, j, Mw_i)
                     *ZneT_ZL(K, i).conjugate()*(b0 + 0.5 + mNj*mNj*f); 
            }
        }
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
    double MLk, mNj;
    for (int k=0; k<6; ++k) {
        MLk = sqrt(mySUSY.getMse2()(k));
        for (int j=0; j<4; ++j) {
            mNj = mySUSY.getMneu()(j);
            b0p = PV.B0p(muIR, 0.0, MLk, mNj);
            b0 = PV.B0(mu, 0.0, MLk, mNj);
            delv += 0.5*L_eLN(intM, k, j, Mw_i)*L_eLN(intJ, k, j, Mw_i).conjugate()
                    *( (MLk*MLk - mNj*mNj)*b0p - b0 );
        }
    }

    /* sneutrino - chargino loops */
    double Msnk, mCi;
    for (int K=0; K<3; ++K) {  /* K=0-3 for left-handed sneutrinos */
        Msnk = sqrt(mySUSY.getMsn2()(K));
        for (int i=0; i<2; ++i) {
            mCi = mySUSY.getMch()(i);
            b0p = PV.B0p(muIR, 0.0, Msnk, mCi);
            b0 = PV.B0(mu, 0.0, Msnk, mCi);
            delv += 0.5*L_esnC(intM, K, i, Mw_i)*L_esnC(intJ, K, i, Mw_i).conjugate()
                    *( (Msnk*Msnk - mCi*mCi)*b0p - b0 );
        }
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
    double MLk, mCi;
    for (int k=0; k<6; ++k) {
        MLk = sqrt(mySUSY.getMse2()(k));
        for (int i=0; i<2; ++i) {
            mCi = mySUSY.getMch()(i);            
            b0p = PV.B0p(muIR, 0.0, MLk, mCi);
            b0 = PV.B0(mu, 0.0, MLk, mCi);
            Sigma += 0.5*L_nLC(intI, k, i, Mw_i)*L_nLC(intJ, k, i, Mw_i).conjugate()
                     *( (MLk*MLk - mCi*mCi)*b0p - b0 );
        }
    }

    /* sneutrino - neutralino loops */
    double Msnk, mNj;
    for (int K=0; K<3; ++K) {  /* K=0-3 for left-handed sneutrinos */
        Msnk = sqrt(mySUSY.getMsn2()(K));
        for (int j=0; j<4; ++j) {
            mNj = mySUSY.getMneu()(j);
            b0p = PV.B0p(muIR, 0.0, Msnk, mNj);
            b0 = PV.B0(mu, 0.0, Msnk, mNj);
            Sigma += 0.5*L_nsnN(intI, K, j, Mw_i)*L_nsnN(intJ, K, j, Mw_i).conjugate()
                     *( (Msnk*Msnk - mNj*mNj)*b0p - b0 );
        }
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

double EWSUSY::DeltaR_SUSY_EW1(const double Mw_i) const
{
    double DeltaR_SM_EW1; 

    
    /* Write codes for DeltaR_SM_EW1! */
    /* How to take into account the discrepancy in the hVV couplings between SM and MSSM? */


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
    
    return ( e/sqrt(2.0)/sW/cW*Zne(J,K).conjugate()*(ZN(0,j)*sW - ZN(1,j)*cW) );
}

complex EWSUSY::L_eLN(const int N, const int k, const int j, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();
    double sW = sqrt(1.0 - cW*cW);

    return ( e/sqrt(2.0)/sW/cW*ZL(N,k)*(ZN(0,j)*sW + ZN(1,j)*cW)
             + Yl(N,N)*ZL(N+3,k)*ZN(2,j) );
}

complex EWSUSY::R_eLN(const int N, const int k, const int j, const double Mw_i) const
{
    double e = sqrt(4.0*M_PI*mySUSY.getAle());
    double cW = Mw_i/mySUSY.getMz();

    return ( - e*sqrt(2.0)/cW*ZL(N+3,k)*ZN(0,j).conjugate()
             + Yl(N,N)*ZL(N,k)*ZN(2,j).conjugate() );
}



