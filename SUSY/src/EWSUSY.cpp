/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSUSY.h"

EWSUSY::EWSUSY(const SUSY& SUSY_in)
: mySUSY(SUSY_in),
        Yu(3,3,0.0), Yd(3,3,0.0), Yl(3,3,0.0),
        Au(3,3,0.0), Ad(3,3,0.0), Al(3,3,0.0),
        Zm(2,2,0.), Zp(2,2,0.), Zn(4,4,0.),
        ZU(6,6,0.), ZD(6,6,0.), ZL(6,6,0.)
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
    Zn = mySUSY.getN().hconjugate();
    ZU = mySUSY.getRu().hconjugate();
    ZD = mySUSY.getRd().transpose();
    ZL = mySUSY.getRl().transpose();
}

complex EWSUSY::FA(const double mu, const double p2,
                   const double mi, const double mj,
                   const double cV_aij, const double cV_bji,
                   const double cA_aij, const double cA_bji) const
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

complex EWSUSY::PiT_Z(const double mu, const double p2, const double Mw) const
{
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double cW2 = Mw*Mw/Mz2;
    double sW2 = 1.0 - cW2;
    double e2 = 4.0*M_PI*mySUSY.getAle();
    double g2squared = e2/sW2;

    complex PiT = complex(0.0, 0.0, false);




    return PiT;
}


