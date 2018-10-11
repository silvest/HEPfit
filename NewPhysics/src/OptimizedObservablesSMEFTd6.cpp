/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "OptimizedObservablesSMEFTd6.h"
#include "StandardModel.h"
#include "HeffDB1.h"

eeWW::eeWW(const StandardModel& SM_i)
: ThObservable(SM_i)
{
};

double eeWW::computeThValue()
{
//    double energy = 1500;
//    double weight;

    // Get phase space point
//    std::vector<double*> p = get_momenta(ninitial, energy, getMasses(), weight);
    std::vector<double*> pp(1, new double[4]);
    pp[0][0] = 7.500000e+02;
    pp[0][1] = 0.;
    pp[0][2] = 0.;
    pp[0][3] = 7.500000e+02;
    pp.push_back(new double[4]);
    pp[1][0] = 7.500000e+02;
    pp[1][1] = 0.;
    pp[1][2] = 0.;
    pp[1][3] = -7.500000e+02;
    pp.push_back(new double[4]);
    pp[2][0] = 1.328270e+02;
    pp[2][1] = -3.315104e+01;
    pp[2][2] = 6.012053e+01;
    pp[2][3] = -1.137081e+02;
    pp.push_back(new double[4]);
    pp[3][0] = 4.924941e+02;
    pp[3][1] = -1.557744e+02;
    pp[3][2] = -4.529006e+02;
    pp[3][3] = 1.147424e+02;
    pp.push_back(new double[4]);
    pp[4][0] = 2.285372e+02;
    pp[4][1] = -1.588214e+02;
    pp[4][2] = -1.465645e+02;
    pp[4][3] = 7.432258e+01;
    pp.push_back(new double[4]);
    pp[5][0] = 6.461417e+02;
    pp[5][1] = 3.477469e+02;
    pp[5][2] = 5.393446e+02;
    pp[5][3] = -7.535681e+01;
    updateParameters();
    
    for (int i = 0; i < 1000000; i++) {
    

    // Set momenta for this event
    setMomenta(pp);

    // Evaluate matrix element
    sigmaKin();
    }
    return *getMatrixElements();
}

void eeWW::updateParameters()
{

    zero = 0;
    ZERO = 0;

    std::vector<int> indices(2, 0);
    Gamma_H = 6.382339e-03;
    Gamma_W = 2.047600e+00;
    Gamma_Z = 2.441404e+00;
    Gamma_T = 1.491500e+00;
    ymtau = 1.777000e+00;
    ymt = 1.730000e+02;
    ymb = 4.700000e+00;
    aS = 1.180000e-01;
    Gf = 1.166390e-05;
    aEWM1 = 1.325070e+02;
    MH = 1.250000e+02;
    MZ = 9.118800e+01;
    MTA = 1.777000e+00;
    MT = 1.730000e+02;
    MB = 4.700000e+00;
    conjg__CKM3x3 = 1.;
    CKM3x3 = 1.;
    conjg__CKM1x1 = 1.;
    complexi = std::complex<double> (0., 1.);
    MZ2 = MZ*MZ;
    MZ4 = MZ2*MZ2;
    MH2 = MH*MH;
    aEW = 1. / aEWM1;
    MW = sqrt(MZ2 / 2. + sqrt(MZ4 / 4. - (aEW * M_PI * MZ2) / (Gf * M_SQRT2)));
    sqrt__aEW = sqrt(aEW);
    ee = 2. * sqrt__aEW * sqrt(M_PI);
    MW2 = MW*MW;
    sw2 = 1. - MW2 / MZ2;
    cw = sqrt(1. - sw2);
    sqrt__sw2 = sqrt(sw2);
    sw = sqrt__sw2;
    g1 = ee / cw;
    gw = ee / sw;
    vev = (2. * MW * sw) / ee;
    vev2 = vev * vev;
    lam = MH2 / (2. * vev2);
    yb = (ymb * M_SQRT2) / vev;
    yt = (ymt * M_SQRT2) / vev;
    ytau = (ymtau * M_SQRT2) / vev;
    muH = sqrt(lam * vev2);
    I1x33 = yb * conjg__CKM3x3;
    I2x33 = yt * conjg__CKM3x3;
    I3x33 = CKM3x3 * yt;
    I4x33 = CKM3x3 * yb;
    ee2 = ee*ee;
    cw2 = 1. - sw2;


    GC_3 = -(ee * complexi);
    GC_51 = (cw * ee * complexi) / (2. * sw);
    GC_53 = (cw * ee * complexi) / sw;
    GC_59 = (ee * complexi * sw) / (2. * cw);
    GC_100 = (ee * complexi * conjg__CKM1x1) / (sw * M_SQRT2);

    sqrt__aS = sqrt(aS);
    G = 2. * sqrt__aS * sqrt(M_PI);
    G2 = G*G;

    mME.assign(6, 0.);
    jamp2[0] = new double[1];

}

void eeWW::ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double> fi[6])
{
    std::complex<double> chi[2];
    double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2];
    int ip, im, nh;
    fi[0] = std::complex<double> (-p[0] * nsf, -p[3] * nsf);
    fi[1] = std::complex<double> (-p[1] * nsf, -p[2] * nsf);
    nh = nhel * nsf;
    if (fmass != 0.0) {
        pp = std::min(p[0], sqrt(p[1] * p[1] + p[2] * p[2] + p[3] * p[3]));
        if (pp == 0.0) {
            sqm[0] = sqrt(std::abs(fmass));
            sqm[1] = Sgn(sqm[0], fmass);
            ip = (1 + nh) / 2;
            im = (1 - nh) / 2;
            fi[2] = ip * sqm[ip];
            fi[3] = im * nsf * sqm[ip];
            fi[4] = ip * nsf * sqm[im];
            fi[5] = im * sqm[im];
        } else {
            sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5;
            sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5;
            omega[0] = sqrt(p[0] + pp);
            omega[1] = fmass / omega[0];
            ip = (1 + nh) / 2;
            im = (1 - nh) / 2;
            sfomega[0] = sf[0] * omega[ip];
            sfomega[1] = sf[1] * omega[im];
            pp3 = std::max(pp + p[3], 0.0);
            chi[0] = std::complex<double> (sqrt(pp3 * 0.5 / pp), 0);
            if (pp3 == 0.0) {
                chi[1] = std::complex<double> (-nh, 0);
            } else {
                chi[1] = std::complex<double> (nh * p[1], p[2]) / sqrt(2.0 * pp * pp3);
            }
            fi[2] = sfomega[0] * chi[im];
            fi[3] = sfomega[0] * chi[ip];
            fi[4] = sfomega[1] * chi[im];
            fi[5] = sfomega[1] * chi[ip];
        }
    } else {
        if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0) {
            sqp0p3 = 0.0;
        } else {
            sqp0p3 = sqrt(std::max(p[0] + p[3], 0.0)) * nsf;
        }
        chi[0] = std::complex<double> (sqp0p3, 0.0);
        if (sqp0p3 == 0.0) {
            chi[1] = std::complex<double> (-nhel * sqrt(2.0 * p[0]), 0.0);
        } else {
            chi[1] = std::complex<double> (nh * p[1], p[2]) / sqp0p3;
        }
        if (nh == 1) {
            fi[2] = std::complex<double> (0.0, 0.0);
            fi[3] = std::complex<double> (0.0, 0.0);
            fi[4] = chi[0];
            fi[5] = chi[1];
        } else {
            fi[2] = chi[1];
            fi[3] = chi[0];
            fi[4] = std::complex<double> (0.0, 0.0);
            fi[5] = std::complex<double> (0.0, 0.0);
        }
    }
    return;
}

double eeWW::Sgn(double a, double b)
{
    return (b < 0) ? - abs(a) : abs(a);
}

void eeWW::txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double>
        tc[18])
{
    std::complex<double> ft[6][4], ep[4], em[4], e0[4];
    double pt, pt2, pp, pzpt, emp, sqh, sqs;
    int i, j;

    sqh = sqrt(0.5);
    sqs = sqrt(0.5 / 3);

    pt2 = p[1] * p[1] + p[2] * p[2];
    pp = std::min(p[0], sqrt(pt2 + p[3] * p[3]));
    pt = std::min(pp, sqrt(pt2));

    ft[4][0] = std::complex<double> (p[0] * nst, p[3] * nst);
    ft[5][0] = std::complex<double> (p[1] * nst, p[2] * nst);

    // construct eps+
    if (nhel >= 0) {
        if (pp == 0) {
            ep[0] = std::complex<double> (0, 0);
            ep[1] = std::complex<double> (-sqh, 0);
            ep[2] = std::complex<double> (0, nst * sqh);
            ep[3] = std::complex<double> (0, 0);
        } else {
            ep[0] = std::complex<double> (0, 0);
            ep[3] = std::complex<double> (pt / pp * sqh, 0);

            if (pt != 0) {
                pzpt = p[3] / (pp * pt) * sqh;
                ep[1] = std::complex<double> (-p[1] * pzpt, -nst * p[2] / pt * sqh);
                ep[2] = std::complex<double> (-p[2] * pzpt, nst * p[1] / pt * sqh);
            } else {
                ep[1] = std::complex<double> (-sqh, 0);
                ep[2] = std::complex<double> (0, nst * Sgn(sqh, p[3]));
            }
        }

    }

    // construct eps-
    if (nhel <= 0) {
        if (pp == 0) {
            em[0] = std::complex<double> (0, 0);
            em[1] = std::complex<double> (sqh, 0);
            em[2] = std::complex<double> (0, nst * sqh);
            em[3] = std::complex<double> (0, 0);
        } else {
            em[0] = std::complex<double> (0, 0);
            em[3] = std::complex<double> (-pt / pp * sqh, 0);

            if (pt != 0) {
                pzpt = -p[3] / (pp * pt) * sqh;
                em[1] = std::complex<double> (-p[1] * pzpt, -nst * p[2] / pt * sqh);
                em[2] = std::complex<double> (-p[2] * pzpt, nst * p[1] / pt * sqh);
            } else {
                em[1] = std::complex<double> (sqh, 0);
                em[2] = std::complex<double> (0, nst * Sgn(sqh, p[3]));
            }
        }
    }

    // construct eps0
    if (std::labs(nhel) <= 1) {
        if (pp == 0) {
            e0[0] = std::complex<double> (0, 0);
            e0[1] = std::complex<double> (0, 0);
            e0[2] = std::complex<double> (0, 0);
            e0[3] = std::complex<double> (1, 0);
        } else {
            emp = p[0] / (tmass * pp);
            e0[0] = std::complex<double> (pp / tmass, 0);
            e0[3] = std::complex<double> (p[3] * emp, 0);

            if (pt != 0) {
                e0[1] = std::complex<double> (p[1] * emp, 0);
                e0[2] = std::complex<double> (p[2] * emp, 0);
            } else {
                e0[1] = std::complex<double> (0, 0);
                e0[2] = std::complex<double> (0, 0);
            }
        }
    }

    if (nhel == 2) {
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++)
                ft[i][j] = ep[i] * ep[j];
        }
    } else if (nhel == -2) {
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++)
                ft[i][j] = em[i] * em[j];
        }
    } else if (tmass == 0) {
        for (j = 0; j < 4; j++) {
            for (i = 0; i < 4; i++)
                ft[i][j] = 0;
        }
    } else if (tmass != 0) {
        if (nhel == 1) {
            for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++)
                    ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]);
            }
        } else if (nhel == 0) {
            for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++)
                    ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
                        + 2.0 * e0[i] * e0[j]);
            }
        } else if (nhel == -1) {
            for (j = 0; j < 4; j++) {
                for (i = 0; i < 4; i++)
                    ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]);
            }
        } else {
            std::cerr << "Invalid helicity in txxxxx.\n";
            std::exit(1);
        }
    }

    tc[0] = ft[4][0];
    tc[1] = ft[5][0];

    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++)
            tc[j * 4 + i + 2] = ft[j][i];
    }
}

void eeWW::vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double> vc[6])
{
    double hel, hel0, pt, pt2, pp, pzpt, emp, sqh;
    int nsvahl;
    sqh = sqrt(0.5);
    hel = double(nhel);
    nsvahl = nsv * std::abs(hel);
    pt2 = (p[1] * p[1]) + (p[2] * p[2]);
    pp = std::min(p[0], sqrt(pt2 + (p[3] * p[3])));
    pt = std::min(pp, sqrt(pt2));
    vc[0] = std::complex<double> (p[0] * nsv, p[3] * nsv);
    vc[1] = std::complex<double> (p[1] * nsv, p[2] * nsv);
    if (vmass != 0.0) {
        hel0 = 1.0 - std::abs(hel);
        if (pp == 0.0) {
            vc[2] = std::complex<double> (0.0, 0.0);
            vc[3] = std::complex<double> (-hel * sqh, 0.0);
            vc[4] = std::complex<double> (0.0, nsvahl * sqh);
            vc[5] = std::complex<double> (hel0, 0.0);
        } else {
            emp = p[0] / (vmass * pp);
            vc[2] = std::complex<double> (hel0 * pp / vmass, 0.0);
            vc[5] = std::complex<double> (hel0 * p[3] * emp + hel * pt / pp * sqh, 0.0);
            if (pt != 0.0) {
                pzpt = p[3] / (pp * pt) * sqh * hel;
                vc[3] = std::complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
                        p[2] / pt * sqh);
                vc[4] = std::complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
                        p[1] / pt * sqh);
            } else {
                vc[3] = std::complex<double> (-hel * sqh, 0.0);
                vc[4] = std::complex<double> (0.0, nsvahl * Sgn(sqh, p[3]));
            }
        }
    } else {
        pp = p[0];
        pt = sqrt((p[1] * p[1]) + (p[2] * p[2]));
        vc[2] = std::complex<double> (0.0, 0.0);
        vc[5] = std::complex<double> (hel * pt / pp * sqh, 0.0);
        if (pt != 0.0) {
            pzpt = p[3] / (pp * pt) * sqh * hel;
            vc[3] = std::complex<double> (-p[1] * pzpt, -nsv * p[2] / pt * sqh);
            vc[4] = std::complex<double> (-p[2] * pzpt, nsv * p[1] / pt * sqh);
        } else {
            vc[3] = std::complex<double> (-hel * sqh, 0.0);
            vc[4] = std::complex<double> (0.0, nsv * Sgn(sqh, p[3]));
        }
    }
    return;
}

void eeWW::sxxxxx(double p[4], int nss, std::complex<double> sc[3])
{
    sc[2] = std::complex<double> (1.00, 0.00);
    sc[0] = std::complex<double> (p[0] * nss, p[3] * nss);
    sc[1] = std::complex<double> (p[1] * nss, p[2] * nss);
    return;
}

void eeWW::oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double> fo[6])
{
    std::complex<double> chi[2];
    double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2];
    int nh, ip, im;
    fo[0] = std::complex<double> (p[0] * nsf, p[3] * nsf);
    fo[1] = std::complex<double> (p[1] * nsf, p[2] * nsf);
    nh = nhel * nsf;
    if (fmass != 0.000) {
        pp = std::min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3])));
        if (pp == 0.000) {
            sqm[0] = sqrt(std::abs(fmass));
            sqm[1] = Sgn(sqm[0], fmass);
            ip = -((1 - nh) / 2) * nhel;
            im = (1 + nh) / 2 * nhel;
            fo[2] = im * sqm[std::abs(ip)];
            fo[3] = ip * nsf * sqm[std::abs(ip)];
            fo[4] = im * nsf * sqm[std::abs(im)];
            fo[5] = ip * sqm[std::abs(im)];
        } else {
            pp = std::min(p[0], sqrt((p[1] * p[1]) + (p[2] * p[2]) + (p[3] * p[3])));
            sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5;
            sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5;
            omega[0] = sqrt(p[0] + pp);
            omega[1] = fmass / omega[0];
            ip = (1 + nh) / 2;
            im = (1 - nh) / 2;
            sfomeg[0] = sf[0] * omega[ip];
            sfomeg[1] = sf[1] * omega[im];
            pp3 = std::max(pp + p[3], 0.00);
            chi[0] = std::complex<double> (sqrt(pp3 * 0.5 / pp), 0.00);
            if (pp3 == 0.00) {
                chi[1] = std::complex<double> (-nh, 0.00);
            } else {
                chi[1] = std::complex<double> (nh * p[1], -p[2]) / sqrt(2.0 * pp * pp3);
            }
            fo[2] = sfomeg[1] * chi[im];
            fo[3] = sfomeg[1] * chi[ip];
            fo[4] = sfomeg[0] * chi[im];
            fo[5] = sfomeg[0] * chi[ip];
        }
    } else {
        if ((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00)) {
            sqp0p3 = 0.00;
        } else {
            sqp0p3 = sqrt(std::max(p[0] + p[3], 0.00)) * nsf;
        }
        chi[0] = std::complex<double> (sqp0p3, 0.00);
        if (sqp0p3 == 0.000) {
            chi[1] = std::complex<double> (-nhel, 0.00) * sqrt(2.0 * p[0]);
        } else {
            chi[1] = std::complex<double> (nh * p[1], -p[2]) / sqp0p3;
        }
        if (nh == 1) {
            fo[2] = chi[0];
            fo[3] = chi[1];
            fo[4] = std::complex<double> (0.00, 0.00);
            fo[5] = std::complex<double> (0.00, 0.00);
        } else {
            fo[2] = std::complex<double> (0.00, 0.00);
            fo[3] = std::complex<double> (0.00, 0.00);
            fo[4] = chi[1];
            fo[5] = chi[0];
        }
    }
    return;
}

void eeWW::VVV1_0(std::complex<double> V1[], std::complex<double> V2[],
        std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
        & vertex)
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    std::complex<double> TMP2;
    double P1[4];
    std::complex<double> TMP10;
    double P2[4];
    std::complex<double> TMP7;
    double P3[4];
    std::complex<double> TMP6;
    std::complex<double> TMP5;
    std::complex<double> TMP4;
    std::complex<double> TMP9;
    std::complex<double> TMP3;
    std::complex<double> TMP8;
    P1[0] = V1[0].real();
    P1[1] = V1[1].real();
    P1[2] = V1[1].imag();
    P1[3] = V1[0].imag();
    P2[0] = V2[0].real();
    P2[1] = V2[1].real();
    P2[2] = V2[1].imag();
    P2[3] = V2[0].imag();
    P3[0] = V3[0].real();
    P3[1] = V3[1].real();
    P3[2] = V3[1].imag();
    P3[3] = V3[0].imag();
    TMP9 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]);
    TMP8 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]);
    TMP5 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]);
    TMP4 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]);
    TMP7 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]);
    TMP6 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]);
    TMP10 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]);
    TMP3 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]);
    TMP2 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]);
    vertex = COUP * (TMP2 * (-cI * (TMP3) + cI * (TMP4)) + (TMP6 * (-cI * (TMP7)
            + cI * (TMP5)) + TMP8 * (-cI * (TMP9) + cI * (TMP10))));
}

void eeWW::FFV4_3(std::complex<double> F1[], std::complex<double> F2[],
        std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    std::complex<double> denom;
    std::complex<double> TMP11;
    double P3[4];
    double OM3;
    std::complex<double> TMP1;
    OM3 = 0.;
    if (M3 != 0.)
        OM3 = 1. / (M3 * M3);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
            F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    TMP11 = (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
            F1[5] * (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])));
    denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
            P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-2. * cI) * (OM3 * -1. / 2. * P3[0] * (TMP1 + 2. * (TMP11)) +
            (+1. / 2. * (F1[2] * F2[4] + F1[3] * F2[5]) + F1[4] * F2[2] + F1[5] *
            F2[3]));
    V3[3] = denom * (-2. * cI) * (OM3 * -1. / 2. * P3[1] * (TMP1 + 2. * (TMP11)) +
            (-1. / 2. * (F1[2] * F2[5] + F1[3] * F2[4]) + F1[4] * F2[3] + F1[5] *
            F2[2]));
    V3[4] = denom * 2. * cI * (OM3 * 1. / 2. * P3[2] * (TMP1 + 2. * (TMP11)) +
            (+1. / 2. * cI * (F1[2] * F2[5]) - 1. / 2. * cI * (F1[3] * F2[4]) - cI *
            (F1[4] * F2[3]) + cI * (F1[5] * F2[2])));
    V3[5] = denom * 2. * cI * (OM3 * 1. / 2. * P3[3] * (TMP1 + 2. * (TMP11)) +
            (+1. / 2. * (F1[2] * F2[4]) - 1. / 2. * (F1[3] * F2[5]) - F1[4] * F2[2] +
            F1[5] * F2[3]));
}

void eeWW::FFV2_3(std::complex<double> F1[], std::complex<double> F2[],
        std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    std::complex<double> denom;
    std::complex<double> TMP1;
    double P3[4];
    double OM3;
    OM3 = 0.;
    if (M3 != 0.)
        OM3 = 1. / (M3 * M3);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP1 = (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
            F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])));
    denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
            P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] - P3[0] * OM3 * TMP1);
    V3[3] = denom * (-cI) * (-F1[2] * F2[5] - F1[3] * F2[4] - P3[1] * OM3 *
            TMP1);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5]) + cI * (F1[3] * F2[4]) - P3[2]
            * OM3 * TMP1);
    V3[5] = denom * (-cI) * (F1[3] * F2[5] - F1[2] * F2[4] - P3[3] * OM3 * TMP1);
}

void eeWW::FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[],
        std::complex<double> COUP1, std::complex<double> COUP2, double M3, double
        W3, std::complex<double> V3[])
{
    int i;
    std::complex<double> Vtmp[6];
    FFV2_3(F1, F2, COUP1, M3, W3, V3);
    FFV4_3(F1, F2, COUP2, M3, W3, Vtmp);
    i = 2;
    while (i < 6) {
        V3[i] = V3[i] + Vtmp[i];
        i++;
    }
}

void eeWW::FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
        std::complex<double> COUP, double M3, double W3, std::complex<double> V3[])
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    double P3[4];
    std::complex<double> denom;
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    denom = COUP / ((P3[0] * P3[0]) - (P3[1] * P3[1]) - (P3[2] * P3[2]) - (P3[3] *
            P3[3]) - M3 * (M3 - cI * W3));
    V3[2] = denom * (-cI) * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] +
            F1[5] * F2[3]);
    V3[3] = denom * (-cI) * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] -
            F1[3] * F2[4]);
    V3[4] = denom * (-cI) * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3]
            * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * (-cI) * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] -
            F1[5] * F2[3]);
}

void eeWW::FFV2_1(std::complex<double> F2[], std::complex<double> V3[],
        std::complex<double> COUP, double M1, double W1, std::complex<double> F1[])
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    double P1[4];
    std::complex<double> denom;
    F1[0] = +F2[0] + V3[0];
    F1[1] = +F2[1] + V3[1];
    P1[0] = -F1[0].real();
    P1[1] = -F1[1].real();
    P1[2] = -F1[1].imag();
    P1[3] = -F1[0].imag();
    denom = COUP / ((P1[0] * P1[0]) - (P1[1] * P1[1]) - (P1[2] * P1[2]) - (P1[3] *
            P1[3]) - M1 * (M1 - cI * W1));
    F1[2] = denom * cI * M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
            (V3[4])));
    F1[3] = denom * - cI * M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] * (V3[5]
            - V3[2]));
    F1[4] = denom * (-cI) * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
            (V3[4]) - V3[3]) + (P1[2] * (-1.) * (V3[4] + cI * (V3[3])) - P1[3] *
            (V3[2] + V3[5])))) + F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
            (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
            (V3[3] + cI * (V3[4]))))));
    F1[5] = denom * (-cI) * (F2[4] * (P1[0] * (V3[3] - cI * (V3[4])) + (P1[1] *
            (-1.) * (V3[2] + V3[5]) + (P1[2] * (+cI * (V3[2] + V3[5])) + P1[3] *
            (V3[3] - cI * (V3[4]))))) + F2[5] * (P1[0] * (V3[2] - V3[5]) + (P1[1] *
            (-1.) * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) + P1[3]
            * (V3[2] - V3[5])))));
}

void eeWW::FFV2_0(std::complex<double> F1[], std::complex<double> F2[],
        std::complex<double> V3[], std::complex<double> COUP, std::complex<double>
        & vertex)
{
    static std::complex<double> cI = std::complex<double> (0., 1.);
    std::complex<double> TMP0;
    TMP0 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
            F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
    vertex = COUP * - cI * TMP0;
}

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > w+ w- WEIGHTED<=4 @1
// *   Decay: w+ > mu+ vm WEIGHTED<=2
// *   Decay: w- > u~ d WEIGHTED<=2

//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void eeWW::sigmaKin()
{
    // Set the parameters which change event by event

    // Reset color flows
    for (int i = 0; i < 1; i++)
        jamp2[0][i] = 0.;
    
    // Local variables and constants
    const int ncomb = 64;
    static bool goodhel[ncomb] = {ncomb * false};
    static int ntry = 0, sum_hel = 0, ngood = 0;
    static int igood[ncomb];
    static int jhel;
    //  std::complex<double> * * wfs; 
    double t[nprocesses];
    
    // Helicities for the process
    static const int helicities[ncomb][nexternal] = {
        {-1, -1, -1, -1, -1, -1},
        {-1, -1, -1, -1, -1, 1},
        {-1, -1, -1, -1, 1, -1},
        {-1, -1, -1, -1, 1, 1},
        {-1, -1, -1, 1, -1, -1},
        {-1, -1, -1, 1, -1, 1},
        {-1, -1, -1, 1, 1, -1},
        {-1, -1, -1, 1, 1, 1},
        {-1, -1, 1, -1, -1, -1},
        {-1, -1, 1, -1, -1, 1},
        {-1, -1, 1, -1, 1, -1},
        {-1, -1, 1, -1, 1, 1},
        {-1, -1, 1, 1, -1, -1},
        {-1, -1, 1, 1, -1, 1},
        {-1, -1, 1, 1, 1, -1},
        {-1, -1, 1, 1, 1, 1},
        {-1,1, -1, -1, -1, -1},
        {-1, 1, -1, -1, -1, 1},
        {-1, 1, -1, -1, 1, -1},
        {-1, 1, -1, -1, 1, 1},
        {-1, 1, -1, 1, -1, -1},
        {-1, 1, -1, 1, -1, 1},
        {-1, 1, -1, 1, 1, -1},
        {-1, 1, -1, 1, 1, 1},
        {-1, 1, 1, -1, -1, -1},
        {-1, 1, 1, -1, -1, 1},
        {-1, 1, 1, -1, 1, -1},
        {-1, 1, 1, -1, 1, 1},
        {-1, 1, 1, 1, -1, -1},
        {-1, 1, 1, 1, -1, 1},
        {-1, 1, 1, 1, 1, -1},
        {-1, 1, 1, 1, 1, 1},
        {1, -1, -1, -1, -1, -1},
        {1, -1, -1, -1, -1, 1},
        {1, -1, -1, -1, 1, -1},
        {1, -1, -1, -1, 1, 1},
        {1, -1, -1, 1, -1, -1},
        {1, -1, -1, 1, -1, 1},
        {1, -1, -1, 1, 1, -1},
        {1, -1, -1, 1, 1, 1},
        {1, -1, 1, -1, -1, -1},
        {1, -1, 1, -1, -1, 1},
        {1, -1, 1, -1, 1, -1},
        {1, -1, 1, -1, 1, 1},
        {1, -1, 1, 1, -1, -1},
        {1, -1, 1, 1, -1, 1},
        {1, -1, 1, 1, 1, -1},
        {1, -1, 1, 1, 1, 1},
        {1, 1, -1, -1, -1, -1},
        {1, 1, -1, -1, -1, 1},
        {1, 1, -1, -1, 1, -1},
        {1, 1, -1, -1, 1, 1},
        {1, 1, -1, 1, -1, -1},
        {1, 1, -1, 1, -1, 1},
        {1, 1, -1, 1, 1, -1},
        {1, 1, -1, 1, 1, 1},
        {1, 1, 1, -1, -1, -1},
        {1, 1, 1, -1, -1, 1},
        {1, 1, 1, -1, 1, -1},
        {1, 1, 1, -1, 1, 1},
        {1, 1, 1, 1, -1, -1},
        {1,1, 1, 1, -1, 1},
        {1, 1, 1, 1, 1, -1},
        {1, 1, 1, 1, 1, 1}};
    // Denominators: spins, colors and identical particles
    const int denominators[nprocesses] = {4};
    
    ntry = ntry + 1;
    
    // Reset the matrix elements
    for (int i = 0; i < nprocesses; i++) {
        matrix_element[i] = 0.;
    }
    
    // Define permutation   
    for (int i = 0; i < nexternal; i++) {
        perm[i] = i;
    }
    
    if (sum_hel == 0 || ntry < 10) {
        // Calculate the matrix element for all helicities
        for (int ihel = 0; ihel < ncomb; ihel++) {
            if (goodhel[ihel] || ntry < 2) {
                calculate_wavefunctions(perm, helicities[ihel]);
                t[0] = matrix_1_epem_wpwm_wp_mupvm_wm_uxd();
                double tsum = 0;
                for (int iproc = 0; iproc < nprocesses; iproc++) {
                    matrix_element[iproc] += t[iproc];
                    tsum += t[iproc];
                }
                // Store which helicities give non-zero result
                if (tsum != 0. && !goodhel[ihel]) {
                    goodhel[ihel] = true;
                    ngood++;
                    igood[ngood] = ihel;
                }
            }
        }
        jhel = 0;
        sum_hel = std::min(sum_hel, ngood);
    } else {
        // Only use the "good" helicities
        for (int j = 0; j < sum_hel; j++) {
            jhel++;
            if (jhel >= ngood) jhel = 0;
            double hwgt = double(ngood) / double(sum_hel);
            int ihel = igood[jhel];
            calculate_wavefunctions(perm, helicities[ihel]);
            t[0] = matrix_1_epem_wpwm_wp_mupvm_wm_uxd();

            for (int iproc = 0; iproc < nprocesses; iproc++) {
                matrix_element[iproc] += t[iproc] * hwgt;
            }
        }
    }
    
    for (int i = 0; i < nprocesses; i++)
        matrix_element[i] /= denominators[i];
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double eeWW::sigmaHat()
{
    // Select between the different processes
    if (id1 == -11 && id2 == 11) {
        // Add matrix elements for processes with beams (-11, 11)
        return matrix_element[0];
    } else {
        // Return 0 if not correct initial state assignment
        return 0.;
    }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void eeWW::calculate_wavefunctions(const int perm[], const int hel[])
{
    // Calculate wavefunctions for all processes

    // Calculate all wavefunctions
    oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
    ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]);
    ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
    FFV2_3(w[2], w[3], GC_100, MW, Gamma_W, w[4]);
    ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[5]);
    oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[6]);
    FFV2_3(w[5], w[6], GC_100, MW, Gamma_W, w[7]);
    FFV1P0_3(w[1], w[0], GC_3, ZERO, ZERO, w[8]);
    FFV2_4_3(w[1], w[0], -GC_51, GC_59, MZ, Gamma_Z,
            w[9]);
    FFV2_1(w[0], w[4], GC_100, ZERO, ZERO, w[10]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0
    VVV1_0(w[8], w[7], w[4], -GC_3, amp[0]);
    VVV1_0(w[7], w[4], w[9], GC_53, amp[1]);
    FFV2_0(w[1], w[10], w[7], GC_100, amp[2]);

}

double eeWW::matrix_1_epem_wpwm_wp_mupvm_wm_uxd()
{
    int i, j;
    // Local variables
//    const int ngraphs = 3;
    const int ncolor = 1;
    std::complex<double> ztemp;
    std::complex<double> jamp[ncolor];
    // The color matrix;
    static const double denom[ncolor] = {1};
    static const double cf[ncolor][ncolor] = {
        {3}
    };

    // Calculate color flows
    jamp[0] = -amp[0] - amp[1] - amp[2];

    // Sum and square the color flows to get the matrix element
    double matrix = 0;
    for (i = 0; i < ncolor; i++) {
        ztemp = 0.;
        for (j = 0; j < ncolor; j++)
            ztemp = ztemp + cf[i][j] * jamp[j];
        matrix = matrix + real(ztemp * conj(jamp[i])) / denom[i];
    }

    // Store the leading color flows for choice of color
    for (i = 0; i < ncolor; i++)
        jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

    return matrix;
}