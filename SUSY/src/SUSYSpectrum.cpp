/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdexcept>
#include "SUSYSpectrum.h"


SUSYSpectrum::SUSYSpectrum(SUSY& SUSY_in)
: mySUSY(SUSY_in), Mchargino(2,2,0.), Mneutralino(4,4,0.),
        Msup2(6,0.), Msdown2(6,0.), Msneutrino2(6,0.), Mselectron2(6,0.),
        mch(2,0.), mneu(4,0.), m_su2(6,0.), m_sd2(6,0.), m_sn2(6,0.), m_se2(6,0.),
        U(2,2,0.), V(2,2,0.), N(4,4,0.),
        Ru(6,6,0.), Rd(6,6,0.), Rn(6,6,0.), Rl(6,6,0.)
{
}

void SUSYSpectrum::CalcHiggs()
{
    double Mw = mySUSY.Mw_tree();
    double Mz = mySUSY.getMz();
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;

    /* charged Higgs mass */
    mh[3] = mySUSY.mHptree;

    /* pseudo-scalar Higgs mass */
    mh[2] = sqrt(mh[3] * mh[3] - Mw * Mw);

    double temp = mh[2] * mh[2] + Mz * Mz;
    double temp1 = 2.0 * mh[2] * Mz * cos2b;
    double temp2 = sqrt(temp * temp - temp1 * temp1);

    /* light Higgs mass */
    mh[0] = sqrt((temp - temp2)/2.0);
    
    /* heavy Higgs mass */
    mh[1] = sqrt((temp + temp2)/2.0);
}

void SUSYSpectrum::CalcChargino()
{
    double Mw = mySUSY.Mw_tree();

    Mchargino.assign(0, 0, mySUSY.getM2());
    Mchargino.assign(0, 1, sqrt(2) * Mw * mySUSY.getSinb());
    //
    Mchargino.assign(1, 0, sqrt(2) * Mw * mySUSY.getCosb());
    Mchargino.assign(1, 1, mySUSY.getMuH());

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> Utmp(2, 2, 0.), Vtmp(2, 2, 0.);
    Mchargino.singularvalue(Utmp, Vtmp, mch);

    /**
     * SLHA: diag(mChi) = U^* M V^+.
     */
    U = Utmp.transpose();
    V = Vtmp.hconjugate();
}

void SUSYSpectrum::CalcNeutralino()
{
    double Mw = mySUSY.Mw_tree();
    double Mz = mySUSY.getMz();
    double cW2 = Mw*Mw/Mz/Mz;
    double cW = sqrt(cW2);
    double sW = sqrt(1.0 - cW2);
    double sb = mySUSY.getSinb();
    double cb = mySUSY.getCosb();

    Mneutralino.assign(0, 0, mySUSY.getM1());
    Mneutralino.assign(0, 2, - Mz * sW * cb);
    Mneutralino.assign(0, 3, Mz * sW * sb);
    //
    Mneutralino.assign(1, 1, mySUSY.getM2());
    Mneutralino.assign(1, 2, Mz * cW * cb);
    Mneutralino.assign(1, 3, - Mz * cW * sb);
    //
    Mneutralino.assign(2, 0, Mneutralino(0,2));
    Mneutralino.assign(2, 1, Mneutralino(1,2));
    Mneutralino.assign(2, 3, - mySUSY.getMuH());
    //
    Mneutralino.assign(3, 0, Mneutralino(0,3));
    Mneutralino.assign(3, 1, Mneutralino(1,3));
    Mneutralino.assign(3, 2, Mneutralino(2,3));

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> Ntemp1(4, 4, 0.), Ntemp2(4, 4, 0.);
    Mneutralino.singularvalue(Ntemp1, Ntemp2, mneu);

    //matrix<complex> Nleft = Ntemp1.transpose();
    matrix<complex> Nright = Ntemp2.hconjugate();
    //std::cout << "Nleft = " << Nleft << std::endl;
    //std::cout << "Nright = " << Nright << std::endl;

    /* adopt N=Nright as N^* M N^+ = diag, not as N M N^+.
     * As a result, a phase rotation is required for N. */
    matrix<complex> Mdiag_tmp(4, 4, 0.);
    Mdiag_tmp = Nright.hconjugate().transpose() * Mneutralino * Nright.hconjugate();
    //std::cout << "Mdiag_tmp = " << Mdiag_tmp << std::endl;
    vector<complex> v1(4, 0.);
    for(int i = 0; i < 4; i++)
        v1.assign(i, complex(1., Mdiag_tmp(i,i).arg()/2.0, true));
    Nright = matrix<complex>(v1) * Nright;

    N = Nright;
}

void SUSYSpectrum::CalcSup()
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    matrix<complex> CKM(mySUSY.getVCKM());
    matrix<complex> Id3 = matrix<complex>::Id(3);
    
    /* DRbar up-type quark masses at scale Q */
    matrix<double> Mu(3, 3, 0.);
    Mu(0, 0) = mySUSY.Mq_Q(mySUSY.UP);
    Mu(1, 1) = mySUSY.Mq_Q(mySUSY.CHARM);
    Mu(2, 2) = mySUSY.Mq_Q(mySUSY.TOP);

    matrix<complex> uLL( CKM * mySUSY.msQhat2 * CKM.hconjugate() + Mu * Mu
                         + cos2b * Mz2 * (1.0/2.0 - 2.0/3.0 * sW2) * Id3 );
    matrix<complex> uLR( mySUSY.v2()/sqrt(2.0) * mySUSY.getTUhat().hconjugate()
                         - mySUSY.getMuH() * Mu / mySUSY.getTanb() );
    matrix<complex> uRR( mySUSY.msUhat2 + Mu * Mu + cos2b * Mz2 * 2.0/3.0 * sW2 * Id3 );
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) {
            Msup2.assign(i,   j,   uLL(i,j));
            Msup2.assign(i,   j+3, uLR(i,j));
            Msup2.assign(i+3, j,   uLR(j,i).conjugate());
            Msup2.assign(i+3, j+3, uRR(i,j));
        }

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> RuTmp(6, 6, 0.);
    Msup2.eigensystem(RuTmp, m_su2);

    Ru = RuTmp.hconjugate();
}

void SUSYSpectrum::CalcSdown()
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    matrix<complex> Id3 = matrix<complex>::Id(3);

    /* DRbar down-type quark masses at scale Q */
    matrix<double> Md(3, 3, 0.);
    Md(0, 0) = mySUSY.Mq_Q(mySUSY.DOWN);
    Md(1, 1) = mySUSY.Mq_Q(mySUSY.STRANGE);
    Md(2, 2) = mySUSY.Mq_Q(mySUSY.BOTTOM);

    matrix<complex> dLL( mySUSY.msQhat2 + Md * Md
                         + cos2b * Mz2 * (- 1.0/2.0 + 1.0/3.0 * sW2) * Id3 );
    matrix<complex> dLR( mySUSY.v1()/sqrt(2.0) * mySUSY.getTDhat().hconjugate()
                         - mySUSY.getMuH() * Md * mySUSY.getTanb() );
    matrix<complex> dRR( mySUSY.msDhat2 + Md * Md - cos2b * Mz2 /3.0 * sW2 * Id3 );
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) {
            Msdown2.assign(i,   j,   dLL(i,j));
            Msdown2.assign(i,   j+3, dLR(i,j));
            Msdown2.assign(i+3, j,   dLR(j,i).conjugate());
            Msdown2.assign(i+3, j+3, dRR(i,j));
        }

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> RdTmp(6, 6, 0.);
    Msdown2.eigensystem(RdTmp, m_sd2);

    Rd = RdTmp.hconjugate();
}

void SUSYSpectrum::CalcSneutrino()
{
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    matrix<complex> Id3 = matrix<complex>::Id(3);

    matrix<complex> nLL( mySUSY.msLhat2 + cos2b * Mz2 /2.0 * Id3 );
    matrix<complex> nLR( mySUSY.v2()/sqrt(2.0) * mySUSY.getTNhat().hconjugate() );
    matrix<complex> nRR( mySUSY.msNhat2 );
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) {
            Msneutrino2.assign(i,   j,   nLL(i,j));
            Msneutrino2.assign(i,   j+3, nLR(i,j));
            Msneutrino2.assign(i+3, j,   nLR(j,i).conjugate());
            Msneutrino2.assign(i+3, j+3, nRR(i,j));
        }

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> RnTmp(6, 6, 0.);
    Msneutrino2.eigensystem(RnTmp, m_sn2);

    Rn = RnTmp.hconjugate();
}

void SUSYSpectrum::CalcSelectron()
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    matrix<complex> Id3 = matrix<complex>::Id(3);

    matrix<double> Me(3, 3, 0.);
    Me(0, 0) = mySUSY.Ml_Q(mySUSY.ELECTRON);
    Me(1, 1) = mySUSY.Ml_Q(mySUSY.MU);
    Me(2, 2) = mySUSY.Ml_Q(mySUSY.TAU);

    matrix<complex> eLL( mySUSY.msLhat2 + Me * Me
                         + cos2b * Mz2 * (- 1.0/2.0 + sW2) * Id3 );
    matrix<complex> eLR( mySUSY.v1()/sqrt(2.0) * mySUSY.getTEhat().hconjugate()
                         - mySUSY.getMuH() * Me * mySUSY.getTanb() );
    matrix<complex> eRR( mySUSY.msEhat2 + Me * Me - cos2b * Mz2 * sW2 * Id3 );
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) {
            Mselectron2.assign(i,   j,   eLL(i,j));
            Mselectron2.assign(i,   j+3, eLR(i,j));
            Mselectron2.assign(i+3, j,   eLR(j,i).conjugate());
            Mselectron2.assign(i+3, j+3, eRR(i,j));
        }

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    matrix<complex> RlTmp(6, 6, 0.);
    Mselectron2.eigensystem(RlTmp, m_se2);

    Rl = RlTmp.hconjugate();
}



