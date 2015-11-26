/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdexcept>
#include <limits>
#include "SUSYSpectrum.h"
#include "SUSY.h"


SUSYSpectrum::SUSYSpectrum(const SUSY & SUSY_in)
: mySUSY(SUSY_in), Mchargino(2,2,0.), Mneutralino(4,4,0.),
        Msup2(6,0.), Msdown2(6,0.), Msneutrino2(6,0.), Mselectron2(6,0.),
        mch(2,0.), mneu(4,0.), m_su2(6,0.), m_sd2(6,0.), m_sn2(6,0.), m_se2(6,0.),
        U(2,2,0.), V(2,2,0.), N(4,4,0.),
        Ru(6,6,0.), Rd(6,6,0.), Rn(6,6,0.), Rl(6,6,0.)
{
}

bool SUSYSpectrum::CalcHiggs(double mh[4], gslpp::complex& saeff_i)
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
    
    /* CP-even Higgs mixing angle*/
    saeff_i = sin(atan((mh[2]*mh[2]+Mz*Mz)*sqrt(1-cos2b*cos2b)/((mh[2]*mh[2]-Mz*Mz)*cos2b))/2.0);

    return true;
}

bool SUSYSpectrum::CalcChargino(gslpp::matrix<gslpp::complex>& U_i, gslpp::matrix<gslpp::complex>& V_i, gslpp::vector<double>& mch_i)
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
    gslpp::matrix<gslpp::complex> Utmp(2, 2, 0.), Vtmp(2, 2, 0.);
    Mchargino.singularvalue(Utmp, Vtmp, mch_i);

    /**
     * SLHA: diag(mChi) = U^* M V^+.
     */
    U_i = Utmp.transpose();
    V_i = Vtmp.hconjugate();

    return true;
}

bool SUSYSpectrum::CalcNeutralino(gslpp::matrix<gslpp::complex>& N_i, gslpp::vector<double>& mneu_i)
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
    gslpp::matrix<gslpp::complex> Ntemp1(4, 4, 0.), Ntemp2(4, 4, 0.);
    Mneutralino.singularvalue(Ntemp1, Ntemp2, mneu_i);

    //gslpp::matrix<gslpp::complex> Nleft = Ntemp1.transpose();
    gslpp::matrix<gslpp::complex> Nright = Ntemp2.hconjugate();
    //std::cout << "Nleft = " << Nleft << std::endl;
    //std::cout << "Nright = " << Nright << std::endl;

    /* adopt N=Nright as N^* M N^+ = diag, not as N M N^+.
     * As a result, a phase rotation is required for N. */
    gslpp::matrix<gslpp::complex> Mdiag_tmp(4, 4, 0.);
    Mdiag_tmp = Nright.hconjugate().transpose() * Mneutralino * Nright.hconjugate();
    //std::cout << "Mdiag_tmp = " << Mdiag_tmp << std::endl;
    gslpp::vector<gslpp::complex> v1(4, 0.);
    for(int i = 0; i < 4; i++)
        v1.assign(i, gslpp::complex(1., Mdiag_tmp(i,i).arg()/2.0, true));
    Nright = gslpp::matrix<gslpp::complex>(v1) * Nright;

    N_i = Nright;

    return true;
}

bool SUSYSpectrum::CalcSup(gslpp::matrix<gslpp::complex>& Ru_i, gslpp::vector<double>& m_su2_i)
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    gslpp::matrix<gslpp::complex> CKM(mySUSY.getVCKM());
    gslpp::matrix<gslpp::complex> Id3 = gslpp::matrix<gslpp::complex>::Id(3);
    
    /* DRbar up-type quark masses at scale Q */
    gslpp::matrix<double> Mu(3, 3, 0.);
    Mu(0, 0) = mySUSY.Mq_Q(mySUSY.UP);
    Mu(1, 1) = mySUSY.Mq_Q(mySUSY.CHARM);
    Mu(2, 2) = mySUSY.Mq_Q(mySUSY.TOP);

    gslpp::matrix<gslpp::complex> uLL( CKM * mySUSY.msQhat2 * CKM.hconjugate() + Mu * Mu
                         + cos2b * Mz2 * (1.0/2.0 - 2.0/3.0 * sW2) * Id3 );
    gslpp::matrix<gslpp::complex> uLR( mySUSY.v2()/sqrt(2.0) * mySUSY.getTUhat().hconjugate()
                         - mySUSY.getMuH() * Mu / mySUSY.getTanb() );
    gslpp::matrix<gslpp::complex> uRR( mySUSY.msUhat2 + Mu * Mu + cos2b * Mz2 * 2.0/3.0 * sW2 * Id3 );
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
    gslpp::matrix<gslpp::complex> RuTmp(6, 6, 0.);
    Msup2.eigensystem(RuTmp, m_su2_i);

    Ru_i = RuTmp.hconjugate();

    return true;

}

bool SUSYSpectrum::CalcSdown(gslpp::matrix<gslpp::complex>& Rd_i, gslpp::vector<double>& m_sd2_i)
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    gslpp::matrix<gslpp::complex> Id3 = gslpp::matrix<gslpp::complex>::Id(3);

    /* DRbar down-type quark masses at scale Q */
    gslpp::matrix<double> Md(3, 3, 0.);
    Md(0, 0) = mySUSY.Mq_Q(mySUSY.DOWN);
    Md(1, 1) = mySUSY.Mq_Q(mySUSY.STRANGE);
    Md(2, 2) = mySUSY.Mq_Q(mySUSY.BOTTOM);

    gslpp::matrix<gslpp::complex> dLL( mySUSY.msQhat2 + Md * Md
                         + cos2b * Mz2 * (- 1.0/2.0 + 1.0/3.0 * sW2) * Id3 );
    gslpp::matrix<gslpp::complex> dLR( mySUSY.v1()/sqrt(2.0) * mySUSY.getTDhat().hconjugate()
                         - mySUSY.getMuH() * Md * mySUSY.getTanb() );
    gslpp::matrix<gslpp::complex> dRR( mySUSY.msDhat2 + Md * Md - cos2b * Mz2 /3.0 * sW2 * Id3 );
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
    gslpp::matrix<gslpp::complex> RdTmp(6, 6, 0.);
    Msdown2.eigensystem(RdTmp, m_sd2_i);

    Rd_i = RdTmp.hconjugate();

    return true;

}

bool SUSYSpectrum::CalcSneutrino(gslpp::matrix<gslpp::complex>& Rn_i, gslpp::vector<double>& m_sn2_i)
{
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    gslpp::matrix<gslpp::complex> Id3 = gslpp::matrix<gslpp::complex>::Id(3);

    //  this section is useful to re-define the sneutrino matrix...
    
//                double delta12=0.;
//                double delta13=0.1;
//                double delta23=0.1;
//                gslpp::complex sLmass=mySUSY.msLhat2(0,0);
//                gslpp::matrix<gslpp::complex> msLhat2modified(3,3,0.);
//                    msLhat2modified.assign(0, 0, sLmass);
//                    msLhat2modified.assign(1, 1, sLmass);
//                    msLhat2modified.assign(2, 2, sLmass);
//                    msLhat2modified.assign(0, 1, 0.);
//                    msLhat2modified.assign(1, 0, 0.);
//                    msLhat2modified.assign(0, 2, delta13*sLmass);
//                    msLhat2modified.assign(2, 0, delta13*sLmass);
//                    msLhat2modified.assign(1, 2, delta23*sLmass);
//                    msLhat2modified.assign(2, 1, delta23*sLmass);
//                    gslpp::matrix<gslpp::complex> nLL( msLhat2modified + cos2b * Mz2 /2.0 * Id3 );
 
    // ... until here

    gslpp::matrix<gslpp::complex> nLL( mySUSY.msLhat2 + cos2b * Mz2 /2.0 * Id3 );
    gslpp::matrix<gslpp::complex> nLR( mySUSY.v2()/sqrt(2.0) * mySUSY.getTNhat().hconjugate() );
    gslpp::matrix<gslpp::complex> nRR( mySUSY.msNhat2 );
    for(int i = 0; i < 3; i++) {
        for(int j = 0; j < 3; j++) {
            Msneutrino2.assign(i,   j,   nLL(i,j));
            Msneutrino2.assign(i,   j+3, nLR(i,j));
            Msneutrino2.assign(i+3, j,   nLR(j,i).conjugate());
            Msneutrino2.assign(i+3, j+3, nRR(i,j));
        }
    }
    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    gslpp::matrix<gslpp::complex> RnTmp(6, 6, 0.);
    Msneutrino2.eigensystem(RnTmp, m_sn2_i);

    m_sn2_i(0)=std::numeric_limits<double>::max();
    m_sn2_i(1)=std::numeric_limits<double>::max();
    m_sn2_i(2)=std::numeric_limits<double>::max();
    Rn_i = RnTmp.hconjugate();

    return true;

}

bool SUSYSpectrum::CalcSelectron(gslpp::matrix<gslpp::complex>& Rl_i, gslpp::vector<double>& m_se2_i)
{
    double Mw = mySUSY.Mw_tree();
    double Mz2 = mySUSY.getMz()*mySUSY.getMz();
    double sW2 = 1.0 - Mw*Mw/Mz2;
    double cos2b = 2.0 * mySUSY.getCosb() * mySUSY.getCosb() - 1.0;
    gslpp::matrix<gslpp::complex> Id3 = gslpp::matrix<gslpp::complex>::Id(3);

    gslpp::matrix<double> Me(3, 3, 0.);
    Me(0, 0) = mySUSY.Ml_Q(mySUSY.ELECTRON);
    Me(1, 1) = mySUSY.Ml_Q(mySUSY.MU);
    Me(2, 2) = mySUSY.Ml_Q(mySUSY.TAU);

    //  this section is useful to re-define the slepton matrix
    //if you want to use deltas, use this...
                double delta12=0.1;
                double delta13=0.;
                double delta23=0.;
                gslpp::complex sLmass=mySUSY.msLhat2(0,0);
                gslpp::matrix<gslpp::complex> msLhat2modified(3,3,0.);
                    msLhat2modified.assign(0, 0, sLmass);
                    msLhat2modified.assign(1, 1, sLmass);
                    msLhat2modified.assign(2, 2, sLmass);
                    msLhat2modified.assign(0, 1, delta12*sLmass);
                    msLhat2modified.assign(1, 0, delta12*sLmass);
                    msLhat2modified.assign(0, 2, 0.);
                    msLhat2modified.assign(2, 0, 0.);
                    msLhat2modified.assign(1, 2, 0.);
                    msLhat2modified.assign(2, 1, 0.);
                gslpp::matrix<gslpp::complex> msEhat2modified(3,3,0.);
                    msEhat2modified.assign(0, 0, sLmass);
                    msEhat2modified.assign(1, 1, sLmass);
                    msEhat2modified.assign(2, 2, sLmass);
                    msEhat2modified.assign(0, 1, 0.);
                    msEhat2modified.assign(1, 0, 0.);
                    msEhat2modified.assign(0, 2, 0.);
                    msEhat2modified.assign(2, 0, 0.);
                    msEhat2modified.assign(1, 2, 0.);
                    msEhat2modified.assign(2, 1, 0.);
                gslpp::matrix<gslpp::complex> eLL( msLhat2modified + Me * Me
                         + cos2b * Mz2 * (- 1.0/2.0 + sW2) * Id3 );
                gslpp::matrix<gslpp::complex> eLR( mySUSY.v1()/sqrt(2.0) * mySUSY.getTEhat().hconjugate()
                         - mySUSY.getMuH() * Me * mySUSY.getTanb() );
                gslpp::matrix<gslpp::complex> eRR( msEhat2modified + Me * Me - cos2b * Mz2 * sW2 * Id3 );

    // ... until here
    //else use the following...
//    gslpp::matrix<gslpp::complex> eLL( mySUSY.msLhat2 + Me * Me
//                         + cos2b * Mz2 * (- 1.0/2.0 + sW2) * Id3 );
//    gslpp::matrix<gslpp::complex> eLR( mySUSY.v1()/sqrt(2.0) * mySUSY.getTEhat().hconjugate()
//                         - mySUSY.getMuH() * Me * mySUSY.getTanb() );
//    gslpp::matrix<gslpp::complex> eRR( mySUSY.msEhat2 + Me * Me - cos2b * Mz2 * sW2 * Id3 );
    // ... until here

        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                Mselectron2.assign(i,   j,   eLL(i,j));
                Mselectron2.assign(i,   j+3, eLR(i,j));
                Mselectron2.assign(i+3, j,   eLR(j,i).conjugate());
                Mselectron2.assign(i+3, j+3, eRR(i,j));
            }
        }

    /*
     * M.singularvalue(&U, &V, &S) decompose M into
     * M = U diag(S) V^+.
     */
    gslpp::matrix<gslpp::complex> RlTmp(6, 6, 0.);
    Mselectron2.eigensystem(RlTmp, m_se2_i);
    
    Rl_i = RlTmp.hconjugate();

    return true;

}

void SUSYSpectrum::SortSfermionMasses(gslpp::vector<double>& m_sf2, gslpp::matrix<gslpp::complex>& Rf) const
{
    int newIndex[6];
    for (int i = 0; i < 6; i++)
        newIndex[i] = i;

    /* sort sfermion masses in increasing order */
    for (int i = 0; i < 5; i++)
        for (int k = i + 1; k < 6; k++)
            if (m_sf2(i) > m_sf2(k)) {
                std::swap(m_sf2(i), m_sf2(k));
                std::swap(newIndex[i], newIndex[k]);
            }

    /* sort the corresponding rotation matrix, where the first(second) index
     * denotes mass(gauge) eigenstates. */
    gslpp::matrix<gslpp::complex> myRf(6, 6, 0.);
    for (int i = 0; i < 6; i++)
        for (int k = 0; k < 6; k++)
            myRf.assign(k, i, Rf(newIndex[k], i));
    Rf = myRf;
}