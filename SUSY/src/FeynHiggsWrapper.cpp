/*
 * Copyright (C) 2013 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <CSLHA.h>
#include "FeynHiggsWrapper.h"

FeynHiggsWrapper::FeynHiggsWrapper(SUSY& SUSY_in)
: mySUSY(SUSY_in)
{
    int err;
    FHSetFlags(&err,
               4, // Full MSSM
               0, // DRbar
               0, // DRbar
               3, // full 3x3 neutral Higgs mixing in cMSSM (complex MSSM)
               4, // full determination of the propagator matrices’s poles with UHiggs at q^2=0
               2, // include various two-loop contributions
               2, // NLL resummation (for large MCha,MNeu,MGlu,MSUSY)
               1, // running top mass is used in the 1-/2-loop corrections (SM MSbar 2L)
               1, // resum tan beta contributions
               3  // interpolation in phases for missing 2-loop corrections.
                  // The cMSSM a_s a_t corrections are combined with the remaining
                  // corrections, whose complex phases are interpolated in
                  // At, Ab, M_3, MUE
              );
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggsWrapper::FeynHiggsWrapper(): FHSetFlags error " << err;
        throw std::runtime_error(ss.str());
    }
}

bool FeynHiggsWrapper::SetFeynHiggsPars()
{
    int err;

    /* FeynHiggs debug flag */
    //FHSetDebug(2);
    //FHSetDebug(3);

    Mw_FHinput = mySUSY.Mw_tree(); /* Tree-level W-boson mass */
    //Mw_FHinput = mySUSY.StandardModel::Mw(); /* SM prediction, which should not be used, since mHl cannot be set before calling FeynHiggs. */
    //std::cout << "Mw = " << Mw_FHinput << " used in FeynHiggsWrapper::SetFeynHiggsPars()" << std::endl;

    /* Set the FeynHiggs SM input parameters */
    FHSetSMPara(&err,
                1.0/mySUSY.alphaMz(),
                mySUSY.getAlsMz(), mySUSY.getGF(),
                mySUSY.getLeptons(StandardModel::ELECTRON).getMass(),
                mySUSY.getQuarks(QCD::UP).getMass(),
                mySUSY.getQuarks(QCD::DOWN).getMass(),
                mySUSY.getLeptons(StandardModel::MU).getMass(),
                mySUSY.getQuarks(QCD::CHARM).getMass(),
                mySUSY.getQuarks(QCD::STRANGE).getMass(),
                mySUSY.getLeptons(StandardModel::TAU).getMass(),
                mySUSY.getQuarks(QCD::BOTTOM).getMass(),
                Mw_FHinput, mySUSY.getMz(),
                mySUSY.getLambda(), mySUSY.getA(), mySUSY.getRhob(), mySUSY.getEtab());
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::SetFeynHiggsPars(): Error has been detected in SetPara.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    /* Parameters for FeynHiggs */
    double Q_S = mySUSY.Q_SUSY;
    gslpp::complex muHFH = mySUSY.muH;
    gslpp::complex M1FH = mySUSY.m1;
    gslpp::complex M2FH = mySUSY.m2;
    gslpp::matrix<gslpp::complex> MsQ2 = mySUSY.msQhat2;
    gslpp::matrix<gslpp::complex> MsU2 = mySUSY.msUhat2;
    gslpp::matrix<gslpp::complex> MsD2 = mySUSY.msDhat2;
    gslpp::matrix<gslpp::complex> MsL2 = mySUSY.msLhat2;
    gslpp::matrix<gslpp::complex> MsE2 = mySUSY.msEhat2;
    gslpp::matrix<gslpp::complex> KU = mySUSY.TUhat.hconjugate() * mySUSY.v2() / sqrt(2.0);
    gslpp::matrix<gslpp::complex> KD = mySUSY.TDhat.hconjugate() * mySUSY.v1() / sqrt(2.0);
    gslpp::matrix<gslpp::complex> KE = mySUSY.TEhat.hconjugate() * mySUSY.v1() / sqrt(2.0);

    /* MFV trilinear couplings */
    gslpp::vector<gslpp::complex> AU(3,0.), AD(3,0.), AE(3,0.);
    for (int i=0; i<3; i++) {
        int p = (int)mySUSY.UP + 2*i;
        AU.assign(i, KU(i,i) / mySUSY.Mq_Q((QCD::quark)p));
        p = (int)mySUSY.DOWN + 2*i;
        AD.assign(i, KD(i,i) / mySUSY.Mq_Q((QCD::quark)p));
        p = (int)mySUSY.ELECTRON + 2*i;
        AE.assign(i, KE(i,i) / mySUSY.Ml_Q((QCD::lepton)p));
    }

    /* Check if non-minimal flavor-violating (NMFV) entries exist in the
     * sfermion mass matrices. See also IniFV() in SetFV.F of FeynHiggs. */
    NMFVu = true; NMFVd = true; NMFVe = true;// NMFVnu = true; 
    double TMPu = 0.0, TMPd = 0.0, TMPe = 0.0; //TMPnu = 0.0
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
           if (i < j) {
               TMPu += MsQ2(i, j).abs2() + MsU2(i, j).abs2();
               TMPd += MsQ2(i, j).abs2() + MsD2(i, j).abs2();
               //TMPnu += MsL2(i, j).abs2(); /* not used */
               TMPe += MsL2(i, j).abs2() + MsE2(i, j).abs2();
           }
           if (i != j) {
               TMPu += KU(i, j).abs2();
               TMPd += KD(i, j).abs2();
               TMPe += KE(i, j).abs2();
           }
        }
    }
    if (!TMPu) NMFVu = false;
    if (!TMPe) NMFVd = false;
    if (!TMPe) NMFVe = false;

    /* NMFV trilinear couplings. In the case of NMFV, the trilinear couplings
     * AU, AD and AE for FHSetPara() as well as KU, KD and KE for FHSetNMFV()
     * and FHSetLFV() have to be rotated. */
    gslpp::complex muHphase(1.0, - 2.0*muHFH.arg(), true);
    if (NMFVu) AU *= muHphase;
    if (NMFVd) AD *= muHphase;
    if (NMFVe) AE *= muHphase;
    KU *= muHphase;
    KD *= muHphase;
    KE *= muHphase;

    /* NMFV parameters for FeynHiggs */
    gslpp::matrix<gslpp::complex> deltaQLL(3,3,0.);
    gslpp::matrix<gslpp::complex> deltaULR(3,3,0.), deltaURL(3,3,0.), deltaURR(3,3,0.);
    gslpp::matrix<gslpp::complex> deltaDLR(3,3,0.), deltaDRL(3,3,0.), deltaDRR(3,3,0.);
    gslpp::matrix<gslpp::complex> deltaLLL(3,3,0.);
    gslpp::matrix<gslpp::complex> deltaELR(3,3,0.), deltaERL(3,3,0.), deltaERR(3,3,0.);
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++) {
            deltaQLL.assign(i, j, MsQ2(i,j) / sqrt(MsQ2(i,i).real() * MsQ2(j,j).real()));
            deltaULR.assign(i, j, KU(i,j) / sqrt(MsQ2(i,i).real() * MsU2(j,j).real()));
            deltaURL.assign(i, j, KU(j,i).conjugate() / sqrt(MsU2(i,i).real() * MsQ2(j,j).real()));
            deltaURR.assign(i, j, MsU2(i,j) / sqrt(MsU2(i,i).real() * MsU2(j,j).real()));
            deltaDLR.assign(i, j, KD(i,j) / sqrt(MsQ2(i,i).real() * MsD2(j,j).real()));
            deltaDRL.assign(i, j, KD(j,i).conjugate() / sqrt(MsD2(i,i).real() * MsQ2(j,j).real()));
            deltaDRR.assign(i, j, MsD2(i,j) / sqrt(MsD2(i,i).real() * MsD2(j,j).real()));
            deltaLLL.assign(i, j, MsL2(i,j) / sqrt(MsL2(i,i).real() * MsL2(j,j).real()));
            deltaELR.assign(i, j, KE(i,j) / sqrt(MsL2(i,i).real() * MsE2(j,j).real()));
            deltaERL.assign(i, j, KE(j,i).conjugate() / sqrt(MsE2(i,i).real() * MsL2(j,j).real()));
            deltaERR.assign(i, j, MsE2(i,j) / sqrt(MsE2(i,i).real() * MsE2(j,j).real()));
        }

    /* Set the FeynHiggs parameters, where the GUT relation is used for M1=0. */
    FHSetPara(&err,
              mySUSY.mut/mySUSY.quarks[QCD::TOP].getMass(),
              mySUSY.mtpole, mySUSY.tanb,
              mySUSY.mHptree,  // as now used, "mHptree" is a name for MA0. We shall be using mA instead of mHptree
	      -1, // this is now not used, the mHptree 
              //
              sqrt(MsL2(2,2).real()), sqrt(MsE2(2,2).real()),
              sqrt(MsQ2(2,2).real()), sqrt(MsU2(2,2).real()),
              sqrt(MsD2(2,2).real()),
              sqrt(MsL2(1,1).real()), sqrt(MsE2(1,1).real()),
              sqrt(MsQ2(1,1).real()), sqrt(MsU2(1,1).real()),
              sqrt(MsD2(1,1).real()),
              sqrt(MsL2(0,0).real()), sqrt(MsE2(0,0).real()),
              sqrt(MsQ2(0,0).real()), sqrt(MsU2(0,0).real()),
              sqrt(MsD2(0,0).real()),
              //
              ToComplex2(muHFH.real(), muHFH.imag()),
              //
              ToComplex2(AE(2).real(), AE(2).imag()),
              ToComplex2(AU(2).real(), AU(2).imag()),
              ToComplex2(AD(2).real(), AD(2).imag()),
              ToComplex2(AE(1).real(), AE(1).imag()),
              ToComplex2(AU(1).real(), AU(1).imag()),
              ToComplex2(AD(1).real(), AD(1).imag()),
              ToComplex2(AE(0).real(), AE(0).imag()),
              ToComplex2(AU(0).real(), AU(0).imag()),
              ToComplex2(AD(0).real(), AD(0).imag()),
              //
              ToComplex2(M1FH.real(), M1FH.imag()),
              ToComplex2(M2FH.real(), M2FH.imag()),
              ToComplex2(mySUSY.m3, 0.),
              //
              Q_S, Q_S, Q_S);
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::SetFeynHiggsPars(): Error has been detected in SetPara.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    /* Set the non-minimal flavor-violating parameters in the squark sector */
    FHSetNMFV(&err,
              // Q_LL
              ToComplex2(deltaQLL(0,1).real(), deltaQLL(0,1).imag()),
              ToComplex2(deltaQLL(1,2).real(), deltaQLL(1,2).imag()),
              ToComplex2(deltaQLL(0,2).real(), deltaQLL(0,2).imag()),
              // U_LR
              ToComplex2(deltaULR(0,1).real(), deltaULR(0,1).imag()),
              ToComplex2(deltaULR(1,2).real(), deltaULR(1,2).imag()),
              ToComplex2(deltaULR(0,2).real(), deltaULR(0,2).imag()),
              // U_RL
              ToComplex2(deltaURL(0,1).real(), deltaURL(0,1).imag()),
              ToComplex2(deltaURL(1,2).real(), deltaURL(1,2).imag()),
              ToComplex2(deltaURL(0,2).real(), deltaURL(0,2).imag()),
              // U_RR
              ToComplex2(deltaURR(0,1).real(), deltaURR(0,1).imag()),
              ToComplex2(deltaURR(1,2).real(), deltaURR(1,2).imag()),
              ToComplex2(deltaURR(0,2).real(), deltaURR(0,2).imag()),
              // D_LR
              ToComplex2(deltaDLR(0,1).real(), deltaDLR(0,1).imag()),
              ToComplex2(deltaDLR(1,2).real(), deltaDLR(1,2).imag()),
              ToComplex2(deltaDLR(0,2).real(), deltaDLR(0,2).imag()),
              // D_RL
              ToComplex2(deltaDRL(0,1).real(), deltaDRL(0,1).imag()),
              ToComplex2(deltaDRL(1,2).real(), deltaDRL(1,2).imag()),
              ToComplex2(deltaDRL(0,2).real(), deltaDRL(0,2).imag()),
              // D_RR
              ToComplex2(deltaDRR(0,1).real(), deltaDRR(0,1).imag()),
              ToComplex2(deltaDRR(1,2).real(), deltaDRR(1,2).imag()),
              ToComplex2(deltaDRR(0,2).real(), deltaDRR(0,2).imag())
              );
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::SetFeynHiggsPars(): Error was detected in SetFV.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    /* Set the non-minimal flavor-violating parameters in the slepton sector,
     * which are not used to compute the sneutrino mass spectrum. */
    FHSetLFV(&err,
              // L_LL
              ToComplex2(deltaLLL(0,1).real(), deltaLLL(0,1).imag()),
              ToComplex2(deltaLLL(1,2).real(), deltaLLL(1,2).imag()),
              ToComplex2(deltaLLL(0,2).real(), deltaLLL(0,2).imag()),
              // E_LR
              ToComplex2(deltaELR(0,1).real(), deltaELR(0,1).imag()),
              ToComplex2(deltaELR(1,2).real(), deltaELR(1,2).imag()),
              ToComplex2(deltaELR(0,2).real(), deltaELR(0,2).imag()),
              // E_RL
              ToComplex2(deltaERL(0,1).real(), deltaERL(0,1).imag()),
              ToComplex2(deltaERL(1,2).real(), deltaERL(1,2).imag()),
              ToComplex2(deltaERL(0,2).real(), deltaERL(0,2).imag()),
              // E_RR
              ToComplex2(deltaERR(0,1).real(), deltaERR(0,1).imag()),
              ToComplex2(deltaERR(1,2).real(), deltaERR(1,2).imag()),
              ToComplex2(deltaERR(0,2).real(), deltaERR(0,2).imag())
              );
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::SetFeynHiggsPars(): Error was detected in SetFV.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    computeHiggsCouplings = true;
    computeHiggsProd = true;
    computeConstraints = true;
    computeFlavour = true;

    return (true);
}

bool FeynHiggsWrapper::CalcHiggsSpectrum()
{
    int err;
    ComplexType SAeff;
    ComplexType UHiggs[3][3];
    ComplexType ZHiggs[3][3];

    /* Compute the Higgs masses and mixings */
    FHHiggsCorr(&err, mySUSY.mh, &SAeff, UHiggs, ZHiggs);
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::CalcHiggsSpectrum(): Error has been detected in HiggsCorr.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    // the sine of the effective Higgs mixing angle, alpha_eff,
    mySUSY.saeff = gslpp::complex(SAeff.real(), SAeff.imag());

    /* Debug */
    //std::cout << "mh[0] = mh = " << mySUSY.mh[0] << std::endl;
    //std::cout << "mh[1] = mH = " << mySUSY.mh[1] << std::endl;
    //std::cout << "mh[2] = mA = " << mySUSY.mh[2] << std::endl;
    //std::cout << "mh[3] = mH+ = " << mySUSY.mh[3] << std::endl;

    /* Check */
    for (int i = 0; i < 4; i++)
        if(std::isnan(mySUSY.mh[i])) {
            std::cout << "FeynHiggsWrapper::CalcHiggsSpectrum(): mh[" << i << "] is undefined"
                      << std::endl;
            return (false);
        }

    return (true);
}

bool FeynHiggsWrapper::CalcSpectrum()
{
    int err, nmfv;
    double MSf[3][5][2], MASf[5][6], MCha[2], MNeu[4];
    ComplexType USf[3][5][2][2], UASf[5][6][6];
    ComplexType UCha[2][2], VCha[2][2], ZNeu[4][4], Deltab;

    /*
     * Note that the order of indices for an array has to be taken care so much.
     * Foe example, the indices in MSf(s,t,g) for the MFV sfermion masses are
     * defined as
     *   s = 1..2  sfermion index
     *   t = 1..5  sfermion type: nu, e, u, d, d-resummed
     *   g = 1..3  generation index
     * according to the manual of FeynHiggs. On the other hand, it is translated
     * from Fortran to C in CFeynHiggs.h as MSf[3][5][2]. Namely, the order
     * of the indices is reversed.
     */

    /* Get some of the MSSM parameters computed with FeynHiggs:
     *   MSf: the MFV sfermion masses
     *   USf: the MFV sfermion mixing matrices
     *   MASf: the NMFV sfermion masses
     *   UASf: the sfermion mixing matrices
     *   MCha: the chargino masses
     *   UCha, VCha: the chargino mixing matrices
     *   MNeu: the neutralino masses
     *   ZNeu: the neutralino mixing matrix
     *   Deltab: the correction to the bottom Yukawa coupling, Delta_b
     *   FHMGl: the gluino mass
     *   FHMHtree: the tree-level Higgs masses
     *   FHSAtree: the tree-level sin(alpha)
     */
    FHGetPara(&err, &nmfv, MSf, USf, MASf, UASf, MCha, UCha, VCha, MNeu, ZNeu,
              &Deltab, &FHMGl, FHMHtree, &FHSAtree);
    if (err != 0) {
        std::cout << "FeynHiggsWrapper::CalcSpectrum(): Error has been detected in GetPara.F:"
                  << err << std::endl;
        return (false);
    }

    /* sfermions in MFV for debug*/
    //gslpp::vector<double> m_sn2_MFV(6,0.), m_se2_MFV(6,0.), m_su2_MFV(6,0.), m_sd2_MFV(6,0.), m_sdresum2_MFV(6,0.);
    //gslpp::matrix<gslpp::complex> Rn_MFV(6,6,0.), Rl_MFV(6,6,0.), Ru_MFV(6,6,0.), Rd_MFV(6,6,0.), Rdresum_MFV(6,6,0.);
    //for (int g = 0; g < 3; g++) { /* generations */
    //    for (int s = 0; s < 2; s++) { /* left or right */
    //        m_sn2_MFV(g + 3*s) = MSf[g][0][s]*MSf[g][0][s];
    //        m_se2_MFV(g + 3*s) = MSf[g][1][s]*MSf[g][1][s];
    //        m_su2_MFV(g + 3*s) = MSf[g][2][s]*MSf[g][2][s];
    //        m_sd2_MFV(g + 3*s) = MSf[g][3][s]*MSf[g][3][s];
    //        m_sdresum2_MFV(g + 3*s) = MSf[g][4][s]*MSf[g][4][s];
    //        for (int t = 0; t < 2; t++) { /* left or right */
    //            Rn_MFV.assign(s,t, gslpp::complex(USf[g][0][t][s].real(), USf[g][0][t][s].imag()));
    //            Rl_MFV.assign(s,t, gslpp::complex(USf[g][1][t][s].real(), USf[g][1][t][s].imag()));
    //            Ru_MFV.assign(s,t, gslpp::complex(USf[g][2][t][s].real(), USf[g][2][t][s].imag()));
    //            Rd_MFV.assign(s,t, gslpp::complex(USf[g][3][t][s].real(), USf[g][3][t][s].imag()));
    //            Rdresum_MFV.assign(s,t, gslpp::complex(USf[g][4][t][s].real(), USf[g][4][t][s].imag()));
    //        }
    //    }
    //}

    /* sfermions in NMFV */
    for (int i = 0; i < 6; i++) {
        mySUSY.m_sn2(i) = MASf[0][i]*MASf[0][i];// heavy decoupled masses for i=3-5
        mySUSY.m_se2(i) = MASf[1][i]*MASf[1][i];
        mySUSY.m_su2(i) = MASf[2][i]*MASf[2][i];
        mySUSY.m_sd2(i) = MASf[3][i]*MASf[3][i];
        mySUSY.m_sdresum2(i) = MASf[4][i]*MASf[4][i];
        for (int j = 0; j < 6; j++) {
            /* R: first (second) index for mass (gauge) eigenstates */
            /* UASf: second (third) index for gauge (mass) eigenstates */
            mySUSY.Rn.assign(i,j, gslpp::complex(UASf[0][j][i].real(), UASf[0][j][i].imag()));
            mySUSY.Rl.assign(i,j, gslpp::complex(UASf[1][j][i].real(), UASf[1][j][i].imag()));
            mySUSY.Ru.assign(i,j, gslpp::complex(UASf[2][j][i].real(), UASf[2][j][i].imag()));
            mySUSY.Rd.assign(i,j, gslpp::complex(UASf[3][j][i].real(), UASf[3][j][i].imag()));
            mySUSY.Rdresum.assign(i,j, gslpp::complex(UASf[4][j][i].real(), UASf[4][j][i].imag()));
        }
    }

    /* Phase rotations due to differences between the SLHA and FeynHiggs
     * notations, which is necessary in the case of NMFV. */
    gslpp::complex muHphase(1.0, mySUSY.muH.arg(), true);
    gslpp::matrix<gslpp::complex> muHphaseMatrix(6,6,0.);
    for (int i = 0; i < 3; i++) {
        muHphaseMatrix.assign(i, i, muHphase.conjugate());
        muHphaseMatrix.assign(i+3, i+3, muHphase);
    }
    if (NMFVe) mySUSY.Rl = mySUSY.Rl * muHphaseMatrix;
    if (NMFVu) mySUSY.Ru = mySUSY.Ru * muHphaseMatrix;
    if (NMFVd) mySUSY.Rd = mySUSY.Rd * muHphaseMatrix;
    if (NMFVd) mySUSY.Rdresum = mySUSY.Rdresum * muHphaseMatrix;

    /* The SLHA2 convention requires increasing eigenvalues of the sfermion masses.
     * This is not done automatically by FeynHiggs. */
    SortSfermionMasses(mySUSY.m_sn2, mySUSY.Rn);
    SortSfermionMasses(mySUSY.m_se2, mySUSY.Rl);
    SortSfermionMasses(mySUSY.m_su2, mySUSY.Ru);
    SortSfermionMasses(mySUSY.m_sd2, mySUSY.Rd);
    SortSfermionMasses(mySUSY.m_sdresum2, mySUSY.Rdresum);

    /* charginos */
    for (int i = 0; i < 2; i++) {
        mySUSY.mch(i) = MCha[i];
        for (int j = 0; j < 2; j++) {
            /* U and V: first (second) index for mass (gauge) eigenstates */
            /* Ucha and VCha: first (second) index for gauge (mass) eigenstates */
            mySUSY.U.assign(i,j, gslpp::complex(UCha[j][i].real(), UCha[j][i].imag()));
            mySUSY.V.assign(i,j, gslpp::complex(VCha[j][i].real(), VCha[j][i].imag()));
        }
    }

    /* neutralinos */
    for (int i = 0; i < 4; i++) {
        mySUSY.mneu(i) = MNeu[i];
        for (int j = 0; j < 4; j++)
            /* N: first (second) index for mass (gauge) eigenstates */
            /* Zneu: first (second) index for gauge (mass) eigenstates */
            mySUSY.N.assign(i,j, gslpp::complex(ZNeu[j][i].real(), ZNeu[j][i].imag()));
    }

    /* the correction to the bottom Yukawa coupling */
    FHDeltab = gslpp::complex(Deltab.real(), Deltab.imag());

    return (true);
}

void FeynHiggsWrapper::SetFeynHiggsParsSLHA(const char *filename) const
{
    int err;
    COMPLEX slhadata[nslhadata];

    SLHARead(&err, slhadata, filename, 1);
    if(err != 0)
        throw std::runtime_error("FeynHiggsWrapper::SetFeynHiggsParsSLHA(): Error in SLHARead");

    FHSetSLHA(&err, slhadata);
    if(err != 0)
        throw std::runtime_error("FeynHiggsWrapper::SetFeynHiggsParsSLHA(): Error in FHSetSLHA");
}

void FeynHiggsWrapper::OutputSLHA(const char* filename) const
{
    int err;
    COMPLEX slhadata[nslhadata];

    SLHAClear(slhadata);

    FHOutputSLHA(&err, slhadata, -1);
    if (err != 0)
        throw std::runtime_error("FeynHiggsWrapper::SetFeynHiggsPars(): Error in FHOutputSLHA");

    SLHAWrite(&err, slhadata, filename);
    if (err != 0)
        throw std::runtime_error("FeynHiggsWrapper::SetFeynHiggsPars(): Error in SLHAWrite");
}

bool FeynHiggsWrapper::CalcHiggsCouplings()
{
    int err;
    ComplexType couplings[ncouplings];
    ComplexType couplingsms[ncouplingsms];
    double gammas[ngammas];
    double gammasms[ngammasms];

    /* Compute the Higgs couplings, decay widths and branching ratios */
    FHCouplings(&err, couplings, couplingsms, gammas, gammasms, 0);
    if (err != 0) {
        std::cout << "FeynHiggsWrapper::CalcHiggsCouplings(): Error has been detected in Couplings.F:"
                  << err << std::endl;
        return (false);
    }

    /* Write codes to store the results into member variables of the current class */

    computeHiggsCouplings = false;

    return (true);
}

bool FeynHiggsWrapper::CalcHiggsProd(const double& sqrts)
{
    int err;
    double prodxs[nprodxs];

    /* Compute Higgs production cross-sections with FeynHiggs */
    FHHiggsProd(&err, sqrts, prodxs);
    if (err != 0) {
        std::cout << "FeynHiggsWrapper::CalcHiggsProd(): Error has been detected in HiggsProd.F:"
                  << err << std::endl;
        return (false);
    }

    /* Write codes to store the results into member variables of the current class */

    computeHiggsProd = false;

    return (true);
}

bool FeynHiggsWrapper::CalcConstraints()
{
    int err, ccb;

    /* Calculate electroweak precision observables */
    FHConstraints(&err, &FHgm2, &FHdeltarho, &FHMWMSSM, &FHMWSM, &FHSW2MSSM,
                  &FHSW2SM, &FHedmeTh, &FHedmn, &FHedmHg, &ccb);
    if (err != 0) {
        std::cout << "FeynHiggsWrapper::CalcConstraints(): Error has been detected in Constraints.F:"
                  << err << std::endl;
        return (false);
    }
    if (ccb != 0) {
        std::cout << "FeynHiggsWrapper::CalcConstraints(): The parameter point corresponds to a colour-breaking minimum"
                  << std::endl;
        return (false);
    }

    computeConstraints = false;

    return (true);
}

bool FeynHiggsWrapper::CalcFlavour()
{
    int err;

    /* Calculate flavour observables */
    FHFlavour(&err, &FHbsgMSSM, &FHbsgSM, &FHdeltaMsMSSM, &FHdeltaMsSM,
              &FHbsmumuMSSM, &FHbsmumuSM);
    if (err != 0) {
        std::cout << "FeynHiggsWrapper::CalcFlavour(): Error has been detected in Flavour.F:"
                  << err << std::endl;
        return (false);
    }

    computeFlavour = false;

    return (true);
}

void FeynHiggsWrapper::SortSfermionMasses(gslpp::vector<double>& m_sf2, gslpp::matrix<gslpp::complex>& Rf) const
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


