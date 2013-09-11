/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
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
               2, // two loops where available
               1, // running top mass is used in the 1-/2-loop corrections
               1, // resum tan beta contributions
               3  // interpolation in phases for missing 2-loop corrections
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
                mySUSY.AlsMz, mySUSY.GF,
                mySUSY.leptons[StandardModel::ELECTRON].getMass(),
                mySUSY.quarks[mySUSY.UP].getMass(),
                mySUSY.quarks[mySUSY.DOWN].getMass(),
                mySUSY.leptons[mySUSY.MU].getMass(),
                mySUSY.quarks[mySUSY.CHARM].getMass(),
                mySUSY.quarks[mySUSY.STRANGE].getMass(),
                mySUSY.leptons[mySUSY.TAU].getMass(),
                mySUSY.quarks[mySUSY.BOTTOM].getMass(),
                Mw_FHinput, mySUSY.Mz,
                mySUSY.lambda, mySUSY.A, mySUSY.rhob, mySUSY.etab);
    if (err != 0) {
#ifdef FHDEBUG
        std::cout << "FeynHiggsWrapper::SetFeynHiggsPars(): Error has been detected in SetPara.F:"
                  << err << std::endl;
#endif
        return (false);
    }

    /* Parameters for FeynHiggs */
    //complex muHFH = mySUSY.muH.conjugate(); /* Incorrect! See the chargino and neutralino mass matrices. */
    complex muHFH = mySUSY.muH;
    matrix<complex> MsQ2 = mySUSY.msQhat2;
    matrix<complex> MsU2 = mySUSY.msUhat2;
    matrix<complex> MsD2 = mySUSY.msDhat2;
    matrix<complex> MsL2 = mySUSY.msLhat2;
    matrix<complex> MsE2 = mySUSY.msEhat2;
    matrix<complex> KU = mySUSY.TUhat.hconjugate() * mySUSY.v2() / sqrt(2.0);
    matrix<complex> KD = mySUSY.TDhat.hconjugate() * mySUSY.v1() / sqrt(2.0);
    matrix<complex> KE = mySUSY.TEhat.hconjugate() * mySUSY.v1() / sqrt(2.0);
    matrix<complex> deltaQLL(3,3,0.);
    matrix<complex> deltaULR(3,3,0.), deltaURL(3,3,0.), deltaURR(3,3,0.);
    matrix<complex> deltaDLR(3,3,0.), deltaDRL(3,3,0.), deltaDRR(3,3,0.);
    matrix<complex> deltaLLL(3,3,0.);
    matrix<complex> deltaELR(3,3,0.), deltaERL(3,3,0.), deltaERR(3,3,0.);
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

    /* Set the FeynHiggs parameters */
    double Q = mySUSY.Q;
    double x1 = mySUSY.v1()/sqrt(2.);
    double x2 = mySUSY.v2()/sqrt(2.);
    FHSetPara(&err,
              mySUSY.mut/mySUSY.quarks[StandardModel::TOP].getMass(),
              mySUSY.mtpole, mySUSY.tanb,
              -1, // using mHptree instead of mA
              mySUSY.mHptree,
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
              ToComplex2(KE(2,2).real(), KE(2,2).imag())/mySUSY.Ml_Q(mySUSY.TAU),
              ToComplex2(KU(2,2).real(), KU(2,2).imag())/mySUSY.Mq_Q(mySUSY.TOP),
              ToComplex2(KD(2,2).real(), KD(2,2).imag())/mySUSY.Mq_Q(mySUSY.BOTTOM),
              ToComplex2(KE(1,1).real(), KE(1,1).imag())/mySUSY.Ml_Q(mySUSY.MU),
              ToComplex2(KU(1,1).real(), KU(1,1).imag())/mySUSY.Mq_Q(mySUSY.CHARM),
              ToComplex2(KD(1,1).real(), KD(1,1).imag())/mySUSY.Mq_Q(mySUSY.STRANGE),
              ToComplex2(KE(0,0).real(), KE(0,0).imag())/mySUSY.Ml_Q(mySUSY.ELECTRON),
              ToComplex2(KU(0,0).real(), KU(0,0).imag())/mySUSY.Mq_Q(mySUSY.UP),
              ToComplex2(KD(0,0).real(), KD(0,0).imag())/mySUSY.Mq_Q(mySUSY.DOWN),
              //
              ToComplex2(mySUSY.m1.real(), mySUSY.m1.imag()),
              ToComplex2(mySUSY.m2.real(), mySUSY.m2.imag()),
              ToComplex2(mySUSY.m3, 0.),
              //
              Q, Q, Q);
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

    /* Set the non-minimal flavor-violating parameters in the slepton sector */
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
    mySUSY.saeff = complex(SAeff.real(), SAeff.imag());

    /* Debug */
    //std::cout << "mh[0] = mh = " << mySUSY.mh[0] << std::endl;
    //std::cout << "mh[1] = mH = " << mySUSY.mh[1] << std::endl;
    //std::cout << "mh[2] = mA = " << mySUSY.mh[2] << std::endl;
    //std::cout << "mh[3] = mH+ = " << mySUSY.mh[3] << std::endl;

    /* Check */
    for(int i = 0; i < 4; i++)
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
    double MSf[3][4][2], MASf[4][6], MCha[2], MNeu[4];
    ComplexType USf[3][4][2][2], UASf[4][6][6];
    ComplexType UCha[2][2], VCha[2][2], ZNeu[4][4], Deltab;

    /*
     * Note that the order of indices for an array has to be taken care so much.
     * Foe example, the indices in MSf(s,t,g) for the MFV sfermion masses are
     * defined as
     *   s = 1..2  sfermion index
     *   t = 1..4  sfermion type: nu, e, u, d
     *   g = 1..3  generation index
     * according to the manual of FeynHiggs. On the other hand, it is traslated
     * from Fortran to C in CFeynHiggs.h as MSf,[3][4][2]. Namely, the order
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

    /* squark and slepton */
    for(int i = 0; i < 6; i++) {
        mySUSY.m_sn2(i) = MASf[0][i]*MASf[0][i];
        mySUSY.m_se2(i) = MASf[1][i]*MASf[1][i];
        mySUSY.m_su2(i) = MASf[2][i]*MASf[2][i];
        mySUSY.m_sd2(i) = MASf[3][i]*MASf[3][i];
        for(int j = 0; j < 6; j++) {
            /* R: first (second) index for mass (gauge) eigenstates */
            /* UASf: second (third) index for gauge (mass) eigenstates */
            mySUSY.Rn.assign(i,j, complex(UASf[0][j][i].real(), UASf[0][j][i].imag()));
            mySUSY.Rl.assign(i,j, complex(UASf[1][j][i].real(), UASf[1][j][i].imag()));
            mySUSY.Ru.assign(i,j, complex(UASf[2][j][i].real(), UASf[2][j][i].imag()));
            mySUSY.Rd.assign(i,j, complex(UASf[3][j][i].real(), UASf[3][j][i].imag()));
        }
    }

    /* The SLHA2 convention requires increasing eigenvalues of the sfermion masses.
     * This is not done automatically by FeynHiggs. */
    SortSfermionMasses(mySUSY.m_sn2, mySUSY.Rn);
    SortSfermionMasses(mySUSY.m_se2, mySUSY.Rl);
    SortSfermionMasses(mySUSY.m_su2, mySUSY.Ru);
    SortSfermionMasses(mySUSY.m_sd2, mySUSY.Rd);

    /* charginos */
    for(int i = 0; i < 2; i++) {
        mySUSY.mch(i) = MCha[i];
        for(int j = 0; j < 2; j++) {
            /* U and V: first (second) index for mass (gauge) eigenstates */
            /* Ucha and VCha: first (second) index for gauge (mass) eigenstates */
            mySUSY.U.assign(i,j, complex(UCha[j][i].real(), UCha[j][i].imag()));
            mySUSY.V.assign(i,j, complex(VCha[j][i].real(), VCha[j][i].imag()));
        }
    }

    /* neutralinos */
    for(int i = 0; i < 4; i++) {
        mySUSY.mneu(i) = MNeu[i];
        for(int j = 0; j < 4; j++)
            /* N: first (second) index for mass (gauge) eigenstates */
            /* Zneu: first (second) index for gauge (mass) eigenstates */
            mySUSY.N.assign(i,j, complex(ZNeu[j][i].real(), ZNeu[j][i].imag()));
    }

    /* the correction to the bottom Yukawa coupling */
    FHDeltab = complex(Deltab.real(), Deltab.imag());

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

void FeynHiggsWrapper::SortSfermionMasses(vector<double>& m_sf2, matrix<complex>& Rf) const
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
    gslpp::matrix<complex> myRf(6, 6, 0.);
    for (int i = 0; i < 6; i++)
        for (int k = 0; k < 6; k++)
            myRf.assign(k, i, Rf(newIndex[k], i));
    Rf = myRf;
}


