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
#include "FeynHiggs.h"

FeynHiggs::FeynHiggs(SUSY& SUSY_in)
: mySUSY(SUSY_in)
{
    int err;
    FHSetFlags(&err,
               4, // Full MSSM
               0, // DRbar
               0, // DRbar
               3, // Higgs mixing in cMSSM
               4, // UHiggs at q^2=0
               2, // two loops where available
               1, // running top mass,
               1, // resum tan beta contributions
               3  // interpolation in phases for missing 2-loop corrections
              );
    if (err != 0) {
        std::stringstream ss;
        ss << "FeynHiggs::FeynHiggs(): FHSetFlags error " << err;
        throw std::runtime_error(ss.str());
    }
}

bool FeynHiggs::SetFeynHiggsPars()
{
    int err;

    /* Echo input parameters in detail,
     * display Higgs mass matrix at p^2 = 0 and CTs */
    //FHSetDebug(2);

    Mw_FHinput = mySUSY.Mw_tree(); /* Tree-level W-boson mass */
    //Mw_FHinput = mySUSY.StandardModel::Mw(); /* SM prediction, which should not be used, since mHl cannot be set before calling FeynHiggs. */
    //std::cout << "Mw = " << Mw_FHinput << " used in FeynHiggs::SetFeynHiggsPars()" << std::endl;

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
        std::cout << "FeynHiggs::SetFeynHiggsPars(): Error has been detected in SetPara.F:"
                  << err << std::endl;
        return (false);
    }

    /* FeynHiggs input: the FeynHiggs notation is related to the SLHA one by
     * hconjugate operation for T matrices and conjugate operation for mu parameter */
    matrix<complex> TUFH = mySUSY.TU.hconjugate();
    matrix<complex> TDFH = mySUSY.TD.hconjugate();
    matrix<complex> TEFH = mySUSY.TE.hconjugate();
    complex muHFH = mySUSY.muH.conjugate();

    /* Set the FeynHiggs parameters */
    double Q = mySUSY.Q;
    double x1 = mySUSY.v1()/sqrt(2.);
    double x2 = mySUSY.v2()/sqrt(2.);
    FHSetPara(&err,
              mySUSY.mut/mySUSY.quarks[StandardModel::TOP].getMass(),
              mySUSY.mtpole, mySUSY.tanb,
              -1, //using MHptree
              mySUSY.mHptree,
              sqrt(mySUSY.MsL2(2,2).real()), sqrt(mySUSY.MsE2(2,2).real()),
              sqrt(mySUSY.MsQ2(2,2).real()), sqrt(mySUSY.MsU2(2,2).real()),
              sqrt(mySUSY.MsD2(2,2).real()),
              sqrt(mySUSY.MsL2(1,1).real()), sqrt(mySUSY.MsE2(1,1).real()),
              sqrt(mySUSY.MsQ2(1,1).real()), sqrt(mySUSY.MsU2(1,1).real()),
              sqrt(mySUSY.MsD2(1,1).real()),
              sqrt(mySUSY.MsL2(0,0).real()), sqrt(mySUSY.MsE2(0,0).real()),
              sqrt(mySUSY.MsQ2(0,0).real()), sqrt(mySUSY.MsU2(0,0).real()),
              sqrt(mySUSY.MsD2(0,0).real()),
              ToComplex2(muHFH.real(), muHFH.imag()),
              ToComplex2(TEFH(2,2).real(), TEFH(2,2).imag())
                *x1/mySUSY.leptons[StandardModel::TAU].getMass(),
              ToComplex2(TUFH(2,2).real(), TUFH(2,2).imag())
                *x2/mySUSY.MS2DRqmass(Q, mySUSY.mu_Q[2]),
              ToComplex2(TDFH(2,2).real(), TDFH(2,2).imag())
                *x1/mySUSY.MS2DRqmass(Q, mySUSY.md_Q[2]),
              ToComplex2(TEFH(1,1).real(), TEFH(1,1).imag())
                *x1/mySUSY.leptons[StandardModel::MU].getMass(),
              ToComplex2(TUFH(1,1).real(), TUFH(1,1).imag())
                *x2/mySUSY.MS2DRqmass(Q, mySUSY.mu_Q[1]),
              ToComplex2(TDFH(1,1).real(), TDFH(1,1).imag())
                *x1/mySUSY.MS2DRqmass(Q, mySUSY.md_Q[1]),
              ToComplex2(TEFH(0,0).real(), TEFH(0,0).imag())
                *x1/mySUSY.leptons[StandardModel::ELECTRON].getMass(),
              ToComplex2(TUFH(0,0).real(), TUFH(0,0).imag())
                *x2/mySUSY.MS2DRqmass(Q, mySUSY.mu_Q[0]),
              ToComplex2(TDFH(0,0).real(), TDFH(0,0).imag())
                *x1/mySUSY.MS2DRqmass(Q, mySUSY.md_Q[0]),
              ToComplex2(mySUSY.m1.real(), mySUSY.m1.imag()),
              ToComplex2(mySUSY.m2.real(), mySUSY.m2.imag()),
              ToComplex2(mySUSY.m3, 0.),
              Q, Q, Q);
    if (err != 0) {
        std::cout << "FeynHiggs::SetFeynHiggsPars(): Error has been detected in SetPara.F:"
                  << err << std::endl;
        return (false);
    }

    /* Set the non-minimal flavor-violating parameters */
    FHSetNMFV(&err, 
              ToComplex2(mySUSY.MsL2(0,1).real(),
                         mySUSY.MsL2(0,1).imag())/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsL2(1,1).real()),
              ToComplex2(mySUSY.MsL2(1,2).real(),
                         mySUSY.MsL2(1,2).imag())/sqrt(mySUSY.MsL2(1,1).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(mySUSY.MsL2(0,2).real(),
                         mySUSY.MsL2(0,2).imag())/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(TUFH(0,1).real(),
                         -TUFH(0,1).imag())*x2/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsU2(1,1).real()),
              ToComplex2(TUFH(1,2).real(),
                         -TUFH(1,2).imag())*x2/sqrt(mySUSY.MsL2(1,1).real()*mySUSY.MsU2(2,2).real()),
              ToComplex2(TUFH(0,2).real(),
                         -TUFH(0,2).imag())*x2/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsU2(2,2).real()),
              ToComplex2(TUFH(1,0).real(),
                         TUFH(1,0).imag())*x2/sqrt(mySUSY.MsU2(0,0).real()*mySUSY.MsL2(1,1).real()),
              ToComplex2(TUFH(2,1).real(),
                         TUFH(2,1).imag())*x2/sqrt(mySUSY.MsU2(1,1).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(TUFH(2,0).real(),
                         TUFH(2,0).imag())*x2/sqrt(mySUSY.MsU2(0,0).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(mySUSY.MsU2(0,1).real(),
                         mySUSY.MsU2(0,1).imag())/sqrt(mySUSY.MsU2(0,0).real()*mySUSY.MsU2(1,1).real()),
              ToComplex2(mySUSY.MsU2(1,2).real(),
                         mySUSY.MsU2(1,2).imag())/sqrt(mySUSY.MsU2(1,1).real()*mySUSY.MsU2(2,2).real()),
              ToComplex2(mySUSY.MsU2(0,2).real(),
                         mySUSY.MsU2(0,2).imag())/sqrt(mySUSY.MsU2(0,0).real()*mySUSY.MsU2(2,2).real()),
              ToComplex2(TDFH(0,1).real(),
                         -TDFH(0,1).imag())*x1/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsD2(1,1).real()),
              ToComplex2(TDFH(1,2).real(),
                         -TDFH(1,2).imag())*x1/sqrt(mySUSY.MsL2(1,1).real()*mySUSY.MsD2(2,2).real()),
              ToComplex2(TDFH(0,2).real(),
                         -TDFH(0,2).imag())*x1/sqrt(mySUSY.MsL2(0,0).real()*mySUSY.MsD2(2,2).real()),
              ToComplex2(TDFH(1,0).real(),
                         TDFH(1,0).imag())*x1/sqrt(mySUSY.MsD2(0,0).real()*mySUSY.MsL2(1,1).real()),
              ToComplex2(TDFH(2,1).real(),
                         TDFH(2,1).imag())*x1/sqrt(mySUSY.MsD2(1,1).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(TDFH(2,0).real(),
                         TDFH(2,0).imag())*x1/sqrt(mySUSY.MsD2(0,0).real()*mySUSY.MsL2(2,2).real()),
              ToComplex2(mySUSY.MsD2(0,1).real(),
                         mySUSY.MsD2(0,1).imag())/sqrt(mySUSY.MsD2(0,0).real()*mySUSY.MsD2(1,1).real()),
              ToComplex2(mySUSY.MsD2(1,2).real(),
                         mySUSY.MsD2(1,2).imag())/sqrt(mySUSY.MsD2(1,1).real()*mySUSY.MsD2(2,2).real()),
              ToComplex2(mySUSY.MsD2(0,2).real(),
                         mySUSY.MsD2(0,2).imag())/sqrt(mySUSY.MsD2(0,0).real()*mySUSY.MsD2(2,2).real()));
    if (err != 0) {
        std::cout << "FeynHiggs::SetFeynHiggsPars(): Error was detected in SetFV.F:"
                  << err << std::endl;
        return (false);
    }

    computeHiggsCouplings = true;
    computeHiggsProd = true;
    computeConstraints = true;
    computeFlavour = true;

    return (true);
}

bool FeynHiggs::CalcHiggsSpectrum()
{
    int err;
    ComplexType SAeff;
    ComplexType UHiggs[3][3];
    ComplexType ZHiggs[3][3];

    /* Compute the Higgs masses and mixings */
    FHHiggsCorr(&err, mySUSY.mh, &SAeff, UHiggs, ZHiggs);
    if (err != 0) {
        std::cout << "FeynHiggs::CalcHiggsSpectrum(): Error has been detected in HiggsCorr.F:"
                  << err << std::endl;
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
            std::cout << "FeynHiggs::CalcHiggsSpectrum(): mh[" << i << "] is undefined"
                      << std::endl;
            return (false);
        }

    return (true);
}

bool FeynHiggs::CalcSpectrum()
{
    int err, nmfv;
    double MSf[3][4][2], MASf[4][6], MCha[2], MNeu[4];
    ComplexType USf[3][4][2][2], UASf[4][6][6];
    ComplexType UCha[2][2], VCha[2][2], ZNeu[4][4], Deltab;

    /*
     * Note that the order of indices for an array has to be take care so much!
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
        std::cout << "FeynHiggs::CalcSpectrum(): Error has been detected in GetPara.F:"
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
            /* UASf: second (third) index for gauge (mass) eigenstates */
            mySUSY.Rn.assign(j,i, complex(UASf[0][i][j].real(), UASf[0][i][j].imag()));
            mySUSY.Rl.assign(j,i, complex(UASf[1][i][j].real(), UASf[1][i][j].imag()));
            mySUSY.Ru.assign(j,i, complex(UASf[2][i][j].real(), UASf[2][i][j].imag()));
            mySUSY.Rd.assign(j,i, complex(UASf[3][i][j].real(), UASf[3][i][j].imag()));
        }
    }

    /* charginos */
    for(int i = 0; i < 2; i++) {
        mySUSY.mch(i) = MCha[i];
        for(int j = 0; j < 2; j++) {
            /* Ucha and VCha: first (second) index for gauge (mass) eigenstates */
            mySUSY.U.assign(j,i, complex(UCha[i][j].real(), UCha[i][j].imag()));
            mySUSY.V.assign(j,i, complex(VCha[i][j].real(), VCha[i][j].imag()));
        }
    }

    /* neutralinos */
    for(int i = 0; i < 4; i++) {
        mySUSY.mneu(i) = MNeu[i];
        for(int j = 0; j < 4; j++)
            /* Zneu: first (second) index for gauge (mass) eigenstates */
            mySUSY.N.assign(j,i, complex(ZNeu[i][j].real(), ZNeu[i][j].imag()));
    }

    /* the correction to the bottom Yukawa coupling */
    FHDeltab = complex(Deltab.real(), Deltab.imag());

    return (true);
}

void FeynHiggs::OutputSLHA(const char* filename) const
{
    int err;
    COMPLEX slhadata[nslhadata];

    FHOutputSLHA(&err, slhadata, 2);
    if (err != 0)
        throw std::runtime_error("FeynHiggs::SetFeynHiggsPars(): Error in FHOutputSLHA");

    SLHAWrite(&err, slhadata, filename);
    if (err != 0)
        throw std::runtime_error("FeynHiggs::SetFeynHiggsPars(): Error in SLHAWrite");

    /* This function does not work correctly! */

}

bool FeynHiggs::CalcHiggsCouplings()
{
    int err;
    ComplexType couplings[ncouplings];
    ComplexType couplingsms[ncouplingsms];
    double gammas[ngammas];
    double gammasms[ngammasms];

    /* Compute the Higgs couplings, decay widths and branching ratios */
    FHCouplings(&err, couplings, couplingsms, gammas, gammasms, 0);
    if (err != 0) {
        std::cout << "FeynHiggs::CalcHiggsCouplings(): Error has been detected in Couplings.F:"
                  << err << std::endl;
        return (false);
    }

    /* Write codes to store the results into member variables of the current class */

    computeHiggsCouplings = false;

    return (true);
}

bool FeynHiggs::CalcHiggsProd(const double& sqrts)
{
    int err;
    double prodxs[nprodxs];

    /* Compute Higgs production cross-sections with FeynHiggs */
    FHHiggsProd(&err, sqrts, prodxs);
    if (err != 0) {
        std::cout << "FeynHiggs::CalcHiggsProd(): Error has been detected in HiggsProd.F:"
                  << err << std::endl;
        return (false);
    }

    /* Write codes to store the results into member variables of the current class */

    computeHiggsProd = false;

    return (true);
}

bool FeynHiggs::CalcConstraints()
{
    int err, ccb;

    /* Calculate electroweak precision observables */
    FHConstraints(&err, &FHgm2, &FHdeltarho, &FHMWMSSM, &FHMWSM, &FHSW2MSSM,
                  &FHSW2SM, &FHedmeTh, &FHedmn, &FHedmHg, &ccb);
    if (err != 0) {
        std::cout << "FeynHiggs::CalcConstraints(): Error has been detected in Constraints.F:"
                  << err << std::endl;
        return (false);
    }
    if (ccb != 0) {
        std::cout << "FeynHiggs::CalcConstraints(): The parameter point corresponds to a colour-breaking minimum"
                  << std::endl;
        return (false);
    }
    
    computeConstraints = false;

    return (true);
}

bool FeynHiggs::CalcFlavour()
{
    int err;

    /* Calculate flavour observables */
    FHFlavour(&err, &FHbsgMSSM, &FHbsgSM, &FHdeltaMsMSSM, &FHdeltaMsSM,
              &FHbsmumuMSSM, &FHbsmumuSM);
    if (err != 0) {
        std::cout << "FeynHiggs::CalcFlavour(): Error has been detected in Flavour.F:"
                  << err << std::endl;
        return (false);
    }

    computeFlavour = false;

    return (true);
}

