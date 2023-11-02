/* 
 * Copyright (C) 2018. HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6GeneralMatching.h"
#include "NPSMEFTd6General.h"
#include <stdexcept>

NPSMEFTd6GeneralMatching::NPSMEFTd6GeneralMatching(const NPSMEFTd6General & NPSMEFTd6General_i) :
StandardModelMatching(NPSMEFTd6General_i),
mySMEFT(NPSMEFTd6General_i), VuL(3, 0.), VuR(3, 0.), VdL(3, 0.), VdR(3, 0.), VeL(3, 0.), VeR(3, 0.),
mcd2(5, NDR, LO), mcd1(10, NDR, LO), mcbd(5, NDR, LO), mcbs(5, NDR, LO), mck2(5, NDR, LO), mculeptonnu(5, NDR, LO) {
}

void NPSMEFTd6GeneralMatching::updateLEFTGeneralParameters() {

    // Dimension 6 operators with no flavour index are assigned directly here

    LambdaNP2 = mySMEFT.getLambda_NP() * mySMEFT.getLambda_NP();
    v = mySMEFT.v(); //This is vtilde in Angelica's notation
    v2 = v*v; //This is vtilde squared

    // The true VEV, corresponding to vbar in Angelica's notation, is equal to v up to corrections
    double vT = v;
    double delta_vT = mySMEFT.getDelta_v();
    double vTosq2 = vT / sqrt(2.);


    //    CG = mySMEFT.getSMEFTCoeffEW("CG")*LambdaNP2;  
    //    CW = mySMEFT.getSMEFTCoeffEW("CW")*LambdaNP2; 
    //    CHG = mySMEFT.getSMEFTCoeffEW("CHG")*LambdaNP2;  
    //    CHW = mySMEFT.getSMEFTCoeffEW("CHW")*LambdaNP2;  
    //    CHB = mySMEFT.getSMEFTCoeffEW("CHB")*LambdaNP2;  
    //    CHWB = mySMEFT.getSMEFTCoeffEW("CHWB")*LambdaNP2;  
    //    CHD = mySMEFT.getSMEFTCoeffEW("CHD")*LambdaNP2;  
    //    CHbox = mySMEFT.getSMEFTCoeffEW("CHbox")*LambdaNP2;  
    //    CH = mySMEFT.getSMEFTCoeffEW("CH")*LambdaNP2;  
    //    CGtilde = mySMEFT.getSMEFTCoeffEW("CGtilde")*LambdaNP2;  
    //    CWtilde = mySMEFT.getSMEFTCoeffEW("CWtilde")*LambdaNP2;  
    //    CHGtilde = mySMEFT.getSMEFTCoeffEW("CHGtilde")*LambdaNP2;  
    //    CHWtilde = mySMEFT.getSMEFTCoeffEW("CHWtilde")*LambdaNP2;  
    //    CHBtilde = mySMEFT.getSMEFTCoeffEW("CHBtilde")*LambdaNP2;  
    //    CHWtildeB = mySMEFT.getSMEFTCoeffEW("CHWtildeB")*LambdaNP2;  
    //    
    //    //Now we do not use the SILH basis anymore, we'll set these operators to zero 
    //    C2B = 0.; 
    //    C2W = 0.; 
    //    C2BS = 0.; 
    //    C2WS = 0.; 
    //    CDHB = 0.; 
    //    CDHW = 0.; 
    //    CDB = 0.; 
    //    CDW = 0.; 
    //    CT = 0.;

    //For operators with fermionic indices we need to switch to the mass eigenstate basis 

    // Let us first define the full mass matrices, including the effect of dimension six operators

    gslpp::matrix<complex> MU(3, 0.), MD(3, 0.), ME(3, 0.);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            MU.assignre(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YuR", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CuHR", i, j) * v2 / 2.));
            MU.assignim(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YuI", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CuHI", i, j) * v2 / 2.));
            MD.assignre(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YdR", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CdHR", i, j) * v2 / 2.));
            MD.assignim(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YdI", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CdHI", i, j) * v2 / 2.));
            ME.assignre(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YeR", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CeHR", i, j) * v2 / 2.));
            ME.assignim(i, j, vTosq2 * (mySMEFT.getSMEFTCoeffEW("YeI", i, j) * (1. + delta_vT) - mySMEFT.getSMEFTCoeffEW("CeHI", i, j) * v2 / 2.));
        }

    gslpp::vector<double> m2(3);

    MU.singularvalue(VuR, VuL, m2);
    mySMEFT.getQuarks(QCD::UP).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::UP).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(0))));
    mySMEFT.getQuarks(QCD::CHARM).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::CHARM).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(1))));
    mySMEFT.getQuarks(QCD::TOP).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::TOP).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(2))));

    MD.singularvalue(VdR, VdL, m2);
    mySMEFT.getQuarks(QCD::DOWN).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::DOWN).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(0))));
    mySMEFT.getQuarks(QCD::STRANGE).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::STRANGE).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(1))));
    mySMEFT.getQuarks(QCD::BOTTOM).setMass(mySMEFT.Mrun(mySMEFT.getQuarks(QCD::BOTTOM).getMass_scale(), mySMEFT.getMuw(), sqrt(m2(2))));

    ME.singularvalue(VeR, VeL, m2);

    mySMEFT.getLeptons(QCD::ELECTRON).setMass(sqrt(m2(0)));
    mySMEFT.getLeptons(QCD::MU).setMass(sqrt(m2(1)));
    mySMEFT.getLeptons(QCD::TAU).setMass(sqrt(m2(2)));

    //Computing the CKM 
    gslpp::matrix<complex> CKMUnphys = (VuL.hconjugate()) * VdL;

    //std::cout << "CKM unphys = " << CKMUnphys << std::endl;

    mySMEFT.getCKM().computeCKM(CKMUnphys(0, 1).abs(), CKMUnphys(1, 2).abs(), CKMUnphys(0, 2).abs(),
            (-CKMUnphys(0, 0) * CKMUnphys(0, 2).conjugate() / (CKMUnphys(1, 0) * CKMUnphys(1, 2).conjugate())).arg());

    //std::cout << "computed CKM = " << mySMEFT.getCKM().getCKM() << std::endl;

    double a11 = remainder(CKMUnphys(0, 0).arg() - mySMEFT.getCKM().getV_ud().arg(), 2. * M_PI);
    double a12 = remainder(CKMUnphys(0, 1).arg() - mySMEFT.getCKM().getV_us().arg(), 2. * M_PI);
    double a13 = remainder(CKMUnphys(0, 2).arg() - mySMEFT.getCKM().getV_ub().arg(), 2. * M_PI);


    //    double a23 = (gslpp::complex(CKMUnphys(1, 0) / CKM(1, 0))).arg() - a11 + a13;
    //    double a33 = (gslpp::complex(CKMUnphys(2, 0) / CKM(2, 0))).arg() - a11 + a13;
    double a23 = remainder(CKMUnphys(1, 0).arg() - mySMEFT.getCKM().getV_cd().arg(), 2. * M_PI) - a11 + a13;
    double a33 = remainder(CKMUnphys(2, 0).arg() - mySMEFT.getCKM().getV_td().arg(), 2. * M_PI) - a11 + a13;

    gslpp::matrix<gslpp::complex> phi1(3, 3, 0.);
    phi1.assign(0, 0, 1.);
    phi1.assign(1, 1, gslpp::complex(1., a23 - a13, true));
    phi1.assign(2, 2, gslpp::complex(1., a33 - a13, true));

    gslpp::matrix<gslpp::complex> phi2dag(3, 3, 0.);
    phi2dag.assign(0, 0, gslpp::complex(1., -a11, true));
    phi2dag.assign(1, 1, gslpp::complex(1., -a12, true));
    phi2dag.assign(2, 2, gslpp::complex(1., -a13, true));

    gslpp::matrix<gslpp::complex> phie(3, 3, 0.);
    phie.assign(0, 0, gslpp::complex(1., -(VeR(0, 0)).arg(), true));
    phie.assign(1, 1, gslpp::complex(1., -(VeR(1, 1)).arg(), true));
    phie.assign(2, 2, gslpp::complex(1., -(VeR(2, 2)).arg(), true));

    VeR = VeR*phie;
    VeL = VeL*phie;
    VuL = VuL*phi1;
    VuR = VuR*phi1;
    VdL = VdL*phi2dag;
    VdR = VdR*phi2dag;

    // to implement Manohar's matching formulae we define the couplings
    // in his notation. Namely, in the formulae below, the barred quantities are 
    // tree level in the theory scheme.

    double cbar = mySMEFT.getXWZ_tree();
    double sbar = -mySMEFT.getXBZ_tree();
    double sbar2 = sbar*sbar;
    //    double delta_cbar = mySMEFT.getDelta_xWZ(); not needed currently
    double delta_sbar = mySMEFT.getDelta_xBZ();
    double g1bar = mySMEFT.getG1_tree();
    //    double delta_g1bar = mySMEFT.getDelta_g1(); not needed currently
    double g2bar = mySMEFT.getG2_tree();
    //    double delta_g2bar = mySMEFT.getDelta_g2(); not needed currently
    double delta_MZ2 = mySMEFT.getDelta_Mz2();
    double ebar = mySMEFT.getEeMz();
    //    double delta_ebar = mySMEFT.getDelta_ale() / 2.; not needed currently
    // the Z coupling and its correction  were not explicit in Angelica's notes, so they need to be checked
    double gZbar = ebar / sbar / cbar;
    double delta_gZbar = (g1bar * g1bar + g2bar * g2bar) / (2. * g1bar * g2bar) * v2 * mySMEFT.getSMEFTCoeffEW("CHWB");
    double gZ2oMZ2 = gZbar / mySMEFT.getMz();
    gZ2oMZ2 *= gZ2oMZ2;
    double delta_gZ2oMZ2 = 2. * delta_gZbar - delta_MZ2;
    double g22oMW2 = 4. / v2;
    double delta_g22oMW2 = -2. * delta_vT;



    //std::cout << "CKM from rotated UfA = " << (VuL.hconjugate()) * VdL << std::endl;

    //std::cout << "has the diagonalization worked? " << VuR.hconjugate()*MU*VuL << std::endl;
    //std::cout << "has the diagonalization worked? " << VdR.hconjugate()*MD*VdL << std::endl;
    //std::cout << "has the diagonalization worked? " << VeR.hconjugate()*ME*VeL << std::endl;

    //match and rotate following Manohar. This is performed AT LINEAR ORDER for the moment

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            Ceg.at(i).at(j) = 0.;
            Cdg.at(i).at(j) = 0.;
            CdG.at(i).at(j) = 0.;

            for (int p = 0; p < 3; p++)
                for (int r = 0; r < 3; r++) {
                    Ceg.at(i).at(j) += vTosq2 * VeL.hconjugate()(i, p) * (-(mySMEFT.getSMEFTCoeffEW("CeWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", p, r))
                            * sbar + (mySMEFT.getSMEFTCoeffEW("CeBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", p, r))
                            * cbar) * VeR(r, j);
                    Cdg.at(i).at(j) += vTosq2 * VdL.hconjugate()(i, p) * (-(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r))
                            * sbar + (mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r))
                            * cbar) * VdR(r, j);
                    CdG.at(i).at(j) += vTosq2 * VdL.hconjugate()(i, p) * (mySMEFT.getSMEFTCoeffEW("CdGR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdGI", p, r))
                            * VdR(r, j);
                }

        }

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            Cug.at(i).at(j) = 0.;
            CuG.at(i).at(j) = 0.;

            for (int p = 0; p < 3; p++)
                for (int r = 0; r < 3; r++) {
                    Cug.at(i).at(j) += vTosq2 * VuL.hconjugate()(i, p) * ((mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))
                            * sbar + (mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r))
                            * cbar) * VuR(r, j);
                    CuG.at(i).at(j) += vTosq2 * VuL.hconjugate()(i, p) * (mySMEFT.getSMEFTCoeffEW("CuGR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuGI", p, r))
                            * VuR(r, j);

                }
        }

    gslpp::matrix<double> KD = gslpp::matrix<double>::Id(3);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++) {

                    //Leptonic LL operators
                    CnunuVLL.at(i).at(j).at(k).at(l) = 0.;
                    CeeVLL.at(i).at(j).at(k).at(l) = 0.;
                    CnueVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LL operators
                    CnudVLL.at(i).at(j).at(k).at(l) = 0.;
                    CedVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LL operators
                    CddVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Leptonic RR operators
                    CeeVRR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic RR operators
                    CedVRR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic RR operators
                    CddVRR.at(i).at(j).at(k).at(l) = 0.;

                    //Leptonic LLRR operators
                    CnueVLR.at(i).at(j).at(k).at(l) = 0.;
                    CeeVLR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LLRR operators
                    CnudVLR.at(i).at(j).at(k).at(l) = 0.;
                    CedVLR.at(i).at(j).at(k).at(l) = 0.;
                    CdeVLR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LLRR operators
                    CddV1LR.at(i).at(j).at(k).at(l) = 0.;
                    CddV8LR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRRL operators
                    CedSRL.at(i).at(j).at(k).at(l) = 0.;

                    //Leptonic LRLR operators
                    CeeSRR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRLR operators
                    CedSRR.at(i).at(j).at(k).at(l) = 0.;
                    CedTRR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LRLR operators
                    CddS1RR.at(i).at(j).at(k).at(l) = 0.;
                    CddS8RR.at(i).at(j).at(k).at(l) = 0.;

                    for (int p = 0; p < 3; p++)
                        for (int r = 0; r < 3; r++)
                            for (int s = 0; s < 3; s++)
                                for (int t = 0; t < 3; t++) {
                                    CnunuVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VeL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CllR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", p, r, s, t)) +
                                            (gZ2oMZ2 * (
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", s, t)) * KD(p, r) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", s, t)) * KD(p, r) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", s, r)) * KD(p, t) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", s, r)) * KD(p, t) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, t)) * KD(r, s) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, t)) * KD(r, s) -
                                            delta_gZ2oMZ2 * KD(p, t) * KD(r, s) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) * KD(s, t) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) * KD(s, t) -
                                            delta_gZ2oMZ2 * KD(p, r) * KD(s, t))) / 16.) * VeL(r, j) * VeL(t, l);
                                    CeeVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VeL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CllR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", p, r, s, t)) +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2)*(
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", s, t)) * KD(p, r) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", s, t)) * KD(p, r) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", s, r)) * KD(p, t) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", s, r)) * KD(p, t) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, t)) * KD(r, s) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, t)) * KD(r, s) +
                                            delta_gZ2oMZ2 * KD(p, t) * KD(r, s) -
                                            2. * delta_gZ2oMZ2 * sbar2 * KD(p, t) * KD(r, s) -
                                            8. * delta_sbar * sbar2 * KD(p, t) * KD(r, s) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) * KD(s, t) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) * KD(s, t) +
                                            delta_gZ2oMZ2 * KD(p, r) * KD(s, t) -
                                            2. * delta_gZ2oMZ2 * sbar2 * KD(p, r) * KD(s, t) -
                                            8. * delta_sbar * sbar2 * KD(p, r) * KD(s, t))) / 16.) * VeL(r, j) * VeL(t, l);
                                    CnueVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VeL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CllR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", p, r, s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("CllR", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", s, t, p, r)) -
                                            (g22oMW2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", r, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", r, s)) * KD(p, t) +
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, t)) +
                                            delta_g22oMW2 * KD(p, t)) *
                                            KD(r, s))) / 2. +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2)*(
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) -
                                            delta_gZ2oMZ2 * KD(p, r)) * KD(s, t)) / 4. +
                                            (gZ2oMZ2 * KD(p, r)*
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", s, t)) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", s, t)) -
                                            4. * delta_sbar * sbar2 * KD(s, t))) /
                                            4.) * VeL(r, j) * VeL(t, l);
                                    CnudVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", p, r, s, t)) -
                                            (mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t)) +
                                            (gZ2oMZ2 * (-3. + 2. * sbar2)*(
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) -
                                            delta_gZ2oMZ2 * KD(p, r)) * KD(s, t)) / 12. +
                                            (gZ2oMZ2 * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) +
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)) -
                                            4. * delta_sbar * sbar2 * KD(s, t))) / 12.) * VeL(r, j) * VdL(t, l);
                                    CedVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", p, r, s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", p, r, s, t)) -
                                            (gZ2oMZ2 * (-3. + 2. * sbar2)*(-0.5 * (v2 * (
                                            (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)))) +
                                            delta_gZ2oMZ2 * (-0.5 + sbar2) * KD(p, r) +
                                            2. * delta_sbar * sbar2 * KD(p, r)) * KD(s, t)) / 6. -
                                            gZ2oMZ2 * (-0.5 + sbar2) * KD(p, r)*
                                            (-0.5 * (v2 * (
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)))) +
                                            (2. * delta_sbar * sbar2 * KD(s, t)) / 3.)) * VeL(r, j) * VdL(t, l);
                                    CddVLL.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t)) -
                                            (gZ2oMZ2 * (-3. + 2. * sbar2)*(
                                            -3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) * KD(p, r) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)) * KD(p, r) +
                                            (-3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 + 2. * delta_gZ2oMZ2 * sbar2 +
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            72.) * VdL(r, j) * VdL(t, l);
                                    CeeVRR.at(i).at(j).at(k).at(l) += VeR.hconjugate()(i, p) * VeR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CeeR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeeI", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 * (
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, t)) * KD(p, r) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, r)) * KD(p, t) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", p, t)) * KD(r, s) -
                                            2. * delta_gZ2oMZ2 * sbar2 * KD(p, t) * KD(r, s) -
                                            8. * delta_sbar * sbar2 * KD(p, t) * KD(r, s) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", p, r)) * KD(s, t) -
                                            2. * delta_gZ2oMZ2 * sbar2 * KD(p, r) * KD(s, t) -
                                            8. * delta_sbar * sbar2 * KD(p, r) * KD(s, t))) / 8.) * VeR(r, j) * VeR(t, l);
                                    CedVRR.at(i).at(j).at(k).at(l) += VeR.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CedR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CedI", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 * (
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r) +
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", p, r)) -
                                            2. * (delta_gZ2oMZ2 + 4. * delta_sbar) * sbar2 *
                                            KD(p, r)) * KD(s, t))) / 6.) * VeR(r, j) * VdR(t, l);
                                    CddVRR.at(i).at(j).at(k).at(l) += VdR.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CddR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CddI", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 *
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r) +
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r)) -
                                            2. * (delta_gZ2oMZ2 + 4. * delta_sbar) * sbar2 *
                                            KD(p, r)) * KD(s, t))) / 36.) * VdR(r, j) * VdR(t, l);
                                    CnueVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VeR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CleR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CleI", p, r, s, t)) +
                                            (gZ2oMZ2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, t)) * KD(p, r) -
                                            2. * sbar2 * (-(v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r))) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            (delta_gZ2oMZ2 + 2. * delta_sbar) * KD(p, r)) *
                                            KD(s, t))) / 4.) * VeL(r, j) * VeR(t, l);
                                    CeeVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VeR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CleR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CleI", p, r, s, t)) +
                                            (gZ2oMZ2 * ((-1. + 2. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHeR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, t)) * KD(p, r) +
                                            2. * sbar2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            (delta_gZ2oMZ2 + 2. * delta_sbar - 2. * delta_gZ2oMZ2 * sbar2 -
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) / 4.) * VeL(r, j) * VeR(t, l);
                                    CnudVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CldR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CldI", p, r, s, t)) +
                                            (gZ2oMZ2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r)) / 4. -
                                            (gZ2oMZ2 * sbar2 * (-(v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r))) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            (delta_gZ2oMZ2 + 2. * delta_sbar) * KD(p, r)) *
                                            KD(s, t)) / 6.) * VeL(r, j) * VdR(t, l);
                                    CedVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CldR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CldI", p, r, s, t)) +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r)) / 4. -
                                            (gZ2oMZ2 * sbar2 * (-(v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r))) - v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            (-delta_gZ2oMZ2 - 2. * delta_sbar + 2. * delta_gZ2oMZ2 * sbar2 +
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t)) / 6.) * VeL(r, j) * VdR(t, l);
                                    CdeVLR.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VeR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CqeR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CqeI", p, r, s, t)) +
                                            (gZ2oMZ2 * ((-3. + 2. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHeR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, t)) * KD(p, r) +
                                            2. * sbar2 * (3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + 3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (3. * delta_gZ2oMZ2 + 6. * delta_sbar - 2. * delta_gZ2oMZ2 * sbar2 -
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            12.) * VdL(r, j) * VeR(t, l);
                                    CddV1LR.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd1I", p, r, s, t)) +
                                            (gZ2oMZ2 * (3. * (-3. + 2. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r) +
                                            2. * sbar2 * (3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + 3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (3. * delta_gZ2oMZ2 + 6. * delta_sbar - 2. * delta_gZ2oMZ2 * sbar2 -
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            36.) * VdL(r, j) * VdR(t, l);
                                    CddV8LR.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd8I", p, r, s, t))) * VdL(r, j) * VdR(t, l);
                                    CedSRL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CledqR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CledqI", p, r, s, t))) * VeR(r, j) * VdL(t, l);


                                }

                }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++) {

                    //Semileptonic LL operators
                    CnuuVLL.at(i).at(j).at(k).at(l) = 0.;
                    CeuVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic RR operators
                    CeuVRR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LLRR operators
                    CnuuVLR.at(i).at(j).at(k).at(l) = 0.;
                    CeuVLR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LLRR operators
                    CduV1LR.at(i).at(j).at(k).at(l) = 0.;
                    CduV8LR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRRL operators
                    CeuSRL.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRLR operators
                    CeuSRR.at(i).at(j).at(k).at(l) = 0.;
                    CeuTRR.at(i).at(j).at(k).at(l) = 0.;

                    for (int p = 0; p < 3; p++)
                        for (int r = 0; r < 3; r++)
                            for (int s = 0; s < 3; s++)
                                for (int t = 0; t < 3; t++) {
                                    CnuuVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", p, r, s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", p, r, s, t)) -
                                            (gZ2oMZ2 * (-3. + 4. * sbar2)*(
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) -
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) -
                                            delta_gZ2oMZ2 * KD(p, r)) * KD(s, t)) / 12. +
                                            (gZ2oMZ2 * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 12.) * VeL(r, j) * VuL(t, l);
                                    CeuVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", p, r, s, t)) -
                                            (mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", p, r, s, t)) +
                                            (gZ2oMZ2 * (-3. + 4. * sbar2)*(-0.5 * (v2 * (
                                            (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)))) +
                                            delta_gZ2oMZ2 * (-0.5 + sbar2) * KD(p, r) +
                                            2. * delta_sbar * sbar2 * KD(p, r)) * KD(s, t)) / 6. +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2) * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 12.) * VeL(r, j) * VuL(t, l);
                                    CeuVRR.at(i).at(j).at(k).at(l) += VeR.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CeuR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeuI", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 * (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) * KD(p, r) +
                                            2. * (-(v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", p, r))) +
                                            2. * (delta_gZ2oMZ2 + 4. * delta_sbar) * sbar2 *
                                            KD(p, r)) * KD(s, t))) / 6.) * VeR(r, j) * VuR(t, l);
                                    CnuuVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CluR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CluI", p, r, s, t)) +
                                            (gZ2oMZ2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) * KD(p, r)) / 4. +
                                            (gZ2oMZ2 * sbar2 * (-(v2 * (mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r))) +
                                            v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            (delta_gZ2oMZ2 + 2. * delta_sbar) * KD(p, r)) *
                                            KD(s, t)) / 3.) * VeL(r, j) * VuR(t, l);
                                    CeuVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CluR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CluI", p, r, s, t)) +
                                            (2. * gZ2oMZ2 * sbar2 * (-0.5 * (v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)))) +
                                            delta_gZ2oMZ2 * (-0.5 + sbar2) * KD(p, r) +
                                            2. * delta_sbar * sbar2 * KD(p, r)) * KD(s, t)) / 3. +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2) * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 12.) * VeL(r, j) * VuR(t, l);
                                    CduV1LR.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqu1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu1I", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 * (-3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - 3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 + 2. * delta_gZ2oMZ2 * sbar2 + 4. * delta_sbar * sbar2) *
                                            KD(p, r)) * KD(s, t)) / 9. +
                                            (gZ2oMZ2 * (-3. + 2. * sbar2) * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 36.) * VdL(r, j) * VuR(t, l);
                                    CduV8LR.at(i).at(j).at(k).at(l) += VdL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqu8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu8I", p, r, s, t))) * VdL(r, j) * VuR(t, l);
                                    CeuSRL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CluR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CluI", p, r, s, t)) +
                                            (2. * gZ2oMZ2 * sbar2 * (-0.5 * (v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", p, r)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)))) +
                                            delta_gZ2oMZ2 * (-0.5 + sbar2) * KD(p, r) +
                                            2. * delta_sbar * sbar2 * KD(p, r)) * KD(s, t)) / 3. +
                                            (gZ2oMZ2 * (-1. + 2. * sbar2) * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 12.) * VeR(r, j) * VuL(t, l);
                                    CeuSRR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuL.hconjugate()(k, s) *
                                            (-(mySMEFT.getSMEFTCoeffEW("Clequ1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ1I", p, r, s, t))) * VeR(r, j) * VuR(t, l);
                                    CeuTRR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VuL.hconjugate()(k, s) *
                                            (-(mySMEFT.getSMEFTCoeffEW("Clequ3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ3I", p, r, s, t))) * VeR(r, j) * VuR(t, l);
                                }

                }

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++) {

                    //Nonleptonic LL operators
                    CudV1LL.at(i).at(j).at(k).at(l) = 0.;
                    CudV8LL.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic RR operators
                    CudV1RR.at(i).at(j).at(k).at(l) = 0.;
                    CudV8RR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LLRR operators
                    CueVLR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LLRR operators
                    CudV1LR.at(i).at(j).at(k).at(l) = 0.;
                    CudV8LR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LRLR operators
                    CudS1RR.at(i).at(j).at(k).at(l) = 0.;
                    CudS8RR.at(i).at(j).at(k).at(l) = 0.;

                    for (int p = 0; p < 3; p++)
                        for (int r = 0; r < 3; r++)
                            for (int s = 0; s < 3; s++)
                                for (int t = 0; t < 3; t++) {
                                    CudV1LL.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((6. * (mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) +
                                            6. * (mySMEFT.getSMEFTCoeffEW("Cqq1R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", s, t, p, r)) -
                                            6. * (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t)) +
                                            4. * ((mySMEFT.getSMEFTCoeffEW("Cqq3R", p, t, s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, t, s, r)) +
                                            (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, r, p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, r, p, t))) -
                                            6. * (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, t, p, r)) -
                                            g22oMW2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", r, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, s)) * KD(p, t) +
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, t)) +
                                            delta_g22oMW2 * KD(p, t)) * KD(r, s)) +
                                            (gZ2oMZ2 * (-3. + 2. * sbar2)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 + 4. * delta_gZ2oMZ2 * sbar2 +
                                            8. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t)) / 6. + gZ2oMZ2 * (-3. + 4. * sbar2) * KD(p, r)*
                                            (-0.5 * (v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)))) +
                                            (2. * delta_sbar * sbar2 * KD(s, t)) / 3.)) / 6.) * VuL(r, j) * VdL(t, l);
                                    CudV8LL.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            (4. * (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, t, s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, t, s, r)) +
                                            4. * (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, r, p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, r, p, t)) -
                                            g22oMW2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", r, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, s)) * KD(p, t) +
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, t)) +
                                            delta_g22oMW2 * KD(p, t)) * KD(r, s))) * VuL(r, j) * VdL(t, l);
                                    CudV1RR.at(i).at(j).at(k).at(l) += VuR.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cud1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cud1I", p, r, s, t)) +
                                            (gZ2oMZ2 * sbar2 * (
                                            -6. * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r) +
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r)) +
                                            4. * (delta_gZ2oMZ2 + 4. * delta_sbar) * sbar2 *
                                            KD(p, r)) * KD(s, t))) / 18.) * VuR(r, j) * VdR(t, l);
                                    CudV8RR.at(i).at(j).at(k).at(l) += VuR.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cud8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cud8I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                    CueVLR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VeR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CqeR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CqeI", p, r, s, t)) +
                                            (gZ2oMZ2 * ((3. - 4. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHeR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", s, t)) * KD(p, r) +
                                            2. * sbar2 * (3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - 3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 - 6. * delta_sbar + 4. * delta_gZ2oMZ2 * sbar2 +
                                            16. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            12.) * VuL(r, j) * VeR(t, l);
                                    CudV1LR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd1I", p, r, s, t)) +
                                            (gZ2oMZ2 * (3. * (3. - 4. * sbar2) * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHdR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", s, t)) * KD(p, r) +
                                            2. * sbar2 * (3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - 3. * v2 *
                                            (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 - 6. * delta_sbar + 4. * delta_gZ2oMZ2 * sbar2 +
                                            16. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            36.) * VuL(r, j) * VdR(t, l);
                                    CudV8LR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd8I", p, r, s, t))) * VuL(r, j) * VdR(t, l);
                                    CudS1RR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cquqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd1I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                    CudS8RR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cquqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd8I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                }

                }


    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int l = 0; l < 2; l++) {

                    //Nonleptonic LL operators
                    CuuVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic RR operators
                    CuuVRR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LR operators
                    CuuV1LR.at(i).at(j).at(k).at(l) = 0.;
                    CuuV8LR.at(i).at(j).at(k).at(l) = 0.;
                    CudduV1LR.at(i).at(j).at(k).at(l) = 0.;
                    CudduV8LR.at(i).at(j).at(k).at(l) = 0.;

                    //Nonleptonic LRLR operators
                    CuuS1RR.at(i).at(j).at(k).at(l) = 0.;
                    CuuS8RR.at(i).at(j).at(k).at(l) = 0.;
                    CudduS1RR.at(i).at(j).at(k).at(l) = 0.;
                    CudduS8RR.at(i).at(j).at(k).at(l) = 0.;

                    for (int p = 0; p < 3; p++)
                        for (int r = 0; r < 3; r++)
                            for (int s = 0; s < 3; s++)
                                for (int t = 0; t < 3; t++) {
                                    CuuVLL.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VuL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) +
                                            (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t)) -
                                            (gZ2oMZ2 * (-3. + 4. * sbar2)*(
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", s, t)) * KD(p, r) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", s, t)) * KD(p, r) +
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 + 4. * delta_gZ2oMZ2 * sbar2 +
                                            16. * delta_sbar * sbar2) * KD(p, r)) * KD(s, t))) /
                                            72.) * VuL(r, j) * VuL(t, l);
                                    CuuVRR.at(i).at(j).at(k).at(l) += VuR.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CuuR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuuI", p, r, s, t)) -
                                            (gZ2oMZ2 * sbar2 * (
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) * KD(p, r) +
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r)) +
                                            4. * (delta_gZ2oMZ2 + 4. * delta_sbar) * sbar2 *
                                            KD(p, r)) * KD(s, t))) / 18.) * VuR(r, j) * VuR(t, l);
                                    CuuV1LR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqu1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu1I", p, r, s, t)) -
                                            (gZ2oMZ2 * sbar2 * (
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) -
                                            3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)) +
                                            (-3. * delta_gZ2oMZ2 + 4. * delta_gZ2oMZ2 * sbar2 + 8. * delta_sbar * sbar2) *
                                            KD(p, r)) * KD(s, t)) / 9. -
                                            (gZ2oMZ2 * (-3. + 4. * sbar2) * KD(p, r)*
                                            (3. * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", s, t)) +
                                            8. * delta_sbar * sbar2 * KD(s, t))) / 36.) * VuL(r, j) * VuR(t, l);
                                    CuuV8LR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VuR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Cqu8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu8I", p, r, s, t))) * VuL(r, j) * VuR(t, l);
                                    CudduV1LR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            (-0.25 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHudR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", t, s)) * KD(p, r))) *
                                            VdL(r, j) * VuR(t, l);
                                    CudduS1RR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            (-(mySMEFT.getSMEFTCoeffEW("Cquqd1R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd1I", s, t, p, r))) *
                                            VdR(r, j) * VuR(t, l);
                                    CudduS8RR.at(i).at(j).at(k).at(l) += VuL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            (-(mySMEFT.getSMEFTCoeffEW("Cquqd8R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd8I", s, t, p, r))) *
                                            VdR(r, j) * VuR(t, l);

                                }

                }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 2; l++) {

                    //Semileptonic LL operators
                    CnueduVLL.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LR operators
                    CnueduVLR.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRRL operators
                    CnueduSRL.at(i).at(j).at(k).at(l) = 0.;

                    //Semileptonic LRLR operators
                    CnueduSRR.at(i).at(j).at(k).at(l) = 0.;
                    CnueduTRR.at(i).at(j).at(k).at(l) = 0.;

                    for (int p = 0; p < 3; p++)
                        for (int r = 0; r < 3; r++)
                            for (int s = 0; s < 3; s++)
                                for (int t = 0; t < 3; t++) {
                                    CnueduVLL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            (2. * (mySMEFT.getSMEFTCoeffEW("Clq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", p, r, s, t)) -
                                            (g22oMW2 * (v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", t, s)) * KD(p, r) +
                                            (v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", p, r)) +
                                            delta_g22oMW2 * KD(p, r)) * KD(s, t))) / 2.) * VeL(r, j) * VuL(t, l);
                                    CnueduVLR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            (-0.25 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHudR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", t, s)) * KD(p, r)))
                                            * VeL(r, j) * VuR(t, l);
                                    CnueduSRL.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdR.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("CledqR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CledqI", p, r, s, t)))
                                            * VeR(r, j) * VuL(t, l);
                                    CnueduSRR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clequ1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ1I", p, r, s, t)))
                                            * VeR(r, j) * VuR(t, l);
                                    CnueduTRR.at(i).at(j).at(k).at(l) += VeL.hconjugate()(i, p) * VdL.hconjugate()(k, s) *
                                            ((mySMEFT.getSMEFTCoeffEW("Clequ3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ3I", p, r, s, t)))
                                            * VeR(r, j) * VuR(t, l);

                                }

                }




    StandardModelMatching::updateSMParameters();
}

NPSMEFTd6GeneralMatching::~NPSMEFTd6GeneralMatching() {
}

// Matching to the Delta F=2 Hamiltonian in the SUSY Basis, using eq. (199) of 2009.07276

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdk2() {

    vmck2.clear();
    vmck2 = StandardModelMatching::CMdk2();

    mck2.setMu(mySMEFT.getMuw());
    
    switch (mck2.getOrder()) {
        case LO:
            mck2.setCoeff(0, CddVLL.at(0).at(1).at(0).at(1), LO);
            mck2.setCoeff(1, CddS1RR.at(1).at(0).at(1).at(0).conjugate() - CddS8RR.at(1).at(0).at(1).at(0).conjugate()/6., LO);
            mck2.setCoeff(2, CddS8RR.at(1).at(0).at(1).at(0).conjugate()/2., LO);
            mck2.setCoeff(3, CddV1LR.at(0).at(1).at(0).at(1), LO);
            mck2.setCoeff(4, 2.*CddV1LR.at(0).at(1).at(0).at(1) - CddV8LR.at(0).at(1).at(0).at(1)/3., LO);
            break;
        default:
            std::stringstream out;
            out << mck2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented"); 
    }
    
    vmck2.push_back(mck2);

    switch (mck2.getOrder()) {
        case LO:
            mck2.setCoeff(0, CddVRR.at(0).at(1).at(0).at(1), LO);
            mck2.setCoeff(1, CddS1RR.at(0).at(1).at(0).at(1) - CddS8RR.at(0).at(1).at(0).at(1)/6., LO);
            mck2.setCoeff(2, CddS8RR.at(0).at(1).at(0).at(1)/2., LO);
            mck2.setCoeff(3, 0., LO);
            mck2.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mck2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented"); 
    }
    
    vmck2.push_back(mck2);

    return (vmck2);

}

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdd2() {

    vmcd2.clear();
    vmcd2 = StandardModelMatching::CMdd2();

    mcd2.setMu(mySMEFT.getMuw());
    
    switch (mcd2.getOrder()) {
        case LO:
            mcd2.setCoeff(0, CuuVLL.at(0).at(1).at(0).at(1), LO);
            mcd2.setCoeff(1, CuuS1RR.at(1).at(0).at(1).at(0).conjugate() - CuuS8RR.at(1).at(0).at(1).at(0).conjugate()/6., LO);
            mcd2.setCoeff(2, CuuS8RR.at(1).at(0).at(1).at(0).conjugate()/2., LO);
            mcd2.setCoeff(3, CuuV1LR.at(0).at(1).at(0).at(1), LO);
            mcd2.setCoeff(4, 2.*CuuV1LR.at(0).at(1).at(0).at(1) - CuuV8LR.at(0).at(1).at(0).at(1)/3., LO);
            break;
        default:
            std::stringstream out;
            out << mcd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdd2(): order " + out.str() + "not implemented"); 
    }
    
    vmcd2.push_back(mcd2);

    switch (mcd2.getOrder()) {
        case LO:
            mcd2.setCoeff(0, CuuVRR.at(0).at(1).at(0).at(1), LO);
            mcd2.setCoeff(1, CuuS1RR.at(0).at(1).at(0).at(1) - CuuS8RR.at(0).at(1).at(0).at(1)/6., LO);
            mcd2.setCoeff(2, CuuS8RR.at(0).at(1).at(0).at(1)/2., LO);
            mcd2.setCoeff(3, 0., LO);
            mcd2.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdd2(): order " + out.str() + "not implemented"); 
    }
    
    vmcd2.push_back(mcd2);

    return (vmcd2);

}

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdbd2() {

    vmcdb.clear();
    vmcdb = StandardModelMatching::CMdbd2();

    mcbd.setMu(mySMEFT.getMuw());
    
    switch (mcbd.getOrder()) {
        case LO:
            mcbd.setCoeff(0, CddVLL.at(0).at(2).at(0).at(2), LO);
            mcbd.setCoeff(1, CddS1RR.at(2).at(0).at(2).at(0).conjugate() - CddS8RR.at(2).at(0).at(2).at(0).conjugate()/6., LO);
            mcbd.setCoeff(2, CddS8RR.at(2).at(0).at(2).at(0).conjugate()/2., LO);
            mcbd.setCoeff(3, CddV1LR.at(0).at(2).at(0).at(2), LO);
            mcbd.setCoeff(4, 2.*CddV1LR.at(0).at(2).at(0).at(2) - CddV8LR.at(0).at(2).at(0).at(2)/3., LO);
            break;
        default:
            std::stringstream out;
            out << mcbd.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented"); 
    }
    
    vmcdb.push_back(mcbd);

    switch (mcbd.getOrder()) {
        case LO:
            mcbd.setCoeff(0, CddVRR.at(0).at(2).at(0).at(2), LO);
            mcbd.setCoeff(1, CddS1RR.at(0).at(2).at(0).at(2) - CddS8RR.at(0).at(2).at(0).at(2)/6., LO);
            mcbd.setCoeff(2, CddS8RR.at(0).at(2).at(0).at(2)/2., LO);
            mcbd.setCoeff(3, 0., LO);
            mcbd.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbd.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented"); 
    }
    
    vmcdb.push_back(mcbd);

    return (vmcdb);

}

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdbs2() {

    vmcds.clear();
    vmcds = StandardModelMatching::CMdbs2();

    mcbs.setMu(mySMEFT.getMuw());
    
    switch (mcbs.getOrder()) {
        case LO:
            mcbs.setCoeff(0, CddVLL.at(1).at(2).at(1).at(2), LO);
            mcbs.setCoeff(1, CddS1RR.at(2).at(1).at(2).at(1).conjugate() - CddS8RR.at(2).at(1).at(2).at(1).conjugate()/6., LO);
            mcbs.setCoeff(2, CddS8RR.at(2).at(1).at(2).at(1).conjugate()/2., LO);
            mcbs.setCoeff(3, CddV1LR.at(1).at(2).at(1).at(2), LO);
            mcbs.setCoeff(4, 2.*CddV1LR.at(1).at(2).at(1).at(2) - CddV8LR.at(1).at(2).at(1).at(2)/3., LO);
            break;
        default:
            std::stringstream out;
            out << mcbs.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented"); 
    }
    
    vmcds.push_back(mcbs);

    switch (mcbs.getOrder()) {
        case LO:
            mcbs.setCoeff(0, CddVRR.at(1).at(2).at(1).at(2), LO);
            mcbs.setCoeff(1, CddS1RR.at(1).at(2).at(1).at(2) - CddS8RR.at(1).at(2).at(1).at(2)/6., LO);
            mcbs.setCoeff(2, CddS8RR.at(1).at(2).at(1).at(2)/2., LO);
            mcbs.setCoeff(3, 0., LO);
            mcbs.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbs.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented"); 
    }
    
    vmcds.push_back(mcbs);

    return (vmcds);

}

/*******************************************************************************
 * Wilson coefficients matching, LEFT basis [1709.04486] ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij             *
 * ****************************************************************************/
std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdiujleptonknu(int i, int j, int k) {

    vmculeptonnu.clear();
    vmculeptonnu = StandardModelMatching::CMdiujleptonknu(i,j,k);

    mculeptonnu.setMu(mySMEFT.getMuw());

    switch (mculeptonnu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mculeptonnu.setCoeff(0, CnueduVLL.at(k).at(k).at(i).at(j), LO);
            mculeptonnu.setCoeff(1, CnueduVLR.at(k).at(k).at(i).at(j), LO);
            mculeptonnu.setCoeff(2, CnueduSRR.at(k).at(k).at(i).at(j), LO);
            mculeptonnu.setCoeff(3, CnueduSRL.at(k).at(k).at(i).at(j), LO);
            mculeptonnu.setCoeff(4, CnueduTRR.at(k).at(k).at(i).at(j), LO);
            break;
        default:
            std::stringstream out;
            out << mculeptonnu.getOrder();
            throw std::runtime_error("StandardModelMatching::CMuleptonnu(): order " + out.str() + "not implemented");
    }

    vmculeptonnu.push_back(mculeptonnu);
    return (vmculeptonnu);

}
