/*
 * Copyright (C) 2018. HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6GeneralMatching.h"
#include "NPSMEFTd6General.h"
#include <stdexcept>

#define NOLEPTONFLAVOURVIOLATION

NPSMEFTd6GeneralMatching::NPSMEFTd6GeneralMatching(const NPSMEFTd6General &NPSMEFTd6General_i) : StandardModelMatching(NPSMEFTd6General_i),
                                                                                                 mySMEFT(NPSMEFTd6General_i), 
                                                                                                 VuL(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VuLd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VuR(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VuRd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VdL(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VdLd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VdR(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VdRd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VeL(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VeLd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VeR(gslpp::matrix<complex>::Id(3)), 
                                                                                                 VeRd(gslpp::matrix<complex>::Id(3)), 
                                                                                                 MU(3, 0.), 
                                                                                                 MD(3, 0.),
                                                                                                 mcd2(5, NDR, NLO), 
                                                                                                 mcd1(10, NDR, NLO), 
                                                                                                 mcbd(5, NDR, NLO), 
                                                                                                 mcbs(5, NDR, NLO), 
                                                                                                 mck2(5, NDR, NLO), 
                                                                                                 mculeptonnu(5, NDR, LO), 
                                                                                                 mckpnn(2, NDR, NLO, NLO_QED11),
                                                                                                 mcbsg(8, NDR, NNLO),
                                                                                                 mcprimebsg(8, NDR, NNLO)
                                                                                                 
{
}

void NPSMEFTd6GeneralMatching::updateLEFTGeneralParameters()
{

        // Dimension 6 operators with no flavour index are assigned directly here

        LambdaNP2 = mySMEFT.getLambda_NP() * mySMEFT.getLambda_NP();
        v = mySMEFT.v(); // This is vtilde in Angelica's notation
        v2 = v * v;      // This is vtilde squared

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

        // For operators with quark indices we need to switch to the mass eigenstate basis; leptons are already in the mass eigenstate basis since we do not have any lepton flavour violation

        VuL = mySMEFT.getVuL();
        VdL = mySMEFT.getVdL();
        VeL = mySMEFT.getVeL();
        VuLd = mySMEFT.getVuLd();
        VdLd = mySMEFT.getVdLd();
        VeLd = mySMEFT.getVeLd();
        VuR = mySMEFT.getVuR();
        VdR = mySMEFT.getVdR();
        VeR = mySMEFT.getVeR();
        VuRd = mySMEFT.getVuRd();
        VdRd = mySMEFT.getVdRd();
        VeRd = mySMEFT.getVeRd();
        
        // to implement Manohar's matching formulae we define the couplings
        // in his notation. Namely, in the formulae below, the barred quantities are
        // tree level in the theory scheme.

        double cbar = mySMEFT.getXWZ_tree();
        double sbar = -mySMEFT.getXBZ_tree();
        double sbar2 = sbar * sbar;
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
        // indeed a piece was missing, I add it here
        delta_gZbar += -(g1bar * g1bar * g1bar * g1bar + g2bar * g2bar * g2bar * g2bar) / (2. * (g1bar * g1bar + g2bar * g2bar) * g1bar * g2bar) * v2 * mySMEFT.getSMEFTCoeffEW("CHWB");
        double gZ2oMZ2 = gZbar / mySMEFT.getMz();
        gZ2oMZ2 *= gZ2oMZ2;
        double delta_gZ2oMZ2 = 2. * delta_gZbar - delta_MZ2;
        double g22oMW2 = 4. / v2;
        double delta_g22oMW2 = -2. * delta_vT;
        //new terms from now on
        //Z boson
        double gZbar2oMZ3 = gZ2oMZ2 / mySMEFT.getMz();
        double gZbar2oMZ4 = gZbar2oMZ3 / mySMEFT.getMz();
        double deltagZbar2oMZ4 = 2. * delta_gZbar - 2. * delta_MZ2;
        //W boson
        double g2bar2oMW3 = g22oMW2 / mySMEFT.getMw();
        double g2bar2oMW4 = g2bar2oMW3 / mySMEFT.getMw();
        double deltag2bar2oMW4 = -4. * delta_vT;
        //h boson
        double lambda = mySMEFT.getLambdaH_tree() * 2.;
        double oneoMh2 = 1. / (lambda * v2);
        double deltaoneoMh2 = (- 2. * mySMEFT.getSMEFTCoeffEW("CHbox") + 0.5 * mySMEFT.getSMEFTCoeffEW("CHD") + 3 * mySMEFT.getSMEFTCoeffEW("CH") / lambda) * v2;
        std::array<double, 3> Me = mySMEFT.getMe_LEW();
        std::array<double, 3> Mu = mySMEFT.getMu_LEW();
        std::array<double, 3> Md = mySMEFT.getMd_LEW();

        // std::cout << "CKM from rotated UfA = " << (VuL.hconjugate()) * VdL << std::endl;

        // std::cout << "has the diagonalization worked? " << VuR.hconjugate()*MU*VuL << std::endl;
        // std::cout << "has the diagonalization worked? " << VdR.hconjugate()*MD*VdL << std::endl;

        // match and rotate following Manohar. This is performed AT LINEAR ORDER for the moment

        // fill all coefficients with zeroes first
        gslpp::matrix<complex> VCKM = mySMEFT.getCKM().getCKM();
        gslpp::matrix<complex> VCKMd = VCKM.hconjugate();

        Ceg = zero33;
        Ceg = zero33;
        if (Ceg.at(0).at(1) != 0. || Ceg.at(0).at(2) != 0. || Ceg.at(1).at(0) != 0. || Ceg.at(1).at(2) != 0. || Ceg.at(2).at(0) != 0. || Ceg.at(2).at(1) != 0.)
                throw("Compiler is not putting to zero correctly the 2-d arrays of Wilson coefficients");
        Cdg = zero33;
        CdG = zero33;
        Cug = zero22;
        CuG = zero22;

        CnunuVLL = zero3333;
        CeeVLL = zero3333;
        CnueVLL = zero3333;
        CnuuVLL = zero3322;
        CnudVLL = zero3333;
        CeuVLL = zero3322;
        CedVLL = zero3333;
        CnueduVLL = zero3332;
        CuuVLL = zero2222;
        CddVLL = zero3333;
        CudV1LL = zero2233;
        CudV8LL = zero2233;
        CeeVRR = zero3333;
        CeuVRR = zero3322;
        CedVRR = zero3333;
        CuuVRR = zero2222;
        CddVRR = zero3333;
        CudV1RR = zero2233;
        CudV8RR = zero2233;
        CnueVLR = zero3333;
        CeeVLR = zero3333;
        CnuuVLR = zero3322;
        CnudVLR = zero3333;
        CeuVLR = zero3322;
        CueVLR = zero2233;
        CedVLR = zero3333;
        CdeVLR = zero3333;
        CnueduVLR = zero3332;
        CuuV1LR = zero2222;
        CuuV8LR = zero2222;
        CudV1LR = zero2233;
        CudV8LR = zero2233;
        CduV1LR = zero3322;
        CduV8LR = zero3322;
        CddV1LR = zero3333;
        CddV8LR = zero3333;
        CudduV1LR = zero2332;
        CudduV8LR = zero2332;
        CeuSRL = zero3322;
        CedSRL = zero3333;
        CnueduSRL = zero3332;
        CeeSRR = zero3333;
        CeuSRR = zero3322;
        CeuTRR = zero3322;
        CedSRR = zero3333;
        CedTRR = zero3333;
        CnueduSRR = zero3332;
        CnueduTRR = zero3332;
        CuuS1RR = zero2222;
        CuuS8RR = zero2222;
        CudS1RR = zero2233;
        CudS8RR = zero2233;
        CddS1RR = zero3333;
        CddS8RR = zero3333;
        CudduS1RR = zero2332;
        CudduS8RR = zero2332;

#ifdef NOLEPTONFLAVOURVIOLATION

        // matching of operators with two external indices and zero internal indices

        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                        Ceg.at(i).at(j) += vTosq2 * (-(mySMEFT.getSMEFTCoeffEW("CeWR", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", i, j)) * sbar + (mySMEFT.getSMEFTCoeffEW("CeBR", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", i, j)) * cbar);

                        CnunuVLL.at(i).at(i).at(j).at(j) += -0.0625 * (delta_gZ2oMZ2 * gZ2oMZ2);
                        CnunuVLL.at(i).at(j).at(j).at(i) += -0.0625 * (delta_gZ2oMZ2 * gZ2oMZ2);
                        CeeVLL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (-1 + 2 * sbar2) * (delta_gZ2oMZ2 - 2 * delta_gZ2oMZ2 * sbar2 - 8 * delta_sbar * sbar2)) / 16.;
                        CeeVLL.at(i).at(j).at(j).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * (delta_gZ2oMZ2 - 2 * delta_gZ2oMZ2 * sbar2 - 8 * delta_sbar * sbar2)) / 16.;
                        CnueVLL.at(i).at(i).at(j).at(j) += -0.25 * (gZ2oMZ2 * (4 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-1 + 2 * sbar2)));
                        CnueVLL.at(i).at(j).at(j).at(i) += -0.5 * (delta_g22oMW2 * g22oMW2);
                        CnudVLL.at(i).at(i).at(j).at(j) += -0.08333333333333333 * (gZ2oMZ2 * (4 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 2 * sbar2)));
                        CedVLL.at(i).at(i).at(j).at(j) += -0.08333333333333333 * (gZ2oMZ2 * (16 * delta_sbar * (-1 + sbar2) * sbar2 + delta_gZ2oMZ2 * (3 - 8 * sbar2 + 4 * sbar2 * sbar2)));
                        CddVLL.at(i).at(i).at(j).at(j) += -0.013888888888888888 * (gZ2oMZ2 * (-3 + 2 * sbar2) * (8 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 2 * sbar2)));
                        CeeVRR.at(i).at(i).at(j).at(j) += -0.25 * ((delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2);
                        CeeVRR.at(i).at(j).at(j).at(i) += -0.25 * ((delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2);
                        CedVRR.at(i).at(i).at(j).at(j) += -0.3333333333333333 * ((delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2);
                        CddVRR.at(i).at(i).at(j).at(j) += -0.05555555555555555 * ((delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2);
                        CnueVLR.at(i).at(i).at(j).at(j) += -0.5 * ((delta_gZ2oMZ2 + 2 * delta_sbar) * gZ2oMZ2 * sbar2);
                        CeeVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (delta_gZ2oMZ2 + delta_sbar * (2 - 8 * sbar2) - 2 * delta_gZ2oMZ2 * sbar2)) / 2.;
                        CnudVLR.at(i).at(i).at(j).at(j) += -0.16666666666666666 * ((delta_gZ2oMZ2 + 2 * delta_sbar) * gZ2oMZ2 * sbar2);
                        CedVLR.at(i).at(i).at(j).at(j) += -0.16666666666666666 * (gZ2oMZ2 * sbar2 * (-delta_gZ2oMZ2 - 2 * delta_sbar + 2 * delta_gZ2oMZ2 * sbar2 + 8 * delta_sbar * sbar2));
                        CdeVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (delta_sbar * (6 - 8 * sbar2) + delta_gZ2oMZ2 * (3 - 2 * sbar2)) * sbar2) / 6.;
                        CddV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (delta_sbar * (6 - 8 * sbar2) + delta_gZ2oMZ2 * (3 - 2 * sbar2)) * sbar2) / 18.;
                        if (mySMEFT.FlagNewTerms) {
                            CnueVLR.at(i).at(j).at(j).at(i) += -0.25*((1 + deltag2bar2oMW4)*g2bar2oMW4*Me[i]*Me[j]);
                            CeeVLR.at(i).at(j).at(j).at(i) += -0.125*(gZbar2oMZ4*(-1 + 2*sbar2)*(-1 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2))*Me[i]*Me[j]);
                            CeeVLR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.;
                            CeeVLR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.;
                            CeeVLR.at(i).at(j).at(j).at(i) += -0.5*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Me[j]);
                            CddV1LR.at(i).at(j).at(j).at(i) += -0.004629629629629629*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                            CddV1LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/108.;
                            CddV1LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/108.;
                            CddV1LR.at(i).at(j).at(j).at(i) += -0.018518518518518517*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                            CddV8LR.at(i).at(j).at(j).at(i) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                            CddV8LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.;
                            CddV8LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.;
                            CddV8LR.at(i).at(j).at(j).at(i) += -0.1111111111111111*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                            CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 - 8*(1 + 2*delta_sbar)*sbar2 + 4*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 8*sbar2 + 4*pow(sbar2,2)))*Md[j]*Me[i])/12.;
                            CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(1 + deltagZbar2oMZ4 + delta_sbar*(2 - 8*sbar2) - 2*sbar2 - 2*deltagZbar2oMZ4*sbar2)*Md[j]*Me[i])/6.;
                            CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 + delta_sbar*(6 - 8*sbar2) + deltagZbar2oMZ4*(3 - 2*sbar2) - 2*sbar2)*sbar2*Md[j]*Me[i])/6.;
                            CedSRL.at(i).at(i).at(j).at(j) += ((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Me[i])/3.;
                            CeeSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.;
                            CeeSRR.at(i).at(i).at(j).at(j) += -0.125*(gZbar2oMZ4*(-1 + 2*sbar2)*(-1 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2))*Me[i]*Me[j]);
                            CeeSRR.at(i).at(i).at(j).at(j) += -0.5*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Me[j]);
                            CeeSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.;
                            CedSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Md[j]*Me[i])/6.;
                            CedSRR.at(i).at(i).at(j).at(j) += -0.08333333333333333*(gZbar2oMZ4*(3 - 8*(1 + 2*delta_sbar)*sbar2 + 4*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 8*sbar2 + 4*pow(sbar2,2)))*Md[j]*Me[i]);
                            CedSRR.at(i).at(i).at(j).at(j) += -0.3333333333333333*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Me[i]);
                            CedSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[j]*Me[i])/6.;
                            CddS1RR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.;
                            CddS1RR.at(i).at(i).at(j).at(j) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                            CddS1RR.at(i).at(i).at(j).at(j) += -0.1111111111111111*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                            CddS1RR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.;
                        }
                        for (int k = 0; k < 3; k++)
                        {
                                CnunuVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CnunuVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CnunuVLL.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CnunuVLL.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", k, j)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", k, j)))) / 16.;
                                CeeVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CeeVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CeeVLL.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 16.;
                                CeeVLL.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", k, j)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", k, j)))) / 16.;
                                CnueVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 4.;
                                CnueVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 4.;
                                CnueVLL.at(j).at(i).at(i).at(k) += -0.5 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)));
                                CnueVLL.at(i).at(j).at(k).at(i) += -0.5 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)));
                                CnudVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 12.;
                                CedVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 12.;
                                CeeVRR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 8.;
                                CeeVRR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 8.);
                                CeeVRR.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 8.;
                                CeeVRR.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", k, j))) / 8.;
                                CedVRR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 6.;
                                CnueVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 4.;
                                CnueVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 2.;
                                CeeVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 4.;
                                CeeVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 2.;
                                CnudVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 6.;
                                CedVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 6.;
                                CdeVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 12.;
                                if (mySMEFT.FlagNewTerms) {
                                    CeeVLL.at(i).at(i).at(j).at(k) += (gZbar2oMZ3*(0.5 - sbar2)*v2*((sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/(4.*sqrt(2));
                                    CeeVLL.at(j).at(k).at(i).at(i) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CeeVLL.at(j).at(i).at(i).at(k) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CeeVLL.at(i).at(j).at(k).at(i) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CnueVLL.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*v2*((sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CnueVLL.at(j).at(i).at(i).at(k) += (g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[k])/sqrt(2);
                                    CnueVLL.at(i).at(j).at(k).at(i) += -((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[k])/sqrt(2));
                                    CedVLL.at(j).at(k).at(i).at(i) += -0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CeeVRR.at(i).at(i).at(j).at(k) += -0.25*(gZbar2oMZ3*sbar2*v2*((sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CeeVRR.at(j).at(k).at(i).at(i) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CeeVRR.at(j).at(i).at(i).at(k) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CeeVRR.at(i).at(j).at(k).at(i) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CedVRR.at(j).at(k).at(i).at(i) += -0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CnueVLR.at(i).at(j).at(k).at(i) += -0.25*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k))*Me[i]*Me[k]);
                                    CnueVLR.at(j).at(i).at(i).at(k) += -0.25*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k))*Me[i]*Me[k]);
                                    CnueVLR.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CeeVLR.at(i).at(j).at(k).at(i) += (gZbar2oMZ4*v2*Me[i]*((mySMEFT.getSMEFTCoeffEW("CHeR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", k, j))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", k, j)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", k, j)))*Me[k]))/8.;
                                    CeeVLR.at(j).at(i).at(i).at(k) += (gZbar2oMZ4*v2*Me[i]*((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k]))/8.;
                                    CeeVLR.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CeeVLR.at(j).at(k).at(i).at(i) += -((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2));
                                    CedVLR.at(j).at(k).at(i).at(i) += -0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/sqrt(2);
                                    CdeVLR.at(i).at(i).at(j).at(k) += -0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*((sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/sqrt(2);
                                    CedSRL.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Md[i]*(-((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j]) + ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k]))/4.;
                                    CeeSRR.at(i).at(i).at(j).at(k) += (gZbar2oMZ4*v2*Me[i]*((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k]))/8.;
                                    CeeSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Me[i]*((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k]))/8.;
                                    CedSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Md[i]*((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k]))/4.;
                                }
                                for (int p = 0; p < 3; p++)
                                        for (int r = 0; r < 3; r++)
                                        {
                                                CnudVLL.at(i).at(i).at(j).at(k) += VdLd(j, p) * ((gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 4.) * VdL(r, k);
                                                CedVLL.at(i).at(i).at(j).at(k) += VdLd(j, p) * ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 4.) * VdL(r, k);
                                                CddVLL.at(i).at(i).at(j).at(k) += VdLd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 24.) * VdL(r, k);
                                                CddVLL.at(j).at(k).at(i).at(i) += VdLd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 24.) * VdL(r, k);
                                                CedVRR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 2.) * VdR(r, k);
                                                CddVRR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 12.) * VdR(r, k);
                                                CddVRR.at(j).at(k).at(i).at(i) += VdRd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 12.) * VdR(r, k);
                                                CnudVLR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 4.) * VdR(r, k);
                                                CnudVLR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2)) * VdR(r,k);
                                                CedVLR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 4.) * VdR(r, k);
                                                CdeVLR.at(j).at(k).at(i).at(i) += VdLd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 2.) * VdL(r, k);
                                                CddV1LR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 12.) * VdR(r, k);
                                                CddV1LR.at(j).at(k).at(i).at(i) += VdLd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 6.) * VdL(r, k);
                                                if (mySMEFT.FlagNewTerms) {
                                                    CnudVLL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2)) * VdL(r,k);
                                                    CnudVLL.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2)) * VdR(r,k);
                                                    CedVLL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2))) * VdL(r,k);
                                                    CedVLL.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2))) * VdR(r,k);
                                                    CddVLL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2)) * VdL(r,k);
                                                    CddVLL.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2)) * VdR(r,k);
                                                    CddVLL.at(j).at(k).at(i).at(i) += VdRd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2)) * VdL(r,k);
                                                    CddVLL.at(j).at(k).at(i).at(i) += VdLd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2)) * VdR(r,k);
                                                    CedVRR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2))) * VdL(r,k);
                                                    CedVRR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2))) * VdR(r,k);
                                                    CddVRR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2)) * VdL(r,k);
                                                    CddVRR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2)) * VdR(r,k);
                                                    CddVRR.at(j).at(k).at(i).at(i) += VdRd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2)) * VdL(r,k);
                                                    CddVRR.at(j).at(k).at(i).at(i) += VdLd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2)) * VdR(r,k);
                                                    CnudVLR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2)) * VdL(r,k);
                                                    CnudVLR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2)) * VdR(r,k);
                                                    CedVLR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2))) * VdL(r,k);
                                                    CedVLR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2))) * VdR(r,k);
                                                    CdeVLR.at(j).at(k).at(i).at(i) += VdRd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2))) * VdL(r,k);
                                                    CdeVLR.at(j).at(k).at(i).at(i) += VdLd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2))) * VdR(r,k);
                                                    CddV1LR.at(i).at(j).at(k).at(i) += VdLd(k,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,j);
                                                    CddV1LR.at(i).at(j).at(k).at(i) += VdRd(k,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/24.) * VdR(r,j);
                                                    CddV1LR.at(j).at(i).at(i).at(k) += VdLd(j,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,k);
                                                    CddV1LR.at(j).at(i).at(i).at(k) += VdRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/24.) * VdR(r,k);
                                                    CddV1LR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/sqrt(2)) * VdL(r,k);
                                                    CddV1LR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/sqrt(2)) * VdR(r,k);
                                                    CddV1LR.at(j).at(k).at(i).at(i) += VdRd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/sqrt(2)) * VdL(r,k);
                                                    CddV1LR.at(j).at(k).at(i).at(i) += VdLd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/sqrt(2)) * VdR(r,k);
                                                    CddV8LR.at(i).at(j).at(k).at(i) += VdLd(k,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,j);
                                                    CddV8LR.at(i).at(j).at(k).at(i) += VdRd(k,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/4.) * VdR(r,j);
                                                    CddV8LR.at(j).at(i).at(i).at(k) += VdLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,k);
                                                    CddV8LR.at(j).at(i).at(i).at(k) += VdRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/4.) * VdR(r,k);
                                                    CedSRL.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[j]*Me[i])/4.) * VdL(r,k);
                                                    CedSRL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[k]*Me[i])) * VdR(r,k);
                                                    CedSRR.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[j]*Me[i])/4.) * VdR(r,k);
                                                    CedSRR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[k]*Me[i])) * VdL(r,k);
                                                    CddS1RR.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/4.) * VdR(r,k);
                                                    CddS1RR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,k);
                                                    CddS1RR.at(j).at(k).at(i).at(i) += VdLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Md[k])) * VdL(r,k);
                                                    CddS1RR.at(j).at(k).at(i).at(i) += VdRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[i]*Md[j])/4.) * VdR(r,k);
                                                }
                                        }

                                for (int l = 0; l < 3; l++)
                                {
                                        CnunuVLL.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CllR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", i, j, k, l));
                                        CeeVLL.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CllR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", i, j, k, l));
                                        CnueVLL.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CllR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", i, j, k, l)) + (mySMEFT.getSMEFTCoeffEW("CllR", k, l, i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CllI", k, l, i, j));
                                        CeeVRR.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CeeR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeeI", i, j, k, l));
                                        CnueVLR.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CleR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CleI", i, j, k, l));
                                        CeeVLR.at(i).at(j).at(k).at(l) += (mySMEFT.getSMEFTCoeffEW("CleR", i, j, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CleI", i, j, k, l));
                                        if (mySMEFT.FlagNewTerms) {
                                            CeeVLR.at(i).at(j).at(k).at(l) += (oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", l, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", l, i))*(mySMEFT.getSMEFTCoeffEW("YeR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", k, j)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", k, j)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YeR", k, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", k, j)))*(mySMEFT.getSMEFTCoeffEW("YeR", l, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", l, i))))/4.;
                                            CedSRL.at(i).at(j).at(k).at(l) += VCKMd(k,l) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHl3R", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", i, j))*Md[k]*Me[j])/2.);
                                            CeeSRR.at(i).at(j).at(k).at(l) += -0.25*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", l, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", l, k))*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", j, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", j, i)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)))*(mySMEFT.getSMEFTCoeffEW("YeR", l, k) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", l, k))));
                                        }
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                {
                                                        CnudVLL.at(i).at(j).at(k).at(l) += VdLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clq1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", i, j, p, r)) - (mySMEFT.getSMEFTCoeffEW("Clq3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", i, j, p, r))) * VdL(r, l);
                                                        CedVLL.at(i).at(j).at(k).at(l) += VdLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clq1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", i, j, p, r)) + (mySMEFT.getSMEFTCoeffEW("Clq3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", i, j, p, r))) * VdL(r, l);
                                                        CedVRR.at(i).at(j).at(k).at(l) += VdRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CedR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CedI", i, j, p, r))) * VdR(r, l);
                                                        CnudVLR.at(i).at(j).at(k).at(l) += VdRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CldR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CldI", i, j, p, r))) * VdR(r, l);
                                                        CedVLR.at(i).at(j).at(k).at(l) += VdRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CldR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CldI", i, j, p, r))) * VdR(r, l);
                                                        CdeVLR.at(i).at(j).at(k).at(l) += VdLd(i, p) * ((mySMEFT.getSMEFTCoeffEW("CqeR", p, r, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CqeI", p, r, k, l))) * VdL(r, j);
                                                        CedSRL.at(i).at(j).at(k).at(l) += VdRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CledqR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CledqI", i, j, p, r))) * VdL(r, l);
                                                        if (mySMEFT.FlagNewTerms) {
                                                            CedSRL.at(i).at(j).at(k).at(l) += VdRd(k,p) * (-0.5*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", j, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", j, i))*(mySMEFT.getSMEFTCoeffEW("YdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", p, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", p, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", p, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", p, r)))*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)))))* VdL(r,l);
                                                            CedSRR.at(i).at(j).at(k).at(l) += VdLd(k,p) * (-0.5*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", j, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", j, i))*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", r, p)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)))*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i))))) * VdR(r,l);

                                                        }
                                                        for (int s = 0; s < 3; s++)
                                                                for (int t = 0; t < 3; t++)
                                                                {
                                                                        CddVLL.at(i).at(j).at(k).at(l) += VdLd(i, p) * VdLd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) + (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t))) * VdL(r, j) * VdL(t, l);
                                                                        CddVRR.at(i).at(j).at(k).at(l) += VdRd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("CddR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CddI", p, r, s, t))) * VdR(r, j) * VdR(t, l);
                                                                        CddV1LR.at(i).at(j).at(k).at(l) += VdLd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd1I", p, r, s, t))) * VdL(r, j) * VdR(t, l);
                                                                        CddV8LR.at(i).at(j).at(k).at(l) += VdLd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd8I", p, r, s, t))) * VdL(r, j) * VdR(t, l);
                                                                        if (mySMEFT.FlagNewTerms) {
                                                                            CddV1LR.at(i).at(j).at(k).at(l) += VdLd(i,p) * VdRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YdR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, p))))/12.) * VdL(r,j) * VdR(t,l);
                                                                            CddV8LR.at(i).at(j).at(k).at(l) += VdLd(i,p) * VdRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YdR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, p))))/2.) * VdL(r,j) * VdR(t,l);
                                                                            CddS1RR.at(i).at(j).at(k).at(l) += VdLd(i,p) * VdLd(k,s) * (-0.25*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", t, s) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", t, s))*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", r, p)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)))*(mySMEFT.getSMEFTCoeffEW("YdR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, s))))) * VdR(r,j) * VdR(t,l);
                                                                        }
                                                                }
                                                }
                                }
                                for (int l = 0; l < 2; l++)
                                {
                                        CnueduVLL.at(i).at(j).at(k).at(l) += VCKMd(k, l) * (-0.5 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHl3R", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", i, j))));
                                        if (mySMEFT.FlagNewTerms) {
                                            CnueduVLL.at(i).at(j).at(k).at(l) += VCKMd(k,l) *((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CeWR", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", i, j))*Me[j])/sqrt(2));
                                            CnueduSRR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHl3R", i, j) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", i, j))*Me[j]*Mu[l]));
                                        }
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                {
                                                        CnueduVLL.at(i).at(j).at(k).at(l) += VdLd(k, p) * (2 * (mySMEFT.getSMEFTCoeffEW("Clq3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", i, j, p, r))) * VuL(r, l);
                                                        CnueduSRL.at(i).at(j).at(k).at(l) += VdRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CledqR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CledqI", i, j, p, r))) * VuL(r, l);
                                                        CnueduSRR.at(i).at(j).at(k).at(l) += VdLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clequ1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ1I", i, j, p, r))) * VuR(r, l);
                                                        CnueduTRR.at(i).at(j).at(k).at(l) += VdLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clequ3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ3I", i, j, p, r))) * VuR(r, l);
                                                }
                                }
                        }
                        for (int k = 0; k < 2; k++)
                        {
                                CnueduVLL.at(i).at(i).at(j).at(k) += VCKMd(j, k) * (-0.5 * (delta_g22oMW2 * g22oMW2));
                                if (mySMEFT.FlagNewTerms) {
                              		CnueduSRR.at(i).at(i).at(j).at(k) += VCKMd(j,k) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Me[i]*Mu[k]));
                              	}
                                for (int p = 0; p < 3; p++)
                                        for (int r = 0; r < 3; r++)
                                        {
                                                CnueduVLL.at(i).at(i).at(j).at(k) += VdLd(j, p) * (-0.5 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p)))) * VuL(r, k);
                                                CnueduVLR.at(i).at(i).at(j).at(k) += VdRd(j, p) * (-0.25 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p)))) * VuR(r, k);
                                                if (mySMEFT.FlagNewTerms) {
                                                    CnueduVLL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p))*Md[j])/sqrt(2))) * VuL(r,k);
                                                    CnueduVLL.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))*Mu[k])/sqrt(2))) * VuR(r,k);
                                                    CnueduVLR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p))*Mu[k])/sqrt(2))) * VuL(r,k);
                                                    CnueduVLR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))*Md[j])/sqrt(2))) * VuR(r,k);
                                                    CnueduSRL.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Md[j]*Me[i])/2.) * VuL(r,k);
                                                    CnueduSRL.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.25*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Me[i]*Mu[k])) * VuR(r,k);
                                                    CnueduSRR.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[j]*Me[i])/4.) * VuR(r,k);
                                                    CnueduSRR.at(i).at(i).at(j).at(k) += VdLd(j,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Me[i]*Mu[k])) * VuL(r,k);
                                                }
                                        }
                                for (int l = 0; l < 2; l++)
                                {
                                if (mySMEFT.FlagNewTerms) {
                                    CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VCKM(k,j) * (-0.08333333333333333*((1 + deltag2bar2oMW4)*g2bar2oMW4*Mu[k]*Mu[l]));
                                    CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VCKM(k,j) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Mu[k]*Mu[l]));
                                }
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                {
                                                        CnuuVLL.at(i).at(j).at(k).at(l) += VuLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clq1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", i, j, p, r)) + (mySMEFT.getSMEFTCoeffEW("Clq3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", i, j, p, r))) * VuL(r, l);
                                                        CeuVLL.at(i).at(j).at(k).at(l) += VuLd(k, p) * ((mySMEFT.getSMEFTCoeffEW("Clq1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq1I", i, j, p, r)) - (mySMEFT.getSMEFTCoeffEW("Clq3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clq3I", i, j, p, r))) * VuL(r, l);
                                                        CeuVRR.at(i).at(j).at(k).at(l) += VuRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CeuR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeuI", i, j, p, r))) * VuR(r, l);
                                                        CnuuVLR.at(i).at(j).at(k).at(l) += VuRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CluR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CluI", i, j, p, r))) * VuR(r, l);
                                                        CeuVLR.at(i).at(j).at(k).at(l) += VuRd(k, p) * ((mySMEFT.getSMEFTCoeffEW("CluR", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CluI", i, j, p, r))) * VuR(r, l);
                                                        CeuSRR.at(i).at(j).at(k).at(l) += VuLd(k, p) * (-(mySMEFT.getSMEFTCoeffEW("Clequ1R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ1I", i, j, p, r))) * VuR(r, l);
                                                        CeuTRR.at(i).at(j).at(k).at(l) += VuLd(k, p) * (-(mySMEFT.getSMEFTCoeffEW("Clequ3R", i, j, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Clequ3I", i, j, p, r))) * VuR(r, l);
                                                        if (mySMEFT.FlagNewTerms) {
                                                            CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VuLd(k,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))*Mu[k]*Mu[l])) * VdL(r,j);
                                                            CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VuRd(k,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", p, r))*Md[j]*Mu[l])/24.) * VdR(r,j);
                                                            CduV1LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VdLd(i,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Mu[k]*Mu[l])) * VuL(r,l);
                                                            CduV1LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VdRd(i,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[i]*Mu[k])/24.) * VuR(r,l);
                                                            CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VuLd(k,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))*Mu[k]*Mu[l])) * VdL(r,j);
                                                            CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VuRd(k,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", p, r))*Md[j]*Mu[l])/4.) * VdR(r,j);
                                                            CduV8LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VdLd(i,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Mu[k]*Mu[l])) * VuL(r,l);
                                                            CduV8LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VdRd(i,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[i]*Mu[k])/4.) * VuR(r,l);
                                                            CeuSRL.at(i).at(j).at(k).at(l) += VuRd(k,p) * (-0.5*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", p, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", p, r))*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", j, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", j, i)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)))*(mySMEFT.getSMEFTCoeffEW("YuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", p, r)))))* VuL(r,l);
                                                            CeuSRR.at(i).at(j).at(k).at(l) += VuLd(k,p) * (-0.5*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", r, p))*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CeHR", j, i) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeHI", j, i)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YeR", j, i) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YeI", j, i)))*(mySMEFT.getSMEFTCoeffEW("YuR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", r, p))))) * VuR(r,l);
                                                        }
                                                        for (int s = 0; s < 3; s++)
                                                                for (int t = 0; t < 3; t++)
                                                                {
                                                                        CduV1LR.at(i).at(j).at(k).at(l) += VdLd(i, p) * VuRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqu1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu1I", p, r, s, t))) * VdL(r, j) * VuR(t, l);
                                                                        CduV8LR.at(i).at(j).at(k).at(l) += VdLd(i, p) * VuRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqu8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu8I", p, r, s, t))) * VdL(r, j) * VuR(t, l);
                                                                }
                                                }
                                }
                        }
                        for (int p = 0; p < 3; p++)
                                for (int r = 0; r < 3; r++)
                                {
                                        Cdg.at(i).at(j) += vTosq2 * VdLd(i, p) * (-(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)) * sbar + (mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) * cbar) * VdR(r, j);
                                        CdG.at(i).at(j) += vTosq2 * VdLd(i, p) * (mySMEFT.getSMEFTCoeffEW("CdGR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdGI", p, r)) * VdR(r, j);
                                }
                }

        for (int i = 0; i < 3; i++)
                for (int j = 0; j < 2; j++)
                {
                        CnuuVLL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (8 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 4 * sbar2))) / 12.;
                        CeuVLL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (4 * delta_sbar * sbar2 * (-5 + 8 * sbar2) + delta_gZ2oMZ2 * (3 - 10 * sbar2 + 8 * sbar2 * sbar2))) / 12.;
                        CeuVRR.at(i).at(i).at(j).at(j) += (2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 3.;
                        CnuuVLR.at(i).at(i).at(j).at(j) += ((delta_gZ2oMZ2 + 2 * delta_sbar) * gZ2oMZ2 * sbar2) / 3.;
                        CeuVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (-delta_gZ2oMZ2 - 2 * delta_sbar + 2 * delta_gZ2oMZ2 * sbar2 + 8 * delta_sbar * sbar2)) / 3.;
                        CduV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (-3 * delta_gZ2oMZ2 - 6 * delta_sbar + 2 * delta_gZ2oMZ2 * sbar2 + 8 * delta_sbar * sbar2)) / 9.;
                        if (mySMEFT.FlagNewTerms) {
                            CeuSRL.at(i).at(i).at(j).at(j) += -0.08333333333333333*(gZbar2oMZ4*(3 - 10*(1 + 2*delta_sbar)*sbar2 + 8*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 10*sbar2 + 8*pow(sbar2,2)))*Me[i]*Mu[j]);
                            CeuSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Mu[j])/3.;
                            CeuSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Me[i]*Mu[j])/6.;
                            CeuSRL.at(i).at(i).at(j).at(j) += (-2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Mu[j])/3.;
                            CeuSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(1 + deltagZbar2oMZ4 + delta_sbar*(2 - 8*sbar2) - 2*sbar2 - 2*deltagZbar2oMZ4*sbar2)*Me[i]*Mu[j])/3.;
                            CeuSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 - 10*(1 + 2*delta_sbar)*sbar2 + 8*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 10*sbar2 + 8*pow(sbar2,2)))*Me[i]*Mu[j])/12.;
                            CeuSRR.at(i).at(i).at(j).at(j) += (2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Mu[j])/3.;
                            CeuSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 + delta_sbar*(6 - 16*sbar2) + deltagZbar2oMZ4*(3 - 4*sbar2) - 4*sbar2)*sbar2*Me[i]*Mu[j])/6.;
                        }
                        for (int k = 0; k < 2; k++)
                                for (int p = 0; p < 3; p++)
                                        for (int r = 0; r < 3; r++)
                                        {
                                                CudV1LL.at(j).at(k).at(i).at(i) += VuLd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 12.) * VuL(r, k);
                                                CnuuVLL.at(i).at(i).at(j).at(k) += VuLd(j, p) * ((gZ2oMZ2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 4.) * VuL(r, k);
                                                CeuVLL.at(i).at(i).at(j).at(k) += VuLd(j, p) * ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 4.) * VuL(r, k);
                                                CudV1RR.at(j).at(k).at(i).at(i) += VuRd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 6.) * VuR(r, k);
                                                CeuVRR.at(i).at(i).at(j).at(k) += VuRd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 2.) * VuR(r, k);
                                                CnuuVLR.at(i).at(i).at(j).at(k) += VuRd(j, p) * ((gZ2oMZ2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 4.) * VuR(r, k);
                                                CeuVLR.at(i).at(i).at(j).at(k) += VuRd(j, p) * ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 4.) * VuR(r, k);
                                                CduV1LR.at(i).at(i).at(j).at(k) += VuRd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 12.) * VuR(r, k);
                                                CudV1LR.at(j).at(k).at(i).at(i) += VuLd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 6.) * VuL(r, k);
                                                CueVLR.at(j).at(k).at(i).at(i) += VuLd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 2.) * VuL(r, k);
                                                if (mySMEFT.FlagNewTerms) {
                                                    CnuuVLL.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/sqrt(2)) * VuL(r,k);
                                                    CnuuVLL.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/sqrt(2)) * VuR(r,k);
                                                    CeuVLL.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/sqrt(2))) * VuL(r,k);
                                                    CeuVLL.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/sqrt(2))) * VuR(r,k);
                                                    CudV1LL.at(j).at(k).at(i).at(i) += VuRd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/sqrt(2)) * VuL(r,k);
                                                    CudV1LL.at(j).at(k).at(i).at(i) += VuLd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/sqrt(2)) * VuR(r,k);
                                                    CeuVRR.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/sqrt(2))) * VuL(r,k);
                                                    CeuVRR.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/sqrt(2))) * VuR(r,k);
                                                    CudV1RR.at(j).at(k).at(i).at(i) += VuRd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/sqrt(2)) * VuL(r,k);
                                                    CudV1RR.at(j).at(k).at(i).at(i) += VuLd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/sqrt(2)) * VuR(r,k);
                                                    CnuuVLR.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/sqrt(2)) * VuL(r,k);
                                                    CnuuVLR.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/sqrt(2)) * VuR(r,k);
                                                    CeuVLR.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/sqrt(2))) * VuL(r,k);
                                                    CeuVLR.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/sqrt(2))) * VuR(r,k);
                                                    CueVLR.at(j).at(k).at(i).at(i) += VuRd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/sqrt(2))) * VuL(r,k);
                                                    CueVLR.at(j).at(k).at(i).at(i) += VuLd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/sqrt(2))) * VuR(r,k);
                                                    CudV1LR.at(j).at(k).at(i).at(i) += VuRd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/sqrt(2)) * VuL(r,k);
                                                    CudV1LR.at(j).at(k).at(i).at(i) += VuLd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/sqrt(2)) * VuR(r,k);
                                                    CduV1LR.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/sqrt(2)) * VuL(r,k);
                                                    CduV1LR.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/sqrt(2)) * VuR(r,k);
                                                    CudduV1LR.at(j).at(i).at(i).at(k) += VuLd(j,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Mu[k])) * VuL(r,k);
                                                    CudduV1LR.at(j).at(i).at(i).at(k) += VuRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Md[i]*Mu[j])/24.) * VuR(r,k);
                                                    CudduV8LR.at(j).at(i).at(i).at(k) += VuLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Mu[k])) * VuL(r,k);
                                                    CudduV8LR.at(j).at(i).at(i).at(k) += VuRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Md[i]*Mu[j])/4.) * VuR(r,k);
                                                    CeuSRL.at(i).at(i).at(j).at(k) += VuLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Me[i]*Mu[j])/4.) * VuL(r,k);
                                                    CeuSRL.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Me[i]*Mu[k])) * VuR(r,k);
                                                    CeuSRR.at(i).at(i).at(j).at(k) += VuRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Me[i]*Mu[j])/4.) * VuR(r,k);
                                                    CeuSRR.at(i).at(i).at(j).at(k) += VuLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Me[i]*Mu[k])) * VuL(r,k);
                                                    CudS1RR.at(j).at(k).at(i).at(i) += VuLd(j,p) * (-0.25*(gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[i]*Mu[k])) * VuL(r,k);
                                                    CudS1RR.at(j).at(k).at(i).at(i) += VuRd(j,p) * ((gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Md[i]*Mu[j])/4.) * VuR(r,k);
                                                }
                                        }
                }

        for (int i = 0; i < 2; i++)
                for (int j = 0; j < 3; j++)
                {
                        CudV1LL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (4 * delta_sbar * sbar2 * (-9 + 8 * sbar2) + delta_gZ2oMZ2 * (9 - 18 * sbar2 + 8 * sbar2 * sbar2))) / 36.;
                        CudV1RR.at(i).at(i).at(j).at(j) += (2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 9.;
                        CudV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (delta_gZ2oMZ2 * (-3 + 4 * sbar2) + 2 * delta_sbar * (-3 + 8 * sbar2))) / 18.;
                        CueVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (delta_gZ2oMZ2 * (-3 + 4 * sbar2) + 2 * delta_sbar * (-3 + 8 * sbar2))) / 6.;
                        if (mySMEFT.FlagNewTerms) {
                            CudduV1LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*(9 - 18*(1 + 2*delta_sbar)*sbar2 + 8*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(9 - 18*sbar2 + 8*pow(sbar2,2)))*Md[j]*Mu[i])/216.;
                            CudduV1LR.at(i).at(j).at(j).at(i) += -0.009259259259259259*(gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Md[j]*Mu[i]);
                            CudduV1LR.at(i).at(j).at(j).at(i) += -0.018518518518518517*(gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[j]*Mu[i]);
                            CudduV1LR.at(i).at(j).at(j).at(i) += ((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Mu[i])/27.;
                            CudduV8LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*(9 - 18*(1 + 2*delta_sbar)*sbar2 + 8*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(9 - 18*sbar2 + 8*pow(sbar2,2)))*Md[j]*Mu[i])/36.;
                            CudduV8LR.at(i).at(j).at(j).at(i) += -0.05555555555555555*(gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Md[j]*Mu[i]);
                            CudduV8LR.at(i).at(j).at(j).at(i) += -0.1111111111111111*(gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[j]*Mu[i]);
                            CudduV8LR.at(i).at(j).at(j).at(i) += (2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Mu[i])/9.;
                            CudS1RR.at(i).at(i).at(j).at(j) += -0.05555555555555555*(gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Md[j]*Mu[i]);
                            CudS1RR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(9 - 18*(1 + 2*delta_sbar)*sbar2 + 8*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(9 - 18*sbar2 + 8*pow(sbar2,2)))*Md[j]*Mu[i])/36.;
                            CudS1RR.at(i).at(i).at(j).at(j) += (2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Mu[i])/9.;
                            CudS1RR.at(i).at(i).at(j).at(j) += -0.1111111111111111*(gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[j]*Mu[i]);
                        }
                        for (int k = 0; k < 3; k++)
                        {
                                CeuVLL.at(j).at(k).at(i).at(i) += -0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k))));
                                CnuuVLL.at(j).at(k).at(i).at(i) += -0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) - (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k))));
                                CeuVRR.at(j).at(k).at(i).at(i) += -0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k)));
                                CnuuVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * (-(mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))) / 3.;
                                CeuVLR.at(j).at(k).at(i).at(i) += -0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k))));
                                CueVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))) / 12.;
                                if (mySMEFT.FlagNewTerms) {
                                    CeuVLL.at(j).at(k).at(i).at(i) += (gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/(6.*sqrt(2));
                                    CeuVRR.at(j).at(k).at(i).at(i) += (sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/3.;
                                    CeuVLR.at(j).at(k).at(i).at(i) += (sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k)))*Me[k]))/3.;
                                    CueVLR.at(i).at(i).at(j).at(k) += (gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", j, k))*Me[j] + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", j, k))*Me[j] + (sbar*(mySMEFT.getSMEFTCoeffEW("CeBR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeBI", k, j)) + cbar*(mySMEFT.getSMEFTCoeffEW("CeWR", k, j) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CeWI", k, j)))*Me[k]))/(6.*sqrt(2));
                                    CeuSRL.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j] - ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k])*Mu[i])/4.;
                                    CeuSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*(-((mySMEFT.getSMEFTCoeffEW("CHeR", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHeI", j, k))*Me[j]) + ((mySMEFT.getSMEFTCoeffEW("CHl1R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl1I", j, k)) + (mySMEFT.getSMEFTCoeffEW("CHl3R", j, k) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHl3I", j, k)))*Me[k])*Mu[i])/4.;
                                }
                                for (int p = 0; p < 3; p++)
                                        for (int r = 0; r < 3; r++)
                                        {
                                                CudV1LL.at(i).at(i).at(j).at(k) += VdLd(j, p) * (-0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))))) * VdL(r, k);
                                                CudV1RR.at(i).at(i).at(j).at(k) += VdRd(j, p) * (-0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r)))) * VdR(r, k);
                                                CudV1LR.at(i).at(i).at(j).at(k) += VdRd(j, p) * ((gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))) / 12.) * VdR(r, k);
                                                CduV1LR.at(j).at(k).at(i).at(i) += VdLd(j, p) * (-0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))))) * VdL(r, k);
                                                if (mySMEFT.FlagNewTerms) {
                                                    CudV1LL.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/(6.*sqrt(2))) * VdL(r,k);
                                                    CudV1LL.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/(6.*sqrt(2))) * VdR(r,k);
                                                    CudV1RR.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/3.) * VdL(r,k);
                                                    CudV1RR.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/3.) * VdR(r,k);
                                                    CudV1LR.at(i).at(i).at(j).at(k) += VdRd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[k])/(6.*sqrt(2))) * VdL(r,k);
                                                    CudV1LR.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[j])/(6.*sqrt(2))) * VdR(r,k);
                                                    CduV1LR.at(j).at(k).at(i).at(i) += VdRd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p)))*Md[j])/3.) * VdL(r,k);
                                                    CduV1LR.at(j).at(k).at(i).at(i) += VdLd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CdBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r)))*Md[k])/3.) * VdR(r,k);
                                                    CudduV1LR.at(i).at(j).at(k).at(i) += VdLd(k,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[k]*Mu[i])/24.) * VdL(r,j);
                                                    CudduV1LR.at(i).at(j).at(k).at(i) += VdRd(k,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[j]*Mu[i])) * VdR(r,j);
                                                    CudduV8LR.at(i).at(j).at(k).at(i) += VdLd(k,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[k]*Mu[i])/4.) * VdL(r,j);
                                                    CudduV8LR.at(i).at(j).at(k).at(i) += VdRd(k,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[j]*Mu[i])) * VdR(r,j);
                                                    CudS1RR.at(i).at(i).at(j).at(k) += VdRd(j,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHdR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHdI", p, r))*Md[j]*Mu[i])) * VdR(r,k);
                                                    CudS1RR.at(i).at(i).at(j).at(k) += VdLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Md[k]*Mu[i])/4.) * VdL(r,k);
                                                }
                                        }
                                for (int l = 0; l < 2; l++)
                                {
                                	if (mySMEFT.FlagNewTerms) {	
                                 		CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VCKMd(k,l) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[j]*Mu[l]));
                                 	}
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                {
                                                        CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i, j) * VdRd(k, p) * (-0.25 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p)))) * VuR(r, l);
                                                        if (mySMEFT.FlagNewTerms) {
                                                            CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VdRd(k,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p))*Mu[l])/sqrt(2))) * VuL(r,l);
                                                            CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VdLd(k,p) * (-((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))*Md[k])/sqrt(2))) * VuR(r,l);
                                                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VdRd(k,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[j]*Md[k])/4.) * VuR(r,l);
                                                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VdLd(k,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Md[j]*Mu[l])) * VuL(r,l);
                                                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * VuLd(i,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))*Md[j]*Mu[l])) * VdL(r,j);
                                                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * VuRd(i,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", p, r))*Mu[i]*Mu[l])/4.) * VdR(r,j);
                                                        }
                                                        for (int s = 0; s < 3; s++)
                                                                for (int t = 0; t < 3; t++)
                                                                {
                                                                        CudduS1RR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * (-(mySMEFT.getSMEFTCoeffEW("Cquqd1R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd1I", s, t, p, r))) * VdR(r, j) * VuR(t, l);
                                                                        CudduS8RR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * (-(mySMEFT.getSMEFTCoeffEW("Cquqd8R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd8I", s, t, p, r))) * VdR(r, j) * VuR(t, l);
                                                                        if (mySMEFT.FlagNewTerms) {
                                                                            CudduV1LR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VdRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YuR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", t, p))))/12.) * VdL(r,j) * VuR(t,l);
                                                                            CudduV8LR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VdRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YuR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", t, p))))/2.) * VdL(r,j) * VuR(t,l);
                                                                            CudduS1RR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VdLd(k,s) * (-0.25*(oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", t, s) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", t, s))*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", r, p)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", r, p)))*(mySMEFT.getSMEFTCoeffEW("YdR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, s))))) * VdR(r,j) * VuR(t,l);
                                                                        }
                                                                }
                                                }
                                }
                        }
                }

        for (int i = 0; i < 2; i++)
                for (int j = 0; j < 2; j++)
                {
                        CuuVLL.at(i).at(i).at(j).at(j) += -0.013888888888888888 * (gZ2oMZ2 * (-3 + 4 * sbar2) * (16 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 4 * sbar2)));
                        CuuVRR.at(i).at(i).at(j).at(j) += (-2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 9.;
                        CuuV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (2 * delta_sbar * (3 - 8 * sbar2) + delta_gZ2oMZ2 * (3 - 4 * sbar2)) * sbar2) / 9.;
                        if (mySMEFT.FlagNewTerms) {
                            CuuV1LR.at(i).at(j).at(j).at(i) += -0.004629629629629629*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                            CuuV1LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/54.;
                            CuuV1LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/54.;
                            CuuV1LR.at(i).at(j).at(j).at(i) += (-2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/27.;
                            CuuV8LR.at(i).at(j).at(j).at(i) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                            CuuV8LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/9.;
                            CuuV8LR.at(i).at(j).at(j).at(i) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/9.;
                            CuuV8LR.at(i).at(j).at(j).at(i) += (-4*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/9.;
                            CuuS1RR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/18.;
                            CuuS1RR.at(i).at(i).at(j).at(j) += -0.013888888888888888*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                            CuuS1RR.at(i).at(i).at(j).at(j) += (-2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/9.;
                            CuuS1RR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/18.;
                        }
                        for (int k = 0; k < 3; k++)
                        {
                                for (int l = 0; l < 3; l++)
                                {
                                        CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, j) * (-0.16666666666666666 * (delta_g22oMW2 * g22oMW2));
                                        CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, j) * (-(delta_g22oMW2 * g22oMW2));
                                        if (mySMEFT.FlagNewTerms) {
                                            CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,j) * (-0.08333333333333333*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[k]*Md[l]));
                                            CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,j) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[k]*Md[l]));
                                        }
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                {
                                                        CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k, j) * VuLd(i, p) * (-0.16666666666666666 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) * VdL(r, l);
                                                        CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VdLd(k, p) * (-0.16666666666666666 * (g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p)))) * VuL(r, j);
                                                        CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k, j) * VuLd(i, p) * (-(g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) * VdL(r, l);
                                                        CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VdLd(k, p) * (-(g22oMW2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p)))) * VuL(r, j);
                                                        CueVLR.at(i).at(j).at(k).at(l) += VuLd(i, p) * ((mySMEFT.getSMEFTCoeffEW("CqeR", p, r, k, l) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CqeI", p, r, k, l))) * VuL(r, j);
                                                        if (mySMEFT.FlagNewTerms) {
                                                            CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuRd(i,p) * ((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p))*Mu[i])/(3.*sqrt(2))) * VdL(r,l);
                                                            CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuLd(i,p) * ((g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r))*Md[l])/(3.*sqrt(2))) * VdR(r,l);
                                                            CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdRd(k,p) * (-0.3333333333333333*(g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p))*Md[k])/sqrt(2)) * VuL(r,j);
                                                            CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdLd(k,p) * (-0.3333333333333333*(g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))*Mu[j])/sqrt(2)) * VuR(r,j);
                                                            CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuRd(i,p) * (sqrt(2)*g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p))*Mu[i]) * VdL(r,l);
                                                            CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuLd(i,p) * (sqrt(2)*g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", p, r))*Md[l]) * VdR(r,l);
                                                            CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdRd(k,p) * (-(sqrt(2)*g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CdWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdWI", r, p))*Md[k])) * VuL(r,j);
                                                            CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdLd(k,p) * (-(sqrt(2)*g2bar2oMW3*v2*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r))*Mu[j])) * VuR(r,j);
                                                            CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdLd(k,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Md[k]*Md[l])) * VuL(r,j);
                                                            CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdRd(k,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[l]*Mu[j])/24.) * VuR(r,j);
                                                            CudV1LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuLd(i,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))*Md[k]*Md[l])) * VdL(r,l);
                                                            CudV1LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuRd(i,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", p, r))*Md[k]*Mu[i])/24.) * VdR(r,l);
                                                            CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdLd(k,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", r, p))*Md[k]*Md[l])) * VuL(r,j);
                                                            CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VdRd(k,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", r, p))*Md[l]*Mu[j])/4.) * VuR(r,j);
                                                            CudV8LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuLd(i,p) * (-0.5*(g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))*Md[k]*Md[l])) * VdL(r,l);
                                                            CudV8LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * VuRd(i,p) * ((g2bar2oMW4*v2*(mySMEFT.getSMEFTCoeffEW("CHudR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHudI", p, r))*Md[k]*Mu[i])/4.) * VdR(r,l);
                                                        }
                                                        for (int s = 0; s < 3; s++)
                                                                for (int t = 0; t < 3; t++)
                                                                {
                                                                        CudV1LL.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) + (mySMEFT.getSMEFTCoeffEW("Cqq1R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", s, t, p, r)) - (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t)) + (2 * (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, t, s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, t, s, r))) / 3. + (2 * (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, r, p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, r, p, t))) / 3. - (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, t, p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, t, p, r))) * VuL(r, j) * VdL(t, l);
                                                                        CudV8LL.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * (4 * ((mySMEFT.getSMEFTCoeffEW("Cqq3R", p, t, s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, t, s, r)) + (mySMEFT.getSMEFTCoeffEW("Cqq3R", s, r, p, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", s, r, p, t)))) * VuL(r, j) * VdL(t, l);
                                                                        CudV1RR.at(i).at(j).at(k).at(l) += VuRd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cud1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cud1I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                                                        CudV8RR.at(i).at(j).at(k).at(l) += VuRd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cud8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cud8I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                                                        CudV1LR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd1I", p, r, s, t))) * VuL(r, j) * VdR(t, l);
                                                                        CudV8LR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqd8I", p, r, s, t))) * VuL(r, j) * VdR(t, l);
                                                                        CudS1RR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cquqd1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd1I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                                                        CudS8RR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VdLd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cquqd8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cquqd8I", p, r, s, t))) * VuR(r, j) * VdR(t, l);
                                                                        if (mySMEFT.FlagNewTerms) {
                                                                            CudS1RR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VdLd(k,s) * ((oneoMh2*(-3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", r, p))*(mySMEFT.getSMEFTCoeffEW("YdR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, s)) + (-3*v2*(mySMEFT.getSMEFTCoeffEW("CdHR", t, s) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CdHI", t, s)) + 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YdR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YdI", t, s)))*(mySMEFT.getSMEFTCoeffEW("YuR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", r, p))))/2.) * VuR(r,j) * VdR(t,l);
                                                                        }
                                                                }
                                                }
                                }
                        }
                        for (int k = 0; k < 2; k++)
                        {
                                for (int p = 0; p < 3; p++)
                                        for (int r = 0; r < 3; r++)
                                        {
                                                CuuVLL.at(i).at(i).at(j).at(k) += VuLd(j, p) * (-0.041666666666666664 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))))) * VuL(r, k);
                                                CuuVLL.at(j).at(k).at(i).at(i) += VuLd(j, p) * (-0.041666666666666664 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r))))) * VuL(r, k);
                                                CuuVRR.at(i).at(i).at(j).at(k) += VuRd(j, p) * (-0.16666666666666666 * (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r)))) * VuR(r, k);
                                                CuuVRR.at(j).at(k).at(i).at(i) += VuRd(j, p) * (-0.16666666666666666 * (gZ2oMZ2 * sbar2 * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r)))) * VuR(r, k);
                                                CuuV1LR.at(i).at(i).at(j).at(k) += VuRd(j, p) * ((gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))) / 12.) * VuR(r, k);
                                                CuuV1LR.at(j).at(k).at(i).at(i) += VuLd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * (-(mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) + (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))) / 3.) * VuL(r, k);
                                                if (mySMEFT.FlagNewTerms) {
                                                CuuVLL.at(i).at(i).at(j).at(k) += VuRd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/(12.*sqrt(2))) * VuL(r,k);
                                                CuuVLL.at(i).at(i).at(j).at(k) += VuLd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/(12.*sqrt(2))) * VuR(r,k);
                                                CuuVLL.at(j).at(k).at(i).at(i) += VuRd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/(12.*sqrt(2))) * VuL(r,k);
                                                CuuVLL.at(j).at(k).at(i).at(i) += VuLd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/(12.*sqrt(2))) * VuR(r,k);
                                                CuuVRR.at(i).at(i).at(j).at(k) += VuRd(j,p) * ((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/(3.*sqrt(2))) * VuL(r,k);
                                                CuuVRR.at(i).at(i).at(j).at(k) += VuLd(j,p) * ((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/(3.*sqrt(2))) * VuR(r,k);
                                                CuuVRR.at(j).at(k).at(i).at(i) += VuRd(j,p) * ((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/(3.*sqrt(2))) * VuL(r,k);
                                                CuuVRR.at(j).at(k).at(i).at(i) += VuLd(j,p) * ((gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/(3.*sqrt(2))) * VuR(r,k);
                                                CuuV1LR.at(i).at(j).at(k).at(i) += VuLd(k,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/24.) * VuL(r,j);
                                                CuuV1LR.at(i).at(j).at(k).at(i) += VuRd(k,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,j);
                                                CuuV1LR.at(j).at(i).at(i).at(k) += VuLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/24.) * VuL(r,k);
                                                CuuV1LR.at(j).at(i).at(i).at(k) += VuRd(j,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,k);
                                                CuuV1LR.at(i).at(i).at(j).at(k) += VuRd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[k])/(6.*sqrt(2))) * VuL(r,k);
                                                CuuV1LR.at(i).at(i).at(j).at(k) += VuLd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[j])/(6.*sqrt(2))) * VuR(r,k);
                                                CuuV1LR.at(j).at(k).at(i).at(i) += VuRd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", r, p)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", r, p)))*Mu[j])/3.) * VuL(r,k);
                                                CuuV1LR.at(j).at(k).at(i).at(i) += VuLd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) + cbar*(mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)))*Mu[k])/3.) * VuR(r,k);
                                                CuuV8LR.at(i).at(j).at(k).at(i) += VuLd(k,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/4.) * VuL(r,j);
                                                CuuV8LR.at(i).at(j).at(k).at(i) += VuRd(k,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,j);
                                                CuuV8LR.at(j).at(i).at(i).at(k) += VuLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/4.) * VuL(r,k);
                                                CuuV8LR.at(j).at(i).at(i).at(k) += VuRd(j,p) * (-0.25*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,k);
                                                CuuS1RR.at(i).at(i).at(j).at(k) += VuRd(j,p) * (-0.125*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,k);
                                                CuuS1RR.at(i).at(i).at(j).at(k) += VuLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/8.) * VuL(r,k);
                                                CuuS1RR.at(j).at(k).at(i).at(i) += VuLd(j,p) * ((gZbar2oMZ4*v2*((mySMEFT.getSMEFTCoeffEW("CHq1R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq1I", p, r)) - (mySMEFT.getSMEFTCoeffEW("CHq3R", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHq3I", p, r)))*Mu[i]*Mu[k])/8.) * VuL(r,k);
                                                CuuS1RR.at(j).at(k).at(i).at(i) += VuRd(j,p) * (-0.125*(gZbar2oMZ4*v2*(mySMEFT.getSMEFTCoeffEW("CHuR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CHuI", p, r))*Mu[i]*Mu[j])) * VuR(r,k);
                                                }
                                        }
                                for (int l = 0; l < 2; l++)
                                        for (int p = 0; p < 3; p++)
                                                for (int r = 0; r < 3; r++)
                                                        for (int s = 0; s < 3; s++)
                                                                for (int t = 0; t < 3; t++)
                                                                {
                                                                        CuuVLL.at(i).at(j).at(k).at(l) += VuLd(i, p) * VuLd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqq1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq1I", p, r, s, t)) + (mySMEFT.getSMEFTCoeffEW("Cqq3R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqq3I", p, r, s, t))) * VuL(r, j) * VuL(t, l);
                                                                        CuuVRR.at(i).at(j).at(k).at(l) += VuRd(i, p) * VuRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("CuuR", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuuI", p, r, s, t))) * VuR(r, j) * VuR(t, l);
                                                                        CuuV1LR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VuRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqu1R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu1I", p, r, s, t))) * VuL(r, j) * VuR(t, l);
                                                                        CuuV8LR.at(i).at(j).at(k).at(l) += VuLd(i, p) * VuRd(k, s) * ((mySMEFT.getSMEFTCoeffEW("Cqu8R", p, r, s, t) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("Cqu8I", p, r, s, t))) * VuL(r, j) * VuR(t, l);
                                                                        if (mySMEFT.FlagNewTerms) {
                                                                            CuuV1LR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VuRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YuR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YuR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YuR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", t, p))))/12.) * VuL(r,j) * VuR(t,l);
                                                                            CuuV8LR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VuRd(k,s) * ((oneoMh2*(3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", t, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", t, p))*(mySMEFT.getSMEFTCoeffEW("YuR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", s, r)) + (3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", s, r) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", s, r)) - 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YuR", s, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", s, r)))*(mySMEFT.getSMEFTCoeffEW("YuR", t, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", t, p))))/2.) * VuL(r,j) * VuR(t,l);
                                                                            CuuS1RR.at(i).at(j).at(k).at(l) += VuLd(i,p) * VuLd(k,s) * ((oneoMh2*(-3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", t, s) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", t, s))*(mySMEFT.getSMEFTCoeffEW("YuR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", r, p)) + (-3*v2*(mySMEFT.getSMEFTCoeffEW("CuHR", r, p) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuHI", r, p)) + 2*(1 + 2*(mySMEFT.getSMEFTCoeffEW("CHbox") - 0.25 * mySMEFT.getSMEFTCoeffEW("CHD")) * v2 + deltaoneoMh2)*(mySMEFT.getSMEFTCoeffEW("YuR", r, p) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", r, p)))*(mySMEFT.getSMEFTCoeffEW("YuR", t, s) - gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("YuI", t, s))))/4.) * VuR(r,j) * VuR(t,l);
                                                                        }
                                                                }
                        }

                        for (int p = 0; p < 3; p++)
                                for (int r = 0; r < 3; r++)
                                {
                                        Cug.at(i).at(j) += vTosq2 * VuLd(i, p) * ((mySMEFT.getSMEFTCoeffEW("CuWR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuWI", p, r)) * sbar + (mySMEFT.getSMEFTCoeffEW("CuBR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuBI", p, r)) * cbar) * VuR(r, j);
                                        CuG.at(i).at(j) += vTosq2 * VuLd(i, p) * (mySMEFT.getSMEFTCoeffEW("CuGR", p, r) + gslpp::complex::i() * mySMEFT.getSMEFTCoeffEW("CuGI", p, r)) * VuR(r, j);
                                }
                }

#endif

        StandardModelMatching::updateSMParameters();
}

NPSMEFTd6GeneralMatching::~NPSMEFTd6GeneralMatching()
{
}

// Matching to the Delta F=2 Hamiltonian in the SUSY Basis, checked using 1512.02830

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMdk2()
{

        vmck2.clear();
        vmck2 = StandardModelMatching::CMdk2();

        mck2.setMu(mySMEFT.getMuw());

        switch (mck2.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mck2.setCoeff(l, 0., NLO);
        case LO:
                mck2.setCoeff(0, -CddVLL.at(0).at(1).at(0).at(1), LO);
                mck2.setCoeff(1, -(CddS1RR.at(1).at(0).at(1).at(0).conjugate() - CddS8RR.at(1).at(0).at(1).at(0).conjugate() / 6.), LO);
                mck2.setCoeff(2, -CddS8RR.at(1).at(0).at(1).at(0).conjugate() / 2., LO);
                mck2.setCoeff(3, CddV8LR.at(0).at(1).at(0).at(1), LO);
                mck2.setCoeff(4, 2. * CddV1LR.at(0).at(1).at(0).at(1) - CddV8LR.at(0).at(1).at(0).at(1) / 3., LO);
                break;
        default:
                std::stringstream out;
                out << mck2.getOrder();
                throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented");
        }

        // std::cout << "NPSMEFTd6GeneralMatching::CMk2(): Matching to the Delta F=2 Hamiltonian in the SUSY Basis, checked using 1512.02830" << std::endl;
        // std::cout << "C1 = " << (*(mck2.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2 = " << (*(mck2.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3 = " << (*(mck2.getCoeff(LO)))(2) << std::endl;
        // std::cout << "C4 = " << (*(mck2.getCoeff(LO)))(3) << std::endl;
        // std::cout << "C5 = " << (*(mck2.getCoeff(LO)))(4) << std::endl;
        // //mck2.setCoeff(0, 0., LO);


        vmck2.push_back(mck2);

        switch (mck2.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mck2.setCoeff(l, 0., NLO);
        case LO:
                mck2.setCoeff(0, -CddVRR.at(0).at(1).at(0).at(1), LO);
                mck2.setCoeff(1, -(CddS1RR.at(0).at(1).at(0).at(1) - CddS8RR.at(0).at(1).at(0).at(1) / 6.), LO);
                mck2.setCoeff(2, -CddS8RR.at(0).at(1).at(0).at(1) / 2., LO);
                mck2.setCoeff(3, 0., LO);
                mck2.setCoeff(4, 0., LO);
                break;
        default:
                std::stringstream out;
                out << mck2.getOrder();
                throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented");
        }

        vmck2.push_back(mck2);

        // std::cout << "C1t = " << (*(mck2.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2t = " << (*(mck2.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3t = " << (*(mck2.getCoeff(LO)))(2) << std::endl;

        return (vmck2);
}

std::vector<WilsonCoefficient> &NPSMEFTd6GeneralMatching::CMdd2()
{

        vmcd2.clear();
        vmcd2 = StandardModelMatching::CMdd2();

        mcd2.setMu(mySMEFT.getMuw());

        switch (mcd2.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcd2.setCoeff(l, 0., NLO);
        case LO:
                mcd2.setCoeff(0, -CuuVLL.at(0).at(1).at(0).at(1), LO);
                mcd2.setCoeff(1, -(CuuS1RR.at(1).at(0).at(1).at(0).conjugate() - CuuS8RR.at(1).at(0).at(1).at(0).conjugate() / 6.), LO);
                mcd2.setCoeff(2, -CuuS8RR.at(1).at(0).at(1).at(0).conjugate() / 2., LO);
                mcd2.setCoeff(3, CuuV8LR.at(0).at(1).at(0).at(1), LO);
                mcd2.setCoeff(4, 2. * CuuV1LR.at(0).at(1).at(0).at(1) - CuuV8LR.at(0).at(1).at(0).at(1) / 3., LO);
                break;
        default:
                std::stringstream out;
                out << mcd2.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdd2(): order " + out.str() + "not implemented");
        }

        // std::cout << "NPSMEFTd6GeneralMatching::CMdd2(): Matching to the Delta B=2 Hamiltonian in the SUSY Basis, checked using 1512.02830" << std::endl;
        // std::cout << "C1 = " << (*(mcd2.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2 = " << (*(mcd2.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3 = " << (*(mcd2.getCoeff(LO)))(2) << std::endl;
        // std::cout << "C4 = " << (*(mcd2.getCoeff(LO)))(3) << std::endl;
        // std::cout << "C5 = " << (*(mcd2.getCoeff(LO)))(4) << std::endl;
        // //mcd2.setCoeff(0, 0., LO);

        vmcd2.push_back(mcd2);

        switch (mcd2.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcd2.setCoeff(l, 0., NLO);
        case LO:
                mcd2.setCoeff(0, -CuuVRR.at(0).at(1).at(0).at(1), LO);
                mcd2.setCoeff(1, -(CuuS1RR.at(0).at(1).at(0).at(1) - CuuS8RR.at(0).at(1).at(0).at(1) / 6.), LO);
                mcd2.setCoeff(2, -CuuS8RR.at(0).at(1).at(0).at(1) / 2., LO);
                mcd2.setCoeff(3, 0., LO);
                mcd2.setCoeff(4, 0., LO);
                break;
        default:
                std::stringstream out;
                out << mcd2.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdd2(): order " + out.str() + "not implemented");
        }

        vmcd2.push_back(mcd2);

        // std::cout << "C1t = " << (*(mcd2.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2t = " << (*(mcd2.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3t = " << (*(mcd2.getCoeff(LO)))(2) << std::endl;

        return (vmcd2);
}

std::vector<WilsonCoefficient> &NPSMEFTd6GeneralMatching::CMdbd2()
{

        vmcdb.clear();
        vmcdb = StandardModelMatching::CMdbd2();

        mcbd.setMu(mySMEFT.getMuw());

        switch (mcbd.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcbd.setCoeff(l, 0., NLO);
        case LO:
                mcbd.setCoeff(0, -CddVLL.at(0).at(2).at(0).at(2), LO);
                mcbd.setCoeff(1, -(CddS1RR.at(2).at(0).at(2).at(0).conjugate() - CddS8RR.at(2).at(0).at(2).at(0).conjugate() / 6.), LO);
                mcbd.setCoeff(2, -CddS8RR.at(2).at(0).at(2).at(0).conjugate() / 2., LO);
                mcbd.setCoeff(3, CddV8LR.at(0).at(2).at(0).at(2), LO);
                mcbd.setCoeff(4, 2. * CddV1LR.at(0).at(2).at(0).at(2) - CddV8LR.at(0).at(2).at(0).at(2) / 3., LO);
                break;
        default:
                std::stringstream out;
                out << mcbd.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented");
        }


        // std::cout << "NPSMEFTd6GeneralMatching::CMdbd2(): Matching to the Delta B=2 Hamiltonian in the SUSY Basis, checked using 1512.02830" << std::endl;
        // std::cout << "C1 = " << (*(mcbd.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2 = " << (*(mcbd.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3 = " << (*(mcbd.getCoeff(LO)))(2) << std::endl;
        // std::cout << "C4 = " << (*(mcbd.getCoeff(LO)))(3) << std::endl;
        // std::cout << "C5 = " << (*(mcbd.getCoeff(LO)))(4) << std::endl;
        // //mcbd.setCoeff(0, 0., LO);

        vmcdb.push_back(mcbd);

        switch (mcbd.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcbd.setCoeff(l, 0., NLO);
        case LO:
                mcbd.setCoeff(0, -CddVRR.at(0).at(2).at(0).at(2), LO);
                mcbd.setCoeff(1, -(CddS1RR.at(0).at(2).at(0).at(2) - CddS8RR.at(0).at(2).at(0).at(2) / 6.), LO);
                mcbd.setCoeff(2, -CddS8RR.at(0).at(2).at(0).at(2) / 2., LO);
                mcbd.setCoeff(3, 0., LO);
                mcbd.setCoeff(4, 0., LO);
                break;
        default:
                std::stringstream out;
                out << mcbd.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented");
        }

        vmcdb.push_back(mcbd);

        // std::cout << "C1t = " << (*(mcbd.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2t = " << (*(mcbd.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3t = " << (*(mcbd.getCoeff(LO)))(2) << std::endl;

        return (vmcdb);
}

std::vector<WilsonCoefficient> &NPSMEFTd6GeneralMatching::CMdbs2()
{

        vmcds.clear();
        vmcds = StandardModelMatching::CMdbs2();

        mcbs.setMu(mySMEFT.getMuw());

        switch (mcbs.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcbs.setCoeff(l, 0., NLO);
        case LO:
                mcbs.setCoeff(0, -CddVLL.at(1).at(2).at(1).at(2), LO);
                mcbs.setCoeff(1, -(CddS1RR.at(2).at(1).at(2).at(1).conjugate() - CddS8RR.at(2).at(1).at(2).at(1).conjugate() / 6.), LO);
                mcbs.setCoeff(2, -CddS8RR.at(2).at(1).at(2).at(1).conjugate() / 2., LO);
                mcbs.setCoeff(3, CddV8LR.at(1).at(2).at(1).at(2), LO);
                mcbs.setCoeff(4, 2. * CddV1LR.at(1).at(2).at(1).at(2) - CddV8LR.at(1).at(2).at(1).at(2) / 3., LO);
                break;
        default:
                std::stringstream out;
                out << mcbs.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented");
        }


        // std::cout << "NPSMEFTd6GeneralMatching::CMdbs2(): Matching to the Delta BS=2 Hamiltonian in the SUSY Basis, checked using 1512.02830" << std::endl;
        // std::cout << "C1 = " << (*(mcbs.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2 = " << (*(mcbs.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3 = " << (*(mcbs.getCoeff(LO)))(2) << std::endl;
        // std::cout << "C4 = " << (*(mcbs.getCoeff(LO)))(3) << std::endl;
        // std::cout << "C5 = " << (*(mcbs.getCoeff(LO)))(4) << std::endl;
        // //mcbs.setCoeff(0, 0., LO);

        vmcds.push_back(mcbs);

        switch (mcbs.getOrder())
        {
        case NNLO:
        case NLO:
            for (int l = 0; l < 5; l++)
                mcbs.setCoeff(l, 0., NLO);
        case LO:
                mcbs.setCoeff(0, -CddVRR.at(1).at(2).at(1).at(2), LO);
                mcbs.setCoeff(1, -(CddS1RR.at(1).at(2).at(1).at(2) - CddS8RR.at(1).at(2).at(1).at(2) / 6.), LO);
                mcbs.setCoeff(2, - CddS8RR.at(1).at(2).at(1).at(2) / 2., LO);
                mcbs.setCoeff(3, 0., LO);
                mcbs.setCoeff(4, 0., LO);
                break;
        default:
                std::stringstream out;
                out << mcbs.getOrder();
                throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented");
        }

        vmcds.push_back(mcbs);
        // std::cout << "C1t = " << (*(mcbs.getCoeff(LO)))(0) << std::endl;
        // std::cout << "C2t = " << (*(mcbs.getCoeff(LO)))(1) << std::endl;
        // std::cout << "C3t = " << (*(mcbs.getCoeff(LO)))(2) << std::endl;

        return (vmcds);
}

/*******************************************************************************
 * Wilson coefficients matching, LEFT basis [1709.04486] ordered as CnueduVLLkkij^+, CnueduVLRkkij^+, CnueduSRRkkij^+, CnueduSRLkkij^+, CnueduTRRkkij^+             *
 * ****************************************************************************/
std::vector<WilsonCoefficient> &NPSMEFTd6GeneralMatching::CMdiujleptonknu(int i, int j, int k)
{

        vmculeptonnu.clear();
        vmculeptonnu = StandardModelMatching::CMdiujleptonknu(i, j, k);

        mculeptonnu.setMu(mySMEFT.getMuw());

        switch (mculeptonnu.getOrder())
        {
        case NNLO:
        case NLO:
        case LO:
                mculeptonnu.setCoeff(0, -(CnueduVLL.at(k).at(k).at(i).at(j)).conjugate(), LO);
                mculeptonnu.setCoeff(1, -(CnueduVLR.at(k).at(k).at(i).at(j)).conjugate(), LO);
                mculeptonnu.setCoeff(2, -(CnueduSRR.at(k).at(k).at(i).at(j)).conjugate(), LO);
                mculeptonnu.setCoeff(3, -(CnueduSRL.at(k).at(k).at(i).at(j)).conjugate(), LO);
                mculeptonnu.setCoeff(4, -(CnueduTRR.at(k).at(k).at(i).at(j)).conjugate(), LO);
                break;
        default:
                std::stringstream out;
                out << mculeptonnu.getOrder();
                throw std::runtime_error("StandardModelMatching::CMuleptonnu(): order " + out.str() + "not implemented");
        }

        vmculeptonnu.push_back(mculeptonnu);
        return (vmculeptonnu);
}

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMkpnn() {
        
    vmckpnn = StandardModelMatching::CMkpnn();

    mckpnn.setMu(mySMEFT.getMuw());
 
    switch (mckpnn.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
        // assume lepton universality for now
            mckpnn.setCoeff(0, -(CnudVLL.at(0).at(0).at(1).at(0) + CnudVLL.at(1).at(1).at(1).at(0) + CnudVLL.at(2).at(2).at(1).at(0))/3., LO);
            mckpnn.setCoeff(1, -(CnudVLR.at(0).at(0).at(1).at(0) + CnudVLR.at(1).at(1).at(1).at(0) + CnudVLR.at(2).at(2).at(1).at(0))/3., LO);
            break;
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("NPSMEFTd6GeneralMatching::CMkpnn(): order " + out.str() + " not implemented"); 
    }

    switch (mckpnn.getOrder_qed()) {
        case NLO_QED11:
        case LO_QED:
            break; 
        default:
            std::stringstream out;
            out << mckpnn.getOrder_qed();
            throw std::runtime_error("NPSMEFTd6GeneralMatching::CMkpnn(): qed order " + out.str() + " not implemented"); 
    }

    vmckpnn.push_back(mckpnn);
    return (vmckpnn);

}

/*******************************************************************************
 * Wilson coefficients misiak basis for b -> s g                               * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *
 * ****************************************************************************/
std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMbsg() {

    vmcbsg = StandardModelMatching::CMbsg();

    mcbsg.setMu(mySMEFT.getMuw());

    gslpp::complex LEFT_factor = sqrt(2.) / 4. / mySMEFT.getGF() / mySMEFT.getCKM().computelamt_s() ;
    gslpp::complex LEFT_factor_radiative = 16. * M_PI * M_PI / mySMEFT.getQuarks(QCD::BOTTOM).getMass() * LEFT_factor / sqrt(4. * M_PI * mySMEFT.getAle());
    
    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
        //  {O1, O2} = {{-(1/N), (-1 + N^2)/(2 N^2)}, {2, 1/N}} {OV8LLud,OV1LLud}
            mcbsg.setCoeff(0, (-1./mySMEFT.getNc() * getCudV8LL(1,1,1,2) + .5 * (1. - 1./mySMEFT.getNc()/mySMEFT.getNc()) * getCudV1LL(1,1,1,2)) * LEFT_factor, LO);
            mcbsg.setCoeff(1, (2. * getCudV8LL(1,1,1,2) + 1./mySMEFT.getNc() * getCudV1LL(1,1,1,2)) * LEFT_factor, LO);
            // Add penguin operators in the future
            mcbsg.setCoeff(6, getCdg(1,2) * LEFT_factor_radiative, LO);
            mcbsg.setCoeff(7, getCdG(1,2) * LEFT_factor_radiative * (mySMEFT.getAle()/mySMEFT.Als(mySMEFT.getQuarks(QCD::BOTTOM).getMass())), LO);
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented");
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}

std::vector<WilsonCoefficient>& NPSMEFTd6GeneralMatching::CMprimebsg() {

    vmcprimebsg = StandardModelMatching::CMprimebsg();

    mcprimebsg.setMu(mySMEFT.getMuw());

    gslpp::complex LEFT_factor = sqrt(2.) / 4. / mySMEFT.getGF() / mySMEFT.getCKM().computelamt_s() ;
    gslpp::complex LEFT_factor_radiative = 16. * M_PI * M_PI / mySMEFT.getQuarks(QCD::BOTTOM).getMass() * LEFT_factor / sqrt(4. * M_PI * mySMEFT.getAle());

    switch (mcprimebsg.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
        //  {O1prime, O2prime} = {{-(1/N), (-1 + N^2)/(2 N^2)}, {2, 1/N}} {OV8RRud,OV1RRud}
            mcprimebsg.setCoeff(0, (-1./mySMEFT.getNc() * getCudV8RR(1,1,1,2) + .5 * (1. - 1./mySMEFT.getNc()/mySMEFT.getNc()) * getCudV1RR(1,1,1,2)) * LEFT_factor, LO);
            mcprimebsg.setCoeff(1, (2. * getCudV8RR(1,1,1,2) + 1./mySMEFT.getNc() * getCudV1RR(1,1,1,2)) * LEFT_factor, LO);
            // Add penguin operators in the future
            mcprimebsg.setCoeff(6, getCdg(2,1).conjugate() * LEFT_factor_radiative, LO);
            mcprimebsg.setCoeff(7, getCdG(2,1).conjugate() * LEFT_factor_radiative * (mySMEFT.getAle()/mySMEFT.Als(mySMEFT.getQuarks(QCD::BOTTOM).getMass())), LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented");
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

/*******************************************************************************
 * Methods to get the NP contribution to the LEFT basis [1709.04486]          *
 * ****************************************************************************/

//dimension 6 four-fermion operators involving all left-handed fields

const gslpp::complex NPSMEFTd6GeneralMatching::getCnunuVLL(int i, int j, int k, int l) const
{
    return (CnunuVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCeeVLL(int i, int j, int k, int l) const
{
    return (CeeVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnueVLL(int i, int j, int k, int l) const
{
    return (CnueVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnuuVLL(int i, int j, int k, int l) const
{
    return (CnuuVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnudVLL(int i, int j, int k, int l) const
{
    return (CnudVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCeuVLL(int i, int j, int k, int l) const
{
    return (CeuVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCedVLL(int i, int j, int k, int l) const
{
    return (CedVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnueduVLL(int i, int j, int k, int l) const
{
    return (CnueduVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCuuVLL(int i, int j, int k, int l) const
{
    return (CuuVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCddVLL(int i, int j, int k, int l) const
{
    return (CddVLL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV1LL(int i, int j, int k, int l) const
{
    return (CudV1LL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV8LL(int i, int j, int k, int l) const
{
    return (CudV8LL.at(i).at(j).at(k).at(l));
}



//dimension 6 four-fermion operators involving all right-handed fields

const gslpp::complex NPSMEFTd6GeneralMatching::getCeeVRR(int i, int j, int k, int l) const
{
    return (CeeVRR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCeuVRR(int i, int j, int k, int l) const
{
    return (CeuVRR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCedVRR(int i, int j, int k, int l) const
{
    return (CedVRR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCuuVRR(int i, int j, int k, int l) const
{
    return (CuuVRR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCddVRR(int i, int j, int k, int l) const
{
    return (CddVRR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV1RR(int i, int j, int k, int l) const
{
    return (CudV1RR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV8RR(int i, int j, int k, int l) const
{
    return (CudV8RR.at(i).at(j).at(k).at(l));
}



//dimension 6 four-fermion operators involving a left-handed vector current and a right-handed vector current

const gslpp::complex NPSMEFTd6GeneralMatching::getCnueVLR(int i, int j, int k, int l) const
{
    return (CnueVLR.at(i).at(j).at(k).at(l));
}


const gslpp::complex NPSMEFTd6GeneralMatching::getCeeVLR(int i, int j, int k, int l) const
{
    return (CeeVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnuuVLR(int i, int j, int k, int l) const
{
    return (CnuuVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnudVLR(int i, int j, int k, int l) const
{
    return (CnudVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCeuVLR(int i, int j, int k, int l) const
{
    return (CeuVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCedVLR(int i, int j, int k, int l) const
{
    return (CedVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCueVLR(int i, int j, int k, int l) const
{
    return (CueVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCdeVLR(int i, int j, int k, int l) const
{
    return (CdeVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCnueduVLR(int i, int j, int k, int l) const
{
    return (CnueduVLR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCuuV1LR(int i, int j, int k, int l) const
{
    return (CuuV1LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCuuV8LR(int i, int j, int k, int l) const
{
    return (CuuV8LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV1LR(int i, int j, int k, int l) const
{
    return (CudV1LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudV8LR(int i, int j, int k, int l) const
{
    return (CudV8LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCduV1LR(int i, int j, int k, int l) const
{
    return (CduV1LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCduV8LR(int i, int j, int k, int l) const
{
    return (CduV8LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCddV1LR(int i, int j, int k, int l) const
{
    return (CddV1LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCddV8LR(int i, int j, int k, int l) const
{
    return (CddV8LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudduV1LR(int i, int j, int k, int l) const
{
    return (CudduV1LR.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCudduV8LR(int i, int j, int k, int l) const
{
    return (CudduV8LR.at(i).at(j).at(k).at(l));
}

    //dimension 6 four-fermion operators involving two right-handed scalar densities or tensor currents
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCeeSRR(int i, int j, int k, int l) const
{
    return (CeeSRR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCeuSRR(int i, int j, int k, int l) const
{
    return (CeuSRR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCeuTRR(int i, int j, int k, int l) const
{
    return (CeuTRR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCedSRR(int i, int j, int k, int l) const
{
    return (CedSRR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCedTRR(int i, int j, int k, int l) const
{
    return (CedTRR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCnueduSRR(int i, int j, int k, int l) const
{
    return (CnueduSRR.at(i).at(j).at(k).at(l));
}
    
 
const gslpp::complex NPSMEFTd6GeneralMatching::getCnueduTRR(int i, int j, int k, int l) const
{
    return (CnueduTRR.at(i).at(j).at(k).at(l));
}
    
  
const gslpp::complex NPSMEFTd6GeneralMatching::getCuuS1RR(int i, int j, int k, int l) const
{
    return (CuuS1RR.at(i).at(j).at(k).at(l));
}


const gslpp::complex NPSMEFTd6GeneralMatching::getCuuS8RR(int i, int j, int k, int l) const
{
    return (CuuS8RR.at(i).at(j).at(k).at(l));
}
    
   
const gslpp::complex NPSMEFTd6GeneralMatching::getCudS1RR(int i, int j, int k, int l) const
{
    return (CudS1RR.at(i).at(j).at(k).at(l));
}
    
    
const gslpp::complex NPSMEFTd6GeneralMatching::getCudS8RR(int i, int j, int k, int l) const
{
    return (CudS8RR.at(i).at(j).at(k).at(l));
}
    
    
const gslpp::complex NPSMEFTd6GeneralMatching::getCddS1RR(int i, int j, int k, int l) const
{
    return (CddS1RR.at(i).at(j).at(k).at(l));
}
    
   
const gslpp::complex NPSMEFTd6GeneralMatching::getCddS8RR(int i, int j, int k, int l) const
{
    return (CddS8RR.at(i).at(j).at(k).at(l));
}
    
   
const gslpp::complex NPSMEFTd6GeneralMatching::getCudduS1RR(int i, int j, int k, int l) const
{
    return (CudduS1RR.at(i).at(j).at(k).at(l));
}
    

const gslpp::complex NPSMEFTd6GeneralMatching::getCudduS8RR(int i, int j, int k, int l) const
{
    return (CudduS8RR.at(i).at(j).at(k).at(l));
}


//dimension 6 four-fermion operators involving a right-handed and a left-handed scalar density, plus hermitian conjugates
 
const gslpp::complex NPSMEFTd6GeneralMatching::getCeuSRL(int i, int j, int k, int l) const
{
    return (CeuSRL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCedSRL(int i, int j, int k, int l) const
{
    return (CedSRL.at(i).at(j).at(k).at(l));
}


const gslpp::complex NPSMEFTd6GeneralMatching::getCnueduSRL(int i, int j, int k, int l) const
{
    return (CnueduSRL.at(i).at(j).at(k).at(l));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCdG(int i, int j) const
{
    return (CdG.at(i).at(j));
}

const gslpp::complex NPSMEFTd6GeneralMatching::getCdg(int i, int j) const
{
    return (Cdg.at(i).at(j));
}

//Fermion rotation matrices to mass-eigenstate basis
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVuL() const
{
    return VuL;
}   

const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVuR() const
{
    return VuR;
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVdL() const
{
    return VdL;
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVdR() const
{
    return VdL;
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVeL() const
{
    return VeL;
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVeR() const
{
    return VeR;
}   
