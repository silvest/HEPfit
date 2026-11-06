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

namespace {
// EW-scale Wilson coefficients cached locally with callable flavour-index access.
struct Coeff2 {
    gslpp::complex v[9];
    explicit Coeff2(const gslpp::complex (&c)[3][3])
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                v[3 * i + j] = c[i][j];
    }
    const gslpp::complex &operator()(int i, int j) const { return v[3 * i + j]; }
};

struct Coeff4 {
    gslpp::complex v[81];
    explicit Coeff4(const gslpp::complex (&c)[3][3][3][3])
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        v[27 * i + 9 * j + 3 * k + l] = c[i][j][k][l];
    }
    const gslpp::complex &operator()(int i, int j, int k, int l) const { return v[27 * i + 9 * j + 3 * k + l]; }
};
}

NPSMEFTd6GeneralMatching::NPSMEFTd6GeneralMatching(const NPSMEFTd6General &NPSMEFTd6General_i) : StandardModelMatching(NPSMEFTd6General_i),
                                                                                                 mySMEFT(NPSMEFTd6General_i), 
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
        const double c_CHbox = mySMEFT.getCHbox_EW();
        const double c_CHD = mySMEFT.getCHD_EW();
        const double c_CH = mySMEFT.getCH_EW();
        const double c_CHWB = mySMEFT.getCHWB_EW();

        //    CG = mySMEFT.getSMEFTCoeffEW("CG")*LambdaNP2;
        //    CW = mySMEFT.getSMEFTCoeffEW("CW")*LambdaNP2;
        //    CHG = mySMEFT.getSMEFTCoeffEW("CHG")*LambdaNP2;
        //    CHW = mySMEFT.getSMEFTCoeffEW("CHW")*LambdaNP2;
        //    CHB = mySMEFT.getSMEFTCoeffEW("CHB")*LambdaNP2;
        //    CHWB = c_CHWB*LambdaNP2;
        //    CHD = c_CHD*LambdaNP2;
        //    CHbox = c_CHbox*LambdaNP2;
        //    CH = c_CH*LambdaNP2;
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
        double delta_gZbar = (g1bar * g1bar + g2bar * g2bar) / (2. * g1bar * g2bar) * v2 * c_CHWB;
        // indeed a piece was missing, I add it here
        delta_gZbar += -(g1bar * g1bar * g1bar * g1bar + g2bar * g2bar * g2bar * g2bar) / (2. * (g1bar * g1bar + g2bar * g2bar) * g1bar * g2bar) * v2 * c_CHWB;
        double gZ2oMZ2 = gZbar / mySMEFT.getMz();
        gZ2oMZ2 *= gZ2oMZ2;
        double delta_gZ2oMZ2 = 2. * delta_gZbar - delta_MZ2;
        double g22oMW2 = 4. / v2;
        double delta_g22oMW2 = -2. * delta_vT;
        //new terms from now on
        //Z boson
        double gZbar2oMZ3 = mySMEFT.FlagMomVert * gZ2oMZ2 / mySMEFT.getMz();
        double gZbar2oMZ4 = mySMEFT.FlagMomProp * gZbar2oMZ3 / mySMEFT.getMz();
        double deltagZbar2oMZ4 = 2. * delta_gZbar - 2. * delta_MZ2;
        //W boson
        double g2bar2oMW3 = mySMEFT.FlagMomVert * g22oMW2 / mySMEFT.getMw();
        double g2bar2oMW4 = mySMEFT.FlagMomProp * g2bar2oMW3 / mySMEFT.getMw();
        double deltag2bar2oMW4 = -4. * delta_vT;
        //h boson
        double lambda = mySMEFT.getLambdaH_tree() * 2.;
        double oneoMh2 = mySMEFT.FlagHiggsExch * 1. / (lambda * v2);
        double deltaoneoMh2 = (- 2. * c_CHbox + 0.5 * c_CHD + 3 * c_CH / lambda) * v2;
        
        std::array<double, 3> Me = {
    	mySMEFT.getLeptons(StandardModel::ELECTRON).getMass(),
    	mySMEFT.getLeptons(StandardModel::MU).getMass(),
    	mySMEFT.getLeptons(StandardModel::TAU).getMass()
};
        std::array<double, 3> Mu = {
        mySMEFT.getQuarks(QCD::UP).getMass(),   
        mySMEFT.getQuarks(QCD::CHARM).getMass()
};
        std::array<double, 3> Md = {
        mySMEFT.getQuarks(QCD::DOWN).getMass(),
        mySMEFT.getQuarks(QCD::STRANGE).getMass(),
        mySMEFT.getQuarks(QCD::BOTTOM).getMass(),
};

        // std::cout << "CKM from rotated UfA = " << (VuL.hconjugate()) * VdL << std::endl;

        // std::cout << "has the diagonalization worked? " << VuR.hconjugate()*MU*VuL << std::endl;
        // std::cout << "has the diagonalization worked? " << VdR.hconjugate()*MD*VdL << std::endl;

        // match and rotate following Manohar. This is performed AT LINEAR ORDER for the moment

        // Coefficients are read once here instead of at every use inside the loops below

        const Coeff2 c_CHl1(mySMEFT.getCHl1_EW());
        const Coeff2 c_CHl3(mySMEFT.getCHl3_EW());
        const Coeff2 c_CHe(mySMEFT.getCHe_EW());
        const Coeff2 c_CHq1(mySMEFT.getCHq1_EW());
        const Coeff2 c_CHq3(mySMEFT.getCHq3_EW());
        const Coeff2 c_CHu(mySMEFT.getCHu_EW());
        const Coeff2 c_CHd(mySMEFT.getCHd_EW());
        const Coeff2 c_CHud(mySMEFT.getCHud_EW());
        const Coeff2 c_CeH(mySMEFT.getCeH_EW());
        const Coeff2 c_CuH(mySMEFT.getCuH_EW());
        const Coeff2 c_CdH(mySMEFT.getCdH_EW());
        const Coeff2 c_CuG(mySMEFT.getCuG_EW());
        const Coeff2 c_CuW(mySMEFT.getCuW_EW());
        const Coeff2 c_CuB(mySMEFT.getCuB_EW());
        const Coeff2 c_CdG(mySMEFT.getCdG_EW());
        const Coeff2 c_CdW(mySMEFT.getCdW_EW());
        const Coeff2 c_CdB(mySMEFT.getCdB_EW());
        const Coeff2 c_CeW(mySMEFT.getCeW_EW());
        const Coeff2 c_CeB(mySMEFT.getCeB_EW());
        const Coeff4 c_Cll(mySMEFT.getCll_EW());
        const Coeff4 c_Clq1(mySMEFT.getClq1_EW());
        const Coeff4 c_Clq3(mySMEFT.getClq3_EW());
        const Coeff4 c_Cee(mySMEFT.getCee_EW());
        const Coeff4 c_Ceu(mySMEFT.getCeu_EW());
        const Coeff4 c_Ced(mySMEFT.getCed_EW());
        const Coeff4 c_Cle(mySMEFT.getCle_EW());
        const Coeff4 c_Clu(mySMEFT.getClu_EW());
        const Coeff4 c_Cld(mySMEFT.getCld_EW());
        const Coeff4 c_Cqe(mySMEFT.getCqe_EW());
        const Coeff4 c_Cledq(mySMEFT.getCledq_EW());
        const Coeff4 c_Cqq1(mySMEFT.getCqq1_EW());
        const Coeff4 c_Cqq3(mySMEFT.getCqq3_EW());
        const Coeff4 c_Cuu(mySMEFT.getCuu_EW());
        const Coeff4 c_Cdd(mySMEFT.getCdd_EW());
        const Coeff4 c_Cud1(mySMEFT.getCud1_EW());
        const Coeff4 c_Cud8(mySMEFT.getCud8_EW());
        const Coeff4 c_Cqu1(mySMEFT.getCqu1_EW());
        const Coeff4 c_Cqu8(mySMEFT.getCqu8_EW());
        const Coeff4 c_Cqd1(mySMEFT.getCqd1_EW());
        const Coeff4 c_Cqd8(mySMEFT.getCqd8_EW());
        const Coeff4 c_Cquqd1(mySMEFT.getCquqd1_EW());
        const Coeff4 c_Cquqd8(mySMEFT.getCquqd8_EW());
        const Coeff4 c_Clequ1(mySMEFT.getClequ1_EW());
        const Coeff4 c_Clequ3(mySMEFT.getClequ3_EW());

        const gslpp::matrix<gslpp::complex> c_Ye = mySMEFT.getYe();
        const gslpp::matrix<gslpp::complex> c_Yu = mySMEFT.getYu();
        const gslpp::matrix<gslpp::complex> c_Yd = mySMEFT.getYd();


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
        // Work in the up-quark mass basis: VuL = VuR = VdR = 1, VdL = VCKM.
        // Identity rotations have already been contracted; only CKM sums remain.

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Ceg.at(i).at(j) += vTosq2 * (-(c_CeW(i, j)) * sbar + (c_CeB(i, j)) * cbar);
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
                if (mySMEFT.FlagNewTerms)
                {
                    CnueVLR.at(i).at(j).at(j).at(i) += -0.25*((1 + deltag2bar2oMW4)*g2bar2oMW4*Me[i]*Me[j]);
                    CeeVLR.at(i).at(j).at(j).at(i) += -0.125*(gZbar2oMZ4*(-1 + 2*sbar2)*(-1 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2))*Me[i]*Me[j]);
                    CeeVLR.at(i).at(j).at(j).at(i) += 2. * ((gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.);
                    CeeVLR.at(i).at(j).at(j).at(i) += -0.5*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Me[j]);
                    CddV1LR.at(i).at(j).at(j).at(i) += -0.004629629629629629*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                    CddV1LR.at(i).at(j).at(j).at(i) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/108.);
                    CddV1LR.at(i).at(j).at(j).at(i) += -0.018518518518518517*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                    CddV8LR.at(i).at(j).at(j).at(i) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                    CddV8LR.at(i).at(j).at(j).at(i) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.);
                    CddV8LR.at(i).at(j).at(j).at(i) += -0.1111111111111111*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                    CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 - 8*(1 + 2*delta_sbar)*sbar2 + 4*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 8*sbar2 + 4*pow(sbar2,2)))*Md[j]*Me[i])/12.;
                    CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(1 + deltagZbar2oMZ4 + delta_sbar*(2 - 8*sbar2) - 2*sbar2 - 2*deltagZbar2oMZ4*sbar2)*Md[j]*Me[i])/6.;
                    CedSRL.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*(3 + delta_sbar*(6 - 8*sbar2) + deltagZbar2oMZ4*(3 - 2*sbar2) - 2*sbar2)*sbar2*Md[j]*Me[i])/6.;
                    CedSRL.at(i).at(i).at(j).at(j) += ((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Me[i])/3.;
                    CeeSRR.at(i).at(i).at(j).at(j) += 2. * ((gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Me[i]*Me[j])/4.);
                    CeeSRR.at(i).at(i).at(j).at(j) += -0.125*(gZbar2oMZ4*(-1 + 2*sbar2)*(-1 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2))*Me[i]*Me[j]);
                    CeeSRR.at(i).at(i).at(j).at(j) += -0.5*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Me[i]*Me[j]);
                    CedSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-1 + 2*sbar2 + deltagZbar2oMZ4*(-1 + 2*sbar2) + delta_sbar*(-2 + 8*sbar2))*Md[j]*Me[i])/6.;
                    CedSRR.at(i).at(i).at(j).at(j) += -0.08333333333333333*(gZbar2oMZ4*(3 - 8*(1 + 2*delta_sbar)*sbar2 + 4*(1 + 4*delta_sbar)*pow(sbar2,2) + deltagZbar2oMZ4*(3 - 8*sbar2 + 4*pow(sbar2,2)))*Md[j]*Me[i]);
                    CedSRR.at(i).at(i).at(j).at(j) += -0.3333333333333333*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[j]*Me[i]);
                    CedSRR.at(i).at(i).at(j).at(j) += (gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[j]*Me[i])/6.;
                    CddS1RR.at(i).at(i).at(j).at(j) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 2*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2) + delta_sbar*(-6 + 8*sbar2))*Md[i]*Md[j])/18.);
                    CddS1RR.at(i).at(i).at(j).at(j) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 2*sbar2)*(-3 + (2 + 8*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 2*sbar2))*Md[i]*Md[j]);
                    CddS1RR.at(i).at(i).at(j).at(j) += -0.1111111111111111*((1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Md[i]*Md[j]);
                }
                for (int k = 0; k < 3; k++)
                {
                    CnunuVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 16.;
                    CnunuVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 16.;
                    CnunuVLL.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 16.;
                    CnunuVLL.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * v2 * ((c_CHl1(k, j)) - (c_CHl3(k, j)))) / 16.;
                    CeeVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 16.;
                    CeeVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 16.;
                    CeeVLL.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 16.;
                    CeeVLL.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHl1(k, j)) + (c_CHl3(k, j)))) / 16.;
                    CnueVLL.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 4.;
                    CnueVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 4.;
                    CnueVLL.at(j).at(i).at(i).at(k) += -0.5 * (g22oMW2 * v2 * (c_CHl3(j, k)));
                    CnueVLL.at(i).at(j).at(k).at(i) += -0.5 * (g22oMW2 * v2 * (c_CHl3(j, k).conjugate()));
                    CnudVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 12.;
                    CedVLL.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 12.;
                    CeeVRR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * sbar2 * v2 * (c_CHe(j, k))) / 8.;
                    CeeVRR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHe(j, k))) / 8.);
                    CeeVRR.at(j).at(i).at(i).at(k) += (gZ2oMZ2 * sbar2 * v2 * (c_CHe(j, k))) / 8.;
                    CeeVRR.at(i).at(j).at(k).at(i) += (gZ2oMZ2 * sbar2 * v2 * (c_CHe(k, j))) / 8.;
                    CedVRR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * (c_CHe(j, k))) / 6.;
                    CnueVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * v2 * (c_CHe(j, k))) / 4.;
                    CnueVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 2.;
                    CeeVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (c_CHe(j, k))) / 4.;
                    CeeVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 2.;
                    CnudVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k)))) / 6.;
                    CedVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k)))) / 6.;
                    CdeVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (c_CHe(j, k))) / 12.;
                    if (mySMEFT.FlagNewTerms)
                    {
                        CeeVLL.at(i).at(i).at(j).at(k) += (gZbar2oMZ3*(0.5 - sbar2)*v2*((sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/(4.*sqrt(2));
                        CeeVLL.at(j).at(k).at(i).at(i) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CeeVLL.at(j).at(i).at(i).at(k) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CeeVLL.at(i).at(j).at(k).at(i) += -0.125*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(c_CeB(k, j))*Me[j] + cbar*(c_CeW(k, j))*Me[j] + (sbar*(c_CeB(j, k).conjugate()) + cbar*(c_CeW(j, k).conjugate()))*Me[k]))/sqrt(2);
                        CnueVLL.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*v2*((sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CnueVLL.at(j).at(i).at(i).at(k) += (g2bar2oMW3*v2*(c_CeW(j, k))*Me[k])/sqrt(2);
                        CnueVLL.at(i).at(j).at(k).at(i) += -((g2bar2oMW3*v2*(c_CeW(j, k).conjugate())*Me[k])/sqrt(2));
                        CedVLL.at(j).at(k).at(i).at(i) += -0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CeeVRR.at(i).at(i).at(j).at(k) += -0.25*(gZbar2oMZ3*sbar2*v2*((sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CeeVRR.at(j).at(k).at(i).at(i) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CeeVRR.at(j).at(i).at(i).at(k) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CeeVRR.at(i).at(j).at(k).at(i) += -0.25*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CedVRR.at(j).at(k).at(i).at(i) += -0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CnueVLR.at(i).at(j).at(k).at(i) += -0.25*(g2bar2oMW4*v2*(c_CHl3(j, k).conjugate())*Me[i]*Me[k]);
                        CnueVLR.at(j).at(i).at(i).at(k) += -0.25*(g2bar2oMW4*v2*(c_CHl3(j, k))*Me[i]*Me[k]);
                        CnueVLR.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CeeVLR.at(i).at(j).at(k).at(i) += (gZbar2oMZ4*v2*Me[i]*((c_CHe(k, j))*Me[j] - ((c_CHl1(k, j)) + (c_CHl3(k, j)))*Me[k]))/8.;
                        CeeVLR.at(j).at(i).at(i).at(k) += (gZbar2oMZ4*v2*Me[i]*((c_CHe(j, k))*Me[j] - ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k]))/8.;
                        CeeVLR.at(i).at(i).at(j).at(k) += -0.5*(gZbar2oMZ3*(-1 + 2*sbar2)*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CeeVLR.at(j).at(k).at(i).at(i) += -((gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2));
                        CedVLR.at(j).at(k).at(i).at(i) += -0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/sqrt(2);
                        CdeVLR.at(i).at(i).at(j).at(k) += -0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*((sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/sqrt(2);
                        CedSRL.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Md[i]*(-((c_CHe(j, k))*Me[j]) + ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k]))/4.;
                        CeeSRR.at(i).at(i).at(j).at(k) += (gZbar2oMZ4*v2*Me[i]*((c_CHe(j, k))*Me[j] - ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k]))/8.;
                        CeeSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Me[i]*((c_CHe(j, k))*Me[j] - ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k]))/8.;
                        CedSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*Md[i]*((c_CHe(j, k))*Me[j] - ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k]))/4.;
                        for (int p = 0; p < 3; p++)
                        {
                            CnudVLL.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2)) * VCKM(p,k);
                            CnudVLL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2));
                            CedVLL.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2))) * VCKM(p,k);
                            CedVLL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2)));
                            CddVLL.at(i).at(i).at(j).at(k) += (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2)) * VCKM(p,k);
                            CddVLL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2));
                            CddVLL.at(j).at(k).at(i).at(i) += (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2)) * VCKM(p,k);
                            CddVLL.at(j).at(k).at(i).at(i) += VCKMd(j,p) * (-0.041666666666666664*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2));
                            CedVRR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2))) * VCKM(p,k);
                            CedVRR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2)));
                            CddVRR.at(i).at(i).at(j).at(k) += (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2)) * VCKM(p,k);
                            CddVRR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2));
                            CddVRR.at(j).at(k).at(i).at(i) += (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2)) * VCKM(p,k);
                            CddVRR.at(j).at(k).at(i).at(i) += VCKMd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2));
                            CnudVLR.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2)) * VCKM(p,k);
                            CnudVLR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2));
                            CedVLR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2))) * VCKM(p,k);
                            CedVLR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2)));
                            CdeVLR.at(j).at(k).at(i).at(i) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2))) * VCKM(p,k);
                            CdeVLR.at(j).at(k).at(i).at(i) += VCKMd(j,p) * (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2)));
                            for (int r = 0; r < 3; r++)
                            {
                                CddV1LR.at(i).at(j).at(k).at(i) += VCKMd(k,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,j);
                                CddV1LR.at(j).at(i).at(i).at(k) += VCKMd(j,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,k);
                                CddV8LR.at(i).at(j).at(k).at(i) += VCKMd(k,p) * (-0.25*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,j);
                                CddV8LR.at(j).at(i).at(i).at(k) += VCKMd(j,p) * (-0.25*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,k);
                                CedSRL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[j]*Me[i])/4.) * VCKM(r,k);
                                CedSRR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.25*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k]*Me[i])) * VCKM(r,k);
                                CddS1RR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.25*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,k);
                                CddS1RR.at(j).at(k).at(i).at(i) += VCKMd(j,p) * (-0.25*(gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[i]*Md[k])) * VCKM(r,k);
                            }
                            CddV1LR.at(i).at(i).at(j).at(k) += (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/sqrt(2)) * VCKM(p,k);
                            CddV1LR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2));
                            CddV1LR.at(j).at(k).at(i).at(i) += (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/sqrt(2)) * VCKM(p,k);
                            CddV1LR.at(j).at(k).at(i).at(i) += VCKMd(j,p) * (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/sqrt(2));
                        }
                        CddV1LR.at(i).at(j).at(k).at(i) += ((gZbar2oMZ4*v2*(c_CHd(k, j))*Md[i]*Md[j])/24.);
                        CddV1LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*(c_CHd(j, k))*Md[i]*Md[j])/24.);
                        CddV8LR.at(i).at(j).at(k).at(i) += ((gZbar2oMZ4*v2*(c_CHd(k, j))*Md[i]*Md[j])/4.);
                        CddV8LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*(c_CHd(j, k))*Md[i]*Md[j])/4.);
                        CedSRL.at(i).at(i).at(j).at(k) += (-0.25*(gZbar2oMZ4*v2*(c_CHd(j, k))*Md[k]*Me[i]));
                        CedSRR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ4*v2*(c_CHd(j, k))*Md[j]*Me[i])/4.);
                        CddS1RR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ4*v2*(c_CHd(j, k))*Md[i]*Md[j])/4.);
                        CddS1RR.at(j).at(k).at(i).at(i) += ((gZbar2oMZ4*v2*(c_CHd(j, k))*Md[i]*Md[j])/4.);
                    }
                    for (int p = 0; p < 3; p++)
                    {
                        for (int r = 0; r < 3; r++)
                        {
                            CnudVLL.at(i).at(i).at(j).at(k) += VCKMd(j, p) * ((gZ2oMZ2 * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 4.) * VCKM(r, k);
                            CedVLL.at(i).at(i).at(j).at(k) += VCKMd(j, p) * ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 4.) * VCKM(r, k);
                            CddVLL.at(i).at(i).at(j).at(k) += VCKMd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 24.) * VCKM(r, k);
                            CddVLL.at(j).at(k).at(i).at(i) += VCKMd(j, p) * ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 24.) * VCKM(r, k);
                            CdeVLR.at(j).at(k).at(i).at(i) += VCKMd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 2.) * VCKM(r, k);
                            CddV1LR.at(j).at(k).at(i).at(i) += VCKMd(j, p) * ((gZ2oMZ2 * sbar2 * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r)))) / 6.) * VCKM(r, k);
                        }
                        CnudVLR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/sqrt(2));
                    }
                    CedVRR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHd(j, k))) / 2.);
                    CddVRR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHd(j, k))) / 12.);
                    CddVRR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHd(j, k))) / 12.);
                    CnudVLR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * v2 * (c_CHd(j, k))) / 4.);
                    CedVLR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (c_CHd(j, k))) / 4.);
                    CddV1LR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (c_CHd(j, k))) / 12.);
                    for (int l = 0; l < 3; l++)
                    {
                        CnunuVLL.at(i).at(j).at(k).at(l) += (c_Cll(i, j, k, l));
                        CeeVLL.at(i).at(j).at(k).at(l) += (c_Cll(i, j, k, l));
                        CnueVLL.at(i).at(j).at(k).at(l) += (c_Cll(i, j, k, l)) + (c_Cll(k, l, i, j));
                        CeeVRR.at(i).at(j).at(k).at(l) += (c_Cee(i, j, k, l));
                        CnueVLR.at(i).at(j).at(k).at(l) += (c_Cle(i, j, k, l));
                        CeeVLR.at(i).at(j).at(k).at(l) += (c_Cle(i, j, k, l));
                        if (mySMEFT.FlagNewTerms)
                        {
                            CeeVLR.at(i).at(j).at(k).at(l) += (oneoMh2*(3*v2*(c_CeH(i, l))*(c_Ye(k, j)) + (3*v2*(c_CeH(j, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Ye(k, j)))*(c_Ye(l, i).conjugate())))/4.;
                            CedSRL.at(i).at(j).at(k).at(l) += VCKMd(k,l) * ((g2bar2oMW4*v2*(c_CHl3(i, j))*Md[k]*Me[j])/2.);
                            CeeSRR.at(i).at(j).at(k).at(l) += -0.25*(oneoMh2*(3*v2*(c_CeH(k, l))*(c_Ye(j, i).conjugate()) + (3*v2*(c_CeH(i, j)) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Ye(j, i).conjugate()))*(c_Ye(l, k).conjugate())));
                            for (int p = 0; p < 3; p++)
                            {
                                CedSRL.at(i).at(j).at(k).at(l) += (-0.5*(oneoMh2*(3*v2*(c_CeH(i, j))*(c_Yd(k, p)) + (3*v2*(c_CdH(p, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(k, p)))*(c_Ye(j, i).conjugate()))))* VCKM(p,l);
                                CedSRR.at(i).at(j).at(k).at(l) += VCKMd(k,p) * (-0.5*(oneoMh2*(3*v2*(c_CeH(i, j))*(c_Yd(l, p).conjugate()) + (3*v2*(c_CdH(p, l)) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(l, p).conjugate()))*(c_Ye(j, i).conjugate()))));
                                for (int r = 0; r < 3; r++)
                                {
                                    CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((oneoMh2*(3*v2*(c_CdH(p, l))*(c_Yd(k, r)) + (3*v2*(c_CdH(r, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(k, r)))*(c_Yd(l, p).conjugate())))/12.) * VCKM(r,j);
                                    CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((oneoMh2*(3*v2*(c_CdH(p, l))*(c_Yd(k, r)) + (3*v2*(c_CdH(r, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(k, r)))*(c_Yd(l, p).conjugate())))/2.) * VCKM(r,j);
                                    CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * (-0.25*(oneoMh2*(3*v2*(c_CdH(r, l))*(c_Yd(j, p).conjugate()) + (3*v2*(c_CdH(p, j)) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(j, p).conjugate()))*(c_Yd(l, r).conjugate()))));
                                }
                            }
                        }
                        for (int p = 0; p < 3; p++)
                        {
                            for (int r = 0; r < 3; r++)
                            {
                                CnudVLL.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Clq1(i, j, p, r)) - (c_Clq3(i, j, p, r))) * VCKM(r, l);
                                CedVLL.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Clq1(i, j, p, r)) + (c_Clq3(i, j, p, r))) * VCKM(r, l);
                                CdeVLR.at(i).at(j).at(k).at(l) += VCKMd(i, p) * ((c_Cqe(p, r, k, l))) * VCKM(r, j);
                                for (int s = 0; s < 3; s++)
                                {
                                    for (int t = 0; t < 3; t++)
                                    {
                                        CddVLL.at(i).at(j).at(k).at(l) += VCKMd(i, p) * VCKMd(k, s) * ((c_Cqq1(p, r, s, t)) + (c_Cqq3(p, r, s, t))) * VCKM(r, j) * VCKM(t, l);
                                    }
                                }
                                CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i, p) * ((c_Cqd1(p, r, k, l))) * VCKM(r, j);
                                CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i, p) * ((c_Cqd8(p, r, k, l))) * VCKM(r, j);
                            }
                            CedSRL.at(i).at(j).at(k).at(l) += ((c_Cledq(i, j, k, p))) * VCKM(p, l);
                        }
                        CedVRR.at(i).at(j).at(k).at(l) += ((c_Ced(i, j, k, l)));
                        CnudVLR.at(i).at(j).at(k).at(l) += ((c_Cld(i, j, k, l)));
                        CedVLR.at(i).at(j).at(k).at(l) += ((c_Cld(i, j, k, l)));
                        CddVRR.at(i).at(j).at(k).at(l) += ((c_Cdd(i, j, k, l)));
                        if (mySMEFT.FlagQuadraticWC)
                        {
                            for (int p = 0; p < 3; p++)
                            {
                                for (int r = 0; r < 3; r++)
                                {
                                    for (int s = 0; s < 3; s++)
                                    {
                                        for (int t = 0; t < 3; t++)
                                        {
                                            CddVLL.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,s) * (-0.125*(gZ2oMZ2*v2*v2*(c_CHq1(s, t))*(c_CHq3(p, r))) - (gZ2oMZ2*v2*v2*(c_CHq3(p, r))*(c_CHq3(s, t)))/8. - (gZ2oMZ2*v2*v2*(c_CHq1(p, r))*((c_CHq1(s, t)) + (c_CHq3(s, t))))/8.) * VCKM(r,j) * VCKM(t,l);
                                        }
                                    }
                                    CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.25*(gZ2oMZ2*v2*v2*(c_CHd(k, l))*((c_CHq1(p, r)) + (c_CHq3(p, r))))) * VCKM(r,j);
                                }
                            }
                            CddVRR.at(i).at(j).at(k).at(l) += ((gZ2oMZ2*v2*v2*(c_CHd(i, j))*(c_CHd(k, l)))/8.);
                            if (mySMEFT.FlagNewTerms)
                            {
                                for (int p = 0; p < 3; p++)
                                {
                                    for (int r = 0; r < 3; r++)
                                    {
                                        for (int s = 0; s < 3; s++)
                                        {
                                            CddVLL.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(s, k).conjugate()) + cbar*(c_CdW(s, k).conjugate()))*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k])/(8.*sqrt(2))) * VCKM(r,j) * VCKM(s,l);
                                            CddVLL.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,s) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(s, l)) + cbar*(c_CdW(s, l)))*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[l])/(8.*sqrt(2))) * VCKM(r,j);
                                            CddVLL.at(i).at(j).at(k).at(l) += VCKMd(k,r) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, i).conjugate()) + cbar*(c_CdW(p, i).conjugate()))*((c_CHq1(r, s)) + (c_CHq3(r, s)))*Md[i])/(8.*sqrt(2))) * VCKM(p,j) * VCKM(s,l);
                                            CddVLL.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, j)) + cbar*(c_CdW(p, j)))*((c_CHq1(r, s)) + (c_CHq3(r, s)))*Md[j])/(8.*sqrt(2))) * VCKM(s,l);
                                        }
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.041666666666666664*(oneoMh2*(3*v2*(c_CdH(p, l))*(3*v2*(c_CdH(r, k).conjugate()))))) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(l, k).conjugate()) + cbar*(c_CdW(l, k).conjugate()))*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[l])/(2.*sqrt(2))) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(k, l)) + cbar*(c_CdW(k, l)))*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k])/(2.*sqrt(2))) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(r, p).conjugate()) + cbar*(c_CdW(r, p).conjugate()))*(c_CHd(k, l))*Md[i])/(2.*sqrt(2))) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, r)) + cbar*(c_CdW(p, r)))*(c_CHd(k, l))*Md[j])/(2.*sqrt(2))) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*v2*((c_CHq1(k, r)) + (c_CHq3(k, r)))*((c_CHq1(p, l)) + (c_CHq3(p, l)))*Md[k]*Md[l])) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ4*v2*v2*(c_CHd(k, r))*((c_CHq1(p, l)) + (c_CHq3(p, l)))*Md[j]*Md[l])/24.) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ4*v2*v2*(c_CHd(p, l))*((c_CHq1(k, r)) + (c_CHq3(k, r)))*Md[i]*Md[k])/24.) * VCKM(r,j);
                                        CddV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.041666666666666664*(gZbar2oMZ4*v2*v2*(c_CHd(p, l))*(c_CHd(k, r))*Md[i]*Md[j])) * VCKM(r,j);
                                        CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.25*(oneoMh2*(3*v2*(c_CdH(p, l))*(3*v2*(c_CdH(r, k).conjugate()))))) * VCKM(r,j);
                                        CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.25*(gZbar2oMZ4*v2*v2*((c_CHq1(k, r)) + (c_CHq3(k, r)))*((c_CHq1(p, l)) + (c_CHq3(p, l)))*Md[k]*Md[l])) * VCKM(r,j);
                                        CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ4*v2*v2*(c_CHd(k, r))*((c_CHq1(p, l)) + (c_CHq3(p, l)))*Md[j]*Md[l])/4.) * VCKM(r,j);
                                        CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ4*v2*v2*(c_CHd(p, l))*((c_CHq1(k, r)) + (c_CHq3(k, r)))*Md[i]*Md[k])/4.) * VCKM(r,j);
                                        CddV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * (-0.25*(gZbar2oMZ4*v2*v2*(c_CHd(p, l))*(c_CHd(k, r))*Md[i]*Md[j])) * VCKM(r,j);
                                        CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * ((oneoMh2*(3*v2*(c_CdH(p, j))*(3*v2*(c_CdH(r, l)))))/8.);
                                        CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * ((gZbar2oMZ4*v2*v2*(c_CHd(r, l))*((c_CHq1(p, j)) + (c_CHq3(p, j)))*Md[j]*Md[k])/4.);
                                        CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * (-0.25*(gZbar2oMZ4*v2*v2*((c_CHq1(p, j)) + (c_CHq3(p, j)))*((c_CHq1(r, l)) + (c_CHq3(r, l)))*Md[j]*Md[l]));
                                        CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * (-0.25*(gZbar2oMZ4*v2*v2*(c_CHd(p, j))*(c_CHd(r, l))*Md[i]*Md[k]));
                                        CddS1RR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * VCKMd(k,r) * ((gZbar2oMZ4*v2*v2*(c_CHd(p, j))*((c_CHq1(r, l)) + (c_CHq3(r, l)))*Md[i]*Md[l])/4.);
                                    }
                                    CddVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, k).conjugate()) + cbar*(c_CdW(p, k).conjugate()))*(c_CHd(i, j))*Md[l])/(4.*sqrt(2)))* VCKM(p,l);
                                    CddVRR.at(i).at(j).at(k).at(l) += VCKMd(k,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, l)) + cbar*(c_CdW(p, l)))*(c_CHd(i, j))*Md[k])/(4.*sqrt(2)));
                                    CddVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, i).conjugate()) + cbar*(c_CdW(p, i).conjugate()))*(c_CHd(k, l))*Md[j])/(4.*sqrt(2))) * VCKM(p,j);
                                    CddVRR.at(i).at(j).at(k).at(l) += VCKMd(i,p) * ((gZbar2oMZ3*v2*v2*(sbar*(c_CdB(p, j)) + cbar*(c_CdW(p, j)))*(c_CHd(k, l))*Md[i])/(4.*sqrt(2)));
                                }
                            }
                        }
                    }
                    for (int l = 0; l < 2; l++)
                    {
                        CnueduVLL.at(i).at(j).at(k).at(l) += VCKMd(k, l) * (-0.5 * (g22oMW2 * v2 * (c_CHl3(i, j))));
                        if (mySMEFT.FlagNewTerms)
                        {
                            CnueduVLL.at(i).at(j).at(k).at(l) += VCKMd(k,l) *((g2bar2oMW3*v2*(c_CeW(i, j))*Me[j])/sqrt(2));
                            CnueduSRR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * (-0.5*(g2bar2oMW4*v2*(c_CHl3(i, j))*Me[j]*Mu[l]));
                        }
                        for (int p = 0; p < 3; p++)
                        {
                            CnueduVLL.at(i).at(j).at(k).at(l) += VCKMd(k, p) * (2 * (c_Clq3(i, j, p, l)));
                            CnueduSRR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Clequ1(i, j, p, l)));
                            CnueduTRR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Clequ3(i, j, p, l)));
                        }
                        CnueduSRL.at(i).at(j).at(k).at(l) += ((c_Cledq(i, j, k, l)));
                    }
                }
                for (int k = 0; k < 2; k++)
                {
                    CnueduVLL.at(i).at(i).at(j).at(k) += VCKMd(j, k) * (-0.5 * (delta_g22oMW2 * g22oMW2));
                    if (mySMEFT.FlagNewTerms)
                    {
                        CnueduSRR.at(i).at(i).at(j).at(k) += VCKMd(j,k) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Me[i]*Mu[k]));
                        CnueduVLL.at(i).at(i).at(j).at(k) += (-((g2bar2oMW3*v2*(c_CdW(k, j).conjugate())*Md[j])/sqrt(2)));
                        for (int p = 0; p < 3; p++)
                        {
                            CnueduVLL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-((g2bar2oMW3*v2*(c_CuW(p, k))*Mu[k])/sqrt(2)));
                            CnueduVLR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-((g2bar2oMW3*v2*(c_CuW(p, k))*Md[j])/sqrt(2)));
                            CnueduSRL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((g2bar2oMW4*v2*(c_CHq3(k, p).conjugate())*Md[j]*Me[i])/2.);
                            CnueduSRR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(k, p).conjugate())*Me[i]*Mu[k]));
                        }
                        CnueduVLR.at(i).at(i).at(j).at(k) += (-((g2bar2oMW3*v2*(c_CdW(k, j).conjugate())*Mu[k])/sqrt(2)));
                        CnueduSRL.at(i).at(i).at(j).at(k) += (-0.25*(g2bar2oMW4*v2*(c_CHud(k, j).conjugate())*Me[i]*Mu[k]));
                        CnueduSRR.at(i).at(i).at(j).at(k) += ((g2bar2oMW4*v2*(c_CHud(k, j).conjugate())*Md[j]*Me[i])/4.);
                    }
                    for (int p = 0; p < 3; p++)
                    {
                        CnueduVLL.at(i).at(i).at(j).at(k) += VCKMd(j, p) * (-0.5 * (g22oMW2 * v2 * (c_CHq3(k, p).conjugate())));
                    }
                    CnueduVLR.at(i).at(i).at(j).at(k) += (-0.25 * (g22oMW2 * v2 * (c_CHud(k, j).conjugate())));
                    for (int l = 0; l < 2; l++)
                    {
                        if (mySMEFT.FlagNewTerms)
                        {
                            CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VCKM(k,j) * (-0.08333333333333333*((1 + deltag2bar2oMW4)*g2bar2oMW4*Mu[k]*Mu[l]));
                            CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * VCKM(k,j) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Mu[k]*Mu[l]));
                            for (int p = 0; p < 3; p++)
                            {
                                CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * (-0.08333333333333333*(g2bar2oMW4*v2*(c_CHq3(k, p))*Mu[k]*Mu[l])) * VCKM(p,j);
                                CduV1LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VCKMd(i,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(c_CHq3(l, p).conjugate())*Mu[k]*Mu[l]));
                                CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(k, p))*Mu[k]*Mu[l])) * VCKM(p,j);
                                CduV8LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * VCKMd(i,p) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(l, p).conjugate())*Mu[k]*Mu[l]));
                            }
                            CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * ((g2bar2oMW4*v2*(c_CHud(k, j))*Md[j]*Mu[l])/24.);
                            CduV1LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * ((g2bar2oMW4*v2*(c_CHud(l, i).conjugate())*Md[i]*Mu[k])/24.);
                            CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i,l) * ((g2bar2oMW4*v2*(c_CHud(k, j))*Md[j]*Mu[l])/4.);
                            CduV8LR.at(i).at(j).at(k).at(l) += VCKM(k,j) * ((g2bar2oMW4*v2*(c_CHud(l, i).conjugate())*Md[i]*Mu[k])/4.);
                            CeuSRL.at(i).at(j).at(k).at(l) += (-0.5*(oneoMh2*(3*v2*(c_CuH(l, k).conjugate())*(c_Ye(j, i).conjugate()) + (3*v2*(c_CeH(i, j)) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Ye(j, i).conjugate()))*(c_Yu(k, l)))));
                            CeuSRR.at(i).at(j).at(k).at(l) += (-0.5*(oneoMh2*(3*v2*(c_CuH(k, l))*(c_Ye(j, i).conjugate()) + (3*v2*(c_CeH(i, j)) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Ye(j, i).conjugate()))*(c_Yu(l, k).conjugate()))));
                        }
                        CnuuVLL.at(i).at(j).at(k).at(l) += ((c_Clq1(i, j, k, l)) + (c_Clq3(i, j, k, l)));
                        CeuVLL.at(i).at(j).at(k).at(l) += ((c_Clq1(i, j, k, l)) - (c_Clq3(i, j, k, l)));
                        CeuVRR.at(i).at(j).at(k).at(l) += ((c_Ceu(i, j, k, l)));
                        CnuuVLR.at(i).at(j).at(k).at(l) += ((c_Clu(i, j, k, l)));
                        CeuVLR.at(i).at(j).at(k).at(l) += ((c_Clu(i, j, k, l)));
                        CeuSRR.at(i).at(j).at(k).at(l) += (-(c_Clequ1(i, j, k, l)));
                        CeuTRR.at(i).at(j).at(k).at(l) += (-(c_Clequ3(i, j, k, l)));
                        for (int p = 0; p < 3; p++)
                        {
                            for (int r = 0; r < 3; r++)
                            {
                                CduV1LR.at(i).at(j).at(k).at(l) += VCKMd(i, p) * ((c_Cqu1(p, r, k, l))) * VCKM(r, j);
                                CduV8LR.at(i).at(j).at(k).at(l) += VCKMd(i, p) * ((c_Cqu8(p, r, k, l))) * VCKM(r, j);
                            }
                        }
                    }
                }
                for (int p = 0; p < 3; p++)
                {
                    Cdg.at(i).at(j) += vTosq2 * VCKMd(i, p) * (-(c_CdW(p, j)) * sbar + (c_CdB(p, j)) * cbar);
                    CdG.at(i).at(j) += vTosq2 * VCKMd(i, p) * (c_CdG(p, j));
                }
            }
            for (int j = 0; j < 2; j++)
            {
                CnuuVLL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (8 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 4 * sbar2))) / 12.;
                CeuVLL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (4 * delta_sbar * sbar2 * (-5 + 8 * sbar2) + delta_gZ2oMZ2 * (3 - 10 * sbar2 + 8 * sbar2 * sbar2))) / 12.;
                CeuVRR.at(i).at(i).at(j).at(j) += (2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 3.;
                CnuuVLR.at(i).at(i).at(j).at(j) += ((delta_gZ2oMZ2 + 2 * delta_sbar) * gZ2oMZ2 * sbar2) / 3.;
                CeuVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (-delta_gZ2oMZ2 - 2 * delta_sbar + 2 * delta_gZ2oMZ2 * sbar2 + 8 * delta_sbar * sbar2)) / 3.;
                CduV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (-3 * delta_gZ2oMZ2 - 6 * delta_sbar + 2 * delta_gZ2oMZ2 * sbar2 + 8 * delta_sbar * sbar2)) / 9.;
                if (mySMEFT.FlagNewTerms)
                {
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
                {
                    CudV1LL.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))) / 12.);
                    CnuuVLL.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))) / 4.);
                    CeuVLL.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))) / 4.);
                    CudV1RR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHu(j, k))) / 6.);
                    CeuVRR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * sbar2 * v2 * (c_CHu(j, k))) / 2.);
                    CnuuVLR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * v2 * (c_CHu(j, k))) / 4.);
                    CeuVLR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (-1 + 2 * sbar2) * v2 * (c_CHu(j, k))) / 4.);
                    CduV1LR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (-3 + 2 * sbar2) * v2 * (c_CHu(j, k))) / 12.);
                    CudV1LR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))) / 6.);
                    CueVLR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))) / 2.);
                    if (mySMEFT.FlagNewTerms)
                    {
                        CnuuVLL.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/sqrt(2));
                        CnuuVLL.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/sqrt(2));
                        CeuVLL.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/sqrt(2)));
                        CeuVLL.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/sqrt(2)));
                        CudV1LL.at(j).at(k).at(i).at(i) += (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/sqrt(2));
                        CudV1LL.at(j).at(k).at(i).at(i) += (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/sqrt(2));
                        CeuVRR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/sqrt(2)));
                        CeuVRR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/sqrt(2)));
                        CudV1RR.at(j).at(k).at(i).at(i) += (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/sqrt(2));
                        CudV1RR.at(j).at(k).at(i).at(i) += (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/sqrt(2));
                        CnuuVLR.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/sqrt(2));
                        CnuuVLR.at(i).at(i).at(j).at(k) += (-0.5*(gZbar2oMZ3*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/sqrt(2));
                        CeuVLR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/sqrt(2)));
                        CeuVLR.at(i).at(i).at(j).at(k) += (-((gZbar2oMZ3*(-0.5 + sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/sqrt(2)));
                        CueVLR.at(j).at(k).at(i).at(i) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/sqrt(2)));
                        CueVLR.at(j).at(k).at(i).at(i) += (-((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/sqrt(2)));
                        CudV1LR.at(j).at(k).at(i).at(i) += (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/sqrt(2));
                        CudV1LR.at(j).at(k).at(i).at(i) += (-0.3333333333333333*(gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/sqrt(2));
                        CduV1LR.at(i).at(i).at(j).at(k) += (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/sqrt(2));
                        CduV1LR.at(i).at(i).at(j).at(k) += (-0.16666666666666666*(gZbar2oMZ3*(-3 + 2*sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/sqrt(2));
                        CudduV1LR.at(j).at(i).at(i).at(k) += (-0.041666666666666664*(gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Md[i]*Mu[k]));
                        CudduV1LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*(c_CHu(j, k))*Md[i]*Mu[j])/24.);
                        CudduV8LR.at(j).at(i).at(i).at(k) += (-0.25*(gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Md[i]*Mu[k]));
                        CudduV8LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*(c_CHu(j, k))*Md[i]*Mu[j])/4.);
                        CeuSRL.at(i).at(i).at(j).at(k) += ((gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Me[i]*Mu[j])/4.);
                        CeuSRL.at(i).at(i).at(j).at(k) += (-0.25*(gZbar2oMZ4*v2*(c_CHu(j, k))*Me[i]*Mu[k]));
                        CeuSRR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ4*v2*(c_CHu(j, k))*Me[i]*Mu[j])/4.);
                        CeuSRR.at(i).at(i).at(j).at(k) += (-0.25*(gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Me[i]*Mu[k]));
                        CudS1RR.at(j).at(k).at(i).at(i) += (-0.25*(gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Md[i]*Mu[k]));
                        CudS1RR.at(j).at(k).at(i).at(i) += ((gZbar2oMZ4*v2*(c_CHu(j, k))*Md[i]*Mu[j])/4.);
                    }
                }
            }
        }
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                CudV1LL.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (4 * delta_sbar * sbar2 * (-9 + 8 * sbar2) + delta_gZ2oMZ2 * (9 - 18 * sbar2 + 8 * sbar2 * sbar2))) / 36.;
                CudV1RR.at(i).at(i).at(j).at(j) += (2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 9.;
                CudV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (delta_gZ2oMZ2 * (-3 + 4 * sbar2) + 2 * delta_sbar * (-3 + 8 * sbar2))) / 18.;
                CueVLR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * sbar2 * (delta_gZ2oMZ2 * (-3 + 4 * sbar2) + 2 * delta_sbar * (-3 + 8 * sbar2))) / 6.;
                if (mySMEFT.FlagNewTerms)
                {
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
                    CeuVLL.at(j).at(k).at(i).at(i) += -0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k))));
                    CnuuVLL.at(j).at(k).at(i).at(i) += -0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((c_CHl1(j, k)) - (c_CHl3(j, k))));
                    CeuVRR.at(j).at(k).at(i).at(i) += -0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * (c_CHe(j, k)));
                    CnuuVLR.at(j).at(k).at(i).at(i) += (gZ2oMZ2 * sbar2 * v2 * (-(c_CHl1(j, k)) + (c_CHl3(j, k)))) / 3.;
                    CeuVLR.at(j).at(k).at(i).at(i) += -0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * ((c_CHl1(j, k)) + (c_CHl3(j, k))));
                    CueVLR.at(i).at(i).at(j).at(k) += (gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (c_CHe(j, k))) / 12.;
                    if (mySMEFT.FlagNewTerms)
                    {
                        CeuVLL.at(j).at(k).at(i).at(i) += (gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/(6.*sqrt(2));
                        CeuVRR.at(j).at(k).at(i).at(i) += (sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/3.;
                        CeuVLR.at(j).at(k).at(i).at(i) += (sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CeB(k, j).conjugate())*Me[j] + cbar*(c_CeW(k, j).conjugate())*Me[j] + (sbar*(c_CeB(j, k)) + cbar*(c_CeW(j, k)))*Me[k]))/3.;
                        CueVLR.at(i).at(i).at(j).at(k) += (gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CeB(j, k))*Me[j] + cbar*(c_CeW(j, k))*Me[j] + (sbar*(c_CeB(k, j).conjugate()) + cbar*(c_CeW(k, j).conjugate()))*Me[k]))/(6.*sqrt(2));
                        CeuSRL.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*((c_CHe(j, k))*Me[j] - ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k])*Mu[i])/4.;
                        CeuSRR.at(j).at(k).at(i).at(i) += (gZbar2oMZ4*v2*(-((c_CHe(j, k))*Me[j]) + ((c_CHl1(j, k)) + (c_CHl3(j, k)))*Me[k])*Mu[i])/4.;
                        for (int p = 0; p < 3; p++)
                        {
                            CudV1LL.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/(6.*sqrt(2))) * VCKM(p,k);
                            CudV1LL.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/(6.*sqrt(2)));
                            CudV1RR.at(i).at(i).at(j).at(k) += ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/3.) * VCKM(p,k);
                            CudV1RR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/3.);
                            CudV1LR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[k])/(6.*sqrt(2))) * VCKM(p,k);
                            CudV1LR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[j])/(6.*sqrt(2)));
                            CduV1LR.at(j).at(k).at(i).at(i) += ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, j).conjugate()) + cbar*(c_CdW(p, j).conjugate()))*Md[j])/3.) * VCKM(p,k);
                            CduV1LR.at(j).at(k).at(i).at(i) += VCKMd(j,p) * ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CdB(p, k)) + cbar*(c_CdW(p, k)))*Md[k])/3.);
                            for (int r = 0; r < 3; r++)
                            {
                                CudduV1LR.at(i).at(j).at(k).at(i) += VCKMd(k,p) * ((gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k]*Mu[i])/24.) * VCKM(r,j);
                                CudduV8LR.at(i).at(j).at(k).at(i) += VCKMd(k,p) * ((gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k]*Mu[i])/4.) * VCKM(r,j);
                                CudS1RR.at(i).at(i).at(j).at(k) += VCKMd(j,p) * ((gZbar2oMZ4*v2*((c_CHq1(p, r)) + (c_CHq3(p, r)))*Md[k]*Mu[i])/4.) * VCKM(r,k);
                            }
                        }
                        CudduV1LR.at(i).at(j).at(k).at(i) += (-0.041666666666666664*(gZbar2oMZ4*v2*(c_CHd(k, j))*Md[j]*Mu[i]));
                        CudduV8LR.at(i).at(j).at(k).at(i) += (-0.25*(gZbar2oMZ4*v2*(c_CHd(k, j))*Md[j]*Mu[i]));
                        CudS1RR.at(i).at(i).at(j).at(k) += (-0.25*(gZbar2oMZ4*v2*(c_CHd(j, k))*Md[j]*Mu[i]));
                    }
                    for (int p = 0; p < 3; p++)
                    {
                        for (int r = 0; r < 3; r++)
                        {
                            CudV1LL.at(i).at(i).at(j).at(k) += VCKMd(j, p) * (-0.08333333333333333 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r))))) * VCKM(r, k);
                            CduV1LR.at(j).at(k).at(i).at(i) += VCKMd(j, p) * (-0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * ((c_CHq1(p, r)) + (c_CHq3(p, r))))) * VCKM(r, k);
                        }
                    }
                    CudV1RR.at(i).at(i).at(j).at(k) += (-0.3333333333333333 * (gZ2oMZ2 * sbar2 * v2 * (c_CHd(j, k))));
                    CudV1LR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (c_CHd(j, k))) / 12.);
                    for (int l = 0; l < 2; l++)
                    {
                        if (mySMEFT.FlagNewTerms)
                        {
                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VCKMd(k,l) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[j]*Mu[l]));
                            CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i,j) * (-((g2bar2oMW3*v2*(c_CdW(l, k).conjugate())*Mu[l])/sqrt(2)));
                            for (int p = 0; p < 3; p++)
                            {
                                CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VCKMd(k,p) * (-((g2bar2oMW3*v2*(c_CuW(p, l))*Md[k])/sqrt(2)));
                                CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * VCKMd(k,p) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(l, p).conjugate())*Md[j]*Mu[l]));
                                CudduS1RR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(i, p))*Md[j]*Mu[l])) * VCKM(p,j);
                                CudduV1LR.at(i).at(j).at(k).at(l) += ((oneoMh2*(3*v2*(c_CuH(i, l))*(c_Yd(k, p)) + (3*v2*(c_CdH(p, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(k, p)))*(c_Yu(l, i).conjugate())))/12.) * VCKM(p,j);
                                CudduV8LR.at(i).at(j).at(k).at(l) += ((oneoMh2*(3*v2*(c_CuH(i, l))*(c_Yd(k, p)) + (3*v2*(c_CdH(p, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(k, p)))*(c_Yu(l, i).conjugate())))/2.) * VCKM(p,j);
                            }
                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKM(i,j) * ((g2bar2oMW4*v2*(c_CHud(l, k).conjugate())*Md[j]*Md[k])/4.);
                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKMd(k,l) * ((g2bar2oMW4*v2*(c_CHud(i, j))*Mu[i]*Mu[l])/4.);
                        }
                        CudduV1LR.at(i).at(j).at(k).at(l) += VCKM(i, j) * (-0.25 * (g22oMW2 * v2 * (c_CHud(l, k).conjugate())));
                        for (int p = 0; p < 3; p++)
                        {
                            CudduS1RR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * (-(c_Cquqd1(p, l, i, j)));
                            CudduS8RR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * (-(c_Cquqd8(p, l, i, j)));
                        }
                    }
                }
            }
            for (int j = 0; j < 2; j++)
            {
                CuuVLL.at(i).at(i).at(j).at(j) += -0.013888888888888888 * (gZ2oMZ2 * (-3 + 4 * sbar2) * (16 * delta_sbar * sbar2 + delta_gZ2oMZ2 * (-3 + 4 * sbar2)));
                CuuVRR.at(i).at(i).at(j).at(j) += (-2 * (delta_gZ2oMZ2 + 4 * delta_sbar) * gZ2oMZ2 * sbar2 * sbar2) / 9.;
                CuuV1LR.at(i).at(i).at(j).at(j) += (gZ2oMZ2 * (2 * delta_sbar * (3 - 8 * sbar2) + delta_gZ2oMZ2 * (3 - 4 * sbar2)) * sbar2) / 9.;
                if (mySMEFT.FlagNewTerms)
                {
                    CuuV1LR.at(i).at(j).at(j).at(i) += -0.004629629629629629*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                    CuuV1LR.at(i).at(j).at(j).at(i) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/54.);
                    CuuV1LR.at(i).at(j).at(j).at(i) += (-2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/27.;
                    CuuV8LR.at(i).at(j).at(j).at(i) += -0.027777777777777776*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                    CuuV8LR.at(i).at(j).at(j).at(i) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/9.);
                    CuuV8LR.at(i).at(j).at(j).at(i) += (-4*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/9.;
                    CuuS1RR.at(i).at(i).at(j).at(j) += 2. * ((gZbar2oMZ4*sbar2*(-3 + 4*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2) + 2*delta_sbar*(-3 + 8*sbar2))*Mu[i]*Mu[j])/18.);
                    CuuS1RR.at(i).at(i).at(j).at(j) += -0.013888888888888888*(gZbar2oMZ4*(-3 + 4*sbar2)*(-3 + 4*(1 + 4*delta_sbar)*sbar2 + deltagZbar2oMZ4*(-3 + 4*sbar2))*Mu[i]*Mu[j]);
                    CuuS1RR.at(i).at(i).at(j).at(j) += (-2*(1 + deltagZbar2oMZ4 + 4*delta_sbar)*gZbar2oMZ4*pow(sbar2,2)*Mu[i]*Mu[j])/9.;
                }
                for (int k = 0; k < 3; k++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, j) * (-0.16666666666666666 * (delta_g22oMW2 * g22oMW2));
                        CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, j) * (-(delta_g22oMW2 * g22oMW2));
                        if (mySMEFT.FlagNewTerms)
                        {
                            CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,j) * (-0.08333333333333333*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[k]*Md[l]));
                            CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,j) * (-0.5*((1 + deltag2bar2oMW4)*g2bar2oMW4*Md[k]*Md[l]));
                            for (int p = 0; p < 3; p++)
                            {
                                CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * ((g2bar2oMW3*v2*(c_CuW(p, i).conjugate())*Mu[i])/(3.*sqrt(2))) * VCKM(p,l);
                                CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,p) * (-0.3333333333333333*(g2bar2oMW3*v2*(c_CuW(p, j))*Mu[j])/sqrt(2));
                                CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * (sqrt(2)*g2bar2oMW3*v2*(c_CuW(p, i).conjugate())*Mu[i]) * VCKM(p,l);
                                CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,p) * (-(sqrt(2)*g2bar2oMW3*v2*(c_CuW(p, j))*Mu[j]));
                                CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,p) * (-0.08333333333333333*(g2bar2oMW4*v2*(c_CHq3(j, p).conjugate())*Md[k]*Md[l]));
                                CudV1LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * (-0.08333333333333333*(g2bar2oMW4*v2*(c_CHq3(i, p))*Md[k]*Md[l])) * VCKM(p,l);
                                CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * VCKMd(k,p) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(j, p).conjugate())*Md[k]*Md[l]));
                                CudV8LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * (-0.5*(g2bar2oMW4*v2*(c_CHq3(i, p))*Md[k]*Md[l])) * VCKM(p,l);
                                CudS1RR.at(i).at(j).at(k).at(l) += VCKMd(k,p) * ((oneoMh2*(-3*v2*(c_CuH(i, j))*(c_Yd(l, p).conjugate()) + (-3*v2*(c_CdH(p, l)) + 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yd(l, p).conjugate()))*(c_Yu(j, i).conjugate())))/2.);
                            }
                            CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * ((g2bar2oMW3*v2*(c_CdW(i, l))*Md[l])/(3.*sqrt(2)));
                            CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * (-0.3333333333333333*(g2bar2oMW3*v2*(c_CdW(j, k).conjugate())*Md[k])/sqrt(2));
                            CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k,j) * (sqrt(2)*g2bar2oMW3*v2*(c_CdW(i, l))*Md[l]);
                            CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i,l) * (-(sqrt(2)*g2bar2oMW3*v2*(c_CdW(j, k).conjugate())*Md[k]));
                            CudV1LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * ((g2bar2oMW4*v2*(c_CHud(j, k).conjugate())*Md[l]*Mu[j])/24.);
                            CudV1LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * ((g2bar2oMW4*v2*(c_CHud(i, l))*Md[k]*Mu[i])/24.);
                            CudV8LR.at(i).at(j).at(k).at(l) += VCKM(i,l) * ((g2bar2oMW4*v2*(c_CHud(j, k).conjugate())*Md[l]*Mu[j])/4.);
                            CudV8LR.at(i).at(j).at(k).at(l) += VCKMd(k,j) * ((g2bar2oMW4*v2*(c_CHud(i, l))*Md[k]*Mu[i])/4.);
                        }
                        for (int p = 0; p < 3; p++)
                        {
                            CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k, j) * (-0.16666666666666666 * (g22oMW2 * v2 * (c_CHq3(i, p)))) * VCKM(p, l);
                            CudV1LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, p) * (-0.16666666666666666 * (g22oMW2 * v2 * (c_CHq3(j, p).conjugate())));
                            CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k, j) * (-(g22oMW2 * v2 * (c_CHq3(i, p)))) * VCKM(p, l);
                            CudV8LL.at(i).at(j).at(k).at(l) += VCKM(i, l) * VCKMd(k, p) * (-(g22oMW2 * v2 * (c_CHq3(j, p).conjugate())));
                            for (int r = 0; r < 3; r++)
                            {
                                CudV1LL.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Cqq1(i, j, p, r)) + (c_Cqq1(p, r, i, j)) - (c_Cqq3(i, j, p, r)) + (2 * (c_Cqq3(i, r, p, j))) / 3. + (2 * (c_Cqq3(p, j, i, r))) / 3. - (c_Cqq3(p, r, i, j))) * VCKM(r, l);
                                CudV8LL.at(i).at(j).at(k).at(l) += VCKMd(k, p) * (4 * ((c_Cqq3(i, r, p, j)) + (c_Cqq3(p, j, i, r)))) * VCKM(r, l);
                            }
                            CudS1RR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Cquqd1(i, j, p, l)));
                            CudS8RR.at(i).at(j).at(k).at(l) += VCKMd(k, p) * ((c_Cquqd8(i, j, p, l)));
                        }
                        CueVLR.at(i).at(j).at(k).at(l) += ((c_Cqe(i, j, k, l)));
                        CudV1RR.at(i).at(j).at(k).at(l) += ((c_Cud1(i, j, k, l)));
                        CudV8RR.at(i).at(j).at(k).at(l) += ((c_Cud8(i, j, k, l)));
                        CudV1LR.at(i).at(j).at(k).at(l) += ((c_Cqd1(i, j, k, l)));
                        CudV8LR.at(i).at(j).at(k).at(l) += ((c_Cqd8(i, j, k, l)));
                    }
                }
                for (int k = 0; k < 2; k++)
                {
                    CuuVLL.at(i).at(i).at(j).at(k) += (-0.041666666666666664 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))));
                    CuuVLL.at(j).at(k).at(i).at(i) += (-0.041666666666666664 * (gZ2oMZ2 * (-3 + 4 * sbar2) * v2 * ((c_CHq1(j, k)) - (c_CHq3(j, k)))));
                    CuuVRR.at(i).at(i).at(j).at(k) += (-0.16666666666666666 * (gZ2oMZ2 * sbar2 * v2 * (c_CHu(j, k))));
                    CuuVRR.at(j).at(k).at(i).at(i) += (-0.16666666666666666 * (gZ2oMZ2 * sbar2 * v2 * (c_CHu(j, k))));
                    CuuV1LR.at(i).at(i).at(j).at(k) += ((gZ2oMZ2 * (3 - 4 * sbar2) * v2 * (c_CHu(j, k))) / 12.);
                    CuuV1LR.at(j).at(k).at(i).at(i) += ((gZ2oMZ2 * sbar2 * v2 * (-(c_CHq1(j, k)) + (c_CHq3(j, k)))) / 3.);
                    if (mySMEFT.FlagNewTerms)
                    {
                        CuuVLL.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/(12.*sqrt(2)));
                        CuuVLL.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/(12.*sqrt(2)));
                        CuuVLL.at(j).at(k).at(i).at(i) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/(12.*sqrt(2)));
                        CuuVLL.at(j).at(k).at(i).at(i) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/(12.*sqrt(2)));
                        CuuVRR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/(3.*sqrt(2)));
                        CuuVRR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/(3.*sqrt(2)));
                        CuuVRR.at(j).at(k).at(i).at(i) += ((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/(3.*sqrt(2)));
                        CuuVRR.at(j).at(k).at(i).at(i) += ((gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/(3.*sqrt(2)));
                        CuuV1LR.at(i).at(j).at(k).at(i) += ((gZbar2oMZ4*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*Mu[i]*Mu[k])/24.);
                        CuuV1LR.at(i).at(j).at(k).at(i) += (-0.041666666666666664*(gZbar2oMZ4*v2*(c_CHu(k, j))*Mu[i]*Mu[j]));
                        CuuV1LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Mu[i]*Mu[k])/24.);
                        CuuV1LR.at(j).at(i).at(i).at(k) += (-0.041666666666666664*(gZbar2oMZ4*v2*(c_CHu(j, k))*Mu[i]*Mu[j]));
                        CuuV1LR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[k])/(6.*sqrt(2)));
                        CuuV1LR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ3*(-3 + 4*sbar2)*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[j])/(6.*sqrt(2)));
                        CuuV1LR.at(j).at(k).at(i).at(i) += ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(k, j).conjugate()) + cbar*(c_CuW(k, j).conjugate()))*Mu[j])/3.);
                        CuuV1LR.at(j).at(k).at(i).at(i) += ((sqrt(2)*gZbar2oMZ3*sbar2*v2*(sbar*(c_CuB(j, k)) + cbar*(c_CuW(j, k)))*Mu[k])/3.);
                        CuuV8LR.at(i).at(j).at(k).at(i) += ((gZbar2oMZ4*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*Mu[i]*Mu[k])/4.);
                        CuuV8LR.at(i).at(j).at(k).at(i) += (-0.25*(gZbar2oMZ4*v2*(c_CHu(k, j))*Mu[i]*Mu[j]));
                        CuuV8LR.at(j).at(i).at(i).at(k) += ((gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Mu[i]*Mu[k])/4.);
                        CuuV8LR.at(j).at(i).at(i).at(k) += (-0.25*(gZbar2oMZ4*v2*(c_CHu(j, k))*Mu[i]*Mu[j]));
                        CuuS1RR.at(i).at(i).at(j).at(k) += (-0.125*(gZbar2oMZ4*v2*(c_CHu(j, k))*Mu[i]*Mu[j]));
                        CuuS1RR.at(i).at(i).at(j).at(k) += ((gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Mu[i]*Mu[k])/8.);
                        CuuS1RR.at(j).at(k).at(i).at(i) += ((gZbar2oMZ4*v2*((c_CHq1(j, k)) - (c_CHq3(j, k)))*Mu[i]*Mu[k])/8.);
                        CuuS1RR.at(j).at(k).at(i).at(i) += (-0.125*(gZbar2oMZ4*v2*(c_CHu(j, k))*Mu[i]*Mu[j]));
                    }
                    for (int l = 0; l < 2; l++)
                    {
                        CuuVLL.at(i).at(j).at(k).at(l) += ((c_Cqq1(i, j, k, l)) + (c_Cqq3(i, j, k, l)));
                        CuuVRR.at(i).at(j).at(k).at(l) += ((c_Cuu(i, j, k, l)));
                        CuuV1LR.at(i).at(j).at(k).at(l) += ((c_Cqu1(i, j, k, l)));
                        CuuV8LR.at(i).at(j).at(k).at(l) += ((c_Cqu8(i, j, k, l)));
                        if (mySMEFT.FlagQuadraticWC)
                        {
                            CuuVLL.at(i).at(j).at(k).at(l) += ((gZ2oMZ2*v2*v2*(c_CHq1(k, l))*(c_CHq3(i, j)) - gZ2oMZ2*v2*v2*(c_CHq3(i, j))*(c_CHq3(k, l)) + gZ2oMZ2*v2*v2*(c_CHq1(i, j))*(-(c_CHq1(k, l)) + (c_CHq3(k, l))))/8.);
                            CuuVRR.at(i).at(j).at(k).at(l) += (-0.125*(gZ2oMZ2*v2*v2*(c_CHu(i, j))*(c_CHu(k, l))));
                            CuuV1LR.at(i).at(j).at(k).at(l) += (-0.25*(gZ2oMZ2*v2*v2*(c_CHq1(i, j))*(c_CHu(k, l))) + (gZ2oMZ2*v2*v2*(c_CHq3(i, j))*(c_CHu(k, l)))/4.);
                        }
                        if (mySMEFT.FlagNewTerms)
                        {
                            CuuV1LR.at(i).at(j).at(k).at(l) += ((oneoMh2*(3*v2*(c_CuH(i, l))*(c_Yu(k, j)) + (3*v2*(c_CuH(j, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yu(k, j)))*(c_Yu(l, i).conjugate())))/12.);
                            CuuV8LR.at(i).at(j).at(k).at(l) += ((oneoMh2*(3*v2*(c_CuH(i, l))*(c_Yu(k, j)) + (3*v2*(c_CuH(j, k).conjugate()) - 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yu(k, j)))*(c_Yu(l, i).conjugate())))/2.);
                            CuuS1RR.at(i).at(j).at(k).at(l) += ((oneoMh2*(-3*v2*(c_CuH(k, l))*(c_Yu(j, i).conjugate()) + (-3*v2*(c_CuH(i, j)) + 2*(1 + 2*(c_CHbox - 0.25 * c_CHD) * v2 + deltaoneoMh2)*(c_Yu(j, i).conjugate()))*(c_Yu(l, k).conjugate())))/4.);
                            if (mySMEFT.FlagQuadraticWC)
                            {
                                CuuVLL.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*(sbar*(c_CuB(l, k).conjugate()) + cbar*(c_CuW(l, k).conjugate()))*Mu[k])/(4.*sqrt(2)));
                                CuuVLL.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*(sbar*(c_CuB(k, l)) + cbar*(c_CuW(k, l)))*Mu[l])/(4.*sqrt(2)));
                                CuuVLL.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(k, l)) - (c_CHq3(k, l)))*(sbar*(c_CuB(j, i).conjugate()) + cbar*(c_CuW(j, i).conjugate()))*Mu[i])/(4.*sqrt(2)));
                                CuuVLL.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(k, l)) - (c_CHq3(k, l)))*(sbar*(c_CuB(i, j)) + cbar*(c_CuW(i, j)))*Mu[j])/(4.*sqrt(2)));
                                CuuVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(i, j))*(sbar*(c_CuB(l, k).conjugate()) + cbar*(c_CuW(l, k).conjugate()))*Mu[l])/(4.*sqrt(2)));
                                CuuVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(i, j))*(sbar*(c_CuB(k, l)) + cbar*(c_CuW(k, l)))*Mu[k])/(4.*sqrt(2)));
                                CuuVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(k, l))*(sbar*(c_CuB(j, i).conjugate()) + cbar*(c_CuW(j, i).conjugate()))*Mu[j])/(4.*sqrt(2)));
                                CuuVRR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(k, l))*(sbar*(c_CuB(i, j)) + cbar*(c_CuW(i, j)))*Mu[i])/(4.*sqrt(2)));
                                CuuV1LR.at(i).at(j).at(k).at(l) += (-0.041666666666666664*(oneoMh2*(3*v2*(c_CuH(i, l))*(3*v2*(c_CuH(j, k).conjugate())))));
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*(sbar*(c_CuB(l, k).conjugate()) + cbar*(c_CuW(l, k).conjugate()))*Mu[l])/(2.*sqrt(2)));
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*(sbar*(c_CuB(k, l)) + cbar*(c_CuW(k, l)))*Mu[k])/(2.*sqrt(2)));
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(k, l))*(sbar*(c_CuB(j, i).conjugate()) + cbar*(c_CuW(j, i).conjugate()))*Mu[i])/(2.*sqrt(2)));
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ3*v2*v2*(c_CHu(k, l))*(sbar*(c_CuB(i, j)) + cbar*(c_CuW(i, j)))*Mu[j])/(2.*sqrt(2)));
                                CuuV1LR.at(i).at(j).at(k).at(l) += (-0.041666666666666664*(gZbar2oMZ4*v2*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*((c_CHq1(i, l)) - (c_CHq3(i, l)))*Mu[k]*Mu[l]));
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(i, l)) - (c_CHq3(i, l)))*(c_CHu(k, j))*Mu[j]*Mu[l])/24.);
                                CuuV1LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*(c_CHu(i, l))*Mu[i]*Mu[k])/24.);
                                CuuV1LR.at(i).at(j).at(k).at(l) += (-0.041666666666666664*(gZbar2oMZ4*v2*v2*(c_CHu(i, l))*(c_CHu(k, j))*Mu[i]*Mu[j]));
                                CuuV8LR.at(i).at(j).at(k).at(l) += (-0.25*(oneoMh2*(3*v2*(c_CuH(i, l))*(3*v2*(c_CuH(j, k).conjugate())))));
                                CuuV8LR.at(i).at(j).at(k).at(l) += (-0.25*(gZbar2oMZ4*v2*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*((c_CHq1(i, l)) - (c_CHq3(i, l)))*Mu[k]*Mu[l]));
                                CuuV8LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(i, l)) - (c_CHq3(i, l)))*(c_CHu(k, j))*Mu[j]*Mu[l])/4.);
                                CuuV8LR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(k, j)) - (c_CHq3(k, j)))*(c_CHu(i, l))*Mu[i]*Mu[k])/4.);
                                CuuV8LR.at(i).at(j).at(k).at(l) += (-0.25*(gZbar2oMZ4*v2*v2*(c_CHu(i, l))*(c_CHu(k, j))*Mu[i]*Mu[j]));
                                CuuS1RR.at(i).at(j).at(k).at(l) += ((oneoMh2*(3*v2*(c_CuH(i, j))*(3*v2*(c_CuH(k, l)))))/8.);
                                CuuS1RR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*(c_CHu(k, l))*Mu[j]*Mu[k])/8.);
                                CuuS1RR.at(i).at(j).at(k).at(l) += (-0.125*(gZbar2oMZ4*v2*v2*((c_CHq1(i, j)) - (c_CHq3(i, j)))*((c_CHq1(k, l)) - (c_CHq3(k, l)))*Mu[j]*Mu[l]));
                                CuuS1RR.at(i).at(j).at(k).at(l) += (-0.125*(gZbar2oMZ4*v2*v2*(c_CHu(i, j))*(c_CHu(k, l))*Mu[i]*Mu[k]));
                                CuuS1RR.at(i).at(j).at(k).at(l) += ((gZbar2oMZ4*v2*v2*((c_CHq1(k, l)) - (c_CHq3(k, l)))*(c_CHu(i, j))*Mu[i]*Mu[l])/8.);
                            }
                        }
                    }
                }
                Cug.at(i).at(j) += vTosq2 * ((c_CuW(i, j)) * sbar + (c_CuB(i, j)) * cbar);
                CuG.at(i).at(j) += vTosq2 * (c_CuG(i, j));
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
        //SUSY basis for Delta B=2 Hamiltonian: we are using bbar d bbar d operators

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
                mcbd.setCoeff(0, -CddVLL.at(2).at(0).at(2).at(0), LO);
                mcbd.setCoeff(1, -(CddS1RR.at(2).at(0).at(2).at(0) - CddS8RR.at(2).at(0).at(2).at(0) / 6.), LO);
                mcbd.setCoeff(2, -CddS8RR.at(2).at(0).at(2).at(0) / 2., LO);
                mcbd.setCoeff(3, CddV8LR.at(2).at(0).at(2).at(0), LO);
                mcbd.setCoeff(4, 2. * CddV1LR.at(2).at(0).at(2).at(0) - CddV8LR.at(2).at(0).at(2).at(0) / 3., LO);
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
                mcbd.setCoeff(0, -CddVRR.at(2).at(0).at(2).at(0), LO);
                mcbd.setCoeff(1, -(CddS1RR.at(2).at(0).at(2).at(0) - CddS8RR.at(2).at(0).at(2).at(0) / 6.), LO);
                mcbd.setCoeff(2, -CddS8RR.at(2).at(0).at(2).at(0) / 2., LO);
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
        //SUSY basis for Delta B=2 Hamiltonian: we are using bbar s bbar s operators

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
                mcbs.setCoeff(0, -CddVLL.at(2).at(1).at(2).at(1), LO);
                mcbs.setCoeff(1, -(CddS1RR.at(2).at(1).at(2).at(1) - CddS8RR.at(2).at(1).at(2).at(1) / 6.), LO);
                mcbs.setCoeff(2, -CddS8RR.at(2).at(1).at(2).at(1) / 2., LO);
                mcbs.setCoeff(3, CddV8LR.at(2).at(1).at(2).at(1), LO);
                mcbs.setCoeff(4, 2. * CddV1LR.at(2).at(1).at(2).at(1) - CddV8LR.at(2).at(1).at(2).at(1) / 3., LO);
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
                mcbs.setCoeff(0, -CddVRR.at(2).at(1).at(2).at(1), LO);
                mcbs.setCoeff(1, -(CddS1RR.at(2).at(1).at(2).at(1) - CddS8RR.at(2).at(1).at(2).at(1) / 6.), LO);
                mcbs.setCoeff(2, - CddS8RR.at(2).at(1).at(2).at(1) / 2., LO);
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
    gslpp::complex LEFT_factor_radiative = 16. * M_PI * M_PI / mySMEFT.Mrun(mySMEFT.getMuw(),mySMEFT.getQuarks(QCD::BOTTOM).getMass(),QCD::BOTTOM) * LEFT_factor / sqrt(4. * M_PI * mySMEFT.getAle());
    
    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
        //  {O1Misiak, O2Misiak} = {{sbarTac cbarTab}, {sbarc cbar b}} = {{-(1/N), (-1 + N^2)/(2 N^2)}, {2, 1/N}} {OV8LLud_Fierzed,OV1LLud_Fierzed}
        //  C_JMS^T O_JMS = C_JMS^T M O_Misiak = C_Misiak^T O_Misiak
        //  C_Misiak = M^T C_JMS = {{-(1/N), 2}, {(-1 + N^2)/(2 N^2), 1/N}} {CudV8LL, CudV1LL}
            mcbsg.setCoeff(0, (-1./mySMEFT.getNc() * getCudV8LL(1,1,1,2) + 2. * getCudV1LL(1,1,1,2)) * LEFT_factor, LO);
            mcbsg.setCoeff(1, (.5 * (1. - 1./mySMEFT.getNc()/mySMEFT.getNc()) * getCudV8LL(1,1,1,2) + 1./mySMEFT.getNc() * getCudV1LL(1,1,1,2)) * LEFT_factor, LO);
            // Add penguin operators in the future
            mcbsg.setCoeff(6, getCdg(1,2) * LEFT_factor_radiative, LO);
            mcbsg.setCoeff(7, getCdG(1,2) * LEFT_factor_radiative * sqrt(mySMEFT.getAle()/mySMEFT.Als(mySMEFT.getMuw())), LO);
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
    gslpp::complex LEFT_factor_radiative = 16. * M_PI * M_PI / mySMEFT.Mrun(mySMEFT.getMuw(),mySMEFT.getQuarks(QCD::BOTTOM).getMass(),QCD::BOTTOM) * LEFT_factor / sqrt(4. * M_PI * mySMEFT.getAle());

    switch (mcprimebsg.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
        //  {O1Misiak, O2Misiak} = {{sbarTac cbarTab}, {sbarc cbar b}} = {{-(1/N), (-1 + N^2)/(2 N^2)}, {2, 1/N}} {OV8LLud_Fierzed,OV1LLud_Fierzed}
        //  C_JMS^T O_JMS = C_JMS^T M O_Misiak = C_Misiak^T O_Misiak
        //  C_Misiak = M^T C_JMS = {{-(1/N), 2}, {(-1 + N^2)/(2 N^2), 1/N}} {CudV8LL, CudV1LL}
            mcprimebsg.setCoeff(0, (-1./mySMEFT.getNc() * getCudV8RR(1,1,1,2) + 2. * getCudV1RR(1,1,1,2)) * LEFT_factor, LO);
            mcprimebsg.setCoeff(1, (.5 * (1. - 1./mySMEFT.getNc()/mySMEFT.getNc()) * getCudV8RR(1,1,1,2) + 1./mySMEFT.getNc() * getCudV1RR(1,1,1,2)) * LEFT_factor, LO);
            // Add penguin operators in the future
            mcprimebsg.setCoeff(6, getCdg(2,1).conjugate() * LEFT_factor_radiative, LO);
            mcprimebsg.setCoeff(7, getCdG(2,1).conjugate() * LEFT_factor_radiative * (mySMEFT.getAle()/mySMEFT.Als(mySMEFT.getMuw())), LO);
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
    return gslpp::matrix<gslpp::complex>::Id(3);
}   

const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVuR() const
{
    return gslpp::matrix<gslpp::complex>::Id(3);
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVdL() const
{
    return mySMEFT.getCKM().getCKM();
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVdR() const
{
    return gslpp::matrix<gslpp::complex>::Id(3);
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVeL() const
{
    return mySMEFT.getVeL();
}   
    
const gslpp::matrix<gslpp::complex> NPSMEFTd6GeneralMatching::getVeR() const
{
    return mySMEFT.getVeR();
}   
