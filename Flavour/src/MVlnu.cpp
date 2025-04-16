/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "MVlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <boost/bind/bind.hpp>
#include <limits>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_expint.h>
#include <limits>
using namespace boost::placeholders;

MVlnu::MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: mySM(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    CLNflag = false;
    BGLflag = false;
    DMflag = false;
    btocNPpmflag = false;
    NPanalysis = false;
    
    checkcache_int_tau = false;
    checkcache_int_mu = false;
    checkcache_int_el = false;
    
    double max_double = std::numeric_limits<double>::max();
    
    hA1w1_cache = max_double;
    rho2_cache = max_double;
    R1w1_cache = max_double;
    R2w1_cache = max_double;
    N_A_cache = max_double;
    N_1_cache =max_double;
    N_2_cache = max_double;
    j_A_cache= max_double;
    j_0_cache = max_double;
    j_1_cache = max_double;
    j_2_cache = max_double;
    k_A_cache = max_double;
    k_0_cache = max_double;
    k_1_cache = max_double;
    k_2_cache = max_double;
    l_A_cache = max_double;

    af0_cache = max_double;
    af1_cache = max_double;
    af2_cache = max_double;
    ag0_cache = max_double;
    ag1_cache = max_double;
    ag2_cache = max_double;
    aF11_cache = max_double;
    aF12_cache = max_double;
    aF13_cache = max_double;
    aF21_cache = max_double;
    aF22_cache = max_double;
    aF23_cache = max_double;
    
    af_1_cache = max_double;
    ag_1_cache = max_double;
    aF1_1_cache = max_double;
    aF2_1_cache = max_double;
    af_2_cache = max_double;
    ag_2_cache = max_double;
    aF1_2_cache = max_double;
    aF2_2_cache = max_double;
    af_3_cache = max_double;
    ag_3_cache = max_double;
    aF1_3_cache = max_double;
    aF2_3_cache = max_double;
    af_4_cache = max_double;
    ag_4_cache = max_double;
    aF1_4_cache = max_double;
    aF2_4_cache = max_double;
    af_5_cache = max_double;
    ag_5_cache = max_double;
    aF1_5_cache = max_double;
    aF2_5_cache = max_double;
    af_6_cache = max_double;
    ag_6_cache = max_double;
    aF1_6_cache = max_double;
    aF2_6_cache = max_double;
    af_7_cache = max_double;
    ag_7_cache = max_double;
    aF1_7_cache = max_double;
    aF2_7_cache = max_double;
    af_8_cache = max_double;
    ag_8_cache = max_double;
    aF1_8_cache = max_double;
    aF2_8_cache = max_double;
    af_9_cache = max_double;
    ag_9_cache = max_double;
    aF1_9_cache = max_double;
    aF2_9_cache = max_double;
    af_10_cache = max_double;
    ag_10_cache = max_double;
    aF1_10_cache = max_double;
    aF2_10_cache = max_double;
    bf_1_cache = max_double;
    bg_1_cache = max_double;
    bF1_1_cache = max_double;
    bF2_1_cache = max_double;
    bf_2_cache = max_double;
    bg_2_cache = max_double;
    bF1_2_cache = max_double;
    bF2_2_cache = max_double;
    bf_3_cache = max_double;
    bg_3_cache = max_double;
    bF1_3_cache = max_double;
    bF2_3_cache = max_double;
    bf_4_cache = max_double;
    bg_4_cache = max_double;
    bF1_4_cache = max_double;
    bF2_4_cache = max_double;
    bf_5_cache = max_double;
    bg_5_cache = max_double;
    bF1_5_cache = max_double;
    bF2_5_cache = max_double;
    bf_6_cache = max_double;
    bg_6_cache = max_double;
    bF1_6_cache = max_double;
    bF2_6_cache = max_double;
    bf_7_cache = max_double;
    bg_7_cache = max_double;
    bF1_7_cache = max_double;
    bF2_7_cache = max_double;
    bf_8_cache = max_double;
    bg_8_cache = max_double;
    bF1_8_cache = max_double;
    bF2_8_cache = max_double;
    bf_9_cache = max_double;
    bg_9_cache = max_double;
    bF1_9_cache = max_double;
    bF2_9_cache = max_double;
    bf_10_cache = max_double;
    bg_10_cache = max_double;
    bF1_10_cache = max_double;
    bF2_10_cache = max_double;

    CS_cache = max_double;
    CSp_cache = max_double;
    CP_cache = max_double;
    CPp_cache = max_double;
    CV_cache = max_double;
    CVp_cache = max_double;
    CA_cache = max_double;
    CAp_cache = max_double;
    CT_cache = max_double;
    CTp_cache = max_double;
}

MVlnu::~MVlnu() {
}

std::vector<std::string> MVlnu::initializeMVlnuParameters()
{
    CLNflag = mySM.getFlavour().getFlagCLN();
    BGLflag = mySM.getFlavour().getFlagBGL();
    DMflag = mySM.getFlavour().getFlagDM();
    btocNPpmflag = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0);
    NPanalysis = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0 || mySM.getModelName().compare("RealWeakEFTCC") == 0);   
    
    if (CLNflag + BGLflag + DMflag != true) throw std::runtime_error("MVlnu: Set only one among CLNflag, BGLflag, DMflag to true");
    else mvlnuParameters = make_vector<std::string>();
    if (CLNflag) {
        mvlnuParameters.clear();
        if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>()
            << "hA1w1" << "rho2" << "R1w1" << "R2w1"
            << "N_A" << "N_1" << "N_2" << "j_A" << "j_0" << "j_1" << "j_2"
            << "k_A" << "k_0" << "k_1" << "k_2" << "l_A";
    } 
    else if (BGLflag) {
        mvlnuParameters.clear();
        if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>()
            << "af0" << "af1" << "af2" << "ag0" << "ag1" << "ag2"
            << "aF11" << "aF12" << "aF13" << "aF21" << "aF22" << "aF23"
            << "mBcstV1" << "mBcstV2" << "mBcstV3" << "mBcstV4"
            << "mBcstA1" << "mBcstA2" << "mBcstA3" << "mBcstA4"
            << "mBcstP1" << "mBcstP2" << "mBcstP3"
            << "chiTV" << "chiTA" << "chiTP" << "nI";
    }
    else if (DMflag){
        mvlnuParameters.clear();
        if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>() 
            << "af_1" << "ag_1" << "aF1_1" << "aF2_1"
            << "af_2" << "ag_2" << "aF1_2" << "aF2_2"
            << "af_3" << "ag_3" << "aF1_3" << "aF2_3"
            << "af_4" << "ag_4" << "aF1_4" << "aF2_4"
            << "af_5" << "ag_5" << "aF1_5" << "aF2_5"
            << "af_6" << "ag_6" << "aF1_6" << "aF2_6"
            << "af_7" << "ag_7" << "aF1_7" << "aF2_7"
            << "af_8" << "ag_8" << "aF1_8" << "aF2_8"
            << "af_9" << "ag_9" << "aF1_9" << "aF2_9"
            << "af_10" << "ag_10" << "aF1_10" << "aF2_10"
            << "bf_1" << "bg_1" << "bF1_1" << "bF2_1"
            << "bf_2" << "bg_2" << "bF1_2" << "bF2_2"
            << "bf_3" << "bg_3" << "bF1_3" << "bF2_3"
            << "bf_4" << "bg_4" << "bF1_4" << "bF2_4"
            << "bf_5" << "bg_5" << "bF1_5" << "bF2_5"
            << "bf_6" << "bg_6" << "bF1_6" << "bF2_6"
            << "bf_7" << "bg_7" << "bF1_7" << "bF2_7"
            << "bf_8" << "bg_8" << "bF1_8" << "bF2_8"
            << "bf_9" << "bg_9" << "bF1_9" << "bF2_9"
            << "bf_10" << "bg_10" << "bF1_10" << "bF2_10";
    }
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }

    mySM.initializeMeson(meson);
    mySM.initializeMeson(vectorM);
    return mvlnuParameters;
}

void MVlnu::updateParameters() 
{
    if (!mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) return;

    Mlep = mySM.getLeptons(lep).getMass();
    Mnu = 0.; // neutrinos assumed to be massless
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    width = mySM.getMesons(meson).computeWidth();
    w0 = (MM*MM+MV*MV)/(2.*MM*MV);
    RV = 2.*sqrt(MM*MV)/(MM+MV);
    mu_b = MM; // mySM.getMub();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass(); // add the PS b mass
    ale_mub = mySM.Ale(mu_b,FULLNLO);
    /* Amplitude propto 4*GF*Vij/sqrt(2) & kinematics requires 1/(2^9 pi^3 MB^3) */
    amplsq_factor = 1./(64.*M_PI*M_PI*M_PI*MM*MM*MM);
    q2min = Mlep*Mlep;
    q2max = (MM-MV)*(MM-MV);
    
    MV_o_MM = MV / MM;
    sqrtMV_o_MM = sqrt(MV_o_MM);
    
    /* SM + NP Wilson coefficients */
    gslpp::complex norm = 4./sqrt(2.);
    gslpp::vector<gslpp::complex> ** allcoeff_bclnu = mySM.getFlavour().ComputeCoeffdiujlknu(2,1,0,mu_b);
    CV = (*(allcoeff_bclnu[LO]))(0)/norm*(1.+ale_mub/M_PI*log(mySM.getMz()/mu_b))/2.;
    CA = -CV;
    CVp = (*(allcoeff_bclnu[LO]))(1)/norm/2.;
    CAp = -CVp;
    CS = (*(allcoeff_bclnu[LO]))(2)/norm/2.;
    CSp = (*(allcoeff_bclnu[LO]))(3)/norm/2.;
    CP = -CS;
    CPp = -CSp;
    C7 = 0.;
    C7p = 0.;
    CT = (*(allcoeff_bclnu[LO]))(4)/norm/2.;
    CTp = 0.;

    /* SM + NP Wilson coefficients */
    if (NPanalysis) {
        if (lep == StandardModel::TAU) {
            if (btocNPpmflag) {
                CV += (mySM.getCCC3() - mySM.getCCC4()) / 2. / M_SQRT2;
                CVp = (mySM.getCCC3() + mySM.getCCC4()) / 2. / M_SQRT2;
                CA -= (mySM.getCCC3() - mySM.getCCC4()) / 2. / M_SQRT2;
                CAp = -(mySM.getCCC3() + mySM.getCCC4()) / 2. / M_SQRT2;
                CS = (mySM.getCCC1() - mySM.getCCC2()) / 2. / M_SQRT2;
                CSp = (mySM.getCCC1() + mySM.getCCC2()) / 2. / M_SQRT2;
                CP = -(mySM.getCCC1() - mySM.getCCC2()) / 2. / M_SQRT2;
                CPp = -(mySM.getCCC1() + mySM.getCCC2()) / 2. / M_SQRT2;
                CTp = mySM.getOptionalParameter("CT_NP");
            } else {
                CV += mySM.getCCC3() / 2.;
                CVp = mySM.getCCC4() / 2.;
                CA -= mySM.getCCC3() / 2.;
                CAp = -mySM.getCCC4() / 2.;
                CS = mySM.getCCC1() / 2.;
                CSp = mySM.getCCC2() / 2.;
                CP = -mySM.getCCC1() / 2.;
                CPp = -mySM.getCCC2() / 2.;
                CTp = mySM.getCCC5();
            }
        }
    }

    switch (vectorM) {
        case StandardModel::D_star_P:
            if (CLNflag) {
                hA1w1 = mySM.getOptionalParameter("hA1w1");
                rho2 = mySM.getOptionalParameter("rho2");
                R1w1 = mySM.getOptionalParameter("R1w1");
                R2w1 = mySM.getOptionalParameter("R2w1");
                N_A = mySM.getOptionalParameter("N_A");
                N_1 = mySM.getOptionalParameter("N_1");
                N_2 = mySM.getOptionalParameter("N_2"); 
                j_A = mySM.getOptionalParameter("j_A");
                j_0 = mySM.getOptionalParameter("j_0");
                j_1 = mySM.getOptionalParameter("j_1"); 
                j_2 = mySM.getOptionalParameter("j_2");
                k_A = mySM.getOptionalParameter("k_A");
                k_0 = mySM.getOptionalParameter("k_0");
                k_1 = mySM.getOptionalParameter("k_1");
                k_2 = mySM.getOptionalParameter("k_2");
                l_A = mySM.getOptionalParameter("l_A");
            } 
            else if (BGLflag) {
                af0 = mySM.getOptionalParameter("af0");
                af1 = mySM.getOptionalParameter("af1");
                af2 = mySM.getOptionalParameter("af2");
                ag0 = mySM.getOptionalParameter("ag0");
                ag1 = mySM.getOptionalParameter("ag1");
                ag2 = mySM.getOptionalParameter("ag2");
                aF11 = mySM.getOptionalParameter("aF11");
                aF12 = mySM.getOptionalParameter("aF12");
                aF13 = mySM.getOptionalParameter("aF13");
                aF21 = mySM.getOptionalParameter("aF21");
                aF22 = mySM.getOptionalParameter("aF22");
                aF23 = mySM.getOptionalParameter("aF23");
                mBcstV1 = mySM.getOptionalParameter("mBcstV1");
                mBcstV2 = mySM.getOptionalParameter("mBcstV2");
                mBcstV3 = mySM.getOptionalParameter("mBcstV3");
                mBcstV4 = mySM.getOptionalParameter("mBcstV4");
                mBcstA1 = mySM.getOptionalParameter("mBcstA1");
                mBcstA2 = mySM.getOptionalParameter("mBcstA2");
                mBcstA3 = mySM.getOptionalParameter("mBcstA3");
                mBcstA4 = mySM.getOptionalParameter("mBcstA4");
                mBcstP1 = mySM.getOptionalParameter("mBcstP1");
                mBcstP2 = mySM.getOptionalParameter("mBcstP2");
                mBcstP3 = mySM.getOptionalParameter("mBcstP3");
                chiTV = mySM.getOptionalParameter("chiTV");
                chiTA = mySM.getOptionalParameter("chiTA");
                chiTP = mySM.getOptionalParameter("chiTP");
                nI = mySM.getOptionalParameter("nI");
                
                zV1 = sqrt((MM+MV)*(MM+MV)-mBcstV1*mBcstV1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zV1 /= (sqrt((MM+MV)*(MM+MV)-mBcstV1*mBcstV1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zV2 = sqrt((MM+MV)*(MM+MV)-mBcstV2*mBcstV2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zV2 /= (sqrt((MM+MV)*(MM+MV)-mBcstV2*mBcstV2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zV3 = sqrt((MM+MV)*(MM+MV)-mBcstV3*mBcstV3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zV3 /= (sqrt((MM+MV)*(MM+MV)-mBcstV3*mBcstV3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zV4 = sqrt((MM+MV)*(MM+MV)-mBcstV4*mBcstV4)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zV4 /= (sqrt((MM+MV)*(MM+MV)-mBcstV4*mBcstV4)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));

                zA1 = sqrt((MM+MV)*(MM+MV)-mBcstA1*mBcstA1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zA1 /= (sqrt((MM+MV)*(MM+MV)-mBcstA1*mBcstA1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zA2 = sqrt((MM+MV)*(MM+MV)-mBcstA2*mBcstA2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zA2 /= (sqrt((MM+MV)*(MM+MV)-mBcstA2*mBcstA2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zA3 = sqrt((MM+MV)*(MM+MV)-mBcstA3*mBcstA3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zA3 /= (sqrt((MM+MV)*(MM+MV)-mBcstA3*mBcstA3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zA4 = sqrt((MM+MV)*(MM+MV)-mBcstA4*mBcstA4)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zA4 /= (sqrt((MM+MV)*(MM+MV)-mBcstA4*mBcstA4)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));

                zP1 = sqrt((MM+MV)*(MM+MV)-mBcstP1*mBcstP1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zP1 /= (sqrt((MM+MV)*(MM+MV)-mBcstP1*mBcstP1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zP2 = sqrt((MM+MV)*(MM+MV)-mBcstP2*mBcstP2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zP2 /= (sqrt((MM+MV)*(MM+MV)-mBcstP2*mBcstP2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
                zP3 = sqrt((MM+MV)*(MM+MV)-mBcstP3*mBcstP3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
                zP3 /= (sqrt((MM+MV)*(MM+MV)-mBcstP3*mBcstP3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
            }
            else if (DMflag) {
                af_1 = mySM.getOptionalParameter("af_1");
                ag_1 = mySM.getOptionalParameter("ag_1");
                aF1_1 = mySM.getOptionalParameter("aF1_1");
                aF2_1 = mySM.getOptionalParameter("aF2_1");
                af_2 = mySM.getOptionalParameter("af_2");
                ag_2 = mySM.getOptionalParameter("ag_2");
                aF1_2 = mySM.getOptionalParameter("aF1_2");
                aF2_2 = mySM.getOptionalParameter("aF2_2");
                af_3 = mySM.getOptionalParameter("af_3");
                ag_3 = mySM.getOptionalParameter("ag_3");
                aF1_3 = mySM.getOptionalParameter("aF1_3");
                aF2_3 = mySM.getOptionalParameter("aF2_3");
                af_4 = mySM.getOptionalParameter("af_4");
                ag_4 = mySM.getOptionalParameter("ag_4");
                aF1_4 = mySM.getOptionalParameter("aF1_4");
                aF2_4 = mySM.getOptionalParameter("aF2_4");
                af_5 = mySM.getOptionalParameter("af_5");
                ag_5 = mySM.getOptionalParameter("ag_5");
                aF1_5 = mySM.getOptionalParameter("aF1_5");
                aF2_5 = mySM.getOptionalParameter("aF2_5");
                af_6 = mySM.getOptionalParameter("af_6");
                ag_6 = mySM.getOptionalParameter("ag_6");
                aF1_6 = mySM.getOptionalParameter("aF1_6");
                aF2_6 = mySM.getOptionalParameter("aF2_6");
                af_7 = mySM.getOptionalParameter("af_7");
                ag_7 = mySM.getOptionalParameter("ag_7");
                aF1_7 = mySM.getOptionalParameter("aF1_7");
                aF2_7 = mySM.getOptionalParameter("aF2_7");
                af_8 = mySM.getOptionalParameter("af_8");
                ag_8 = mySM.getOptionalParameter("ag_8");
                aF1_8 = mySM.getOptionalParameter("aF1_8");
                aF2_8 = mySM.getOptionalParameter("aF2_8");
                af_9 = mySM.getOptionalParameter("af_9");
                ag_9 = mySM.getOptionalParameter("ag_9");
                aF1_9 = mySM.getOptionalParameter("aF1_9");
                aF2_9 = mySM.getOptionalParameter("aF2_9");
                af_10 = mySM.getOptionalParameter("af_10");
                ag_10 = mySM.getOptionalParameter("ag_10");
                aF1_10 = mySM.getOptionalParameter("aF1_10");
                aF2_10 = mySM.getOptionalParameter("aF2_10");
                bf_1 = mySM.getOptionalParameter("bf_1");
                bg_1 = mySM.getOptionalParameter("bg_1");
                bF1_1 = mySM.getOptionalParameter("bF1_1");
                bF2_1 = mySM.getOptionalParameter("bF2_1");
                bf_2 = mySM.getOptionalParameter("bf_2");
                bg_2 = mySM.getOptionalParameter("bg_2");
                bF1_2 = mySM.getOptionalParameter("bF1_2");
                bF2_2 = mySM.getOptionalParameter("bF2_2");
                bf_3 = mySM.getOptionalParameter("bf_3");
                bg_3 = mySM.getOptionalParameter("bg_3");
                bF1_3 = mySM.getOptionalParameter("bF1_3");
                bF2_3 = mySM.getOptionalParameter("bF2_3");
                bf_4 = mySM.getOptionalParameter("bf_4");
                bg_4 = mySM.getOptionalParameter("bg_4");
                bF1_4 = mySM.getOptionalParameter("bF1_4");
                bF2_4 = mySM.getOptionalParameter("bF2_4");
                bf_5 = mySM.getOptionalParameter("bf_5");
                bg_5 = mySM.getOptionalParameter("bg_5");
                bF1_5 = mySM.getOptionalParameter("bF1_5");
                bF2_5 = mySM.getOptionalParameter("bF2_5");
                bf_6 = mySM.getOptionalParameter("bf_6");
                bg_6 = mySM.getOptionalParameter("bg_6");
                bF1_6 = mySM.getOptionalParameter("bF1_6");
                bF2_6 = mySM.getOptionalParameter("bF2_6");
                bf_7 = mySM.getOptionalParameter("bf_7");
                bg_7 = mySM.getOptionalParameter("bg_7");
                bF1_7 = mySM.getOptionalParameter("bF1_7");
                bF2_7 = mySM.getOptionalParameter("bF2_7");
                bf_8 = mySM.getOptionalParameter("bf_8");
                bg_8 = mySM.getOptionalParameter("bg_8");
                bF1_8 = mySM.getOptionalParameter("bF1_8");
                bF2_8 = mySM.getOptionalParameter("bF2_8");
                bf_9 = mySM.getOptionalParameter("bf_9");
                bg_9 = mySM.getOptionalParameter("bg_9");
                bF1_9 = mySM.getOptionalParameter("bF1_9");
                bF2_9 = mySM.getOptionalParameter("bF2_9");
                bf_10 = mySM.getOptionalParameter("bf_10");
                bg_10 = mySM.getOptionalParameter("bg_10");
                bF1_10 = mySM.getOptionalParameter("bF1_10");
                bF2_10 = mySM.getOptionalParameter("bF2_10");
            }
        else{};
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }

    if ((hA1w1 != hA1w1_cache) || (rho2 != rho2_cache) || (R1w1 != R1w1_cache) || (R2w1 != R2w1_cache)
            || (N_A != N_A_cache) || (N_1 != N_1_cache) || (N_2 != N_2_cache) || (l_A != l_A_cache)
            || (j_A != j_A_cache) || (j_0 != j_0_cache) || (j_1 != j_1_cache) || (j_2 != j_2_cache)
            || (k_A != k_A_cache) || (k_0 != k_0_cache) || (k_1 != k_1_cache) || (k_2 != k_2_cache)
            || (af0 != af0_cache) || (af1 != af1_cache) || (af2 != af2_cache) 
            || (ag0 != ag0_cache) || (ag1 != af1_cache) || (ag2 != af2_cache)
            || (aF11 != aF11_cache) || (aF12 != aF12_cache) || (aF13 != aF13_cache) 
            || (aF21 != aF21_cache) || (aF22 != aF22_cache) || (aF23 != aF23_cache)
            || (af_1 != af_1_cache) || (ag_1 != ag_1_cache) || (aF1_1 != aF1_1_cache) || (aF2_1 != aF2_1_cache)
            || (af_2 != af_2_cache) || (ag_2 != ag_2_cache) || (aF1_2 != aF1_2_cache) || (aF2_2 != aF2_2_cache)
            || (af_3 != af_3_cache) || (ag_3 != ag_3_cache) || (aF1_3 != aF1_3_cache) || (aF2_3 != aF2_3_cache)
            || (af_4 != af_4_cache) || (ag_4 != ag_4_cache) || (aF1_4 != aF1_4_cache) || (aF2_4 != aF2_4_cache)
            || (af_5 != af_5_cache) || (ag_5 != ag_5_cache) || (aF1_5 != aF1_5_cache) || (aF2_5 != aF2_5_cache)
            || (af_6 != af_6_cache) || (ag_6 != ag_6_cache) || (aF1_6 != aF1_6_cache) || (aF2_6 != aF2_6_cache)
            || (af_7 != af_7_cache) || (ag_7 != ag_7_cache) || (aF1_7 != aF1_7_cache) || (aF2_7 != aF2_7_cache)
            || (af_8 != af_8_cache) || (ag_8 != ag_8_cache) || (aF1_8 != aF1_8_cache) || (aF2_8 != aF2_8_cache)
            || (af_9 != af_9_cache) || (ag_9 != ag_9_cache) || (aF1_9 != aF1_9_cache) || (aF2_9 != aF2_9_cache)
            || (af_10 != af_10_cache) || (ag_10 != ag_10_cache) || (aF1_10 != aF1_10_cache) || (aF2_10 != aF2_10_cache)
            || (bf_1 != bf_1_cache) || (bg_1 != bg_1_cache) || (bF1_1 != bF1_1_cache) || (bF2_1 != bF2_1_cache)
            || (bf_2 != bf_2_cache) || (bg_2 != bg_2_cache) || (bF1_2 != bF1_2_cache) || (bF2_2 != bF2_2_cache)
            || (bf_3 != bf_3_cache) || (bg_3 != bg_3_cache) || (bF1_3 != bF1_3_cache) || (bF2_3 != bF2_3_cache)
            || (bf_4 != bf_4_cache) || (bg_4 != bg_4_cache) || (bF1_4 != bF1_4_cache) || (bF2_4 != bF2_4_cache)
            || (bf_5 != bf_5_cache) || (bg_5 != bg_5_cache) || (bF1_5 != bF1_5_cache) || (bF2_5 != bF2_5_cache)
            || (bf_6 != bf_6_cache) || (bg_6 != bg_6_cache) || (bF1_6 != bF1_6_cache) || (bF2_6 != bF2_6_cache)
            || (bf_7 != bf_7_cache) || (bg_7 != bg_7_cache) || (bF1_7 != bF1_7_cache) || (bF2_7 != bF2_7_cache)
            || (bf_8 != bf_8_cache) || (bg_8 != bg_8_cache) || (bF1_8 != bF1_8_cache) || (bF2_8 != bF2_8_cache)
            || (bf_9 != bf_9_cache) || (bg_9 != bg_9_cache) || (bF1_9 != bF1_9_cache) || (bF2_9 != bF2_9_cache)
            || (bf_10 != bf_10_cache) || (bg_10 != bg_10_cache) || (bF1_10 != bF1_10_cache) || (bF2_10 != bF2_10_cache)
            || (CS != CS_cache) || (CSp != CSp_cache)
            || (CP != CP_cache) || (CPp != CPp_cache)
            || (CV != CV_cache) || (CVp != CVp_cache)
            || (CA != CA_cache) || (CAp != CAp_cache)
            || (CT != CT_cache) || (CTp != CTp_cache)) {
        checkcache_int_tau = false;
        checkcache_int_mu = false;
        checkcache_int_el = false;
    }
    
    if (!checkcache_int_tau || !checkcache_int_mu || !checkcache_int_el) {
        if (lep == StandardModel::TAU) {
            cached_intJ1s_tau = integrateJ(1, q2min, q2max);
            cached_intJ1c_tau = integrateJ(2, q2min, q2max);
            cached_intJ2s_tau = integrateJ(3, q2min, q2max);
            cached_intJ2c_tau = integrateJ(4, q2min, q2max);
            cached_intJ3_tau = integrateJ(5, q2min, q2max);
            cached_intJ6s_tau = integrateJ(8, q2min, q2max);
            cached_intJ6c_tau = integrateJ(9, q2min, q2max);
            cached_intJ9_tau = integrateJ(12, q2min, q2max);
            cached_intJ4_tau = 0.;
            cached_intJ5_tau = 0.;
            cached_intJ7_tau = 0.;
            cached_intJ8_tau = 0.;
            // not needed at present
            /*
            cached_intJ4_tau = integrateJ(6,q2min,q2max);
            cached_intJ5_tau = integrateJ(7,q2min,q2max);
            cached_intJ7_tau = integrateJ(10,q2min,q2max);
            cached_intJ8_tau = integrateJ(11,q2min,q2max);
             */
            checkcache_int_tau = true;
        }
        if (lep == StandardModel::MU) {
            cached_intJ1s_mu = integrateJ(1, q2min, q2max);
            cached_intJ1c_mu = integrateJ(2, q2min, q2max);
            cached_intJ2s_mu = integrateJ(3, q2min, q2max);
            cached_intJ2c_mu = integrateJ(4, q2min, q2max);
            cached_intJ3_mu = integrateJ(5, q2min, q2max);
            cached_intJ6s_mu = integrateJ(8, q2min, q2max);
            cached_intJ6c_mu = integrateJ(9, q2min, q2max);
            cached_intJ9_mu = integrateJ(12, q2min, q2max);
            cached_intJ4_mu = 0.;
            cached_intJ5_mu = 0.;
            cached_intJ7_mu = 0.;
            cached_intJ8_mu = 0.;
            // not needed at present
            /*
            cached_intJ4_mu = integrateJ(6,q2min,q2max);
            cached_intJ5_mu = integrateJ(7,q2min,q2max);
            cached_intJ7_mu = integrateJ(10,q2min,q2max);
            cached_intJ8_mu = integrateJ(11,q2min,q2max);
             */
            checkcache_int_mu = true;
        }
        if (lep == StandardModel::ELECTRON) {
            cached_intJ1s_el = integrateJ(1, q2min, q2max);
            cached_intJ1c_el = integrateJ(2, q2min, q2max);
            cached_intJ2s_el = integrateJ(3, q2min, q2max);
            cached_intJ2c_el = integrateJ(4, q2min, q2max);
            cached_intJ3_el = integrateJ(5, q2min, q2max);
            cached_intJ6s_el = integrateJ(8, q2min, q2max);
            cached_intJ6c_el = integrateJ(9, q2min, q2max);
            cached_intJ9_el = integrateJ(12, q2min, q2max);
            cached_intJ4_el = 0.;
            cached_intJ5_el = 0.;
            cached_intJ7_el = 0.;
            cached_intJ8_el = 0.;
            // not needed at present
            /*
            cached_intJ4_el = integrateJ(6,q2min,q2max);
            cached_intJ5_el = integrateJ(7,q2min,q2max);
            cached_intJ7_el = integrateJ(10,q2min,q2max);
            cached_intJ8_el = integrateJ(11,q2min,q2max);
             */
            checkcache_int_el = true;
        }
    }
    if (CLNflag) {
        hA1w1_cache = hA1w1;
        rho2_cache = rho2;
        R1w1_cache = R1w1;
        R2w1_cache = R2w1;
        N_A_cache = N_A;
        N_1_cache =N_1;
        N_2_cache = N_2;
        j_A_cache= j_A;
        j_0_cache = j_0;
        j_1_cache = j_1;
        j_2_cache = j_2;
        k_A_cache = k_A;
        k_0_cache = k_0;
        k_1_cache = k_1;
        k_2_cache = k_2;
        l_A_cache = l_A;
    }
    else if (BGLflag){
        af0_cache = af0;
        af1_cache = af1;
        af2_cache = af2;
        ag0_cache = ag0;
        ag1_cache = ag1;
        ag2_cache = ag2;
        aF11_cache = aF11;
        aF12_cache = aF12;
        aF13_cache = aF13;
        aF21_cache = aF21;
        aF22_cache = aF22;
        aF23_cache = aF23;
    }
    else if (DMflag){
        af_1_cache = af_1;
        ag_1_cache = ag_1;
        aF1_1_cache = aF1_1;
        aF2_1_cache = aF2_1;
        af_2_cache = af_2;
        ag_2_cache = ag_2;
        aF1_2_cache = aF1_2;
        aF2_2_cache = aF2_2;
        af_3_cache = af_3;
        ag_3_cache = ag_3;
        aF1_3_cache = aF1_3;
        aF2_3_cache = aF2_3;
        af_4_cache = af_4;
        ag_4_cache = ag_4;
        aF1_4_cache = aF1_4;
        aF2_4_cache = aF2_4;
        af_5_cache = af_5;
        ag_5_cache = ag_5;
        aF1_5_cache = aF1_5;
        aF2_5_cache = aF2_5;
        af_6_cache = af_6;
        ag_6_cache = ag_6;
        aF1_6_cache = aF1_6;
        aF2_6_cache = aF2_6;
        af_7_cache = af_7;
        ag_7_cache = ag_7;
        aF1_7_cache = aF1_7;
        aF2_7_cache = aF2_7;
        af_8_cache = af_8;
        ag_8_cache = ag_8;
        aF1_8_cache = aF1_8;
        aF2_8_cache = aF2_8;
        af_9_cache = af_9;
        ag_9_cache = ag_9;
        aF1_9_cache = aF1_9;
        aF2_9_cache = aF2_9;
        af_10_cache = af_10;
        ag_10_cache = ag_10;
        aF1_10_cache = aF1_10;
        aF2_10_cache = aF2_10;
        bf_1_cache = bf_1;
        bg_1_cache = bg_1;
        bF1_1_cache = bF1_1;
        bF2_1_cache = bF2_1;
        bf_2_cache = bf_2;
        bg_2_cache = bg_2;
        bF1_2_cache = bF1_2;
        bF2_2_cache = bF2_2;
        bf_3_cache = bf_3;
        bg_3_cache = bg_3;
        bF1_3_cache = bF1_3;
        bF2_3_cache = bF2_3;
        bf_4_cache = bf_4;
        bg_4_cache = bg_4;
        bF1_4_cache = bF1_4;
        bF2_4_cache = bF2_4;
        bf_5_cache = bf_5;
        bg_5_cache = bg_5;
        bF1_5_cache = bF1_5;
        bF2_5_cache = bF2_5;
        bf_6_cache = bf_6;
        bg_6_cache = bg_6;
        bF1_6_cache = bF1_6;
        bF2_6_cache = bF2_6;
        bf_7_cache = bf_7;
        bg_7_cache = bg_7;
        bF1_7_cache = bF1_7;
        bF2_7_cache = bF2_7;
        bf_8_cache = bf_8;
        bg_8_cache = bg_8;
        bF1_8_cache = bF1_8;
        bF2_8_cache = bF2_8;
        bf_9_cache = bf_9;
        bg_9_cache = bg_9;
        bF1_9_cache = bF1_9;
        bF2_9_cache = bF2_9;
        bf_10_cache = bf_10;
        bg_10_cache = bg_10;
        bF1_10_cache = bF1_10;
        bF2_10_cache = bF2_10;
    }
    else{};
    
    CS_cache = CS;
    CSp_cache = CSp;
    CP_cache = CP;
    CPp_cache = CPp;
    CV_cache = CV;
    CVp_cache = CVp;
    CA_cache = CA;
    CAp_cache = CAp;
    CT_cache = CT;
    CTp_cache = CTp;

    mySM.getFlavour().setUpdateFlag(meson, vectorM, lep, false);

    return;
    
}

/*******************************************************************************
 * Kinematic functions                                                          *
 * ****************************************************************************/

double MVlnu::lambda_half(double a, double b, double c) 
{   
    return sqrt(a*a+b*b+c*c-2.*(a*b+a*c+b*c));
}
/*******************************************************************************
 * Form factors                                                                *
 * ****************************************************************************/

double MVlnu::phi_f(double z)
{
    double prefac = 4. * (MV_o_MM) / MM / MM * sqrt(nI / (3. * M_PI * chiTA));
    double num = (1. + z) * sqrt((1. - z)*(1. - z)*(1. - z));
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::f_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double Pfacf = (z - zA1) / (1. - z * zA1)*(z - zA2) / (1. - z * zA2)*(z - zA3) / (1. - z * zA3)*(z - zA4) / (1. - z * zA4);
    double phif = phi_f(z);
    return (af0 + af1 * z + af2 * z * z) / phif / Pfacf;
}

double MVlnu::phi_g(double z)
{
    double prefac = sqrt(nI / (3. * M_PI * chiTV));
    double num = 16. * (MV_o_MM)*(MV_o_MM)*(1. + z)*(1. + z) / sqrt(1. - z);
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::g_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double Pfacg = (z - zV1) / (1. - z * zV1)*(z - zV2) / (1. - z * zV2)*(z - zV3) / (1. - z * zV3)*(z - zV4) / (1. - z * zV4);
    double phig = phi_g(z);
    return (ag0 + ag1 * z + ag2 * z * z) / phig / Pfacg;
}

double MVlnu::phi_F1(double z)
{
    double prefac = 4. * (MV_o_MM) / MM / MM / MM * sqrt(nI / (6. * M_PI * chiTA));
    double num = (1. + z) * sqrt((1. - z)*(1. - z)*(1. - z)*(1. - z)*(1. - z));
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den5 = den * den * den * den*den;
    return prefac * num / den5;
}

double MVlnu::F1_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double PfacF1 = (z - zA1) / (1. - z * zA1)*(z - zA2) / (1. - z * zA2)*(z - zA3) / (1. - z * zA3)*(z - zA4) / (1. - z * zA4);
    double phiF1 = phi_F1(z);
    double aF10 = (MM - MV)*(phi_F1(0.) / phi_f(0.)) * af0; // F1(z=0) = (MM-MV)*f(z=0)
    return (aF10 + aF11 * z + aF12 * z * z + aF13 * z * z * z) / phiF1 / PfacF1;
}

double MVlnu::phi_F2(double z)
{
    double prefac = 8. * sqrt(2.)*(MV_o_MM)*(MV_o_MM) * sqrt(nI / (M_PI * chiTP));
    double num = (1. + z)*(1. + z) / sqrt(1. - z);
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::F2_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double z0 = (sqrt(w0 + 1.) - M_SQRT2) / (sqrt(w0 + 1.) + M_SQRT2);
    double PfacF2 = (z - zP1) / (1. - z * zP1)*(z - zP2) / (1. - z * zP2)*(z - zP3) / (1. - z * zP3);
    double PfacF2z0 = (z0 - zP1) / (1. - z0 * zP1)*(z0 - zP2) / (1. - z0 * zP2)*(z0 - zP3) / (1. - z0 * zP3);
    double phiF2 = phi_F2(z);
    double phiF2z0 = phi_F2(z0);
    double aF20 = PfacF2z0 * phiF2z0 * 2. * F1_BGL(0.) / (MM * MM - MV * MV) - aF21 * z0 - aF22 * z0*z0; // F2(q2=0) = 2.*F1(q2=0)/(MM*MM-MV*MV)
    return (aF20 + aF21 * z + aF22 * z * z + aF23 * z * z * z) / phiF2 / PfacF2;
}

double MVlnu::hA1(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    if (CLNflag) return hA1w1 * N_A * (1. - j_A * 8. * rho2 * z + k_A * (53. * rho2 - 15.) * z * z - l_A * (231. * rho2 - 91.) * z * z * z);
    else if (BGLflag)  return(f_BGL(q2) / sqrt(MM * MV) / (1. + w));
    else if (DMflag) {
        double w = w0 - q2 / (2. * MM * MV);
        double f_fac = 0.;
        if (w<1.05) {f_fac = af_1 + bf_1 * w;}
        else if (w<1.10) {f_fac = af_2 + bf_2 * w;}
        else if (w<1.15) {f_fac = af_3 + bf_3 * w;}
        else if (w<1.20) {f_fac = af_4 + bf_4 * w;}
        else if (w<1.25) {f_fac = af_5 + bf_5 * w;}
        else if (w<1.30) {f_fac = af_6 + bf_6 * w;}
        else if (w<1.35) {f_fac = af_7 + bf_7 * w;}
        else if (w<1.40) {f_fac = af_8 + bf_8 * w;}
        else if (w<1.45) {f_fac = af_9 + bf_9 * w;}
        else {f_fac = af_10 + bf_10 * w;}
        return f_fac / sqrt(MM * MV) / (1. + w);
    }
    else return 0.;
}

double MVlnu::R1(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    return N_1 * R1w1 - j_1 * 0.12 * (w - 1.) + k_1 * 0.05 * (w - 1.)*(w - 1.);
}

double MVlnu::R2(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    return N_2 * R2w1 + j_2 * 0.11 * (w - 1.) - k_2 * 0.06 * (w - 1.)*(w - 1.);
}

double MVlnu::R0(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    /* form factor relation among A0, A1 and A2 at q2=0 */
    double R2q2at0 = R2(0.);
    double R0q2at0 = (MM + MV - (MM - MV) * R2q2at0) / (2. * MV);
    // caveat: HQET rel at the kinematic endpoint, q2 = 0 ...
    double R0w1 = R0q2at0 + j_0 * 0.11 * (w0 - 1.) - k_0 * 0.01 * (w0 - 1.)*(w0 - 1.);
    // one may consider "lattice" R0w1 = 1.14 +- O(10%) + consistency rel at q2 = 0 ...
    return R0w1 - j_0 * 0.11 * (w - 1.) + k_0 * 0.01 * (w - 1.)*(w - 1.);
}

double MVlnu::V(double q2)
{
    if (CLNflag) return R1(q2) / RV * hA1(q2);
    else if (BGLflag)  return (MM + MV) * g_BGL(q2) / 2.;
    else if (DMflag) {
        double w = w0 - q2 / (2. * MM * MV);
        double g_fac = 0.;
        if (w<1.05) {g_fac = ag_1 + bg_1 * w;}
        else if (w<1.10) {g_fac = ag_2 + bg_2 * w;}
        else if (w<1.15) {g_fac = ag_3 + bg_3 * w;}
        else if (w<1.20) {g_fac = ag_4 + bg_4 * w;}
        else if (w<1.25) {g_fac = ag_5 + bg_5 * w;}
        else if (w<1.30) {g_fac = ag_6 + bg_6 * w;}
        else if (w<1.35) {g_fac = ag_7 + bg_7 * w;}
        else if (w<1.40) {g_fac = ag_8 + bg_8 * w;}
        else if (w<1.45) {g_fac = ag_9 + bg_9 * w;}
        else {g_fac = ag_10 + bg_10 * w;}
        return (MM + MV) * g_fac / 2.;
    }
    else return 0.;
}

double MVlnu::A0(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    // A0 = RV * P1 = R0 * A1
    if (CLNflag) return R0(q2) * (w + 1.) * RV / 2. * hA1(q2);
    // F2 = P1 / (RV / 2.)
    else if (BGLflag) return F2_BGL(q2) / 2.;
    else if (DMflag) {
        double w = w0 - q2 / (2. * MM * MV);
        double F2_fac = 0.;
        if (w<1.05) {F2_fac = aF2_1 + bF2_1 * w;} 
        else if (w<1.10) {F2_fac = aF2_2 + bF2_2 * w;} 
        else if (w<1.15) {F2_fac = aF2_3 + bF2_3 * w;} 
        else if (w<1.20) {F2_fac = aF2_4 + bF2_4 * w;} 
        else if (w<1.25) {F2_fac = aF2_5 + bF2_5 * w;} 
        else if (w<1.30) {F2_fac = aF2_6 + bF2_6 * w;} 
        else if (w<1.35) {F2_fac = aF2_7 + bF2_7 * w;} 
        else if (w<1.40) {F2_fac = aF2_8 + bF2_8 * w;} 
        else if (w<1.45) {F2_fac = aF2_9 + bF2_9 * w;} 
        else {F2_fac = aF2_10 + bF2_10 * w;} 
        return F2_fac / 2.;
    }
    else return 0.;
}

double MVlnu::A1(double q2)
{
    /* form factor in 1:1 with hA1 */
    double w = w0 - q2 / (2. * MM * MV);
    return (w + 1.) * RV / 2. * hA1(q2);
}

double MVlnu::A2(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    if (CLNflag) return R2(q2) / RV * hA1(q2);
    else if (BGLflag) return (MM + MV) / 2. / (w * w - 1.) / MM / MV * ((w - MV_o_MM) * f_BGL(q2) - F1_BGL(q2) / MM);
    else if (DMflag) {
        double w = w0 - q2 / (2. * MM * MV);
        double f_fac = 0.;
        double F1_fac = 0.;
        if (w<1.05) {f_fac = af_1 + bf_1 * w; F1_fac = aF1_1 + bF1_1 * w;}
        else if (w<1.10) {f_fac = af_2 + bf_2 * w; F1_fac = aF1_2 + bF1_2 * w;}
        else if (w<1.15) {f_fac = af_3 + bf_3 * w; F1_fac = aF1_3 + bF1_3 * w;}
        else if (w<1.20) {f_fac = af_4 + bf_4 * w; F1_fac = aF1_4 + bF1_4 * w;}
        else if (w<1.25) {f_fac = af_5 + bf_5 * w; F1_fac = aF1_5 + bF1_5 * w;}
        else if (w<1.30) {f_fac = af_6 + bf_6 * w; F1_fac = aF1_6 + bF1_6 * w;}
        else if (w<1.35) {f_fac = af_7 + bf_7 * w; F1_fac = aF1_7 + bF1_7 * w;}
        else if (w<1.40) {f_fac = af_8 + bf_8 * w; F1_fac = aF1_8 + bF1_8 * w;}
        else if (w<1.45) {f_fac = af_9 + bf_9 * w; F1_fac = aF1_9 + bF1_9 * w;}
        else {f_fac = af_10 + bf_10 * w; F1_fac = aF1_10 + bF1_10 * w;}
        return (MM + MV) / 2. / (w * w - 1.) / MM / MV * ((w - MV_o_MM) * f_fac - F1_fac / MM);
    }
    else return 0.;
}

double MVlnu::A12(double q2)
{
    return (A1(q2)*(MM + MV)*(MM + MV)*(MM * MM - MV * MV - q2) -
            A2(q2)*(MM * MM * MM * MM + (MV * MV - q2)*(MV * MV - q2) -
            2. * MM * MM * (MV * MV + q2))) / (16. * MM * MV * MV * (MM + MV));
}

double MVlnu::T1(double q2)
{
    double delta_T1 = 0.;
    return (Mb + Mc) / (MM + MV) * V(q2)*(1. + delta_T1);
}

double MVlnu::T2(double q2)
{
    double delta_T2 = 0.;
    return (Mb - Mc) / (MM - MV) * A1(q2)*(1. + delta_T2);
}

double MVlnu::T23(double q2)
{
    double delta_T23 = 0.;
    return ((Mb - Mc)*((MM - MV)*(MM - MV) - q2)*((MM + MV)*(MM + MV) - q2) * A0(q2) +
            8 * MM * MV * (MV * MV - MM * MM) * A12(q2)) / (4. * MM * (MV - MM) * MV * q2)*(1. + delta_T23);
}
/********************************************************************************
 * Helicity amplitudes  (normalization such that all H \propto (mass scale)^-1) *
 * *****************************************************************************/

gslpp::complex MVlnu::HV0(double q2)
{
    return 4. * gslpp::complex::i() * MM * MV / (sqrt(q2)*(MM + MV))*((CV - CVp)*(MM + MV) * A12(q2) + Mb * (C7 - C7p) * T23(q2));
}

gslpp::complex MVlnu::HVp(double q2)
{
    return gslpp::complex::i()*((((CV + CVp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CV - CVp) * A1(q2))) / (2. * (MM + MV))
            + (Mb / q2)*((C7 + C7p) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(C7 - C7p)*(MM * MM - MV * MV) * T2(q2)));
}

gslpp::complex MVlnu::HVm(double q2)
{
    return gslpp::complex::i()*(((-(CV + CVp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CV - CVp) * A1(q2))) / (2. * (MM + MV))
            + (Mb / q2)*(-(C7 + C7p) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(C7 - C7p)*(MM * MM - MV * MV) * T2(q2)));
}

gslpp::complex MVlnu::HAp(double q2)
{
    return gslpp::complex::i()*((CA + CAp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CA - CAp) * A1(q2)) / (2. * (MM + MV));
}

gslpp::complex MVlnu::HAm(double q2)
{
    return gslpp::complex::i()*(-(CA + CAp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CA - CAp) * A1(q2)) / (2. * (MM + MV));
}

gslpp::complex MVlnu::HA0(double q2)
{
    return 4. * gslpp::complex::i() * MV * MM / (sqrt(q2))*(CA - CAp) * A12(q2);
}

gslpp::complex MVlnu::HP(double q2)
{
    return gslpp::complex::i() * lambda_half(MM*MM, MV*MV, q2) / 2. * ((CP - CPp) / (Mb + Mc)+(Mlep + Mnu) / q2 * (CA - CAp)) * A0(q2);
}

gslpp::complex MVlnu::HS(double q2)
{
    return gslpp::complex::i() * lambda_half(MM*MM, MV*MV, q2) / 2. * ((CS - CSp) / (Mb + Mc)+(Mlep - Mnu) / q2 * (CV - CVp)) * A0(q2);
}

gslpp::complex MVlnu::HT0(double q2)
{
    return 2. * M_SQRT2 * (MM * MV) / (MM + MV)*(CT + CTp) * T23(q2);
}

gslpp::complex MVlnu::HT0t(double q2)
{
    return 2. * (MM * MV) / (MM + MV)*(CT - CTp) * T23(q2);
}

gslpp::complex MVlnu::HTp(double q2)
{
    return ((CT - CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT + CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}

gslpp::complex MVlnu::HTpt(double q2)
{
    return ((CT + CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT - CTp)*(MM * MM - MV * MV) * T2(q2)) / (2. * sqrt(q2));
}

gslpp::complex MVlnu::HTm(double q2)
{
    return (-(CT - CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT + CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}

gslpp::complex MVlnu::HTmt(double q2)
{
    return (-(CT + CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT - CTp)*(MM * MM - MV * MV) * T2(q2)) / (2. * sqrt(q2));
}
/*******************************************************************************
 * Generalized angular coefficients  (see 1506.03970)                          *
 * ****************************************************************************/

gslpp::complex MVlnu::G000(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep * Mlep + lambda_lep2 / (4. * q2));
    double Enu = sqrt(Mnu * Mnu + lambda_lep2 / (4. * q2));
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * (4. / 9. * (3. * Elep * Enu + lambda_lep2 / (4. * q2))*(HVp(q2).abs2() + HVm(q2).abs2() + HV0(q2).abs2() + HAp(q2).abs2() + HAm(q2).abs2() + HA0(q2).abs2()) +
            4. * Mlep * Mnu / 3. * (HVp(q2).abs2() + HVm(q2).abs2() + HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() - HA0(q2).abs2()) +
            4. / 3. * ((Elep * Enu - Mlep * Mnu + lambda_lep2 / (4. * q2)) * HS(q2).abs2()+(Elep * Enu + Mlep * Mnu + lambda_lep2 / (4. * q2)) * HP(q2).abs2()) +
            16. / 9. * (3. * (Elep * Enu + Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTpt(q2).abs2() + HTmt(q2).abs2() + HT0t(q2).abs2()) +
            8. / 9. * (3. * (Elep * Enu - Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTp(q2).abs2() + HTm(q2).abs2() + HT0(q2).abs2()) +
            16. / 3. * (Mlep * Enu + Mnu * Elep)*(HVp(q2) * HTpt(q2).conjugate() + HVm(q2) * HTmt(q2).conjugate() + HV0(q2) * HT0t(q2).conjugate()).imag() +
            8. * M_SQRT2 / 3. * (Mlep * Enu - Mnu * Elep)*(HAp(q2) * HTp(q2).conjugate() + HAm(q2) * HTm(q2).conjugate() + HA0(q2) * HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G010(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / 3. * lambda_lep * ((HVp(q2) * HAp(q2).conjugate() - HVm(q2) * HAm(q2).conjugate()).real()+
            (2. * M_SQRT2) / q2 * (Mlep * Mlep - Mnu * Mnu)*(HTp(q2) * HTpt(q2).conjugate() - HTm(q2) * HTmt(q2).conjugate()).real() +
            2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HTpt(q2).conjugate() - HAm(q2) * HTmt(q2).conjugate()).imag() +
            M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HTp(q2).conjugate() - HVm(q2) * HTm(q2).conjugate()).imag()-
            (Mlep - Mnu) / sqrt(q2)*(HA0(q2) * HP(q2).conjugate()).real()-(Mlep + Mnu) / sqrt(q2)*(HV0(q2) * HS(q2).conjugate()).real()+
            (M_SQRT2 * HT0(q2) * HP(q2).conjugate() + 2. * HT0t(q2) * HS(q2).conjugate()).imag());

}

gslpp::complex MVlnu::G020(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 2. / 9. * lambda_lep2 / q2 * ((-HVp(q2).abs2() - HVm(q2).abs2() + 2. * HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() + 2. * HA0(q2).abs2()) -
            2. * (2. * HT0(q2).abs2() - HTp(q2).abs2() - HTm(q2).abs2()) - 4. * (2. * HT0t(q2).abs2() - HTpt(q2).abs2() - HTmt(q2).abs2()));
}

gslpp::complex MVlnu::G200(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep * Mlep + lambda_lep2 / (4. * q2));
    double Enu = sqrt(Mnu * Mnu + lambda_lep2 / (4. * q2));
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * (-4. / 9. * (3. * Elep * Enu + lambda_lep2 / (4. * q2))*(HVp(q2).abs2() + HVm(q2).abs2() - 2. * HV0(q2).abs2() + HAp(q2).abs2() + HAm(q2).abs2() - 2. * HA0(q2).abs2()) -
            4. / 3. * Mlep * Mnu * (HVp(q2).abs2() + HVm(q2).abs2() - 2. * HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() + 2. * HA0(q2).abs2()) +
            8. / 3. * (Elep * Enu - Mlep * Mnu + lambda_lep2 / (4. * q2)) * HS(q2).abs2() + 8. / 3. * (Elep * Enu + Mlep * Mnu + lambda_lep2 / (4. * q2)) * HP(q2).abs2() -
            16. / 9. * (3. * (Elep * Enu + Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTpt(q2).abs2() + HTmt(q2).abs2() - 2. * HT0t(q2).abs2()) - 8. / 9. * (3. * (Elep * Enu - Mlep * Mnu) - lambda_lep2 / (4. * q2))*
            (HTp(q2).abs2() + HTm(q2).abs2() - 2 * HT0(q2).abs2()) - 16. / 3. * (Mlep * Enu + Mnu * Elep)*(HVp(q2) * HTpt(q2).conjugate() + HVm(q2) * HTmt(q2).conjugate() - 2. * HV0(q2) * HT0t(q2).conjugate()).imag() -
            8. * M_SQRT2 / 3. * (Mlep * Enu - Mnu * Elep)*(HAp(q2) * HTp(q2).conjugate() + HAm(q2) * HTm(q2).conjugate() - 2. * HA0(q2) * HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G210(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 4. * lambda_lep / 3. * ((HVp(q2) * HAp(q2).conjugate() - HVm(q2) * HAm(q2).conjugate()).real() +
            2. * M_SQRT2 * (Mlep * Mlep - Mnu * Mnu) / q2 * (HTp(q2) * HTpt(q2).conjugate() - HTm(q2) * HTmt(q2).conjugate()).real() +
            2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HTpt(q2).conjugate() - HAm(q2) * HTmt(q2).conjugate()).imag() +
            M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HTp(q2).conjugate() - HVm(q2) * HTm(q2).conjugate()).imag() +
            2. * (Mlep - Mnu) / sqrt(q2)*(HA0(q2) * HP(q2).conjugate()).real() + 2. * (Mlep + Mnu) / sqrt(q2)*(HV0(q2) * HS(q2).conjugate()).real() -
            2. * (M_SQRT2 * HT0(q2) * HP(q2).conjugate() + 2. * HT0t(q2) * HS(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G220(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 2. / 9. * lambda_lep2 / q2 * (HVp(q2).abs2() + HVm(q2).abs2() + 4. * HV0(q2).abs2() + HAp(q2).abs2() +
            HAm(q2).abs2() + 4. * HA0(q2).abs2() - 2. * (HTp(q2).abs2() + HTm(q2).abs2() + 4. * HT0(q2).abs2()) -
            4. * (HTpt(q2).abs2() + HTmt(q2).abs2() + 4. * HT0t(q2).abs2()));
}

gslpp::complex MVlnu::G211(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / sqrt(3.) * lambda_lep * (HVp(q2) * HA0(q2).conjugate() + HAp(q2) * HV0(q2).conjugate() -
            HV0(q2) * HAm(q2).conjugate() - HA0(q2) * HVm(q2).conjugate()+(Mlep + Mnu) / sqrt(q2)*(HVp(q2) * HS(q2).conjugate() + HS(q2) * HVm(q2).conjugate()) -
            gslpp::complex::i() * M_SQRT2 * (HP(q2) * HTm(q2).conjugate() - HTp(q2) * HP(q2).conjugate() + M_SQRT2 * (HS(q2) * HTmt(q2).conjugate() - HTpt(q2) * HS(q2).conjugate()))+
            (Mlep - Mnu) / sqrt(q2)*(HAp(q2) * HP(q2).conjugate() + HP(q2) * HAm(q2).conjugate()) -
            gslpp::complex::i()*2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HT0t(q2).conjugate() + HT0t(q2) * HAm(q2).conjugate() - HTpt(q2) * HA0(q2).conjugate() - HA0(q2) * HTmt(q2).conjugate()) -
            gslpp::complex::i() * M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HT0(q2).conjugate() + HT0(q2) * HVm(q2).conjugate() - HTp(q2) * HV0(q2).conjugate() - HV0(q2) * HTm(q2).conjugate()) +
            2. * M_SQRT2 * (Mlep * Mlep - Mnu * Mnu) / q2 * (HTp(q2) * HT0t(q2).conjugate() + HTpt(q2) * HT0(q2).conjugate() - HT0(q2) * HTmt(q2).conjugate() - HT0t(q2) * HTm(q2).conjugate()));
}

gslpp::complex MVlnu::G221(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / 3. * lambda_lep2 / q2 * (HVp(q2) * HV0(q2).conjugate() + HV0(q2) * HVm(q2).conjugate() +
            HAp(q2) * HA0(q2).conjugate() + HA0(q2) * HAm(q2).conjugate() - 2. * (HTp(q2) * HT0(q2).conjugate() +
            HT0(q2) * HTm(q2).conjugate() + 2. * (HTpt(q2) * HT0t(q2).conjugate() + HT0t(q2) * HTmt(q2).conjugate())));
}

gslpp::complex MVlnu::G222(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 8. / 3. * lambda_lep2 / q2 * (HVp(q2) * HVm(q2).conjugate() + HAp(q2) * HAm(q2).conjugate() -
            2. * (HTp(q2) * HTm(q2).conjugate() + 2. * HTpt(q2) * HTmt(q2).conjugate()));
}
/***************************************************************************
 * 12 independent J angular coefficients  (see again 1506.03970)           *
 * ************************************************************************/

double MVlnu::J1s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (8. * G000(q2) + 2. * G020(q2) - 4. * G200(q2) - G220(q2)).real() / 3.;
}

double MVlnu::J1c(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (8. * G000(q2) + 2. * G020(q2) + 8. * G200(q2) + 2. * G220(q2)).real() / 3.;
}

double MVlnu::J2s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * G020(q2) - G220(q2)).real();
}

double MVlnu::J2c(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * (G020(q2) + G220(q2))).real();
}

double MVlnu::J3(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (G222(q2).real());
}

double MVlnu::J4(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (G221(q2).real());
}

double MVlnu::J5(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * G211(q2).real() / sqrt(3.));
}

double MVlnu::J6s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (4. * (2. * G010(q2) - G210(q2)).real() / 3.);
}

double MVlnu::J6c(double q2)
{
    if (q2 < Mlep * Mlep) return 0.;
    if (q2 > (MM - MV)*(MM - MV)) return 0.;
    return -amplsq_factor * (8. * (G010(q2) + G210(q2)).real() / 3.);
}

double MVlnu::J7(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (2. * sqrt(3.)*(G211(q2).imag()) / 3.);
}

double MVlnu::J8(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (G221(q2).imag());
}

double MVlnu::J9(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (G222(q2).imag());
}
/***************************************************************************
 * Integration of angular coefficients Js                                  *
 * ************************************************************************/

double MVlnu::integrateJ(int i, double q2_min, double q2_max)
{
    switch (i) {
        case 1:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J1s);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 2:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J1c);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 3:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J2s);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 4:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J2c);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 5:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J3);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 6:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J4);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 7:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_el;
              wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J5);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 8:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J6s);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 9:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_mu;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J6c);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 10:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J7);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 11:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J8);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 12:
            if (lep == StandardModel::TAU) if (checkcache_int_tau && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_tau;
            if (lep == StandardModel::MU) if (checkcache_int_mu && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_mu;
            if (lep == StandardModel::ELECTRON) if (checkcache_int_el && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_el;
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::J9);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVlnu::integrateJ: index " + out.str() + " not implemented");
    }
}

double MVlnu::dGammadw(double q2)
{
    updateParameters();

    return 3. / 4. * (2. * J1s(q2) + J1c(q2)) - 1. / 4. * (2. * J2s(q2) + J2c(q2));
}

double MVlnu::getDeltaGammaDeltaw(double w_min, double w_max)
{
    updateParameters();

    double q2_min = (2. * MM * MV)*(w0 - w_max); // min is Mlep*Mlep;
    double q2_max = (2. * MM * MV)*(w0 - w_min); // max is (MM-MV)*(MM-MV);

    double intJ1s = integrateJ(1, q2_min, q2_max);
    double intJ1c = integrateJ(2, q2_min, q2_max);
    double intJ2s = integrateJ(3, q2_min, q2_max);
    double intJ2c = integrateJ(4, q2_min, q2_max);

    return 3. / 4. * (2. * intJ1s + intJ1c) - 1. / 4. * (2. * intJ2s + intJ2c);
}

double MVlnu::dGammadcldq2(double q2, double cl)
{
    updateParameters();

    return 3. / 8. * ((J1s(q2) + 2. * J1c(q2)) + cl * (J6s(q2) + 2. * J6c(q2))+(2. * cl * cl - 1.)*(J2s(q2) + 2. * J2c(q2)));
}

double MVlnu::dGammadcl(double cl)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ6s = integrateJ(8, q2min, q2max);
    double intJ6c = integrateJ(9, q2min, q2max);

    return 3. / 8. * ((intJ1s + 2. * intJ1c) + cl * (intJ6s + 2. * intJ6c)+(2. * cl * cl - 1.)*(intJ2s + 2. * intJ2c));
}

double MVlnu::getDeltaGammaDeltacl(double cl_min, double cl_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ6s = integrateJ(8, q2min, q2max);
    double intJ6c = integrateJ(9, q2min, q2max);

    return 3. / 8. * ((cl_max - cl_min)*(intJ1c + 2. * intJ1s)+
            (cl_max * cl_max - cl_min * cl_min) / 2. * (intJ6c + 2. * intJ6s)+
            (2. / 3. * (cl_max * cl_max * cl_max - cl_min * cl_min * cl_min)-(cl_max - cl_min))*(intJ2c + 2. * intJ2s));
}

double MVlnu::dGammadcVdq2(double q2, double cl)
{
    updateParameters();

    return 3. / 8. * ((J1s(q2) + 2. * J1c(q2)) + cl * (J6s(q2) + 2. * J6c(q2))+(2. * cl * cl - 1.)*(J2s(q2) + 2. * J2c(q2)));
}

double MVlnu::dGammadcV(double cV)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);

    return 3. / 8. * (cV * cV * (3. * intJ1c - intJ2c)+(1. - cV * cV)*(3. * intJ1s - intJ2s));
}

double MVlnu::getDeltaGammaDeltacV(double cV_min, double cV_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);

    return 3. / 8. * ((cV_max * cV_max * cV_max - cV_min * cV_min * cV_min) / 3. * (3. * intJ1c - intJ2c)+
            ((cV_max - cV_min)-(cV_max * cV_max * cV_max - cV_min * cV_min * cV_min) / 3.)*(3. * intJ1s - intJ2s));
}

double MVlnu::dGammadchidq2(double q2, double chi)
{
    updateParameters();

    return (3. * J1c(q2) + 6. * J1s(q2) - J2c(q2) - 2. * J2s(q2)) / 8. / M_PI +
            cos(2. * chi) / 2. / M_PI * J3(q2) + sin(2. * chi) / 2. / M_PI * J9(q2);
}

double MVlnu::dGammadchi(double chi)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ3 = integrateJ(5, q2min, q2max);
    double intJ9 = integrateJ(12, q2min, q2max);

    return ((3. * intJ1c + 6. * intJ1s - intJ2c - 2. * intJ2s) / 4. +
            cos(2. * chi) * intJ3 + sin(2. * chi) * intJ9) / 2. / M_PI;
}

double MVlnu::getDeltaGammaDeltachi(double chi_min, double chi_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ3 = integrateJ(5, q2min, q2max);
    double intJ9 = integrateJ(12, q2min, q2max);

    return ((chi_max - chi_min)*(3. * intJ1c + 6. * intJ1s - intJ2c - 2. * intJ2s) / 4. +
            (sin(2. * chi_max) - sin(2. * chi_min)) / 2. * intJ3 -
            (cos(2. * chi_max) - cos(2. * chi_min)) / 2. * intJ9) / (2. * M_PI);
}

double MVlnu::getFL()
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    
    double DeltaJL = (3. * intJ1c - intJ2c) / 4.;
    double DeltaJ = 3. / 4. * (2. * intJ1s + intJ1c) - 1. / 4. * (2. * intJ2s + intJ2c);
    return DeltaJL/DeltaJ;
    
}

double MVlnu::get_unitarity_V_BGL()
{
    updateParameters();

    return ag0 * ag0 + ag1 * ag1 + ag2*ag2;

}

double MVlnu::get_unitarity_A_BGL()
{
    updateParameters();

    double aF10 = (MM - MV)*(phi_F1(0.) / phi_f(0.)) * af0;
    return af0 * af0 + af1 * af1 + af2 * af2 + aF10 * aF10 + aF11 * aF11 + aF12*aF12 + aF13*aF13;
}

double MVlnu::get_unitarity_P_BGL()
{
    updateParameters();

    double z0 = (sqrt(w0 + 1.) - M_SQRT2) / (sqrt(w0 + 1.) + M_SQRT2);
    double PfacF2z0 = (z0 - zP1) / (1. - z0 * zP1)*(z0 - zP2) / (1. - z0 * zP2)*(z0 - zP3) / (1. - z0 * zP3);
    double phiF2z0 = phi_F2(z0);
    double aF20 = PfacF2z0 * phiF2z0 * 2. * F1_BGL(0.) / (MM * MM - MV * MV) - aF21 * z0 - aF22 * z0*z0 - aF23 * z0*z0*z0;
    return aF20 * aF20 + aF21 * aF21 + aF22*aF22 + aF23*aF23;
}

double MVlnu::get_hA1w1()
{
    updateParameters();

    return hA1(q2max);
}

double MVlnu::get_hA1(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return hA1(q2);
}

double MVlnu::get_hA2(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return (hA1(q2) * (1. + w) - RV * (A0(q2) * (1. + MV_o_MM) - A2(q2) * (MV_o_MM - w))) /  (1. + MV_o_MM*MV_o_MM - 2.* MV_o_MM * w);
    // return (hA1(q2) * (1. + w) - A0(q2) * RV * (1. + MV_o_MM) - A2(q2) * RV * (w - MV_o_MM)) / (1. + MV_o_MM*MV_o_MM - 2.*MV_o_MM*w) ;
}

double MVlnu::get_hA3(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return A2(q2) * RV - MV_o_MM * get_hA2(w) ;
}

double MVlnu::get_hV(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return V(q2) * RV;
}

double MVlnu::get_R1(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R1(q2);
    else if (BGLflag) return V(q2) * RV / hA1(q2);
    return 0.;
}

double MVlnu::get_R2(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R2(q2);
    else if (BGLflag) return A2(q2) * RV / hA1(q2);
    return 0.;
}

double MVlnu::get_R0(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R0(q2);
    else if (BGLflag) return A0(q2) / A1(q2);
    return 0.;
}


/***************************************************************************
 * SM computation  ... lep hel asymmetry, see 1203.2654, 1707.09509        *
 * ************************************************************************/

double MVlnu::Hplus(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return (MM+MV)*A1(q2)-2.*MM/(MM+MV)*abs_p*V(q2);
}

double MVlnu::Hminus(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return (MM+MV)*A1(q2)+2.*MM/(MM+MV)*abs_p*V(q2);
}

double MVlnu::H0(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return ((MM*MM-MV*MV-q2)*(MM+MV)*A1(q2)-4.*MM*MM*abs_p*abs_p/(MM+MV)*A2(q2))/(2.*MV*sqrt(q2));
}

double MVlnu::H0t(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return 2.*MM*abs_p*A0(q2)/sqrt(q2);
}

double MVlnu::dGpdq2(double q2){
        
    updateParameters();
        
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    double lep_factor = (1.-Mlep*Mlep/q2)*(1.-Mlep*Mlep/q2)*Mlep*Mlep/(2.*q2);
    return 2./3.*amplsq_factor*abs_p*q2*lep_factor*(Hplus(q2)*Hplus(q2)+Hminus(q2)*Hminus(q2)+H0(q2)*H0(q2)
            +3.*H0t(q2)*H0t(q2));
}

double MVlnu::dGmdq2(double q2){
    
    updateParameters(); 
    
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    double lep_factor = (1.-Mlep*Mlep/q2)*(1.-Mlep*Mlep/q2);
    return 2./3.*amplsq_factor*abs_p*q2*lep_factor*(Hplus(q2)*Hplus(q2)+Hminus(q2)*Hminus(q2)+H0(q2)*H0(q2));
}

double MVlnu::integrateGpm(int i, double q2_min, double q2_max)
{
    switch (i) {
        case 1:
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::dGpdq2);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        case 2:
            wf=ROOT::Math::Functor1D(&(*this),&MVlnu::dGmdq2);
            ig.SetFunction(wf);
            J_res = ig.Integral(q2_min, q2_max);
            return J_res;
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVlnu::integrateGpm: index " + out.str() + " not implemented");
    }
}

double MVlnu::getPlep()
{
    updateParameters();
     
    double DeltaGammaPlus = integrateGpm(1,q2min,q2max);
    double DeltaGammaMinus = integrateGpm(2,q2min,q2max);
    
    return (DeltaGammaPlus-DeltaGammaMinus)/(DeltaGammaPlus+DeltaGammaMinus);
    
}