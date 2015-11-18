/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVll.h"
#include <gslpp_complex.h>
#include <gsl/gsl_sf.h>
#include <boost/bind.hpp>
#include <limits>
#include <TFitResult.h>

MVll::MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i)
: mySM(SM_i),
N_cache(3, 0.),
V_cache(3, 0.),
A0_cache(3, 0.),
A1_cache(3, 0.),
T1_cache(3, 0.),
T2_cache(3, 0.),
k2_cache(2, 0.),
VL0_cache(3, 0.),
TL0_cache(3, 0.),
SL_cache(2, 0.),
Ycache(2, 0.),
H_V0cache(2, 0.),
H_V1cache(2, 0.),
H_V2cache(2, 0.),
H_Scache(2, 0.),
H_Pcache(4, 0.),
T_cache(5, 0.) 
{

    
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    I0_updated = 0;
    I1_updated = 0;
    I2_updated = 0;
    I3_updated = 0;
    I4_updated = 0;
    I5_updated = 0;
    I6_updated = 0;
    I7_updated = 0;
    I8_updated = 0;
    I9_updated = 0;
    I10_updated = 0;
    I11_updated = 0;

    VL1_updated = 0;
    VL2_updated = 0;
    TL1_updated = 0;
    TL2_updated = 0;
    VR1_updated = 0;
    VR2_updated = 0;
    TR1_updated = 0;
    TR2_updated = 0;
    VL0_updated = 0;
    TL0_updated = 0;
    VR0_updated = 0;
    TR0_updated = 0;
    SL_updated = 0;
    SR_updated = 0;
    
    deltaTparpupdated = 0;
    deltaTparmupdated = 0;
    deltaTperpupdated = 0;

    w_sigma = gsl_integration_cquad_workspace_alloc (100);   
    w_DTPPR = gsl_integration_cquad_workspace_alloc (100);
    w_delta = gsl_integration_cquad_workspace_alloc (100);

}

MVll::~MVll() 
{
}

void MVll::updateParameters() 
{
    if (!mySM.getMyFlavour()->getUpdateFlag(meson, vectorM, lep)) return;


    GF = mySM.getGF();
    ale = mySM.getAle();
    Mlep = mySM.getLeptons(lep).getMass();
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    mu_b = mySM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MW = mySM.Mw();
    lambda_t = mySM.computelamt_s();
    width = mySM.getMesons(meson).computeWidth();

    switch (vectorM) {
        case StandardModel::K_star:
            a_0V = mySM.geta_0V();
            a_1V = mySM.geta_1V();
            a_2V = mySM.geta_2V();
            MRV_2 = mySM.getMRV() * mySM.getMRV();

            a_0A0 = mySM.geta_0A0();
            a_1A0 = mySM.geta_1A0();
            a_2A0 = mySM.geta_2A0();
            MRA0_2 = mySM.getMRA0() * mySM.getMRA0();

            a_0A1 = mySM.geta_0A1();
            a_1A1 = mySM.geta_1A1();
            a_2A1 = mySM.geta_2A1();
            MRA1_2 = mySM.getMRA1() * mySM.getMRA1();

            a_0A12 = mySM.geta_0A12();
            a_1A12 = mySM.geta_1A12();
            a_2A12 = mySM.geta_2A12();
            MRA12_2 = mySM.getMRA12() * mySM.getMRA12();

            a_0T1 = mySM.geta_0T1();
            a_1T1 = mySM.geta_1T1();
            a_2T1 = mySM.geta_2T1();
            MRT1_2 = mySM.getMRT1() * mySM.getMRT1();

            a_0T2 = mySM.geta_0T2();
            a_1T2 = mySM.geta_1T2();
            a_2T2 = mySM.geta_2T2();
            MRT2_2 = mySM.getMRT2() * mySM.getMRT2();

            a_0T23 = mySM.geta_0T23();
            a_1T23 = mySM.geta_1T23();
            a_2T23 = mySM.geta_2T23();
            MRT23_2 = mySM.getMRT23() * mySM.getMRT23();

            b = 1;
            break;
        case StandardModel::PHI:
            a_0V = mySM.geta_0Vphi();
            a_1V = mySM.geta_1Vphi();
            a_2V = mySM.geta_2Vphi();
            MRV_2 = mySM.getMRVphi() * mySM.getMRVphi();

            a_0A0 = mySM.geta_0A0phi();
            a_1A0 = mySM.geta_1A0phi();
            a_2A0 = mySM.geta_2A0phi();
            MRA0_2 = mySM.getMRA0phi() * mySM.getMRA0phi();

            a_0A1 = mySM.geta_0A1phi();
            a_1A1 = mySM.geta_1A1phi();
            a_2A1 = mySM.geta_2A1phi();
            MRA1_2 = mySM.getMRA1phi() * mySM.getMRA1phi();

            a_0A12 = mySM.geta_0A12phi();
            a_1A12 = mySM.geta_1A12phi();
            a_2A12 = mySM.geta_2A12phi();
            MRA12_2 = mySM.getMRA12phi() * mySM.getMRA12phi();

            a_0T1 = mySM.geta_0T1phi();
            a_1T1 = mySM.geta_1T1phi();
            a_2T1 = mySM.geta_2T1phi();
            MRT1_2 = mySM.getMRT1phi() * mySM.getMRT1phi();

            a_0T2 = mySM.geta_0T2phi();
            a_1T2 = mySM.geta_1T2phi();
            a_2T2 = mySM.geta_2T2phi();
            MRT2_2 = mySM.getMRT2phi() * mySM.getMRT2phi();

            a_0T23 = mySM.geta_0T23phi();
            a_1T23 = mySM.geta_1T23phi();
            a_2T23 = mySM.geta_2T23phi();
            MRT23_2 = mySM.getMRT23phi() * mySM.getMRT23phi();

            b = 0.489;
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVll: vector " + out.str() + " not implemented");
    }


    h_0[0] = mySM.geth_0();
    h_0[1] = mySM.geth_p();
    h_0[2] = mySM.geth_m();

    h_1[0] = mySM.geth_0_1();
    h_1[1] = mySM.geth_p_1();
    h_1[2] = mySM.geth_m_1();

    h_2[0] = mySM.geth_0_2();
    h_2[1] = mySM.geth_p_2();
    h_2[2] = mySM.geth_m_2();

    allcoeff = mySM.getMyFlavour()->ComputeCoeffBMll(mu_b); //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getMyFlavour()->ComputeCoeffprimeBMll(mu_b); //check the mass scale, scheme fixed to NDR

    C_1 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
    C_1L_bar = (*(allcoeff[LO]))(0)/2.;
    C_2 = ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
    C_2L_bar = (*(allcoeff[LO]))(1) - (*(allcoeff[LO]))(0)/6.;
    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
    C_4 = ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3));
    C_5 = ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));
    C_7 = ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6));
    C_8L = (*(allcoeff[LO]))(7);
    C_9 = ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8));
    C_10 = ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
    C_S = ((*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10));
    C_P = ((*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11));

    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    C_9p = (*(allcoeffprime[LO]))(8) + (*(allcoeffprime[NLO]))(8);
    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    allcoeffh = mySM.getMyFlavour()->ComputeCoeffBMll(mu_h); //check the mass scale, scheme fixed to NDR

    C_1Lh_bar = (*(allcoeffh[LO]))(0)/2.;
    C_2Lh_bar = (*(allcoeffh[LO]))(1) - (*(allcoeff[LO]))(0)/6.;
    C_8Lh = (*(allcoeffh[LO]))(7);
    
    checkCache();

    t_p = pow(MM + MV, 2.);
    t_m = pow(MM - MV, 2.);
    t_0 = t_p * (1. - sqrt(1. - t_m / t_p)); /*Modify it for Lattice*/
    z_0 = (sqrt(t_p) - sqrt(t_p - t_0)) / (sqrt(t_p) + sqrt(t_p - t_0));
    //else t_0 = 12.;
    MMpMV = MM + MV;
    MMpMV2 = MMpMV * MMpMV;
    MMmMV = MM - MV;
    MMmMV2 = MMmMV * MMmMV;
    MM2 = MM*MM;
    MM4 = MM2*MM2;
    MV2 = MV*MV;
    MV4 = MV2*MV2;
    MMMV = MM*MV;
    MM2mMV2 = MM2 - MV2;
    fourMV = 4. * MV;
    twoMM2 = 2. * MM2;
    twoMV2 = 2. * MV2;
    onepMMoMV = (1. + MV / MM);
    MM_MMpMV = MM * MMpMV;
    twoMM_mbpms = 2. * MM * (Mb + Ms);
    fourMM2 = 4. * MM2;
    Mlep2 = Mlep*Mlep;
    twoMlepMb = 2. * Mlep*Mb;
    MboMW = Mb / MW;
    MsoMb = Ms / Mb;
    ninetysixM_PI3MM3 = 96. * M_PI * M_PI * M_PI * MM * MM*MM;
    sixteenM_PI2 = 16. * M_PI*M_PI;
    sixteenM_PI2MM2 = sixteenM_PI2 * MM*MM;
    twoMboMM = 2 * Mb / MM;
    H_0_pre = 8. / 27. + 4. / 9. * gslpp::complex::i() * M_PI;
    H_0_WC = (C_3 + 4. / 3. * C_4 + 16. * C_5 + 64. / 3. * C_6);
    H_c_WC = (4. / 3. * C_1 + C_2 + 6. * C_3 + 60. * C_5);
    H_b_WC = (7. * C_3 + 4. / 3. * C_4 + 76. * C_5 + 64. / 3. * C_6);
    mu_b2 = mu_b*mu_b;
    Mc2 = Mc*Mc;
    Mb2 = Mb*Mb;
    fourMc2 = 4. * Mc2;
    fourMb2 = 4. * Mb2;
    logMc = log(Mc2 / mu_b2);
    logMb = log(Mb2 / mu_b2);
    fournineth = 4. / 9.;
    half = 1. / 2.;
    twothird = 2. / 3.;
    ihalfMPI = gslpp::complex::i() * M_PI / 2.;
    twoMM3 = 2. * MM2 * MM;
    C2_inv = 1. / (2. * C_2.real());
    gtilde_1_pre = -16. * pow(MM, 3.)*(MM + MV) * pow(M_PI, 2.);
    gtilde_2_pre = -16. * pow(MM, 3.) * pow(M_PI, 2.) / MMpMV;
    gtilde_3_pre = 64. * pow(MM, 3.) * pow(M_PI, 2.) * MV*MMpMV;
    S_L_pre = (-2. * MM * (Mb + Ms));
    
    M_PI2osix = M_PI * M_PI / 6.;
    twoMM = 2.*MM;
    m4downcharge = -4. * mySM.getQuarks(QCD::DOWN).getCharge();
    threeGegen0 = mySM.getMesons(vectorM).getGegenalpha(0)*3.;
    threeGegen1otwo = mySM.getMesons(vectorM).getGegenalpha(1)*3./2.;
    twoMc2 = 2.*Mc2;

    sixMMoMb = 6. * MM / Mb;
    CF =4./3.;
    
    deltaT_0 = mySM.Als(mu_b) * CF / 4. / M_PI;
    deltaT_1par = mySM.Als(mu_h) * CF / 4. * M_PI / 3. * mySM.getMesons(meson).getDecayconst() *
            mySM.getMesons(vectorM).getDecayconst() / MM; 
    deltaT_1perp = mySM.Als(mu_h) * CF / 4. * M_PI / 3. * mySM.getMesons(meson).getDecayconst() *
            mySM.getFKstarp() / MM; 
            
    F87_0=-32. / 9. * log(mu_b / Mb) + 8. / 27. * M_PI * M_PI - 44. / 9. - 8. / 9. * gslpp::complex::i() * M_PI;
    F87_1 = (4. / 3. * M_PI * M_PI - 40. / 3.);
    F87_2 = (32. / 9. * M_PI * M_PI - 316. / 9.);
    F87_3 = (200. / 27. * M_PI * M_PI - 658. / 9.);

    F89_0 = 104. / 9. - 32. / 27. * M_PI * M_PI ;
    F89_1 = 1184. / 27. - 40. / 9. * M_PI * M_PI;
    F89_2 = (-32. / 3. * M_PI * M_PI + 14212. / 135.);
    F89_3 = (-560. / 27. * M_PI * M_PI + 193444. / 945.);

    F29_0 = (256./243. - 32./81.*gslpp::complex::i()*M_PI - 128./9.*log(Mc/Mb))*log(mu_b/Mb) + 512./81.*log(mu_b/Mb)*log(mu_b/Mb) + 5.4082 - 1.0934 * gslpp::complex::i();
    F29_L1 = 32./81.*log(mu_b/Mb) + (0.48576 + 0.31119 * gslpp::complex::i());
    F29_1 = (-32./405. + 64./45./Mc2*Mb2)*log(mu_b/Mb) + (1.9061 + 0.80843 * gslpp::complex::i());
    F29_2 = (-8./945. + 16./105./Mc2*Mb2/Mc2*Mb2)*log(mu_b/Mb) + (-1.8286 + 2.8428 * gslpp::complex::i());
    F29_3 = (-32./25515. + 64./2835./Mc2*Mb2/Mc2*Mb2/Mc2*Mb2)*log(mu_b/Mb) + (-12.113 + 8.1251 * gslpp::complex::i());
    F29_L1_1 = (0.21951 - 0.14852 * gslpp::complex::i());
    F29_L1_2 = (0.13015 - 0.22155 * gslpp::complex::i());
    F29_L1_3 = (-0.079692 - 0.31214 * gslpp::complex::i());

    F27_0 = 416./81. *log(mu_b/Mb) + 3.8367 + 0.3531 * gslpp::complex::i();
    F27_1 = (1.3098 + 0.60185 * gslpp::complex::i());
    F27_2 = (0.13507 + 0.89014 * gslpp::complex::i());
    F27_3 = (-1.0271 + 0.77168 * gslpp::complex::i());
    F27_L1_1 = (-0.031936 - 0.10981 * gslpp::complex::i());
    F27_L1_2 = (-0.14169 - 0.035553 * gslpp::complex::i());
    F27_L1_3 = (-0.13592 + 0.093 * gslpp::complex::i()); 

    a_0A12_LCSR = a_0A0 * (MM2mMV2) / (8. * MMMV);
    a_0T2_LCSR = a_0T1;

    NN = ((4. * GF * MM * ale * lambda_t) / (sqrt(2.)*4. * M_PI)).abs2();
    
    if (mySM.getMyFlavour()->getUpdateFlag(meson, vectorM, lep)) {
        switch (lep) {
            case StandardModel::MU:
                fit_DeltaC9_p_mumu();
                fit_DeltaC9_m_mumu();
                fit_DeltaC9_0_mumu();
                break;
            case StandardModel::ELECTRON:
                fit_DeltaC9_p_ee();
                fit_DeltaC9_m_ee();
                fit_DeltaC9_0_ee();
                break;
            default:
                std::stringstream out;
                out << lep;
                throw std::runtime_error("MVll: lepton " + out.str() + " not implemented");
        }
    }

    std::map<std::pair<double, double>, unsigned int >::iterator it;

    if (I0_updated == 0) for (it = sigma0Cached.begin(); it != sigma0Cached.end(); ++it) it->second = 0;
    if (I1_updated == 0) for (it = sigma1Cached.begin(); it != sigma1Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = sigma2Cached.begin(); it != sigma2Cached.end(); ++it) it->second = 0;
    if (I3_updated == 0) for (it = sigma3Cached.begin(); it != sigma3Cached.end(); ++it) it->second = 0;
    if (I4_updated == 0) for (it = sigma4Cached.begin(); it != sigma4Cached.end(); ++it) it->second = 0;
    if (I5_updated == 0) for (it = sigma5Cached.begin(); it != sigma5Cached.end(); ++it) it->second = 0;
    if (I6_updated == 0) for (it = sigma6Cached.begin(); it != sigma6Cached.end(); ++it) it->second = 0;
    if (I7_updated == 0) for (it = sigma7Cached.begin(); it != sigma7Cached.end(); ++it) it->second = 0;
    if (I9_updated == 0) for (it = sigma9Cached.begin(); it != sigma9Cached.end(); ++it) it->second = 0;
    if (I10_updated == 0) for (it = sigma10Cached.begin(); it != sigma10Cached.end(); ++it) it->second = 0;
    if (I11_updated == 0) for (it = sigma11Cached.begin(); it != sigma11Cached.end(); ++it) it->second = 0;

    if (I0_updated == 0) for (it = delta0Cached.begin(); it != delta0Cached.end(); ++it) it->second = 0;
    if (I1_updated == 0) for (it = delta1Cached.begin(); it != delta1Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = delta2Cached.begin(); it != delta2Cached.end(); ++it) it->second = 0;
    if (I3_updated == 0) for (it = delta3Cached.begin(); it != delta3Cached.end(); ++it) it->second = 0;
    if (I11_updated == 0) for (it = delta11Cached.begin(); it != delta11Cached.end(); ++it) it->second = 0;
    
    std::map<double, unsigned int >::iterator iti;
    if (deltaTparpupdated == 0) for (iti = deltaTparpCached.begin(); iti != deltaTparpCached.end(); ++iti) iti->second = 0;
    if (deltaTparmupdated == 0) for (iti = deltaTparmCached.begin(); iti != deltaTparmCached.end(); ++iti) iti->second = 0;
    if (deltaTperpupdated == 0) for (iti = deltaTparpCached.begin(); iti != deltaTparpCached.end(); ++iti) iti->second = 0;
    
    if (deltaTparpupdated*deltaTparmupdated == 0) for (it = I1Cached.begin(); it != I1Cached.end(); ++it) it->second = 0;

    mySM.getMyFlavour()->setUpdateFlag(meson, vectorM, lep, false);
    
    return;
}

void MVll::checkCache() 
{

    if (MM == k2_cache(0) && MV == k2_cache(1)) {
        k2_updated = 1;
        z_updated = 1;
    } else {
        k2_updated = 0;
        z_updated = 0;
        k2_cache(0) = MM;
        k2_cache(1) = MV;
    }

    if (Mlep == beta_cache) {
        beta_updated = 1;
    } else {
        beta_updated = 0;
        beta_cache = Mlep;
    }

    lambda_updated = k2_updated;
    F_updated = lambda_updated * beta_updated;

    if (GF == N_cache(0) && ale == N_cache(1) && MM == N_cache(2) && lambda_t == Nc_cache) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = ale;
        N_cache(2) = MM;
        Nc_cache = lambda_t;
    }

    if (a_0V == V_cache(0) && a_1V == V_cache(1) && a_2V == V_cache(2)) {
        V_updated = V_updated * z_updated;
    } else {
        V_updated = 0;
        V_cache(0) = a_0V;
        V_cache(1) = a_1V;
        V_cache(2) = a_2V;
    }

    if (a_0A0 == A0_cache(0) && a_1A0 == A0_cache(1) && a_2A0 == A0_cache(2)) {
        A0_updated = A0_updated * z_updated;
    } else {
        A0_updated = 0;
        A0_cache(0) = a_0A0;
        A0_cache(1) = a_1A0;
        A0_cache(2) = a_2A0;
    }

    if (a_0A1 == A1_cache(0) && a_1A1 == A1_cache(1) && a_2A1 == A1_cache(2)) {
        A1_updated = A1_updated * z_updated;
    } else {
        A1_updated = 0;
        A1_cache(0) = a_0A1;
        A1_cache(1) = a_1A1;
        A1_cache(2) = a_2A1;
    }

    if (a_0T1 == T1_cache(0) && a_1T1 == T1_cache(1) && a_2T1 == T1_cache(2)) {
        T1_updated = T1_updated * z_updated;
    } else {
        T1_updated = 0;
        T1_cache(0) = a_0T1;
        T1_cache(1) = a_1T1;
        T1_cache(2) = a_2T1;
    }

    if (a_0T2 == T2_cache(0) && a_1T2 == T2_cache(1) && a_2T2 == T2_cache(2)) {
        T2_updated = T2_updated * z_updated;
    } else {
        T2_updated = 0;
        T2_cache(0) = a_0T2;
        T2_cache(1) = a_1T2;
        T2_cache(2) = a_2T2;
    }

    VL1_updated = k2_updated * lambda_updated * A1_updated * V_updated;
    VL2_updated = VL1_updated;

    TL1_updated = k2_updated * lambda_updated * T1_updated * T2_updated;
    TL2_updated = TL1_updated;

    VR1_updated = VL2_updated;
    VR2_updated = VL1_updated;

    TR1_updated = TL2_updated;
    TR2_updated = TL1_updated;

    if (Mb == SL_cache(0) && Ms == SL_cache(1)) {
        Mb_Ms_updated = 1;
        SL_updated = lambda_updated * A0_updated;
        SR_updated = SL_updated;
    } else {
        Mb_Ms_updated = 0;
        SL_updated = 0;
        SR_updated = SL_updated;
        SL_cache(0) = Mb;
        SL_cache(1) = Ms;
    }

    if (a_0A12 == VL0_cache(0) && a_1A12 == VL0_cache(1) && a_2A12 == VL0_cache(2)) {
        VL0_updated = VL0_updated * z_updated;
        VR0_updated = VL0_updated;
    } else {
        VL0_updated = 0;
        VR0_updated = VL0_updated;
        VL0_cache(0) = a_0A12;
        VL0_cache(1) = a_1A12;
        VL0_cache(2) = a_2A12;
    }

    if (a_0T23 == TL0_cache(0) && a_1T23 == TL0_cache(1) && a_2T23 == TL0_cache(2)) {
        TL0_updated = TL0_updated * z_updated;
        TR0_updated = TL0_updated;
    } else {
        TL0_updated = 0;
        TR0_updated = TL0_updated;
        TL0_cache(0) = a_0T23;
        TL0_cache(1) = a_1T23;
        TL0_cache(2) = a_2T23;
    }


    if (C_1 == C_1_cache) {
        C_1_updated = 1;
    } else {
        C_1_updated = 0;
        C_1_cache = C_1;
    }

    if (C_2 == C_2_cache) {
        C_2_updated = 1;
    } else {
        C_2_updated = 0;
        C_2_cache = C_2;
    }

    if (C_3 == C_3_cache) {
        C_3_updated = 1;
    } else {
        C_3_updated = 0;
        C_3_cache = C_3;
    }

    if (C_4 == C_4_cache) {
        C_4_updated = 1;
    } else {
        C_4_updated = 0;
        C_4_cache = C_4;
    }

    if (C_5 == C_5_cache) {
        C_5_updated = 1;
    } else {
        C_5_updated = 0;
        C_5_cache = C_5;
    }

    if (C_6 == C_6_cache) {
        C_6_updated = 1;
    } else {
        C_6_updated = 0;
        C_6_cache = C_6;
    }

    if (C_7 == C_7_cache) {
        C_7_updated = 1;
    } else {
        C_7_updated = 0;
        C_7_cache = C_7;
    }

    if (C_9 == C_9_cache) {
        C_9_updated = 1;
    } else {
        C_9_updated = 0;
        C_9_cache = C_9;
    }

    if (C_10 == C_10_cache) {
        C_10_updated = 1;
    } else {
        C_10_updated = 0;
        C_10_cache = C_10;
    }

    if (C_S == C_S_cache) {
        C_S_updated = 1;
    } else {
        C_S_updated = 0;
        C_S_cache = C_S;
    }

    if (C_P == C_P_cache) {
        C_P_updated = 1;
    } else {
        C_P_updated = 0;
        C_P_cache = C_P;
    }

    if (C_7p == C_7p_cache) {
        C_7p_updated = 1;
    } else {
        C_7p_updated = 0;
        C_7p_cache = C_7p;
    }

    if (C_9p == C_9p_cache) {
        C_9p_updated = 1;
    } else {
        C_9p_updated = 0;
        C_9p_cache = C_9p;
    }

    if (C_10p == C_10p_cache) {
        C_10p_updated = 1;
    } else {
        C_10p_updated = 0;
        C_10p_cache = C_10p;
    }

    if (C_Sp == C_Sp_cache) {
        C_Sp_updated = 1;
    } else {
        C_Sp_updated = 0;
        C_Sp_cache = C_Sp;
    }

    if (C_Pp == C_Pp_cache) {
        C_Pp_updated = 1;
    } else {
        C_Pp_updated = 0;
        C_Pp_cache = C_Pp;
    }

    if (C_2Lh_bar == C_2Lh_cache) {
        C_2Lh_updated = 1;
    } else {
        C_2Lh_updated = 0;
        C_2Lh_cache = C_2Lh_bar;
    }

    if (C_8Lh == C_8Lh_cache) {
        C_8Lh_updated = 1;
    } else {
        C_8Lh_updated = 0;
        C_8Lh_cache = C_8Lh;
    }

    if (Mb == Ycache(0) && Mc == Ycache(1)) {
        Yupdated = C_1_updated * C_2_updated * C_3_updated * C_4_updated * C_5_updated * C_6_updated;
    } else {
        Yupdated = 0;
        Ycache(0) = Mb;
        Ycache(1) = Mc;
    }

    if (h_0[0] == h0Ccache[0] && h_1[0] == h0Ccache[1] && h_2[0] == h0Ccache[2]) {
        h0_updated = 1;
    } else {
        h0_updated = 0;
        h0Ccache[0] = h_0[0];
        h0Ccache[1] = h_1[0];
        h0Ccache[2] = h_2[0];
    }

    if (h_0[1] == h1Ccache[0] && h_1[1] == h1Ccache[1] && h_2[1] == h1Ccache[2]) {
        h1_updated = 1;
    } else {
        h1_updated = 0;
        h1Ccache[0] = h_0[1];
        h1Ccache[1] = h_1[1];
        h1Ccache[2] = h_2[1];
    }

    if (h_0[2] == h2Ccache[0] && h_1[2] == h2Ccache[1] && h_2[2] == h2Ccache[2]) {
        h2_updated = 1;
    } else {
        h2_updated = 0;
        h2Ccache[0] = h_0[2];
        h2Ccache[1] = h_1[2];
        h2Ccache[2] = h_2[2];
    }

    if (MM == H_V0cache(0) && Mb == H_V0cache(1)) {
        H_V0updated = N_updated * C_9_updated * Yupdated * VL0_updated * C_9p_updated * VR0_updated * C_7_updated * TL0_updated * C_7p_updated * TR0_updated * h0_updated;
    } else {
        H_V0updated = 0;
        H_V0cache(0) = MM;
        H_V0cache(1) = Mb;
    }

    if (MM == H_V1cache(0) && Mb == H_V1cache(1)) {
        H_V1updated = N_updated * C_9_updated * Yupdated * VL1_updated * C_9p_updated * VR1_updated * C_7_updated * TL1_updated * C_7p_updated * TR1_updated * h1_updated;
    } else {
        H_V1updated = 0;
        H_V1cache(0) = MM;
        H_V1cache(1) = Mb;
    }

    if (MM == H_V2cache(0) && Mb == H_V2cache(1)) {
        H_V2updated = N_updated * C_9_updated * Yupdated * VL2_updated * C_9p_updated * VR2_updated * C_7_updated * TL2_updated * C_7p_updated * TR2_updated * h2_updated;
    } else {
        H_V2updated = 0;
        H_V2cache(0) = MM;
        H_V2cache(1) = Mb;
    }

    H_A0updated = N_updated * C_10_updated * VL0_updated * C_10p_updated * VR0_updated;
    H_A1updated = N_updated * C_10_updated * VL1_updated * C_10p_updated * VR1_updated;
    H_A2updated = N_updated * C_10_updated * VL2_updated * C_10p_updated * VR2_updated;

    if (Mb == H_Scache(0) && MW == H_Scache(1)) {
        H_Supdated = N_updated * C_S_updated * SL_updated * C_Sp_updated * SR_updated;
    } else {
        H_Supdated = 0;
        H_Scache(0) = Mb;
        H_Scache(1) = MW;
    }

    if (Mb == H_Pcache(0) && MW == H_Pcache(1) && Mlep == H_Pcache(2) && Ms == H_Pcache(3)) {
        H_Pupdated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SR_updated * C_10_updated * C_10p_updated;
    } else {
        H_Pupdated = 0;
        H_Pcache(0) = Mb;
        H_Pcache(1) = MW;
        H_Pcache(2) = Mlep;
        H_Pcache(3) = Ms;

    }
    
        if (MM == T_cache(0) && Mb == T_cache(1) && Mc == T_cache(2) && 
            mySM.getMesons(vectorM).getGegenalpha(0) == T_cache(3) && mySM.getMesons(vectorM).getGegenalpha(1) == T_cache(4) ) {
        T_updated = 1;
    } else {
        T_updated = 0;
        T_cache(0) = MM;
        T_cache(1) = Mb;
        T_cache(2) = Mc;
        T_cache(3) = mySM.getMesons(vectorM).getGegenalpha(0);
        T_cache(4) = mySM.getMesons(vectorM).getGegenalpha(1);
    }

    deltaTparpupdated = C_2Lh_updated * T_updated;
    deltaTparmupdated = C_2Lh_updated * C_8Lh_updated * T_updated;
    deltaTperpupdated = deltaTparpupdated;

    I0_updated = F_updated * H_V0updated * H_A0updated * H_Pupdated * beta_updated * H_Supdated * deltaTparmupdated;
    I1_updated = F_updated * beta_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * deltaTparmupdated;
    I2_updated = F_updated * beta_updated * H_V0updated * H_A0updated  * deltaTparmupdated;
    I3_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * beta_updated * deltaTparmupdated;
    I4_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * deltaTparmupdated;
    I5_updated = F_updated * H_V0updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * beta_updated * deltaTparmupdated;
    I6_updated = F_updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * H_V0updated * beta_updated * H_Supdated * deltaTparmupdated;
    I7_updated = I4_updated * beta_updated ;
    I8_updated = F_updated * beta_updated * H_Supdated * H_V0updated * deltaTparmupdated;
    I9_updated = I6_updated;
    I10_updated = I5_updated;
    I11_updated = I7_updated;

}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/

double MVll::FF_fit(double q2, double a_0, double a_1, double a_2, double MR_2) 
{
    return 1. / (1. - q2 / MR_2) * (a_0 + a_1 * (z(q2) - z_0) + a_2 * (z(q2) - z_0) * (z(q2) - z_0));
}

double MVll::z(double q2) 
{
    return ( sqrt(t_p - q2) - sqrt(t_p - t_0)) / (sqrt(t_p - q2) + sqrt(t_p - t_0));
}

double MVll::V(double q2) 
{
    return FF_fit(q2, a_0V, a_1V, a_2V, MRV_2);
}

double MVll::A_0(double q2) 
{
    return FF_fit(q2, a_0A0, a_1A0, a_2A0, MRA0_2);
}

double MVll::A_1(double q2) 
{
    return FF_fit(q2, a_0A1, a_1A1, a_2A1, MRA1_2);
}

double MVll::A_2(double q2) 
{
    return (MMpMV2 * (MM2mMV2 - q2) * A_1(q2) - 16. * MM * MV2 * MMpMV * FF_fit(q2, a_0A12_LCSR, a_1A12, a_2A12, MRA12_2)) / lambda(q2);
}

double MVll::T_1(double q2) 
{
    return FF_fit(q2, a_0T1, a_1T1, a_2T1, MRT1_2);
}

double MVll::T_2(double q2) 
{
    return FF_fit(q2, a_0T2_LCSR, a_1T2, a_2T2, MRT2_2);
}

double MVll::V_0t(double q2) 
{
    return fourMV / sqrt(q2) * FF_fit(q2, a_0A12_LCSR, a_1A12, a_2A12, MRA12_2);
}

double MVll::V_p(double q2) 
{
    return half * (onepMMoMV * A_1(q2) - sqrt(lambda(q2)) / (MM_MMpMV) * V(q2));
}

double MVll::V_m(double q2) 
{
    return half * (onepMMoMV * A_1(q2) + sqrt(lambda(q2)) / (MM_MMpMV) * V(q2));
}

double MVll::T_0t(double q2) 
{
    return 2 * sqrt(q2) * MV / MM_MMpMV * FF_fit(q2, a_0T23, a_1T23, a_2T23, MRT23_2);
}

double MVll::T_p(double q2) 
{
    return (MM2mMV2 * T_2(q2) - sqrt(lambda(q2)) * T_1(q2)) / twoMM2;
}

double MVll::T_m(double q2) 
{
    return (MM2mMV2 * T_2(q2) + sqrt(lambda(q2)) * T_1(q2)) / twoMM2;
}

double MVll::S_L(double q2) 
{
    return -sqrt(lambda(q2)) / twoMM_mbpms * A_0(q2);
}

/*******************************************************************************
 * QCD factorization perturbative corrections                                  *
 ******************************************************************************/

gslpp::complex MVll::I1(double u, double q2)
{    
    std::pair<double, double > uq2 = std::make_pair(u,q2);

    if (I1Cached[uq2] == 0) {
        ubar = 1. - u;
        xp = .5 + sqrt(0.25 - (Mc2 - gslpp::complex::i()*1.e-10) / (ubar * MM2 + u * q2));
        xm = .5 - sqrt(0.25 - (Mc2 - gslpp::complex::i()*1.e-10) / (ubar * MM2 + u * q2));
        yp = .5 + sqrt(0.25 - (Mc2 - gslpp::complex::i()*1.e-10) / q2);
        ym = .5 - sqrt(0.25 - (Mc2 - gslpp::complex::i()*1.e-10) / q2);
        L1xp = log(1. - 1. / xp) * log(1. - xp) - M_PI2osix + dilog(xp / (xp - 1.));
        L1xm = log(1. - 1. / xm) * log(1. - xm) - M_PI2osix + dilog(xm / (xm - 1.));
        L1yp = log(1. - 1. / yp) * log(1. - yp) - M_PI2osix + dilog(yp / (yp - 1.));
        L1ym = log(1. - 1. / ym) * log(1. - ym) - M_PI2osix + dilog(ym / (ym - 1.));

        cacheI1[uq2] = 1. + twoMc2 / ubar / (MM2 - q2)*(L1xp + L1xm - L1yp - L1ym);
        I1Cached[uq2] = 1;
    }
    
    return cacheI1[uq2];   
}

gslpp::complex MVll::Tperpplus(double u, double q2) 
{
    Ee = (MM2 - q2) / twoMM;
    ubar = 1. - u;
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ (ubar * MM2 + u * q2) - 1.;
    B01 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ q2 - 1.;
    B00 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    
    gslpp::complex tperp = twoMM / Ee / ubar * I1(u, q2) + q2 / Ee / Ee / ubar / ubar * (B01 - B00);
    return m4downcharge * C_8Lh / (u + ubar * q2 / MM2) + MM / 2. / Mb *
            mySM.getQuarks(QCD::UP).getCharge() * tperp * C_2Lh_bar;
}

gslpp::complex MVll::Tparplus(double u, double q2) 
{
    Ee = (MM2 - q2) / twoMM;
    ubar = 1. - u;
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ (ubar * MM2 + u * q2) - 1.;
    B01 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ q2 - 1.;
    B00 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    
    gslpp::complex tpar = twoMM / Ee / ubar * I1(u, q2) + (ubar * MM2 + u * q2) / Ee / Ee / ubar / ubar * (B01 - B00);
    return MM / Mb * mySM.getQuarks(QCD::UP).getCharge() * tpar*C_2Lh_bar;
}

gslpp::complex MVll::Tparminus(double u, double q2) 
{
    double ubar = 1. - u;
    return mySM.getQuarks(QCD::DOWN).getCharge()*(8. * C_8Lh / (ubar + u * q2 / MM2)
            + sixMMoMb * H_c(ubar * MM2 + u * q2,mu_h*mu_h) * C_2Lh_bar);
}
//////////////////////////////////////////////////////////////////
double MVll::Integrand_ReTperpplus(double up) 
{
    double u = up;
    return (Tperpplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.))).real();
}

double MVll::Integrand_ImTperpplus(double up) 
{
    double u = up;
    return (Tperpplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.))).imag();
}

double MVll::Integrand_ReTparplus(double up) 
{   
    double u = up;
    return ((Tparplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/mySM.getMesons(meson).getLambdaM()).real();
}

double MVll::Integrand_ImTparplus(double up) 
{   
    double u = up;
    return ((Tparplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/mySM.getMesons(meson).getLambdaM()).imag();
}

double MVll::Integrand_ReTparminus(double up) 
{
    double Lambdaplus = mySM.getMesons(meson).getLambdaM();
    gslpp::complex Lambdamin = exp(-tmpq2 / MM / Lambdaplus) / Lambdaplus * (-gsl_sf_expint_Ei(tmpq2 / MM / Lambdaplus) + gslpp::complex::i() * M_PI);
    
    double u = up;
    return ((Tparminus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/Lambdamin).real();
}

double MVll::Integrand_ImTparminus(double up) 
{
    double Lambdaplus = mySM.getMesons(meson).getLambdaM();
    gslpp::complex Lambdamin = exp(-tmpq2 / MM / Lambdaplus) / Lambdaplus * (-gsl_sf_expint_Ei(tmpq2 / MM / Lambdaplus) + gslpp::complex::i() * M_PI);
    
    double u = up;
    return ((Tparminus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/Lambdamin).imag();
}

double MVll::Integrand_ImTpar_pm(double up){
    return Integrand_ImTparplus(up) + Integrand_ImTparminus(up);
}

double MVll::Integrand_ReTpar_pm(double up){
    return Integrand_ReTparplus(up) + Integrand_ReTparminus(up);
}

gslpp::complex MVll::F19(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    double Ls = log(s);
    double Lc = log(Mc/Mb);
    double Lm = log(mu_b/Mb);
    gslpp::complex i = gslpp::complex::i();
    return (-1424./729. + 16./243.*i*M_PI + 64./27.*Lc)*Lm - 16./243.*Lm*Ls + (16./1215. - 32./135./Mc2*Mb2)*Lm*s
            +(4./2835. - 8./315./Mc2*Mb2/Mc2*Mb2)*Lm*s2 + (16./76545. - 32./8505/Mc2*Mb2/Mc2*Mb2/Mc2*Mb2)*Lm*s*s2 - 256./243.*Lm*Lm
            +(-11.65 + 0.18223*i + (-24.709 - 0.13474*i) * s + (-43.588 - 0.4738*i) *s2 + (-86.22 - 1.3542 * i) *s*s2 
            + (-0.080959 - 0.051864*i + (-0.036585 + 0.024753*i) * s + (-0.021692 + 0.036925*i) *s2 + (0.013282 + 0.052023 * i) *s*s2)*Ls );
}

gslpp::complex MVll::F27(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    double Ls = log(s);
    gslpp::complex i = gslpp::complex::i();
    return F27_0 + F27_1 * s + F27_2 * s2 + F27_3 * s * s2 + F27_L1_1 * Ls * s + F27_L1_2 * Ls * s2 + F27_L1_3 * Ls * s * s2;
}

gslpp::complex MVll::F29(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    double Ls = log(s);
    gslpp::complex i = gslpp::complex::i();
    return F29_0 + F29_L1 * Ls + F29_1 * s +F29_2 * s2 +F29_3 * s * s2 + F29_L1_1 * Ls * s + F29_L1_2 * Ls * s2 + F29_L1_3 * Ls * s2 *s;
}

gslpp::complex MVll::F87(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    return F87_0 + F87_1 * s + F87_2 * s2 + F87_3 * s * s2 - 0.888889 * log(s)*(1. + s + s2 + s * s2);
}

double MVll::F89(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    return F89_0 + F89_1 * s + F89_2 * s2 + F89_3 * s * s2 + 1.77778 * log(s)*(1. + s + s2 + s * s2);
}

gslpp::complex MVll::Cperp(double q2) 
{
    return 1. / CF * (-C_2L_bar * F27(q2) - C_8L * F87(q2) - q2 / 2. / Mb / MM * 
            (C_2L_bar * F29(q2) + 2.*C_1L_bar*(F19(q2) + F29(q2)/6.) + C_8L * F89(q2)));
}

gslpp::complex MVll::Cpar(double q2) 
{
    return 1. / CF * (C_2L_bar * F27(q2) + C_8L * F87(q2) + MM / 2. / Mb * 
            (C_2L_bar * F29(q2) + 2.*C_1L_bar*(F19(q2) + F29(q2)/6.) + C_8L * F89(q2)));
}

gslpp::complex MVll::deltaTperp(double q2) 
{
    tmpq2 = q2;
    if (deltaTperpCached[q2] == 0) {
        
        DTPPR = convertToGslFunction(boost::bind(&MVll::Integrand_ReTperpplus, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ReTppint = avaDTPPR;
        
        DTPPR = convertToGslFunction(boost::bind(&MVll::Integrand_ImTperpplus, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ImTppint = avaDTPPR;

        cacheDeltaTperp[q2] = ReTppint + gslpp::complex::i() * ImTppint;
        deltaTperpCached[q2] = 1;
    }

    return deltaT_0 * Cperp(q2) + deltaT_1perp / T_1(q2) / mySM.getMesons(meson).getLambdaM() * cacheDeltaTperp[q2];
}

gslpp::complex MVll::deltaTpar(double q2) 
{
    double Lambdaplus = mySM.getMesons(meson).getLambdaM();
    gslpp::complex Lambdamin = exp(-q2 / MM / Lambdaplus) / Lambdaplus * (-gsl_sf_expint_Ei(q2 / MM / Lambdaplus) + gslpp::complex::i() * M_PI);
    double T3q2 = MM2mMV2/lambda(q2) * ( (MM2 + 3.*MV2 - q2) * T_2(q2) - 8.*MM*MV2/MMpMV * FF_fit(q2, a_0T23, a_1T23, a_2T23, MRT23_2) );
    tmpq2 = q2;
    Ee = (MM2 - tmpq2) / twoMM;
    
    if (deltaTparpCached[q2] == 0) {
        
        DTPPR = convertToGslFunction(boost::bind(&MVll::Integrand_ReTpar_pm, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ReTppint = avaDTPPR;

        DTPPR = convertToGslFunction(boost::bind(&MVll::Integrand_ImTpar_pm, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ImTppint = avaDTPPR;        

        cacheDeltaTparp[q2] = (ReTppint + gslpp::complex::i() * ImTppint);
        deltaTparpCached[q2] = 1;
    }

    return deltaT_0 * Cpar(q2) + deltaT_1par * MV/Ee / (T_1(q2) - T3q2) * (cacheDeltaTparp[q2]);
}



gslpp::complex MVll::DeltaC9_p(double q2)
{
    return 1./q2 * Mb/MM * (MM2mMV2 * (MM2 - q2)/MM2 -
            sqrt(lambda(q2))) * deltaTperp(q2);
}

gslpp::complex MVll::DeltaC9_m(double q2)
{
    return 1./q2 * Mb/MM * (MM2mMV2 * (MM2 - q2)/MM2 +
            sqrt(lambda(q2))) * deltaTperp(q2);
}


gslpp::complex MVll::DeltaC9_0(double q2)
{
    return 1. / 2. / MV / MM / sqrt(q2) * ((MM2mMV2 * (MM2mMV2 - q2) - lambda(q2))* (MM2 - q2) * 
            Mb/MM2/q2 * deltaTperp(q2) - lambda(q2) * (deltaTpar(q2) + deltaTperp(q2))* Mb/MM2mMV2);
}

/*******************************************************************************
 * Fitting routines                                                         *
 * ****************************************************************************/

double MVll::reDC9fit(double* x, double* p)
{
    return p[0]/x[0] + p[1] + p[2]*x[0] + p[3]*x[0]*x[0] + p[4]*x[0]*x[0]*x[0] + p[5]*x[0]*x[0]*x[0]*x[0] + p[6]*x[0]*x[0]*x[0]*x[0]*x[0]; 
}

double MVll::imDC9fit(double* x, double* p)
{
    return p[0]/x[0] + p[1] + p[2]*x[0] + p[3]*x[0]*x[0] + p[4]*x[0]*x[0]*x[0] + p[5]*x[0]*x[0]*x[0]*x[0] + p[6]*x[0]*x[0]*x[0]*x[0]*x[0] + p[7]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
    
    //double thr = 4.*Mc2;
}

void MVll::fit_DeltaC9_p_mumu()
{
    int dim = 0;
    for (double i=0.1; i<SWITCH; i+=0.4) {
        double q2tmp = i;        
        myq2.push_back(q2tmp);
        ReDeltaC9_p_mumu.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_p_mumu.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    for (double i=SWITCH; i<8.2; i+=0.4) {
        double q2tmp = i;        
        myq2.push_back(q2tmp);
        ReDeltaC9_p_mumu.push_back(q2tmp * (1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_p_mumu.push_back(q2tmp * (1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    gr1 =TGraph(dim, myq2.data(), ReDeltaC9_p_mumu.data());
    gr2 =TGraph(dim, myq2.data(), ImDeltaC9_p_mumu.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0.1,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0.1,8.1,8,"MVll","imDC9fit");
    
    refres_p_mumu = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_p_mumu = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_p_mumu.clear();
    ImDeltaC9_p_mumu.clear();
    myq2.clear();
}

void MVll::fit_DeltaC9_p_ee()
{
    int dim = 0;
    for (double i=0.002; i<1.12; i+=0.05) {
        double q2tmp = i;        
        myq2.push_back(q2tmp);
        ReDeltaC9_p_ee.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_p_ee.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 -
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    
    gr1 =TGraph(dim, myq2.data(), ReDeltaC9_p_ee.data());
    gr2 =TGraph(dim, myq2.data(), ImDeltaC9_p_ee.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0,8.1,8,"MVll","imDC9fit");
    
    refres_p_ee = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_p_ee = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_p_ee.clear();
    ImDeltaC9_p_ee.clear();
    myq2.clear();
}

void MVll::fit_DeltaC9_m_mumu()
{
    int dim = 0;
    for (double i=0.1; i<SWITCH; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_m_mumu.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_m_mumu.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    for (double i=SWITCH; i<8.2; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_m_mumu.push_back(q2tmp * (1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_m_mumu.push_back(q2tmp * (1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    
    gr1 = TGraph(dim, myq2.data(), ReDeltaC9_m_mumu.data());
    gr2 = TGraph(dim, myq2.data(), ImDeltaC9_m_mumu.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0,8.1,8,"MVll","imDC9fit");
    
    refres_m_mumu = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_m_mumu = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_m_mumu.clear();
    ImDeltaC9_m_mumu.clear();
    myq2.clear();
}

void MVll::fit_DeltaC9_m_ee()
{
    int dim = 0;
    for (double i=0.002; i<1.12; i+=0.05) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_m_ee.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).real());
        ImDeltaC9_m_ee.push_back((1./q2tmp * Mb/MM * (MM2mMV2 * (MM2 - q2tmp)/MM2 +
                                              sqrt(lambda(q2tmp))) * deltaTperp(q2tmp)).imag());
        dim++;
    }
    
    gr1 = TGraph(dim, myq2.data(), ReDeltaC9_m_ee.data());
    gr2 = TGraph(dim, myq2.data(), ImDeltaC9_m_ee.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0,8.1,8,"MVll","imDC9fit");
    
    refres_m_ee = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_m_ee = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_m_ee.clear();
    ImDeltaC9_m_ee.clear();
    myq2.clear();
}

void MVll::fit_DeltaC9_0_mumu()
{
    int dim = 0;
    for (double i=0.1; i<SWITCH; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_0_mumu.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).real());
        ImDeltaC9_0_mumu.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).imag());
        dim++;
    }
    for (double i=SWITCH; i<8.2; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_0_mumu.push_back(q2tmp * (1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).real());
        ImDeltaC9_0_mumu.push_back(q2tmp * (1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).imag());
        dim++;
    }
    
    gr1 = TGraph(dim, myq2.data(), ReDeltaC9_0_mumu.data());
    gr2 = TGraph(dim, myq2.data(), ImDeltaC9_0_mumu.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0,8.1,8,"MVll","imDC9fit");
    
    refres_0_mumu = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_0_mumu = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_0_mumu.clear();
    ImDeltaC9_0_mumu.clear();
    myq2.clear();
}

void MVll::fit_DeltaC9_0_ee()
{
    int dim = 0;
    for (double i=0.002; i<0.05; i+=0.005) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_0_ee.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).real());
        ImDeltaC9_0_ee.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).imag());
        dim++;
    }
    for (double i=0.05; i<1.12; i+=0.05) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9_0_ee.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).real());
        ImDeltaC9_0_ee.push_back((1. / 2. / MV / MM / sqrt(q2tmp) * ((MM2mMV2 * (MM2mMV2 - q2tmp) - lambda(q2tmp))* (MM2 - q2tmp) *
                                                             Mb/MM2/q2tmp * deltaTperp(q2tmp) - lambda(q2tmp) * (deltaTpar(q2tmp) + deltaTperp(q2tmp))* Mb/MM2mMV2)).imag());
        dim++;
    }
    
    gr1 = TGraph(dim, myq2.data(), ReDeltaC9_0_ee.data());
    gr2 = TGraph(dim, myq2.data(), ImDeltaC9_0_ee.data());
    
    reffit = TF1("reffit",this,&MVll::reDC9fit,0,8.1,7,"MVll","reDC9fit");
    imffit = TF1("imffit",this,&MVll::imDC9fit,0,8.1,8,"MVll","imDC9fit");
    
    refres_0_ee = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres_0_ee = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9_0_ee.clear();
    ImDeltaC9_0_ee.clear();
    myq2.clear();
}

gslpp::complex MVll::fDeltaC9_p(double q2)
{
    switch (lep) {
        
        case StandardModel::MU:
            if (q2 < SWITCH) return (reDC9fit(&q2, const_cast<double *>(refres_p_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_p_mumu->GetParams())));
            else return (reDC9fit(&q2, const_cast<double *>(refres_p_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_p_mumu->GetParams())))/q2;
            break;
        case StandardModel::ELECTRON:
            return (reDC9fit(&q2, const_cast<double *>(refres_p_ee->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_p_ee->GetParams())));;
            break;
        default:
            std::stringstream out;
            out << lep;
            throw std::runtime_error("MVll::fDeltaC9_p : " + out.str() + " not implemented");
    }
}

gslpp::complex MVll::fDeltaC9_m(double q2)
{
    switch (lep) {
        
        case StandardModel::MU:
            if (q2 < SWITCH) return (reDC9fit(&q2, const_cast<double *>(refres_m_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_m_mumu->GetParams())));
            else return (reDC9fit(&q2, const_cast<double *>(refres_m_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_m_mumu->GetParams())))/q2;
            break;
        case StandardModel::ELECTRON:
            return (reDC9fit(&q2, const_cast<double *>(refres_m_ee->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_m_ee->GetParams())));
            break;
        default:
            std::stringstream out;
            out << lep;
            throw std::runtime_error("MVll::fDeltaC9_m : " + out.str() + " not implemented");
    }
}


gslpp::complex MVll::fDeltaC9_0(double q2)
{
    switch (lep) {
        
        case StandardModel::MU:
            if (q2 < SWITCH) return (reDC9fit(&q2, const_cast<double *>(refres_0_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_0_mumu->GetParams())));
            else return (reDC9fit(&q2, const_cast<double *>(refres_0_mumu->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_0_mumu->GetParams())))/q2;
            break;
        case StandardModel::ELECTRON:
            return (reDC9fit(&q2, const_cast<double *>(refres_0_ee->GetParams())) + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres_0_ee->GetParams())));
            break;
        default:
            std::stringstream out;
            out << lep;
            throw std::runtime_error("MVll::fDeltaC9_0 : " + out.str() + " not implemented");
    }
}

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MVll::H_c(double q2, double mu2) 
{
    double x = fourMc2 / q2;
    gslpp::complex par;

    if (x > 1.) par = sqrt(x - 1.) * atan(1. / sqrt(x - 1.));
    else par = sqrt(1. - x) * (log((1. + sqrt(1. - x)) / sqrt(x)) - ihalfMPI);

    return -fournineth * (log(Mc2 / mu2) - twothird - x) - fournineth * (2. + x) * par;
}

gslpp::complex MVll::H_b(double q2) 
{
    double x = fourMb2 / q2;
    gslpp::complex par;

    if (x > 1.) par = sqrt(x - 1.) * atan(1. / sqrt(x - 1.));
    else par = sqrt(1. - x) * (log((1. + sqrt(1. - x)) / sqrt(x)) - ihalfMPI);

    return -fournineth * (logMb - twothird - x) - fournineth * (2. + x) * par;
}

gslpp::complex MVll::H_0(double q2) 
{
    return (H_0_pre - fournineth * log(q2 / mu_b2));
}

gslpp::complex MVll::Y(double q2) 
{
    return -half * H_0(q2) * H_0_WC + H_c(q2,mu_b2) * H_c_WC - half * H_b(q2) * H_b_WC;
}

gslpp::complex MVll::H_V_0(double q2) 
{
    return -(((C_9 + fDeltaC9_0(q2) + Y(q2)) - C_9p) * V_0t(q2) + MM2 / q2 * (twoMboMM * (C_7 - C_7p) * T_0t(q2) - sixteenM_PI2 * (h_0[0] + h_1[0] * q2 + h_2[0] * q2 * q2)));
}

gslpp::complex MVll::H_V_p(double q2) 
{
    return -(((C_9 + fDeltaC9_p(q2) + Y(q2)) * V_p(q2) - C_9p * V_m(q2)) + MM2 / q2 * (twoMboMM * (C_7 * T_p(q2) - C_7p * T_m(q2)) - sixteenM_PI2 * (h_0[1] + h_1[1] * q2 + h_2[1] * q2 * q2)));
}

gslpp::complex MVll::H_V_m(double q2) 
{
    return -(((C_9 + fDeltaC9_m(q2) + Y(q2)) * V_m(q2) - C_9p * V_p(q2)) + MM2 / q2 * (twoMboMM * (C_7 * T_m(q2) - C_7p * T_p(q2)) - sixteenM_PI2 * (h_0[2] + h_1[2] * q2 + h_2[2] * q2 * q2)));
}

gslpp::complex MVll::H_A_0(double q2) 
{
    return ( -C_10 + C_10p) * V_0t(q2);
}

gslpp::complex MVll::H_A_p(double q2) 
{
    return ( -C_10 * V_p(q2) + C_10p * V_m(q2));
}

gslpp::complex MVll::H_A_m(double q2) 
{
    return ( -C_10 * V_m(q2) + C_10p * V_p(q2));
}

gslpp::complex MVll::H_S(double q2) 
{
    return MboMW * (C_S - C_Sp) * S_L(q2);
}

gslpp::complex MVll::H_P(double q2) 
{
    return ( MboMW * (C_P - C_Pp) + twoMlepMb / q2 * (C_10 - C_10p) * (1. + MsoMb)) * S_L(q2);
}

/*******************************************************************************
 * Angular coefficients                                                         *
 * ****************************************************************************/
double MVll::k2(double q2) 
{
    return (MM4 + q2 * q2 + MV4 - twoMV2 * q2 - twoMM2 * (q2 + MV2)) / fourMM2;
}

double MVll::beta(double q2) 
{
    return sqrt(1. - 4. * Mlep2 / q2);
}

double MVll::beta2(double q2) 
{
    return 1. - 4. * Mlep2 / q2;
}

double MVll::lambda(double q2) 
{
    return (MM4 + q2 * q2 + MV4 - twoMV2 * q2 - twoMM2 * (q2 + MV2));
}

double MVll::F(double q2, double b_i) 
{
    return sqrt(lambda(q2)) * beta(q2) * q2 * b_i / (ninetysixM_PI3MM3);
}

double MVll::I_1c(double q2) 
{
    return F(q2, b)*((H_V_0(q2).abs2() + H_A_0(q2).abs2()) / 2. + H_P(q2).abs2() + 2. * Mlep2 / q2 * (H_V_0(q2).abs2()
            - H_A_0(q2).abs2()) + beta2(q2) * H_S(q2).abs2());
}

double MVll::I_1s(double q2) 
{
    return F(q2, b)*((beta2(q2) + 2.) / 8. * (H_V_p(q2).abs2() + H_V_m(q2).abs2() + H_A_p(q2).abs2() + H_A_m(q2).abs2()) +
            Mlep2 / q2 * (H_V_p(q2).abs2() + H_V_m(q2).abs2() - H_A_p(q2).abs2() - H_A_m(q2).abs2()));
}

double MVll::I_2c(double q2) 
{
    return -F(q2, b) * beta2(q2) / 2. * (H_V_0(q2).abs2() + H_A_0(q2).abs2());
}

double MVll::I_2s(double q2) 
{
    return F(q2, b) * beta2(q2) / 8. * (H_V_p(q2).abs2() + H_V_m(q2).abs2() + H_A_p(q2).abs2() + H_A_m(q2).abs2());
}

double MVll::I_3(double q2) 
{
    return -F(q2, b) / 2. * ((H_V_p(q2) * H_V_m(q2).conjugate()).real() + (H_A_p(q2) * H_A_m(q2).conjugate()).real());
}

double MVll::I_4(double q2) 
{
    return F(q2, b) * beta2(q2) / 4. * (((H_V_m(q2) + H_V_p(q2)) * H_V_0(q2).conjugate()).real() + ((H_A_m(q2) + H_A_p(q2)) * H_A_0(q2).conjugate()).real());
}

double MVll::I_5(double q2) 
{
    return F(q2, b)*(beta(q2) / 2. * (((H_V_m(q2) - H_V_p(q2)) * H_A_0(q2).conjugate()).real() + ((H_A_m(q2) - H_A_p(q2)) * H_V_0(q2).conjugate()).real()) -
            beta(q2) * Mlep / sqrt(q2)*(H_S(q2).conjugate()*(H_V_p(q2) + H_V_m(q2))).real());
}

double MVll::I_6s(double q2) 
{
    return F(q2, b) * beta(q2)*(H_V_m(q2)*(H_A_m(q2).conjugate()) - H_V_p(q2)*(H_A_p(q2).conjugate())).real();
}

double MVll::I_6c(double q2) 
{
    return 2. * F(q2, b) * beta(q2) * Mlep / sqrt(q2)*(H_S(q2).conjugate() * H_V_0(q2)).real();
}

double MVll::I_7(double q2) 
{
    return F(q2, b)*(beta(q2) / 2. * (((H_V_m(q2) + H_V_p(q2)) * H_A_0(q2).conjugate()).imag() + ((H_A_m(q2) + H_A_p(q2)) * H_V_0(q2).conjugate()).imag()) -
            beta(q2) * Mlep / sqrt(q2)*(H_S(q2).conjugate()*(H_V_m(q2) - H_V_p(q2))).imag());
}

double MVll::I_8(double q2) 
{
    return F(q2, b) * beta2(q2) / 4. * (((H_V_m(q2) - H_V_p(q2)) * H_V_0(q2).conjugate()).imag() + ((H_A_m(q2) - H_A_p(q2)) * H_A_0(q2).conjugate()).imag());
}

double MVll::I_9(double q2) 
{
    return F(q2, b) * beta2(q2) / 2. * ((H_V_p(q2) * H_V_m(q2).conjugate()).imag() + (H_A_p(q2) * H_A_m(q2).conjugate()).imag());
}

double MVll::Delta(int i, double q2) 
{
    return 0; /* FIX CPV */
    //return (I(i, q2,0) - I(i, q2,1))/2;
}


double MVll::integrateSigma(int i, double q_min, double q_max) 
{
    updateParameters();

    std::pair<double, double > qbin = std::make_pair(q_min, q_max);

    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 0:
            if (sigma0Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma1c, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma0[qbin] = NN*avaSigma;
                sigma0Cached[qbin] = 1;
            }
            return cacheSigma0[qbin];
            break;
        case 1:
            if (sigma1Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma1s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma1[qbin] = NN*avaSigma;
                sigma1Cached[qbin] = 1;
            }
            return cacheSigma1[qbin];
            break;
        case 2:
            if (sigma2Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma2c, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma2[qbin] = NN*avaSigma;
                sigma2Cached[qbin] = 1;
            }
            return cacheSigma2[qbin];
            break;
        case 3:
            if (sigma3Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma2s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma3[qbin] = NN*avaSigma;
                sigma3Cached[qbin] = 1;
            }
            return cacheSigma3[qbin];
            break;
        case 4:
            if (sigma4Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma3, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma4[qbin] = NN*avaSigma;
                sigma4Cached[qbin] = 1;
            }
            return cacheSigma4[qbin];
            break;
        case 5:
            if (sigma5Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma4, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma5[qbin] = NN*avaSigma;
                sigma5Cached[qbin] = 1;
            }
            return cacheSigma5[qbin];
            break;
        case 6:
            if (sigma6Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma5, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma6[qbin] = NN*avaSigma;
                sigma6Cached[qbin] = 1;
            }
            return cacheSigma6[qbin];
            break;
        case 7:
            if (sigma7Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma6s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma7[qbin] = NN*avaSigma;
                sigma7Cached[qbin] = 1;
            }
            return cacheSigma7[qbin];
            break;
        case 9:
            if (sigma9Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma7, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma9[qbin] = NN*avaSigma;
                sigma9Cached[qbin] = 1;
            }
            return cacheSigma9[qbin];
            break;
        case 10:
            if (sigma10Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma8, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma10[qbin] = NN*avaSigma;
                sigma10Cached[qbin] = 1;
            }
            return cacheSigma10[qbin];
            break;
        case 11:
            if (sigma11Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma9, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma11[qbin] = NN*avaSigma;
                sigma11Cached[qbin] = 1;
            }
            return cacheSigma11[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);

}

double MVll::getSigma(int i, double q_2) 
{
    updateParameters();

    switch (i) {
        case 0:
            return getSigma1c(q_2);
            break;
        case 1:
            return getSigma1s(q_2);
            break;
        case 2:
            return getSigma2c(q_2);
            break;
        case 3:
            return getSigma2s(q_2);
            break;
        case 4:
            return getSigma3(q_2);
            break;
        case 5:
            return getSigma4(q_2);
            break;
        case 6:
            return getSigma5(q_2);
            break;
        case 7:
            return getSigma6s(q_2);
            break;
        case 9:
            return getSigma7(q_2);
            break;
        case 10:
            return getSigma8(q_2);
            break;
        case 11:
            return getSigma9(q_2);
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateSigma: index " + out.str() + " not implemented");
    }
}

double MVll::integrateDelta(int i, double q_min, double q_max) 
{
    updateParameters();

    std::pair<double, double > qbin = std::make_pair(q_min, q_max);

    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 0:
            if (delta0Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta0, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta0[qbin] = avaDelta;
                delta0Cached[qbin] = 1;
            }
            return cacheDelta0[qbin];
            break;
        case 1:
            if (delta1Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta1, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta1[qbin] = avaDelta;
                delta1Cached[qbin] = 1;
            }
            return cacheDelta1[qbin];
            break;
        case 2:
            if (delta2Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta2, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta2[qbin] = avaDelta;
                delta2Cached[qbin] = 1;
            }
            return cacheDelta2[qbin];
            break;
        case 3:
            if (delta3Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta3, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta3[qbin] = avaDelta;
                delta3Cached[qbin] = 1;
            }
            return cacheDelta3[qbin];
            break;
        case 7:
            if (delta7Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta7, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta7[qbin] = avaDelta;
                delta7Cached[qbin] = 1;
            }
            return cacheDelta7[qbin];
            break;
        case 11:
            if (delta11Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta11, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta11[qbin] = avaDelta;
                delta11Cached[qbin] = 1;
            }
            return cacheDelta11[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("integrateDelta: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);

}