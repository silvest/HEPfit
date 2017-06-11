/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "StandardModel.h"
#include "MPll.h"
#include "gslpp_complex.h"
#include "std_make_vector.h"
#include <gsl/gsl_sf.h>
#include <boost/bind.hpp>
#include <limits>
#include <TFitResult.h>





MPll::MPll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i) 
:       mySM(SM_i),
        fplus_cache(2, 0.),
        fT_cache(2, 0.),
        fplus_lat_cache(3, 0.),
        fT_lat_cache(3, 0.),
        f0_lat_cache(3, 0.),
        k2_cache(2, 0.),
        SL_cache(2, 0.),
        N_cache(3, 0.),
        Ycache(2, 0.),
        H_V0cache(2, 0.),
        H_Scache(2, 0.),
        H_P_cache(4, 0.),
        T_cache(5, 0.) 
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalar = pseudoscalar_i;

#if NFPOLARBASIS_MPLL
    if (pseudoscalar == StandardModel::K_P) mpllParameters = make_vector<std::string>()
        << "r_1_fplus" << "r_2_fplus" << "m_fit2_fplus" << "r_1_fT" << "r_2_fT" << "m_fit2_fT" << "r_2_f0" << "m_fit2_f0"
        << "b_0_fplus" << "b_1_fplus" << "b_2_fplus" << "m_fit_fplus_lat"
        << "b_0_fT" << "b_1_fT" << "b_2_fT" << "m_fit_fT_lat"
        << "b_0_f0" << "b_1_f0" << "b_2_f0" << "m_fit_f0_lat"
        << "absh_0_MP" << "argh_0_MP" << "absh_1_MP" << "argh_1_MP";
#else
    if (pseudoscalar == StandardModel::K_P) mpllParameters = make_vector<std::string>()
        << "r_1_fplus" << "r_2_fplus" << "m_fit2_fplus" << "r_1_fT" << "r_2_fT" << "m_fit2_fT" << "r_2_f0" << "m_fit2_f0"
        << "b_0_fplus" << "b_1_fplus" << "b_2_fplus" << "m_fit_fplus_lat"
        << "b_0_fT" << "b_1_fT" << "b_2_fT" << "m_fit_fT_lat"
        << "b_0_f0" << "b_1_f0" << "b_2_f0" << "m_fit_f0_lat"
        << "reh_0_MP" << "imh_0_MP" << "reh_1_MP" << "imh_1_MP";
#endif
    else {
        std::stringstream out;
        out << pseudoscalar;
        throw std::runtime_error("MPll: pseudoscalar " + out.str() + " not implemented");
    }
    
    I0_updated = 0;
    I2_updated = 0;
    I8_updated = 0;

    VL_updated = 0;
    TL_updated = 0;
    SL_updated = 0;
    
    deltaTparpupdated = 0;
    deltaTparmupdated = 0;
    
    w_sigma = gsl_integration_cquad_workspace_alloc (100);
    w_delta = gsl_integration_cquad_workspace_alloc (100);   
    w_DTPPR = gsl_integration_cquad_workspace_alloc (100);
}


MPll::~MPll() 
{}

std::vector<std::string> MPll::initializeMPllParameters()
{
    mySM.initializeMeson(meson);
    mySM.initializeMeson(pseudoscalar);
    return mpllParameters;
}

void MPll::updateParameters()
{
    if (!mySM.getFlavour().getUpdateFlag(meson, pseudoscalar, lep)) return;

    

    GF = mySM.getGF();
    ale=mySM.getAle();
    Mlep=mySM.getLeptons(lep).getMass();
    MM=mySM.getMesons(meson).getMass();
    MP=mySM.getMesons(pseudoscalar).getMass();
    Mb=mySM.getQuarks(QCD::BOTTOM).getMass();    // add the PS b mass
    Mc=mySM.getQuarks(QCD::CHARM).getMass();
    Ms=mySM.getQuarks(QCD::STRANGE).getMass();
    MW=mySM.Mw();
    lambda_t=mySM.computelamt_s();
    mu_b = mySM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    width = mySM.getMesons(meson).computeWidth();
    
    switch(pseudoscalar){
        case StandardModel::K_P :
            r_1_fplus = mySM.getOptionalParameter("r_1_fplus");
            r_2_fplus = mySM.getOptionalParameter("r_2_fplus");
            m_fit2_fplus = mySM.getOptionalParameter("m_fit2_fplus");
            r_1_fT = mySM.getOptionalParameter("r_1_fT");
            r_2_fT = mySM.getOptionalParameter("r_2_fT");
            m_fit2_fT = mySM.getOptionalParameter("m_fit2_fT");
            r_2_f0 = mySM.getOptionalParameter("r_2_f0");
            m_fit2_f0 = mySM.getOptionalParameter("m_fit2_f0");
            
            b_0_fplus = mySM.getOptionalParameter("b_0_fplus");
            b_1_fplus = mySM.getOptionalParameter("b_1_fplus");
            b_2_fplus = mySM.getOptionalParameter("b_2_fplus");
            m_fit2_fplus_lat = mySM.getOptionalParameter("m_fit_fplus_lat") * mySM.getOptionalParameter("m_fit_fplus_lat");
            b_0_fT = mySM.getOptionalParameter("b_0_fT");
            b_1_fT = mySM.getOptionalParameter("b_1_fT");
            b_2_fT = mySM.getOptionalParameter("b_2_fT");
            m_fit2_fT_lat = mySM.getOptionalParameter("m_fit_fT_lat") * mySM.getOptionalParameter("m_fit_fT_lat");
            b_0_f0 = mySM.getOptionalParameter("b_0_f0");
            b_1_f0 = mySM.getOptionalParameter("b_1_f0");
            b_2_f0 = mySM.getOptionalParameter("b_2_f0");
            m_fit2_f0_lat = mySM.getOptionalParameter("m_fit_f0_lat") * mySM.getOptionalParameter("m_fit_f0_lat");
            
            spectator_charge = mySM.getQuarks(QCD::UP).getCharge();
            
            etaP = -1;
            angmomP = 0.;
    
            break;
        default:
            std::stringstream out;
            out << pseudoscalar;
            throw std::runtime_error("MPll: pseudoscalar " + out.str() + " not implemented");
    }
    
#if NFPOLARBASIS_MPLL
        h_0 = gslpp::complex(mySM.getOptionalParameter("absh_0_MP"), mySM.getOptionalParameter("argh_0_MP"), true);
        h_1 = gslpp::complex(mySM.getOptionalParameter("absh_1_MP"), mySM.getOptionalParameter("argh_1_MP"), true);
#else
        h_0 = gslpp::complex(mySM.getOptionalParameter("reh_0_MP"), mySM.getOptionalParameter("imh_0_MP"), false);
        h_1 = gslpp::complex(mySM.getOptionalParameter("reh_1_MP"), mySM.getOptionalParameter("imh_1_MP"), false);
#endif
    
    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu_b, lep);   //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getFlavour().ComputeCoeffprimeBMll(mu_b, lep);   //check the mass scale, scheme fixed to NDR

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_1L_bar = (*(allcoeff[LO]))(0)/2.;
    C_2 =  ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
    C_2L_bar = (*(allcoeff[LO]))(1) - (*(allcoeff[LO]))(0)/6.;
    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 =  ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_8L = (*(allcoeff[LO]))(7);
    C_9 = (*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8);
    C_10 = (*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9);
    C_S = (*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10);
    C_P = (*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11);
    
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    C_9p = (*(allcoeffprime[LO]))(8) + (*(allcoeffprime[NLO]))(8);
    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    allcoeffh = mySM.getFlavour().ComputeCoeffBMll(mu_h, lep); //check the mass scale, scheme fixed to NDR

    C_1Lh_bar = (*(allcoeffh[LO]))(0)/2.;
    C_2Lh_bar = (*(allcoeffh[LO]))(1) - (*(allcoeff[LO]))(0)/6.;
    C_8Lh = (*(allcoeffh[LO]))(7);
    
    checkCache();
    
    H_0_pre = 8. / 27. + 4. / 9. * gslpp::complex::i() * M_PI;
    H_0_WC = (C_3 + 4. / 3. * C_4 + 16. * C_5 + 64. / 3. * C_6);
    H_c_WC = (4. / 3. * C_1 + C_2 + 6. * C_3 + 60. * C_5);
    H_b_WC = (7. * C_3 + 4. / 3. * C_4 + 76. * C_5 + 64. / 3. * C_6);
    fournineth = 4. / 9.;
    half = 1. / 2.;
    twothird = 2. / 3.;
    ihalfMPI = gslpp::complex::i() * M_PI / 2.;
    Mc2 = Mc*Mc;
    Mb2 = Mb*Mb;
    logMc = log(Mc2 / mu_b2);
    logMb = log(Mb2 / mu_b2);
    mu_b2 = mu_b*mu_b;
    fourMc2 = 4. * Mc2;
    fourMb2 = 4. * Mb2;
    Mlep2 = Mlep*Mlep;
    NN = ((4. * GF * MM * ale * lambda_t) / (sqrt(2.)*4. * M_PI)).abs2();
    MM2 = MM*MM;
    MM4 = MM2*MM2;
    MP2 = MP*MP;
    MP4 = MP2*MP2;
    MM2mMP2 = MM2 - MP2;
    twoMP2 = 2. * MP2;
    twoMM = 2. * MM;
    twoMM2 = 2. * MM2;
    twoMM2_MMpMP = twoMM2 * (MM+MP);
    twoMM_MbpMs = twoMM * (Mb+Ms);
    S_L_pre = -(MM2mMP2 / twoMM_MbpMs) * ( 1 + MsoMb )/( 1 - MsoMb );
    fourMM2 = 4. * MM2;
    twoMboMM = 2 * Mb / MM;
    sixteenM_PI2 = 16. * M_PI*M_PI;
    ninetysixM_PI3MM3 = 96. * M_PI * M_PI * M_PI * MM * MM*MM;
    MboMW = Mb / MW;
    MboMM = Mb / MM;
    MsoMb = Ms / Mb;
    twoMlepMb = 2. * Mlep*Mb;
    DC9pre = - Mb / MP / twoMM  / (MM2 - MP2);
    threeGegen0 = mySM.getMesons(pseudoscalar).getGegenalpha(0)*3.;
    threeGegen1otwo = mySM.getMesons(pseudoscalar).getGegenalpha(1)*3./2.;
    M_PI2osix = M_PI * M_PI / 6.;
    twoMc2 = 2.*Mc2;
    sixMMoMb = 6. * MM / Mb;
    CF =4./3.;
    deltaT_0 = mySM.Als(mu_b) * MboMM / 4. / M_PI;
    deltaT_1par = mySM.Als(mu_h) * CF / 4. * M_PI / 3. * mySM.getMesons(meson).getDecayconst() *
            mySM.getMesons(pseudoscalar).getDecayconst() / MM * MboMM;
            
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
    
    fit_DeltaC9_mumu();
    
    std::map<std::pair<double, double>, unsigned int >::iterator it;
    
    if (I0_updated == 0) for (it = sigma0Cached.begin(); it != sigma0Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = sigma2Cached.begin(); it != sigma2Cached.end(); ++it) it->second = 0;
    
    if (I0_updated == 0) for (it = delta0Cached.begin(); it != delta0Cached.end(); ++it) it->second = 0;
    if (I2_updated == 0) for (it = delta2Cached.begin(); it != delta2Cached.end(); ++it) it->second = 0;
    
    std::map<double, unsigned int >::iterator iti;
    if (deltaTparpupdated == 0) for (iti = deltaTparpCached.begin(); iti != deltaTparpCached.end(); ++iti) iti->second = 0;
    if (deltaTparmupdated == 0) for (iti = deltaTparmCached.begin(); iti != deltaTparmCached.end(); ++iti) iti->second = 0;
    
    if (deltaTparpupdated*deltaTparmupdated == 0) for (it = I1Cached.begin(); it != I1Cached.end(); ++it) it->second = 0;

    mySM.getFlavour().setUpdateFlag(meson, pseudoscalar, lep, false);
    
    return;
    
}

void MPll::checkCache()
{
    
    if(r_1_fplus == fplus_cache(0) && r_2_fplus == fplus_cache(1)) {
        fplus_updated = 1;
    } else {
        fplus_updated = 0;
        fplus_cache(0) = r_1_fplus;
        fplus_cache(1) = r_2_fplus;
    }
    
    if(r_1_fT == fT_cache(0) && r_2_fT == fT_cache(1)) {
        fT_updated = 1;
    } else {
        fT_updated = 0;
        fT_cache(0) = r_1_fT;
        fT_cache(1) = r_2_fT;
    }
    
    if(r_2_f0 == f0_cache) {
        f0_updated = 1;
    } else {
        f0_updated = 0;
        f0_cache = r_2_f0;
    }
    
    if (MM == k2_cache(0) && MP == k2_cache(1) ) {
        k2_updated = 1;
    } else {
        k2_updated = 0;
        k2_cache(0) = MM;
        k2_cache(1) = MP;
    }
    
    if(b_0_fplus == fplus_lat_cache(0) && b_1_fplus == fplus_lat_cache(1) && b_2_fplus == fplus_lat_cache(2)) {
        fplus_lat_updated = 1;
    } else {
        fplus_lat_updated = 0;
        fplus_lat_cache(0) = b_0_fplus;
        fplus_lat_cache(1) = b_1_fplus;
        fplus_lat_cache(2) = b_2_fplus;
    }
    
    if(b_0_fT == fT_lat_cache(0) && b_1_fT == fT_lat_cache(1) && b_2_fT == fT_lat_cache(2)) {
        fT_lat_updated = 1;
    } else {
        fT_lat_updated = 0;
        fT_lat_cache(0) = b_0_fT;
        fT_lat_cache(1) = b_1_fT;
        fT_lat_cache(2) = b_2_fT;
    }
    
    if(b_0_f0 == f0_lat_cache(0) && b_1_f0 == f0_lat_cache(1) && b_2_f0 == f0_lat_cache(2)) {
        f0_lat_updated = 1;
    } else {
        f0_lat_updated = 0;
        f0_lat_cache(0) = b_0_f0;
        f0_lat_cache(1) = b_1_f0;
        f0_lat_cache(2) = b_2_f0;
    }
    
    if (Mlep == beta_cache) {
        beta_updated = 1;
    } else {
        beta_updated = 0;
        beta_cache = Mlep;
    }
    
    lambda_updated = k2_updated;
    F_updated = lambda_updated * beta_updated;
    
    if(LATTICE) 
        VL_updated = k2_updated * fplus_lat_updated;
    else
        VL_updated = k2_updated * fplus_updated;
    
    if(LATTICE) 
        TL_updated = k2_updated * fT_lat_updated;
    else
        TL_updated = k2_updated * fT_updated;
    
    if (Mb == SL_cache(0) && Ms == SL_cache(1) ){
        if(LATTICE) 
            SL_updated = k2_updated * f0_lat_updated;
        else
            SL_updated = k2_updated * f0_updated;
    } else {
        SL_updated = 0;
        SL_cache(0) = Mb;
        SL_cache(1) = Ms;
    }
    
    if (GF == N_cache(0) && ale == N_cache(1) && MM == N_cache(2) && lambda_t == Nc_cache) {
        N_updated = 1;
    } else {
        N_updated = 0;
        N_cache(0) = GF;
        N_cache(1) = ale;
        N_cache(2) = MM;
        Nc_cache = lambda_t;
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
    
    if (Mb == Ycache(0) && Mc == Ycache(1) ) {
        Yupdated = C_1_updated * C_2_updated * C_3_updated * C_4_updated * C_5_updated * C_6_updated;
    } else {
        Yupdated = 0;
        Ycache(0) = Mb;
        Ycache(1) = Mc;
    }
    
    if (MM == H_V0cache(0) && Mb == H_V0cache(1) && h_0 == H_V0Ccache[0] && h_1 == H_V0Ccache[1]) {
        H_V0updated = N_updated * C_9_updated * Yupdated * VL_updated * C_9p_updated  * C_7_updated * TL_updated * C_7p_updated;
    } else {
        H_V0updated = 0;
        H_V0cache(0) = MM;
        H_V0cache(1) = Mb;
        H_V0Ccache[0] = h_0;
        H_V0Ccache[1] = h_1;
    }
    
    H_A0updated = N_updated * C_10_updated * VL_updated * C_10p_updated ;
    
    if (Mb == H_Scache(0) && MW == H_Scache(1)) {
        H_Supdated = N_updated * C_S_updated * SL_updated * C_Sp_updated ;
    } else {
        H_Supdated = 0;
        H_Scache(0) = Mb;
        H_Scache(1) = MW;
    }
    
    if (Mb == H_P_cache(0) && MW == H_P_cache(1) && Mlep == H_P_cache(2) && Ms == H_P_cache(3)) {
        H_P_updated = N_updated * C_P_updated * SL_updated * C_Pp_updated * SL_updated * C_10_updated * C_10p_updated;
    } else {
        H_P_updated = 0;
        H_P_cache(0) = Mb;
        H_P_cache(1) = MW;
        H_P_cache(2) = Mlep;
        H_P_cache(3) = Ms;   
    }
    
    if (MM == T_cache(0) && Mb == T_cache(1) && Mc == T_cache(2) && 
            mySM.getMesons(pseudoscalar).getGegenalpha(0) == T_cache(3) && mySM.getMesons(pseudoscalar).getGegenalpha(1) == T_cache(4) ) {
        T_updated = 1;
    } else {
        T_updated = 0;
        T_cache(0) = MM;
        T_cache(1) = Mb;
        T_cache(2) = Mc;
        T_cache(3) = mySM.getMesons(pseudoscalar).getGegenalpha(0);
        T_cache(4) = mySM.getMesons(pseudoscalar).getGegenalpha(1);
    }

    deltaTparpupdated = C_2Lh_updated * T_updated;
    deltaTparmupdated = C_2Lh_updated * C_8Lh_updated * T_updated;
    
    I0_updated = F_updated * H_V0updated * H_A0updated * H_P_updated * beta_updated * H_Supdated * deltaTparmupdated;
    I2_updated = F_updated * beta_updated * H_V0updated * H_A0updated * deltaTparmupdated;
    I8_updated = F_updated * beta_updated * H_Supdated * H_V0updated * deltaTparmupdated;
    
}

/*******************************************************************************
 * Transverse Form Factors                                                     *
 * ****************************************************************************/
double MPll::LCSR_fit1(double q2, double r_1, double r_2, double m_fit2)
{
    return r_1/( 1. - q2/m_fit2 ) + r_2/pow( ( 1. - q2/m_fit2 ) , 2.) ;

}

double MPll::LCSR_fit2(double q2, double r_2, double m_fit2)
{
    return r_2/( 1. - q2/m_fit2 ) ; 
}

double MPll::zeta(double q2)
{
    double tp, t0;
    
    tp = (MM + MP)*(MM + MP);
    t0 = (MM + MP)*(sqrt(MM) - sqrt(MP))*(sqrt(MM) - sqrt(MP));
    
    return (sqrt(tp - q2) - sqrt(tp - t0)) / (sqrt(tp - q2) + sqrt(tp - t0));
}

double MPll::LATTICE_fit1(double q2, double b_0, double b_1, double b_2, double m_fit2)
{
    double z2 = zeta(q2)*zeta(q2);
    double z3 = zeta(q2)*z2;
    
    return 1./( 1. - q2/m_fit2 ) * ( b_0 + b_1*(zeta(q2) - 1./3.*z3) + b_2*(z2 + 2./3.*z3) ) ;

}

double MPll::LATTICE_fit2(double q2, double b_0, double b_1, double b_2, double m_fit2)
{
    return 1./( 1. - q2/m_fit2 ) * ( b_0 + b_1*zeta(q2) + b_2*zeta(q2)*zeta(q2) ) ; 
}

double MPll::f_plus(double q2)
{
    if (LATTICE)
        return LATTICE_fit1(q2, b_0_fplus, b_1_fplus, b_2_fplus, m_fit2_fplus_lat);
    else
        return LCSR_fit1(q2, r_1_fplus, r_2_fplus, m_fit2_fplus);
}

double MPll::f_T(double q2)
{
    if (LATTICE)
        return LATTICE_fit1(q2, b_0_fT, b_1_fT, b_2_fT, m_fit2_fT_lat);
    else
        return LCSR_fit1(q2, r_1_fT, r_2_fT, m_fit2_fT);
}

double MPll::f_0(double q2)
{
    if (LATTICE)
        return LATTICE_fit2(q2, b_0_f0, b_1_f0, b_2_f0, m_fit2_f0_lat);
    else
        return LCSR_fit2(q2, r_2_f0, m_fit2_f0);
}

gslpp::complex MPll::V_L(double q2)
{
    return gslpp::complex::i() * sqrt(lambda(q2)) / (twoMM*sqrt(q2)) * f_plus(q2);
}

gslpp::complex MPll::T_L(double q2)
{
    return gslpp::complex::i()  * sqrt(lambda(q2)*q2) / (twoMM2_MMpMP) * f_T(q2);
}

double MPll::S_L(double q2)
{
    return S_L_pre * f_0(q2);
}

/*******************************************************************************
 * QCDF                                                                        *
 * ****************************************************************************/

gslpp::complex MPll::I1(double u, double q2)
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

gslpp::complex MPll::Tparplus(double u, double q2) 
{
    Ee = (MM2 - q2) / twoMM;
    ubar = 1. - u;
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ (ubar * MM2 + u * q2) - 1.;
    B01 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    arg1 = (fourMc2 - gslpp::complex::i()*1.e-10)/ q2 - 1.;
    B00 = -2. * sqrt(arg1) * arctan(1. / sqrt(arg1));
    
    gslpp::complex tpar = twoMM / Ee / ubar * I1(u, q2) + (ubar * MM2 + u * q2) / Ee / Ee / ubar / ubar * (B01 - B00);
    return - MM / Mb * mySM.getQuarks(QCD::CHARM).getCharge() * tpar*C_2Lh_bar;
}

gslpp::complex MPll::Tparminus(double u, double q2) 
{
    double ubar = 1. - u;
    return - spectator_charge*(8. * C_8Lh / (ubar + u * q2 / MM2)
            + sixMMoMb * H_c(ubar * MM2 + u * q2,mu_h*mu_h) * C_2Lh_bar);
}

double MPll::Integrand_ReTparplus(double up) 
{   
    double u = up;
    return ((Tparplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/mySM.getMesons(meson).getLambdaM()).real();
}

double MPll::Integrand_ImTparplus(double up) 
{   
    double u = up;
    return ((Tparplus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/mySM.getMesons(meson).getLambdaM()).imag();
}

double MPll::Integrand_ReTparminus(double up) 
{
    double Lambdaplus = mySM.getMesons(meson).getLambdaM();
    gslpp::complex Lambdamin = exp(-tmpq2 / MM / Lambdaplus) / Lambdaplus * (-gsl_sf_expint_Ei(tmpq2 / MM / Lambdaplus) + gslpp::complex::i() * M_PI);
    
    double u = up;
    return ((Tparminus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/Lambdamin).real();
}

double MPll::Integrand_ImTparminus(double up) 
{
    double Lambdaplus = mySM.getMesons(meson).getLambdaM();
    gslpp::complex Lambdamin = exp(-tmpq2 / MM / Lambdaplus) / Lambdaplus * (-gsl_sf_expint_Ei(tmpq2 / MM / Lambdaplus) + gslpp::complex::i() * M_PI);
    
    double u = up;
    return ((Tparminus(u, tmpq2)*6. * u * (1. - u)*
            (1 + threeGegen0 * (2. * u - 1)
            + threeGegen1otwo * ((10. * u - 5.)*(2. * u - 1.) - 1.)))/Lambdamin).imag();
}

double MPll::Integrand_ImTpar_pm(double up){
    return Integrand_ImTparplus(up) + Integrand_ImTparminus(up);
}

double MPll::Integrand_ReTpar_pm(double up){
    return Integrand_ReTparplus(up) + Integrand_ReTparminus(up);
}

gslpp::complex MPll::F19(double q2) 
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

gslpp::complex MPll::F27(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    double Ls = log(s);
    gslpp::complex i = gslpp::complex::i();
    return F27_0 + F27_1 * s + F27_2 * s2 + F27_3 * s * s2 + F27_L1_1 * Ls * s + F27_L1_2 * Ls * s2 + F27_L1_3 * Ls * s * s2;
}

gslpp::complex MPll::F29(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    double Ls = log(s);
    gslpp::complex i = gslpp::complex::i();
    return F29_0 + F29_L1 * Ls + F29_1 * s +F29_2 * s2 +F29_3 * s * s2 + F29_L1_1 * Ls * s + F29_L1_2 * Ls * s2 + F29_L1_3 * Ls * s2 *s;
}

gslpp::complex MPll::F87(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    return F87_0 + F87_1 * s + F87_2 * s2 + F87_3 * s * s2 - 0.888889 * log(s)*(1. + s + s2 + s * s2);
}

double MPll::F89(double q2) 
{
    double s = q2 / Mb2;
    double s2 = s*s;
    return F89_0 + F89_1 * s + F89_2 * s2 + F89_3 * s * s2 + 1.77778 * log(s)*(1. + s + s2 + s * s2);
}

gslpp::complex MPll::Cpar(double q2) 
{
    return - (C_2L_bar * F27(q2) + C_8L * F87(q2) + MM / 2. / Mb * 
            (C_2L_bar * F29(q2) + 2.*C_1L_bar*(F19(q2) + F29(q2)/6.) + C_8L * F89(q2)));
}

gslpp::complex MPll::deltaTpar(double q2) 
{
    tmpq2 = q2;

    //old_handler = gsl_set_error_handler_off();
    
    if (deltaTparpCached[q2] == 0) {
        
        DTPPR = convertToGslFunction(boost::bind(&MPll::Integrand_ReTpar_pm, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ReTppint = avaDTPPR;
    
        DTPPR = convertToGslFunction(boost::bind(&MPll::Integrand_ImTpar_pm, &(*this), _1));
        if (gsl_integration_cquad(&DTPPR, 0., 1., 1.e-2, 1.e-1, w_DTPPR, &avaDTPPR, &errDTPPR, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double ImTppint = avaDTPPR;   
    
        cacheDeltaTparp[q2] = (ReTppint + gslpp::complex::i() * ImTppint);
        deltaTparpCached[q2] = 1;
    }

    //gsl_set_error_handler(old_handler);

    return deltaT_0 * Cpar(q2) + deltaT_1par * cacheDeltaTparp[q2] / f_plus(q2);
}

double MPll::reDC9fit(double* x, double* p)
{
    return p[0]/x[0] + p[1] + p[2]*x[0] + p[3]*x[0]*x[0] + p[4]*x[0]*x[0]*x[0] + p[5]*x[0]*x[0]*x[0]*x[0] + p[6]*x[0]*x[0]*x[0]*x[0]*x[0]; 
}

double MPll::imDC9fit(double* x, double* p)
{
    return p[0]/x[0] + p[1] + p[2]*x[0] + p[3]*x[0]*x[0] + p[4]*x[0]*x[0]*x[0] + p[5]*x[0]*x[0]*x[0]*x[0] + p[6]*x[0]*x[0]*x[0]*x[0]*x[0] + p[7]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
    
    //double thr = 4.*Mc2;
}

void MPll::fit_DeltaC9_mumu()
{
    int dim = 0;
    for (double i=0.1; i<MPllSWITCH; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9.push_back( (DC9pre * sqrt(q2tmp) * lambda(q2tmp) * deltaTpar(q2tmp)).real() );
        ImDeltaC9.push_back( (DC9pre * sqrt(q2tmp) * lambda(q2tmp) * deltaTpar(q2tmp)).imag() );
        dim++;
    }
    for (double i=MPllSWITCH; i<8.2; i+=0.4) {
        double q2tmp = i;
        myq2.push_back(q2tmp);
        ReDeltaC9.push_back(q2tmp * (DC9pre * sqrt(q2tmp) * lambda(q2tmp) * deltaTpar(q2tmp)).real());
        ImDeltaC9.push_back(q2tmp * (DC9pre * sqrt(q2tmp) * lambda(q2tmp) * deltaTpar(q2tmp)).imag());
        dim++;
    }
    
    gr1 = TGraph(dim, myq2.data(), ReDeltaC9.data());
    gr2 = TGraph(dim, myq2.data(), ImDeltaC9.data());
    
    reffit = TF1("reffit",this,&MPll::reDC9fit,0,8.1,7,"MPll","reDC9fit");
    imffit = TF1("imffit",this,&MPll::imDC9fit,0,8.1,8,"MPll","imDC9fit");
    
    refres = gr1.Fit(&reffit, "SQN0+rob=0.99");
    imfres = gr2.Fit(&imffit, "SQN0+rob=0.99");
    
    ReDeltaC9.clear();
    ImDeltaC9.clear();
    myq2.clear();
}

gslpp::complex MPll::fDeltaC9(double q2)
{
    if (q2 < MPllSWITCH) return (reDC9fit(&q2, const_cast<double *>(refres->GetParams())) 
            + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres->GetParams())));
    else return (reDC9fit(&q2, const_cast<double *>(refres->GetParams())) 
            + gslpp::complex::i()*imDC9fit(&q2, const_cast<double *>(imfres->GetParams())))/q2;
            
}


gslpp::complex MPll::DeltaC9(double q2)
{
    return DC9pre * sqrt(q2) * lambda(q2) * deltaTpar(q2);
}

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MPll::H_c(double q2, double mu2) 
{
    double x = fourMc2 / q2;
    gslpp::complex par;

    if (x > 1.) par = sqrt(x - 1.) * atan(1. / sqrt(x - 1.));
    else par = sqrt(1. - x) * (log((1. + sqrt(1. - x)) / sqrt(x)) - ihalfMPI);
    return -fournineth * (log(Mc2/mu2) - twothird - x) - fournineth * (2. + x) * par;
}

gslpp::complex MPll::H_b(double q2, double mu2) 
{
    double x = fourMb2 / q2;
    gslpp::complex par;

    if (x > 1.) par = sqrt(x - 1.) * atan(1. / sqrt(x - 1.));
    else par = sqrt(1. - x) * (log((1. + sqrt(1. - x)) / sqrt(x)) - ihalfMPI);

    return -fournineth * (log(Mb2/mu2) - twothird - x) - fournineth * (2. + x) * par;
}

gslpp::complex MPll::H_0(double q2) 
{
    return (H_0_pre - fournineth * log(q2 / mu_b2));
}

gslpp::complex MPll::Y(double q2) 
{
    return -half * H_0(q2) * H_0_WC + H_c(q2,mu_b2) * H_c_WC - half * H_b(q2,mu_b2) * H_b_WC;
}

gslpp::complex MPll::H_V(double q2) 
{
    return -( (C_9 + Y(q2) + fDeltaC9(q2) - etaP*pow(-1,angmomP)*C_9p)*V_L(q2) + MM2/q2*( twoMboMM*(C_7 - etaP*pow(-1,angmomP)*C_7p)*T_L(q2) - sixteenM_PI2*(h_0 + h_1 * q2)) );
}

gslpp::complex MPll::H_A(double q2) 
{
    return (- C_10 + etaP*pow(-1,angmomP)*C_10p) *V_L(q2);
}

gslpp::complex MPll::H_S(double q2) 
{
    return MboMW*(C_S - etaP*pow(-1,angmomP)*C_Sp)*S_L(q2) ;
}

gslpp::complex MPll::H_P(double q2) 
{
    return ( MboMW * (C_P - etaP*pow(-1,angmomP)*C_Pp) + twoMlepMb / q2 * ( C_10*(1. + etaP*pow(-1,angmomP)*MsoMb) - C_10p*(etaP*pow(-1,angmomP) + MsoMb) )) * S_L(q2);
}

/*******************************************************************************
 * Angular coefficients                                                         *
 * ****************************************************************************/
double MPll::k2(double q2) 
{
    return (MM4 + q2 * q2 + MP4 - twoMP2 * q2 - twoMM2 * (q2 + MP2)) / fourMM2;
}

double MPll::beta(double q2) 
{
    return sqrt(1. - 4.*Mlep2/q2);
}

double MPll::beta2(double q2) 
{
    return 1. - 4. * Mlep2 / q2;
}
double MPll::lambda(double q2) 
{
    return 4.*MM2*k2(q2);
}

double MPll::F(double q2) 
{
    return sqrt(lambda(q2))*beta(q2)*q2/(ninetysixM_PI3MM3);
}

double MPll::I_1c(double q2) 
{
    return F(q2)*((H_V(q2).abs2() + H_A(q2).abs2()) / 2. + H_P(q2).abs2() + 2. * Mlep2 / q2 * (H_V(q2).abs2()
            - H_A(q2).abs2()) + beta2(q2) * H_S(q2).abs2());
}

double MPll::I_2c(double q2) 
{
    return -F(q2) * beta2(q2) / 2. * (H_V(q2).abs2() + H_A(q2).abs2());
}

double MPll::I_6c(double q2) 
{
    return 4. * F(q2) * beta(q2) * Mlep / sqrt(q2)*(H_S(q2).conjugate() * H_V(q2)).real();
}

double MPll::Delta(int i, double q2) 
{
    return 0; /* FIX CPV */
    //return (I(i, q2,0) - I(i, q2,1))/2;
}

double MPll::integrateSigma(int i, double q_min, double q_max)
{
    updateParameters();

    std::pair<double, double > qbin = std::make_pair(q_min, q_max);
    
    old_handler = gsl_set_error_handler_off();

    switch(i){
        case 0:
            if (sigma0Cached[qbin] == 0) {
                FS = convertToGslFunction( boost::bind( &MPll::getSigma1c, &(*this), _1 ) );
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma0[qbin] = NN*avaSigma;
                sigma0Cached[qbin] = 1;
            }
            return cacheSigma0[qbin];
            break;
        case 2:
            if (sigma2Cached[qbin] == 0) {
                FS = convertToGslFunction( boost::bind( &MPll::getSigma2c, &(*this), _1 ) );
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma2[qbin] = NN*avaSigma;
                sigma2Cached[qbin] = 1;
            }
            return cacheSigma2[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MPll::integrateSigma: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);
    
}

double MPll::getSigma(int i, double q_2) 
{
    updateParameters();

    switch (i) {
        case 0:
            return getSigma1c(q_2);
            break;
        case 2:
            return getSigma2c(q_2);
            break;
        case 8:
            return getSigma6c(q_2);
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MPll::getSigma: index " + out.str() + " not implemented");
    }
}

double MPll::integrateDelta(int i, double q_min, double q_max)
{
    updateParameters();

    std::pair<double, double > qbin = std::make_pair(q_min, q_max);

    old_handler = gsl_set_error_handler_off();

    switch(i){
        case 0:
            if (delta0Cached[qbin] == 0) {
                FD = convertToGslFunction( boost::bind( &MPll::getDelta1c, &(*this), _1 ) );
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta0[qbin] = NN*avaDelta;
                delta0Cached[qbin] = 1;
            }
            return cacheDelta0[qbin];
            break;
        case 2:
            if (delta2Cached[qbin] == 0) {
                FD = convertToGslFunction( boost::bind( &MPll::getDelta2c, &(*this), _1 ) );
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta2[qbin] = NN*avaDelta;
                delta2Cached[qbin] = 1;
            }
            return cacheDelta2[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MPll::integrateDelta: index " + out.str() + " not implemented"); 
    }

    gsl_set_error_handler(old_handler);
 
}