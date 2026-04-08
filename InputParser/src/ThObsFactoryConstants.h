/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THOBSFACTORYCONSTANTS_H
#define THOBSFACTORYCONSTANTS_H

/**
 * @brief Namespace containing collider energy and angle constants
 * used across ThObsFactory registration functions.
 */
namespace ThObsConst {

    //-----  Energies of different colliders  -----
    const double sqrt_s_LEP2_WWcos1 = 0.18266; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos2 = 0.18909; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos3 = 0.19838; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_WWcos4 = 0.20592; ///< the center-of-mass energy in TeV
    //
    const double cos1_LEP2_WW = -0.95;
    const double cos2_LEP2_WW = -0.76;
    const double cos3_LEP2_WW = -0.57;
    const double cos4_LEP2_WW = -0.38;
    const double cos5_LEP2_WW = -0.19;
    const double cos6_LEP2_WW = 0.;
    const double cos7_LEP2_WW = 0.19;
    const double cos8_LEP2_WW = 0.38;
    const double cos9_LEP2_WW = 0.57;
    const double cos10_LEP2_WW = 0.76;
    const double cos11_LEP2_WW = 0.95;
    //
    const double cos1_ee_WW = -1.0;
    const double cos2_ee_WW = -0.8;
    const double cos3_ee_WW = -0.6;
    const double cos4_ee_WW = -0.4;
    const double cos5_ee_WW = -0.2;
    const double cos6_ee_WW = 0.;
    const double cos7_ee_WW = 0.2;
    const double cos8_ee_WW = 0.4;
    const double cos9_ee_WW = 0.6;
    const double cos10_ee_WW = 0.8;
    const double cos11_ee_WW = 1.0;
    //
    const double sqrt_s_LEP2_161 = 0.1613; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_172 = 0.1721; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_183 = 0.1827; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_189 = 0.1886; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_192 = 0.1916; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_196 = 0.1955; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_200 = 0.1995; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_202 = 0.2016; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_205 = 0.2049; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_206 = 0.2059; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_207 = 0.2066; ///< the center-of-mass energy in TeV
    const double sqrt_s_LEP2_208 = 0.208; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_LHC7 = 7.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC8 = 8.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC13 = 13.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC14 = 14.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHC27 = 27.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCC50 = 50.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCC84 = 84.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCC100 = 100.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_TeV = 1.96; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_leptcoll_125 = .125; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_161 = .161; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_230 = .230; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_240 = .240; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_250 = .250; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_350 = .350; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_365 = .365; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_380 = .380; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_500 = .500; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_550 = .550; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1000 = 1.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1400 = 1.4; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_1500 = 1.5; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_3000 = 3.0; ///< the center-of-mass energy in TeV
    const double sqrt_s_leptcoll_10000 = 10.0; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_LHeC_1_2 = 1.2; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHeC_1_3 = 1.3; ///< the center-of-mass energy in TeV
    const double sqrt_s_LHeC_1_8 = 1.8; ///< the center-of-mass energy in TeV
    //
    const double sqrt_s_FCCep_3_5 = 3.5; ///< the center-of-mass energy in TeV
    const double sqrt_s_FCCep_5 = 5.0; ///< the center-of-mass energy in TeV
    // Polarizations at lepton colliders
    const double pol_0 = 0.0;
    const double pol_20 = 20.0;
    const double pol_30 = 30.0;
    const double pol_80 = 80.0;
    // Lists with values of energies/angles for energy/angle dependent definitions
    //
    // Parameters for LEP 2 inclusive observables
    const double sqrt_s[12] = {130., 136., 161., 172., 183., 189.,
        192., 196., 200., 202., 205., 207.};
    const double sqrt_s_LEP2eeff[12] = {130.2, 136.2, 161.3, 172.1, 182.7, 188.6,
        191.6, 195.5, 199.5, 201.6, 204.9, 206.7};
    //
    // Parameters for LEP2 differential observables
    const double sqrt_sDiffll[8] = {183., 189., 192., 196., 200., 202., 205., 207.};
    const double cos_Diffll[10] = {-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9};
    const double cosmin_Diffll[10] = {-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8};
    const double cosmax_Diffll[10] = {-0.8, -0.6, -0.4, -0.2,  0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    //
    const double cos_DiffeeInp[15] = {-0.8,-0.6,-0.5,-0.3,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
    const double cos_Diffee[15] = {-0.81,-0.63,-0.45,-0.27,-0.09,0.045,0.135,0.225,0.315,0.405,0.495,0.585,0.675,0.765,0.855};
    const double cosmin_Diffee[15] = {-0.90,-0.72,-0.54,-0.36,-0.18,0.0 ,0.09,0.18,0.27,0.36,0.45,0.54,0.63,0.72,0.81};
    const double cosmax_Diffee[15] = {-0.72,-0.54,-0.36,-0.18, 0.0 ,0.09,0.18,0.27,0.36,0.45,0.54,0.63,0.72,0.81,0.90};
    //
    // Parameters for future e+e- observables
    const double sqrt_see[16] = {158., 161., 163., 230., 240., 250., 345., 350., 360., 365., 380., 500., 550., 1000., 1500., 3000.};
    const double sqrt_s_eeff[16] = {157.5, 161., 162.5, 230., 240., 250., 345., 350., 360., 365., 380., 500., 550., 1000., 1500., 3000.};
    // Approximate electroweak scale, taken as the Z mass, and Higgs mass, taken as 125 GeV
    const double muEW = 91.188;
    const double muMH = 125.;
    //
    // Parameters for LEP 2 heavy flavour observables
    const double sqrt_s_HF[10] = {133., 167., 183., 189., 192.,
        196., 200., 202., 205., 207.};
    const double sqrt_s_LEP2_HF[10] = {133.2, 166.6, 182.7, 188.6, 191.6,
        195.5, 199.5, 201.6, 204.9, 206.7};

} // namespace ThObsConst

#endif /* THOBSFACTORYCONSTANTS_H */
