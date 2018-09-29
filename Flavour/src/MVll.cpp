/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "MVll.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include "F_1.h"
#include "F_2.h"
#include <gsl/gsl_sf_zeta.h>
#include <boost/bind.hpp>
#include <limits>
#include <TFitResult.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_expint.h>

MVll::MVll(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: mySM(SM_i), myF_1(new F_1()), myF_2(new F_2()),
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
    dispersion = false;
    FixedWCbtos = false;
    mJ2 = 3.096 * 3.096;

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

    w_sigma = gsl_integration_cquad_workspace_alloc(100);
    //    w_DTPPR = gsl_integration_cquad_workspace_alloc (100);
    w_delta = gsl_integration_cquad_workspace_alloc(100);

    acc_Re_T_perp = gsl_interp_accel_alloc();
    acc_Im_T_perp = gsl_interp_accel_alloc();
    acc_Re_T_para = gsl_interp_accel_alloc();
    acc_Im_T_para = gsl_interp_accel_alloc();

    spline_Re_T_perp = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Im_T_perp = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Re_T_para = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Im_T_para = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);

#if COMPUTECP
    acc_Re_T_perp_conj = gsl_interp_accel_alloc();
    acc_Im_T_perp_conj = gsl_interp_accel_alloc();
    acc_Re_T_para_conj = gsl_interp_accel_alloc();
    acc_Im_T_para_conj = gsl_interp_accel_alloc();

    spline_Re_T_perp_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Im_T_perp_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Re_T_para_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
    spline_Im_T_para_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM);
#endif

    acc_Re_deltaC7_QCDF = gsl_interp_accel_alloc();
    acc_Im_deltaC7_QCDF = gsl_interp_accel_alloc();
    acc_Re_deltaC9_QCDF = gsl_interp_accel_alloc();
    acc_Im_deltaC9_QCDF = gsl_interp_accel_alloc();

    spline_Re_deltaC7_QCDF = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Im_deltaC7_QCDF = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Re_deltaC9_QCDF = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Im_deltaC9_QCDF = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);

#if COMPUTECP
    acc_Re_deltaC7_QCDF_conj = gsl_interp_accel_alloc();
    acc_Im_deltaC7_QCDF_conj = gsl_interp_accel_alloc();
    acc_Re_deltaC9_QCDF_conj = gsl_interp_accel_alloc();
    acc_Im_deltaC9_QCDF_conj = gsl_interp_accel_alloc();

    spline_Re_deltaC7_QCDF_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Im_deltaC7_QCDF_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Re_deltaC9_QCDF_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
    spline_Im_deltaC9_QCDF_conj = gsl_spline_alloc(gsl_interp_cspline, GSL_INTERP_DIM_DC);
#endif

    h_pole = false;

    M_PI2 = M_PI*M_PI;

    F87_1 = (4. / 3. * M_PI2 - 40. / 3.);
    F87_2 = (32. / 9. * M_PI2 - 316. / 9.);
    F87_3 = (200. / 27. * M_PI2 - 658. / 9.);

    F89_0 = (104. / 9. - 32. / 27. * M_PI2);
    F89_1 = (1184. / 27. - 40. / 9. * M_PI2);
    F89_2 = (-32. / 3. * M_PI2 + 14212. / 135.);
    F89_3 = (-560. / 27. * M_PI2 + 193444. / 945.);

    CF = 4. / 3.;

}

MVll::~MVll()
{
}

std::vector<std::string> MVll::initializeMVllParameters()
{
    dispersion = mySM.getFlavour().getFlagUseDispersionRelation();
    FixedWCbtos = mySM.getFlavour().getFlagFixedWCbtos();
    
#if NFPOLARBASIS_MVLL
    if (vectorM == StandardModel::PHI) mvllParameters = make_vector<std::string>()
        << "a_0Vphi" << "a_1Vphi" << "a_2Vphi" << "MRV" << "a_0A0phi" << "a_1A0phi" << "a_2A0phi" << "MRA0"
        << "a_0A1phi" << "a_1A1phi" << "a_2A1phi" << "MRA1" << "a_1A12phi" << "a_2A12phi" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
        << "a_0T1phi" << "a_1T1phi" << "a_2T1phi" << "MRT1" << "a_1T2phi" << "a_2T2phi" << "MRT2"
        << "a_0T23phi" << "a_1T23phi" << "a_2T23phi" << "MRT23"
        << "absh_0" << "absh_p" << "absh_m" << "argh_0" << "argh_p" << "argh_m"
        << "absh_0_1" << "absh_p_1" << "absh_m_1" << "argh_0_1" << "argh_p_1" << "argh_m_1"
        << "absh_p_2" << "absh_m_2" << "argh_p_2" << "argh_m_2" << "xs_phi" << "SU3_breaking_abs" << "SU3_breaking_arg";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) mvllParameters = make_vector<std::string>()
        << "a_0V" << "a_1V" << "a_2V" << "MRV" << "a_0A0" << "a_1A0" << "a_2A0" << "MRA0"
        << "a_0A1" << "a_1A1" << "a_2A1" << "MRA1" << "a_1A12" << "a_2A12" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
        << "a_0T1" << "a_1T1" << "a_2T1" << "MRT1" << "a_1T2" << "a_2T2" << "MRT2"
        << "a_0T23" << "a_1T23" << "a_2T23" << "MRT23"
        << "absh_0" << "absh_p" << "absh_m" << "argh_0" << "argh_p" << "argh_m"
        << "absh_0_1" << "absh_p_1" << "absh_m_1" << "argh_0_1" << "argh_p_1" << "argh_m_1"
        << "absh_p_2" << "absh_m_2" << "argh_p_2" << "argh_m_2";
#else 
    if (vectorM == StandardModel::PHI) mvllParameters = make_vector<std::string>()
        << "a_0Vphi" << "a_1Vphi" << "a_2Vphi" << "MRV" << "a_0A0phi" << "a_1A0phi" << "a_2A0phi" << "MRA0"
        << "a_0A1phi" << "a_1A1phi" << "a_2A1phi" << "MRA1" << "a_1A12phi" << "a_2A12phi" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
        << "a_0T1phi" << "a_1T1phi" << "a_2T1phi" << "MRT1" << "a_1T2phi" << "a_2T2phi" << "MRT2"
        << "a_0T23phi" << "a_1T23phi" << "a_2T23phi" << "MRT23"
        << "reh_0" << "reh_p" << "reh_m" << "imh_0" << "imh_p" << "imh_m"
        << "reh_0_1" << "reh_p_1" << "reh_m_1" << "imh_0_1" << "imh_p_1" << "imh_m_1"
        << "reh_p_2" << "reh_m_2" << "imh_p_2" << "imh_m_2" << "xs_phi" << "SU3_breaking_abs" << "SU3_breaking_arg";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) mvllParameters = make_vector<std::string>()
        << "a_0V" << "a_1V" << "a_2V" << "MRV" << "a_0A0" << "a_1A0" << "a_2A0" << "MRA0"
        << "a_0A1" << "a_1A1" << "a_2A1" << "MRA1" << "a_1A12" << "a_2A12" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
        << "a_0T1" << "a_1T1" << "a_2T1" << "MRT1" << "a_1T2" << "a_2T2" << "MRT2"
        << "a_0T23" << "a_1T23" << "a_2T23" << "MRT23"
        << "reh_0" << "reh_p" << "reh_m" << "imh_0" << "imh_p" << "imh_m"
        << "reh_0_1" << "reh_p_1" << "reh_m_1" << "imh_0_1" << "imh_p_1" << "imh_m_1"
        << "reh_p_2" << "reh_m_2" << "imh_p_2" << "imh_m_2";
#endif
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVll: vector " + out.str() + " not implemented");
    }

    if (dispersion) {
        mvllParameters.clear();
        if (vectorM == StandardModel::PHI) mvllParameters = make_vector<std::string>()
            << "a_0Vphi" << "a_1Vphi" << "a_2Vphi" << "MRV" << "a_0A0phi" << "a_1A0phi" << "a_2A0phi" << "MRA0"
            << "a_0A1phi" << "a_1A1phi" << "a_2A1phi" << "MRA1" << "a_1A12phi" << "a_2A12phi" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
            << "a_0T1phi" << "a_1T1phi" << "a_2T1phi" << "MRT1" << "a_1T2phi" << "a_2T2phi" << "MRT2"
            << "a_0T23phi" << "a_1T23phi" << "a_2T23phi" << "MRT23"
            << "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1"
            << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2"
            << "r1_3" << "r2_3" << "deltaC9_3" << "phDC9_3" << "xs_phi" << "SU3_breaking_abs" << "SU3_breaking_arg";
        else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) mvllParameters = make_vector<std::string>()
            << "a_0V" << "a_1V" << "a_2V" << "MRV" << "a_0A0" << "a_1A0" << "a_2A0" << "MRA0"
            << "a_0A1" << "a_1A1" << "a_2A1" << "MRA1" << "a_1A12" << "a_2A12" << "MRA12" /*a_0A12 and a_0T2 are not independent*/
            << "a_0T1" << "a_1T1" << "a_2T1" << "MRT1" << "a_1T2" << "a_2T2" << "MRT2"
            << "a_0T23" << "a_1T23" << "a_2T23" << "MRT23"
            << "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1"
            << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2"
            << "r1_3" << "r2_3" << "deltaC9_3" << "phDC9_3";
    }

    if (FixedWCbtos) mvllParameters.insert(mvllParameters.end(), { "C7_SM", "C9_SM", "C10_SM" });
    
    mySM.initializeMeson(meson);
    mySM.initializeMeson(vectorM);
    return mvllParameters;
}

void MVll::updateParameters()
{
    if (!mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) return;


    GF = mySM.getGF();
    ale = mySM.getAle();
    Mlep = mySM.getLeptons(lep).getMass();
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    mu_b = mySM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    mb_pole = mySM.Mbar2Mp(Mb); /* Conversion to pole mass*/
    mc_pole = mySM.Mbar2Mp(Mc); /* Conversion to pole mass*/
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MW = mySM.Mw();
    lambda_t = mySM.getCKM().computelamt_s();
    lambda_u = mySM.getCKM().computelamu_s();
    width = mySM.getMesons(meson).computeWidth();
    alpha_s_mub = mySM.Als(mu_b);
    fB = mySM.getMesons(meson).getDecayconst();
    fpara = mySM.getMesons(vectorM).getDecayconst();
    fperp = mySM.getMesons(vectorM).getDecayconst_p();

    switch (vectorM) {
        case StandardModel::K_star:
        case StandardModel::K_star_P:
            a_0V = mySM.getOptionalParameter("a_0V");
            a_1V = mySM.getOptionalParameter("a_1V");
            a_2V = mySM.getOptionalParameter("a_2V");
            MRV_2 = mySM.getOptionalParameter("MRV") * mySM.getOptionalParameter("MRV");

            a_0A0 = mySM.getOptionalParameter("a_0A0");
            a_1A0 = mySM.getOptionalParameter("a_1A0");
            a_2A0 = mySM.getOptionalParameter("a_2A0");
            MRA0_2 = mySM.getOptionalParameter("MRA0") * mySM.getOptionalParameter("MRA0");

            a_0A1 = mySM.getOptionalParameter("a_0A1");
            a_1A1 = mySM.getOptionalParameter("a_1A1");
            a_2A1 = mySM.getOptionalParameter("a_2A1");
            MRA1_2 = mySM.getOptionalParameter("MRA1") * mySM.getOptionalParameter("MRA1");

            a_0A12 = a_0A0 * (MM * MM - MV * MV) / (8. * MM * MV);
            a_1A12 = mySM.getOptionalParameter("a_1A12");
            a_2A12 = mySM.getOptionalParameter("a_2A12");
            MRA12_2 = mySM.getOptionalParameter("MRA12") * mySM.getOptionalParameter("MRA12");

            a_0T1 = mySM.getOptionalParameter("a_0T1");
            a_1T1 = mySM.getOptionalParameter("a_1T1");
            a_2T1 = mySM.getOptionalParameter("a_2T1");
            MRT1_2 = mySM.getOptionalParameter("MRT1") * mySM.getOptionalParameter("MRT1");

            a_0T2 = a_0T1;
            a_1T2 = mySM.getOptionalParameter("a_1T2");
            a_2T2 = mySM.getOptionalParameter("a_2T2");
            MRT2_2 = mySM.getOptionalParameter("MRT2") * mySM.getOptionalParameter("MRT2");

            a_0T23 = mySM.getOptionalParameter("a_0T23");
            a_1T23 = mySM.getOptionalParameter("a_1T23");
            a_2T23 = mySM.getOptionalParameter("a_2T23");
            MRT23_2 = mySM.getOptionalParameter("MRT23") * mySM.getOptionalParameter("MRT23");


            if (vectorM == StandardModel::K_star) spectator_charge = mySM.getQuarks(QCD::DOWN).getCharge();
            else spectator_charge = mySM.getQuarks(QCD::UP).getCharge();

            etaV = -1;
            angmomV = 1.;

            b = 1;

            SU3_breaking = 1.;

            break;
        case StandardModel::PHI:
            a_0V = mySM.getOptionalParameter("a_0Vphi");
            a_1V = mySM.getOptionalParameter("a_1Vphi");
            a_2V = mySM.getOptionalParameter("a_2Vphi");
            MRV_2 = mySM.getOptionalParameter("MRV") * mySM.getOptionalParameter("MRV");

            a_0A0 = mySM.getOptionalParameter("a_0A0phi");
            a_1A0 = mySM.getOptionalParameter("a_1A0phi");
            a_2A0 = mySM.getOptionalParameter("a_2A0phi");
            MRA0_2 = mySM.getOptionalParameter("MRA0") * mySM.getOptionalParameter("MRA0");

            a_0A1 = mySM.getOptionalParameter("a_0A1phi");
            a_1A1 = mySM.getOptionalParameter("a_1A1phi");
            a_2A1 = mySM.getOptionalParameter("a_2A1phi");
            MRA1_2 = mySM.getOptionalParameter("MRA1") * mySM.getOptionalParameter("MRA1");

            a_0A12 = a_0A0 * (MM * MM - MV * MV) / (8. * MM * MV);
            a_1A12 = mySM.getOptionalParameter("a_1A12phi");
            a_2A12 = mySM.getOptionalParameter("a_2A12phi");
            MRA12_2 = mySM.getOptionalParameter("MRA12") * mySM.getOptionalParameter("MRA12");

            a_0T1 = mySM.getOptionalParameter("a_0T1phi");
            a_1T1 = mySM.getOptionalParameter("a_1T1phi");
            a_2T1 = mySM.getOptionalParameter("a_2T1phi");
            MRT1_2 = mySM.getOptionalParameter("MRT1") * mySM.getOptionalParameter("MRT1");

            a_0T2 = a_0T1;
            a_1T2 = mySM.getOptionalParameter("a_1T2phi");
            a_2T2 = mySM.getOptionalParameter("a_2T2phi");
            MRT2_2 = mySM.getOptionalParameter("MRT2") * mySM.getOptionalParameter("MRT2");

            a_0T23 = mySM.getOptionalParameter("a_0T23phi");
            a_1T23 = mySM.getOptionalParameter("a_1T23phi");
            a_2T23 = mySM.getOptionalParameter("a_2T23phi");
            MRT23_2 = mySM.getOptionalParameter("MRT23") * mySM.getOptionalParameter("MRT23");

            spectator_charge = mySM.getQuarks(QCD::STRANGE).getCharge();

            ys = mySM.getMesons(QCD::B_S).getDgamma_gamma() / 2.;
            xs = mySM.getOptionalParameter("xs_phi");

            etaV = -1;
            angmomV = 1.;

            b = 0.489;

            SU3_breaking = 1. + gslpp::complex(mySM.getOptionalParameter("SU3_breaking_abs"),
                    mySM.getOptionalParameter("SU3_breaking_arg"), true);

            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVll: vector " + out.str() + " not implemented");
    }

    if (!dispersion) {
#if NFPOLARBASIS_MVLL
        h_0[0] = gslpp::complex(mySM.getOptionalParameter("absh_0"), mySM.getOptionalParameter("argh_0"), true);
        h_0[1] = gslpp::complex(mySM.getOptionalParameter("absh_p"), mySM.getOptionalParameter("argh_p"), true);
        h_0[2] = gslpp::complex(mySM.getOptionalParameter("absh_m"), mySM.getOptionalParameter("argh_m"), true);

        h_1[0] = gslpp::complex(mySM.getOptionalParameter("absh_0_1"), mySM.getOptionalParameter("argh_0_1"), true);
        h_1[1] = gslpp::complex(mySM.getOptionalParameter("absh_p_1"), mySM.getOptionalParameter("argh_p_1"), true);
        h_1[2] = gslpp::complex(mySM.getOptionalParameter("absh_m_1"), mySM.getOptionalParameter("argh_m_1"), true);

        h_2[0] = 0.;
        h_2[1] = gslpp::complex(mySM.getOptionalParameter("absh_p_2"), mySM.getOptionalParameter("argh_p_2"), true);
        h_2[2] = gslpp::complex(mySM.getOptionalParameter("absh_m_2"), mySM.getOptionalParameter("argh_m_2"), true);
#else
        h_0[0] = gslpp::complex(mySM.getOptionalParameter("reh_0"), mySM.getOptionalParameter("imh_0"), false);
        h_0[1] = gslpp::complex(mySM.getOptionalParameter("reh_p"), mySM.getOptionalParameter("imh_p"), false);
        h_0[2] = gslpp::complex(mySM.getOptionalParameter("reh_m"), mySM.getOptionalParameter("imh_m"), false);

        h_1[0] = gslpp::complex(mySM.getOptionalParameter("reh_0_1"), mySM.getOptionalParameter("imh_0_1"), false);
        h_1[1] = gslpp::complex(mySM.getOptionalParameter("reh_p_1"), mySM.getOptionalParameter("imh_p_1"), false);
        h_1[2] = gslpp::complex(mySM.getOptionalParameter("reh_m_1"), mySM.getOptionalParameter("imh_m_1"), false);

        h_2[0] = 0.;
        h_2[1] = gslpp::complex(mySM.getOptionalParameter("reh_p_2"), mySM.getOptionalParameter("imh_p_2"), false);
        h_2[2] = gslpp::complex(mySM.getOptionalParameter("reh_m_2"), mySM.getOptionalParameter("imh_m_2"), false);
#endif
    } else {
        h_0[0] = gslpp::complex(mySM.getOptionalParameter("r1_1"));
        h_0[1] = gslpp::complex(mySM.getOptionalParameter("r1_2"));
        h_0[2] = gslpp::complex(mySM.getOptionalParameter("r1_3"));

        h_1[0] = gslpp::complex(mySM.getOptionalParameter("r2_1"));
        h_1[1] = gslpp::complex(mySM.getOptionalParameter("r2_2"));
        h_1[2] = gslpp::complex(mySM.getOptionalParameter("r2_3"));

        h_2[0] = gslpp::complex(mySM.getOptionalParameter("deltaC9_1"));
        h_2[1] = gslpp::complex(mySM.getOptionalParameter("deltaC9_2"));
        h_2[2] = gslpp::complex(mySM.getOptionalParameter("deltaC9_3"));
        exp_Phase[0] = exp(gslpp::complex::i() * mySM.getOptionalParameter("phDC9_1"));
        exp_Phase[1] = exp(gslpp::complex::i() * mySM.getOptionalParameter("phDC9_2"));
        exp_Phase[2] = exp(gslpp::complex::i() * mySM.getOptionalParameter("phDC9_3"));
    }

    allcoeff = mySM.getFlavour().ComputeCoeffBMll(mu_b, lep); //check the mass scale, scheme fixed to NDR
    allcoeffprime = mySM.getFlavour().ComputeCoeffprimeBMll(mu_b, lep); //check the mass scale, scheme fixed to NDR

    C_1 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
    C_1L_bar = (*(allcoeff[LO]))(0) / 2.;
    C_2 = ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
    C_2L_bar = (*(allcoeff[LO]))(1) - (*(allcoeff[LO]))(0) / 6.;
    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
    C_4 = ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3));
    C_5 = ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));
    C_8 = ((*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7));
    C_8L = (*(allcoeff[LO]))(7);
    C_S = MW / Mb * (((*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10)));
    C_P = MW / Mb * (((*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11)));
    C_9p = (*(allcoeffprime[LO]))(8) + (*(allcoeffprime[NLO]))(8);
    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_Sp = MW / Mb * ((*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10));
    C_Pp = MW / Mb * ((*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11));
    
    if (FixedWCbtos) { /** NOTE: ComputeCoeff with different argumetns cannot be mixed. They have to be called sequentially. **/
        allcoeff_noSM = mySM.getFlavour().ComputeCoeffBMll(mu_b, lep, true); //check the mass scale, scheme fixed to NDR
        C_7 = mySM.getOptionalParameter("C7_SM") + ((*(allcoeff_noSM[LO]))(6) + (*(allcoeff_noSM[NLO]))(6));
        C_9 = mySM.getOptionalParameter("C9_SM") + ((*(allcoeff_noSM[LO]))(8) + (*(allcoeff_noSM[NLO]))(8));
        C_10 = mySM.getOptionalParameter("C10_SM") + ((*(allcoeff_noSM[LO]))(9) + (*(allcoeff_noSM[NLO]))(9));
    } else {
        C_7 = ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6));
        C_9 = ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8));
        C_10 = ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
    }
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    C_7p += MsoMb * C_7 + 1. / 3. * C_3 + 4 / 9 * C_4 + 20. / 3. * C_5 + 80. / 9. * C_6;
    
    allcoeffh = mySM.getFlavour().ComputeCoeffBMll(mu_h, lep); //check the mass scale, scheme fixed to NDR

    C_1Lh_bar = (*(allcoeffh[LO]))(0) / 2.;
    C_2Lh_bar = (*(allcoeffh[LO]))(1) - (*(allcoeff[LO]))(0) / 6.;
    C_8Lh = (*(allcoeffh[LO]))(7);

    checkCache();

    t_p = pow(MM + MV, 2.);
    t_m = pow(MM - MV, 2.);
    t_0 = t_p * (1. - sqrt(1. - t_m / t_p)); /*Modify it for Lattice*/
    z_0 = (sqrt(t_p) - sqrt(t_p - t_0)) / (sqrt(t_p) + sqrt(t_p - t_0));
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
    MM2pMV2 = MM2 + MV2;
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
    sixteenM_PI2 = 16. * M_PI2;
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

    M_PI2osix = M_PI2 / 6.;
    twoMM = 2. * MM;

    N_QCDF = M_PI2 / 3. * fB * fperp / MM;

    deltaT_0 = alpha_s_mub * CF / 4. / M_PI;
    deltaT_1par = mySM.Als(mu_h) * CF / 4. * M_PI / 3. * mySM.getMesons(meson).getDecayconst() *
            mySM.getMesons(vectorM).getDecayconst() / MM;
    deltaT_1perp = mySM.Als(mu_h) * CF / 4. * M_PI / 3. * mySM.getMesons(meson).getDecayconst() *
            mySM.getMesons(vectorM).getDecayconst_p() / MM;

    F87_0 = -32. / 9. * log(mu_b / Mb) + 8. / 27. * M_PI2 - 44. / 9. - 8. / 9. * gslpp::complex::i() * M_PI;

    NN = -(4. * GF * MM * ale * lambda_t) / (sqrt(2.)*4. * M_PI);
    NN_conjugate = -(4. * GF * MM * ale * lambda_t.conjugate()) / (sqrt(2.)*4. * M_PI);

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

    if (deltaTparpupdated * deltaTparmupdated == 0) for (it = I1Cached.begin(); it != I1Cached.end(); ++it) it->second = 0;

#if SPLINE    
    if (mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) spline_QCDF_func();
#else
    if (mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) fit_QCDF_func();
#endif

    mySM.getFlavour().setUpdateFlag(meson, vectorM, lep, false);
    
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

    if (h_0[0] == h0Ccache[0] && h_1[0] == h0Ccache[1] && h_2[0] == h0Ccache[2] && SU3_breaking == h0Ccache[3]) {
        h0_updated = 1;
    } else {
        h0_updated = 0;
        h0Ccache[0] = h_0[0];
        h0Ccache[1] = h_1[0];
        h0Ccache[2] = h_2[0];
        h0Ccache[3] = SU3_breaking;
    }

    if (h_0[1] == h1Ccache[0] && h_1[1] == h1Ccache[1] && h_2[1] == h1Ccache[2] && SU3_breaking == h1Ccache[3]) {
        h1_updated = 1;
    } else {
        h1_updated = 0;
        h1Ccache[0] = h_0[1];
        h1Ccache[1] = h_1[1];
        h1Ccache[2] = h_2[1];
        h1Ccache[3] = SU3_breaking;
    }

    if (h_0[2] == h2Ccache[0] && h_1[2] == h2Ccache[1] && h_2[2] == h2Ccache[2] && SU3_breaking == h2Ccache[3]) {
        h2_updated = 1;
    } else {
        h2_updated = 0;
        h2Ccache[0] = h_0[2];
        h2Ccache[1] = h_1[2];
        h2Ccache[2] = h_2[2];
        h2Ccache[3] = SU3_breaking;
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
            mySM.getMesons(vectorM).getGegenalpha(0) == T_cache(3) && mySM.getMesons(vectorM).getGegenalpha(1) == T_cache(4)) {
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
    I2_updated = F_updated * beta_updated * H_V0updated * H_A0updated * deltaTparmupdated;
    I3_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * beta_updated * deltaTparmupdated;
    I4_updated = F_updated * H_V1updated * H_V2updated * H_A1updated * H_A2updated * deltaTparmupdated;
    I5_updated = F_updated * H_V0updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * beta_updated * deltaTparmupdated;
    I6_updated = F_updated * H_V1updated * H_V2updated * H_A0updated * H_A1updated * H_A2updated * H_V0updated * beta_updated * H_Supdated * deltaTparmupdated;
    I7_updated = I4_updated * beta_updated;
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
    return (MMpMV2 * (MM2mMV2 - q2) * A_1(q2) - 16. * MM * MV2 * MMpMV * FF_fit(q2, a_0A12, a_1A12, a_2A12, MRA12_2)) / lambda(q2);
}

double MVll::T_1(double q2)
{
    return FF_fit(q2, a_0T1, a_1T1, a_2T1, MRT1_2);
}

double MVll::T_2(double q2)
{
    return FF_fit(q2, a_0T2, a_1T2, a_2T2, MRT2_2);
}

double MVll::V_0t(double q2)
{
    return fourMV / sqrt(q2) * FF_fit(q2, a_0A12, a_1A12, a_2A12, MRA12_2);
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
 * QCDF NLO                                                                    *
 * ****************************************************************************/

gslpp::complex MVll::A_Seidel(double q2, double mb2)
{
    double sh = q2 / mb2;
    double z = (4. * mb2) / q2;
    double lsh = log(sh);
    gslpp::complex acsq = arccot((gslpp::complex)sqrt(z - 1.));
    double sh2 = sh*sh;
    double osh2 = (1. - sh)*(1. - sh);
    return (-(104.) / (243.) * log((mb2) / (mu_b2)) + (4. * sh) / (27. * (1. - sh)) * (dilog((gslpp::complex)sh) + lsh * log(1. - sh))
            + (1.) / (729. * osh2) * (6. * sh * (29. - 47. * sh) * lsh + 785. - 1600. * sh + 833. * sh * sh + 6. * M_PI * gslpp::complex::i() * (20. - 49. * sh + 47. * sh2))
            - (2.) / (243. * osh2 * (1. - sh)) * (2. * sqrt(z - 1.) * (-4. + 9. * sh - 15. * sh2 + 4. * sh2 * sh) * acsq + 9. * sh2 * sh * lsh * lsh + 18. * M_PI * gslpp::complex::i() * sh * (1. - 2. * sh) * lsh)
            + (2. * sh) / (243. * osh2 * osh2) * (36. * acsq * acsq + M_PI2 * (-4. + 9. * sh - 9. * sh2 + 3. * sh2 * sh)));
}

gslpp::complex MVll::B_Seidel(double q2, double mb2)
{
    double sh = q2 / mb2;
    double z = (4. * mb2) / q2;
    double sqrt_z_m_1 = sqrt(z - 1.);
    gslpp::complex x1 = 0.5 + gslpp::complex::i() / 2. * sqrt_z_m_1;
    gslpp::complex x2 = 0.5 - gslpp::complex::i() / 2. * sqrt_z_m_1;
    gslpp::complex x3 = 0.5 + gslpp::complex::i() / (2. * sqrt_z_m_1);
    gslpp::complex x4 = 0.5 - gslpp::complex::i() / (2. * sqrt_z_m_1);
    gslpp::complex lx1 = log(x1);
    gslpp::complex lx2 = log(x2);
    gslpp::complex lx3 = log(x3);
    gslpp::complex lx4 = log(x4);
    gslpp::complex lx2_x1 = lx2 - lx1;
    gslpp::complex lzm1 = log(z - 1.);
    gslpp::complex acsq = arccot((gslpp::complex)sqrt_z_m_1);
    double sh2 = sh*sh;
    double lsh = log(sh);
    double osh2 = (1. - sh)*(1. - sh);
    double lmb_mu = log(mb2 / mu_b2);
    return (8. / (243. * sh) * ((4. - 34. * sh - 17. * M_PI * gslpp::complex::i() * sh) * lmb_mu + 8. * sh * lmb_mu * lmb_mu + 17. * sh * lsh * lmb_mu)
            + ((2. + sh) * sqrt_z_m_1) / (729. * sh) * (-48. * lmb_mu * acsq - 18. * M_PI * log(z - 1.) + 3. * gslpp::complex::i() * lzm1 * lzm1
            - 24. * gslpp::complex::i() * dilog(-x2 / x1) - 5. * M_PI2 * gslpp::complex::i()
            + 6. * gslpp::complex::i() * (-9. * lx1 * lx1 + lx2 * lx2 - 2. * lx4 * lx4 + 6. * lx1 * lx2 - 4. * lx1 * lx3 + 8. * lx1 * lx4)
            - 12. * M_PI * (2. * lx1 + lx3 + lx4)) - 2. / (243. * sh * (1 - sh)) * (4. * sh * (-8. + 17. * sh) * (dilog((gslpp::complex)sh) + lsh * log(1. - sh))
            + 3. * (2. + sh) * (3. - sh) * lx2_x1 * lx2_x1 + 12. * M_PI * (-6. - sh + sh2) * acsq) + 2. / (2187. * sh * osh2) * (-18. * sh * (120. - 211. * sh + 73. * sh2) * lsh
            - 288. - 8. * sh + 934. * sh2 - 692. * sh2 * sh + 18. * M_PI * gslpp::complex::i() * sh * (82. - 173. * sh + 73. * sh2))
            - 4. / (243. * sh * osh2 * (1 - sh)) * (-2. * sqrt_z_m_1 * (4. - 3. * sh - 18. * sh2 + 16. * sh2 * sh - 5. * sh2 * sh2) * acsq - 9. * sh * sh2 * lsh * lsh
            + 2. * M_PI * gslpp::complex::i() * sh * (8. - 33. * sh + 51. * sh2 - 17. * sh * sh2) * lsh) + 2. / (729. * sh * osh2 * osh2) * (72. * (3. - 8. * sh + 2. * sh2) * acsq * acsq
            - M_PI2 * (54. - 53. * sh - 286. * sh2 + 612. * sh * sh2 - 446. * sh2 * sh2 + 113. * sh2 * sh2 * sh)));
}

gslpp::complex MVll::C_Seidel(double q2)
{
    return -(16.) / (81.) * log((q2) / (mu_b2)) + (428.) / (243.) - (64.) / (27.) * gsl_sf_zeta_int(3) + (16.) / (81.) * M_PI * gslpp::complex::i();
    /* gsl_sf_zeta_int returns a double */
}

gslpp::complex MVll::deltaC7_QCDF(double q2, bool conjugate, bool spline)
{
#if COMPUTECP && SPLINE
    if (spline && !conjugate) return gsl_spline_eval(spline_Re_deltaC7_QCDF, q2, acc_Re_deltaC7_QCDF);
    else if (spline && conjugate) return gsl_spline_eval(spline_Re_deltaC7_QCDF_conj, q2, acc_Re_deltaC7_QCDF_conj);
#elif SPLINE
    if (spline) return gsl_spline_eval(spline_Re_deltaC7_QCDF, q2, acc_Re_deltaC7_QCDF);
#endif

    double muh = mu_b / mb_pole;
    double z = mc_pole * mc_pole / mb_pole / mb_pole;
    double sh = q2 / mb_pole / mb_pole;
    double sh2 = sh*sh;

#if FULLNLOQCDF    
    gslpp::complex A_Sdl = A_Seidel(q2, mb_pole*mb_pole); /* hep-ph/0403185v2.*/
    gslpp::complex Fu_17 = -A_Sdl; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
    gslpp::complex Fu_27 = 6. * A_Sdl; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
#endif    
    gslpp::complex F_17 = myF_1->F_17re(muh, z, sh, 20) + gslpp::complex::i() * myF_1->F_17im(muh, z, sh, 20); /*arXiv:0810.4077*/
    gslpp::complex F_27 = myF_2->F_27re(muh, z, sh, 20) + gslpp::complex::i() * myF_2->F_27im(muh, z, sh, 20); /*arXiv:0810.4077*/
    gslpp::complex F_87 = F87_0 + F87_1 * sh + F87_2 * sh2 + F87_3 * sh * sh2 - 8. / 9. * log(sh) * (sh + sh2 + sh * sh2);

    if (!conjugate) {
        gslpp::complex delta = C_1 * F_17 + C_2 * F_27;
        gslpp::complex delta_t = C_8 * F_87 + delta;
#if FULLNLOQCDF        
        gslpp::complex delta_u = delta + C_1 * Fu_17 + C_2 * Fu_27;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - lambda_u / lambda_t * delta_u);
#else
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    } else {
        gslpp::complex delta = C_1.conjugate() * F_17 + C_2.conjugate() * F_27;
        gslpp::complex delta_t = C_8.conjugate() * F_87 + delta;
#if FULLNLOQCDF
        gslpp::complex delta_u = delta + C_1.conjugate() * Fu_17 + C_2.conjugate() * Fu_27;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - (lambda_u / lambda_t).conjugate() * delta_u);
#else
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    }
}

gslpp::complex MVll::deltaC9_QCDF(double q2, bool conjugate, bool spline)
{
#if COMPUTECP && SPLINE
    if (spline && !conjugate) return gsl_spline_eval(spline_Re_deltaC9_QCDF, q2, acc_Re_deltaC9_QCDF);
    else if (spline && conjugate) return gsl_spline_eval(spline_Re_deltaC9_QCDF_conj, q2, acc_Re_deltaC9_QCDF_conj);
#elif SPLINE
    if (spline) return gsl_spline_eval(spline_Re_deltaC9_QCDF, q2, acc_Re_deltaC9_QCDF);
#endif

    double muh = mu_b / mb_pole;
    double z = mc_pole * mc_pole / mb_pole / mb_pole;
    double sh = q2 / mb_pole / mb_pole;
    double sh2 = sh*sh;

#if FULLNLOQCDF    
    gslpp::complex B_Sdl = B_Seidel(q2, mb_pole*mb_pole); /* hep-ph/0403185v2.*/
    gslpp::complex C_Sdl = C_Seidel(q2); /* hep-ph/0403185v2.*/
    gslpp::complex Fu_19 = -(B_Sdl + 4. * C_Sdl); /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
    gslpp::complex Fu_29 = -(-6. * B_Sdl + 3. * C_Sdl); /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
#endif    
    gslpp::complex F_19 = myF_1->F_19re(muh, z, sh, 20) + gslpp::complex::i() * myF_1->F_19im(muh, z, sh, 20); /*arXiv:0810.4077*/
    gslpp::complex F_29 = myF_2->F_29re(muh, z, sh, 20) + gslpp::complex::i() * myF_2->F_29im(muh, z, sh, 20); /*arXiv:0810.4077*/
    gslpp::complex F_89 = (F89_0 + F89_1 * sh + F89_2 * sh2 + F89_3 * sh * sh2 + 16. / 9. * log(sh) * (1. + sh + sh2 + sh * sh2));

    if (!conjugate) {
        gslpp::complex delta = C_1 * F_19 + C_2 * F_29;
        gslpp::complex delta_t = C_8 * F_89 + delta;
#if FULLNLOQCDF
        gslpp::complex delta_u = delta + C_1 * Fu_19 + C_2 * Fu_29;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - lambda_u / lambda_t * delta_u);
#else        
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    } else {
        gslpp::complex delta = C_1.conjugate() * F_19 + C_2.conjugate() * F_29;
        gslpp::complex delta_t = C_8.conjugate() * F_89 + delta;
#if FULLNLOQCDF        
        gslpp::complex delta_u = delta + C_1.conjugate() * Fu_19 + C_2.conjugate() * Fu_29;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - (lambda_u / lambda_t).conjugate() * delta_u);
#else        
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    }
}

gslpp::complex MVll::Cq34(bool conjugate)
{
    gslpp::complex T_t = -C_3 + 4. / 3. * (C_4 + 12. * C_5 + 16. * C_6);
    gslpp::complex T_u = 0.; /* 0 for K*0, phi*/
    if (meson == QCD::B_P) T_u = -3. * C_2;
    else if (vectorM == QCD::PHI) T_t = T_t + 6. * (C_3 + 10. * C_5);
    if (!conjugate) return T_t + lambda_u / lambda_t * T_u;
    else return T_t + (lambda_u / lambda_t).conjugate() * T_u;
}

gslpp::complex MVll::T_para_minus_WA(bool conjugate)
{
    return -spectator_charge * 4. * MM / mb_pole * Cq34(conjugate);
}

gslpp::complex MVll::T_perp_WA_1()
{
    return -spectator_charge * 4. / mb_pole * (C_3 + 4. / 3. * (C_4 + 3. * C_5 + 4. * C_6));
}

gslpp::complex MVll::T_perp_WA_2(bool conjugate)
{
    return spectator_charge * 2. / mb_pole * Cq34(conjugate);
}

gslpp::complex MVll::T_perp_plus_O8(double q2, double u)
{
    double ubar = 1. - u;
    double ed = -1. / 3.;

    return -(alpha_s_mub / (3. * M_PI))*4. * ed * C_8 / (u + ubar * q2 / MM2);
}

gslpp::complex MVll::T_para_minus_O8(double q2, double u)
{
    double ubar = 1. - u;

    return (alpha_s_mub / (3. * M_PI))*spectator_charge * 8. * C_8 / (ubar + u * q2 / MM2);
}

gslpp::complex MVll::t_perp(double q2, double u, double m2)
{
    double EV = (MM2 - q2 + MV2) / (2. * MM);
    double ubar = 1. - u;

    return (2. * MM) / (ubar * EV) * I1(q2, u, m2) + q2 / (ubar * ubar * EV * EV) * B0diff(q2, u, m2);

}

gslpp::complex MVll::t_para(double q2, double u, double m2)
{
    double EV = (MM2pMV2 - q2) / (2. * MM);
    double ubar = 1. - u;
    return (2. * MM) / (ubar * EV) * I1(q2, u, m2) + (ubar * MM2 + u * q2) / (ubar * ubar * EV * EV) * B0diff(q2, u, m2);
}

gslpp::complex MVll::I1(double q2, double u, double m2)
{
    if (m2 == 0.) return 1.;

    ubar = 1. - u;
    xp = 0.5 + sqrt(0.25 - ((gslpp::complex) m2) / (ubar * MM2 + u * q2));
    xm = 0.5 - sqrt(0.25 - ((gslpp::complex) m2) / (ubar * MM2 + u * q2));
    yp = 0.5 + sqrt(0.25 - ((gslpp::complex) m2) / q2);
    ym = 0.5 - sqrt(0.25 - ((gslpp::complex) m2) / q2);
    L1xp = log(1. - 1. / xp) * log(1. - xp) - M_PI2osix + dilog(xp / (xp - 1.));
    L1xm = log(1. - 1. / xm) * log(1. - xm) - M_PI2osix + dilog(xm / (xm - 1.));
    L1yp = log(1. - 1. / yp) * log(1. - yp) - M_PI2osix + dilog(yp / (yp - 1.));
    L1ym = log(1. - 1. / ym) * log(1. - ym) - M_PI2osix + dilog(ym / (ym - 1.));

    return 1. + 2. * m2 / ubar / (MM2 - q2)*(L1xp + L1xm - L1yp - L1ym);
}

gslpp::complex MVll::B0diff(double q2, double u, double m2)
{
    double ubar = 1. - u;

    if (m2 == 0.) return -log((gslpp::complex)(-(2. / q2))) + log((gslpp::complex)(-(2. / (q2 * u + MM2 * ubar))));
    else return B0(ubar * MM2 + u * q2, m2) - B0(q2, m2);
}

gslpp::complex MVll::B0(double s, double m2)
{
    if (4. * m2 / s == 1.) return gslpp::complex(0.);
    else return -2. * sqrt(4. * (m2 - gslpp::complex::i()*1.e-10) / s - 1.) * arctan(1. / sqrt(4. * (m2 - gslpp::complex::i()*1.e-10) / s - 1.));
}

gslpp::complex MVll::h_func(double s, double m2)
{
    if (m2 == 0.) return 8. / 27. + 4. * gslpp::complex::i() * M_PI / 9. + 8. * log(mu_b) / 9. - 4. * log(s) / 9.;
    if (s == 0.) return -4. / 9. * (1. + log(m2 / mu_b / mu_b));

    double z = 4 * m2 / s;
    gslpp::complex term;
    if (z > 1) term = atan(1. / sqrt(z - 1.));
    else term = log((1. + sqrt(1. - z)) / sqrt(z)) - ihalfMPI;

    return -4. / 9. * log(m2 / mu_b / mu_b) + 8. / 27. + 4. / 9. * z - 4. / 9. * (2. + z) * sqrt(std::abs(z - 1.)) * term;

}

gslpp::complex MVll::T_perp_plus_QSS(double q2, double u, bool conjugate)
{
    gslpp::complex t_perp_mc = t_perp(q2, u, mc_pole * mc_pole);
    double eu = 0.666666667;
#if FULLNLOQCDF    
    gslpp::complex t_perp_mb = t_perp(q2, u, mb_pole*mb_pole);
    gslpp::complex t_perp_0 = t_perp(q2, u, 0.);
    double ed = -0.333333333;   

    gslpp::complex T_t = (eu * t_perp_mc * (C_1 / 6. + C_2 + 6. * C_6)
            + ed * t_perp_mb * (C_3 - C_4/6. + 16. * C_5 + 10. * C_6/3. + mb_pole / MM * (-C_3 + C_4/6. - 4. * C_5 + 2. * C_6/3.))
            + ed * t_perp_0 * (C_3 - C_4/6. + 16. * C_5 - 8. * C_6/3.));

    gslpp::complex T_u = eu * (t_perp_mc - t_perp_0)*(C_2 - C_1 / 6.);

    if (!conjugate) return alpha_s_mub / (3. * M_PI) * MM / (2. * mb_pole)*(T_t + lambda_u / lambda_t * T_u);
    else return alpha_s_mub / (3. * M_PI) * MM / (2. * mb_pole)*(T_t + (lambda_u / lambda_t).conjugate() * T_u);
#else
    return alpha_s_mub / (3. * M_PI) * MM / (2. * mb_pole)*(eu * t_perp_mc * (C_1 / 6. + C_2 + 6. * C_6));
#endif    
}

gslpp::complex MVll::T_para_plus_QSS(double q2, double u, bool conjugate)
{
    gslpp::complex t_para_mc = t_para(q2, u, mc_pole * mc_pole);
    double eu = 0.666666667;
#if FULLNLOQCDF    
    gslpp::complex t_para_mb = t_para(q2, u, mb_pole*mb_pole);    
    gslpp::complex t_para_0 = t_para(q2, u, 0.);    
    double ed = -0.333333333;

    gslpp::complex T_t = (eu * t_para_mc * (-C_1 / 6. + C_2 + 6. * C_6)
                       + ed * t_para_mb * (C_3 - C_4/6. + 16.*C_5 + 10.*C_6/3.)
                       + ed * t_para_0 * (C_3 - C_4/6. + 16.*C_5 - 8.*C_6/3.));

    gslpp::complex T_u = eu * (t_para_mc - t_para_0) * (C_2 - C_1/6.);

    if (!conjugate) return alpha_s_mub / (3. * M_PI) * MM / mb_pole * (T_t + lambda_u / lambda_t * T_u);
    else return alpha_s_mub / (3. * M_PI) * MM / mb_pole * (T_t + (lambda_u / lambda_t).conjugate() * T_u);
#else
    return alpha_s_mub / (3. * M_PI) * MM / mb_pole * (eu * t_para_mc * (-C_1 / 6. + C_2 + 6. * C_6));
#endif    
}

gslpp::complex MVll::T_para_minus_QSS(double q2, double u, bool conjugate)
{
    double ubar = 1. - u;
    gslpp::complex h_mc = h_func(ubar * MM2 + u*q2, mc_pole * mc_pole);
#if FULLNLOQCDF    
    gslpp::complex h_mb = h_func(ubar*MM2 + u*q2, mb_pole*mb_pole);
    gslpp::complex h_0  = h_func(ubar*MM2 + u*q2, 0);

    gslpp::complex T_t = (h_mc * (-C_1 / 6. + C_2 + C_4 + 10. * C_6)
                       + h_mb * (C_3 + 5.*C_4/6. + 16.*C_5 + 22.*C_6/3.)
                       + h_0 * (C_3 + 17.*C_4/6. + 16.*C_5 + 82.*C_6/3.)
                       - 8./27. * (-15.*C_4/2. + 12.*C_5 - 32.*C_6));

    gslpp::complex T_u = (h_mc - h_0)*(C_2 - C_1/6.); 

    if (!conjugate) return alpha_s_mub / (3. * M_PI) * spectator_charge * 6. * MM / mb_pole * (T_t + lambda_u / lambda_t * T_u);
    else return alpha_s_mub / (3. * M_PI) * spectator_charge * 6. * MM / mb_pole * (T_t + (lambda_u / lambda_t).conjugate() * T_u);
#else
    return alpha_s_mub / (3. * M_PI) * spectator_charge * 6. * MM / mb_pole * (h_mc * (-C_1 / 6. + C_2 + C_4 + 10. * C_6));
#endif    
}

double MVll::phi_V(double u)
{
    return 6. * u * (1. - u) * (1. + mySM.getMesons(vectorM).getGegenalpha(0) * gsl_sf_gegenpoly_1(3. / 2., (2. * u - 1.)) + mySM.getMesons(vectorM).getGegenalpha(1) * gsl_sf_gegenpoly_2(3. / 2., (2. * u - 1.)));
}

gslpp::complex MVll::lambda_B_minus(double q2)
{
    double w0 = mySM.getMesons(meson).getLambdaM();
    return 1. / (exp(-q2 / MM / w0) / w0 * (-gsl_sf_expint_Ei(q2 / MM / w0) + gslpp::complex::i() * M_PI));
}

double MVll::T_perp_real(double q2, double u, bool conjugate)
{
    gslpp::complex T_amp = N_QCDF / mySM.getMesons(meson).getLambdaM() * phi_V(u) * (T_perp_plus_O8(q2, u) + T_perp_plus_QSS(q2, u, conjugate));
#if FULLNLOQCDF 
    double ubar = 1. - u;

    T_amp += N_QCDF/(ubar + u*q2/MM2) * phi_V(u) * T_perp_WA_1() 
           + N_QCDF/mySM.getMesons(meson).getLambdaM() * fpara/fperp * MV/(1. - q2/MM2) * T_perp_WA_2(conjugate);
    /*last term proportional to T_perp_WA_2 is a constant but is included in the integral because u is integrated over the range [0,1]*/
#endif
    return T_amp.real();
}

double MVll::T_perp_imag(double q2, double u, bool conjugate)
{
    gslpp::complex T_amp = N_QCDF / mySM.getMesons(meson).getLambdaM() * phi_V(u) * (T_perp_plus_O8(q2, u) + T_perp_plus_QSS(q2, u, conjugate));
#if FULLNLOQCDF            
    double ubar = 1. - u;
    
    T_amp += N_QCDF/(ubar + u*q2/MM2) * phi_V(u) * T_perp_WA_1() 
           + N_QCDF/mySM.getMesons(meson).getLambdaM() * fpara/fperp * MV/(1. - q2/MM2) * T_perp_WA_2(conjugate);
    /*last term proportional to T_perp_WA_2 is a constant but is included in the integral because u is integrated over the range [0,1]*/
#endif
    return T_amp.imag();
}

double MVll::T_para_real(double q2, double u, bool conjugate)
{
    double N = N_QCDF * (MV / ((MM2pMV2 - q2) / (2. * MM)));

    gslpp::complex T_amp = (N / lambda_B_minus(q2) * (T_para_minus_O8(q2, u) + T_para_minus_QSS(q2, u, conjugate))
            + N / mySM.getMesons(meson).getLambdaM() * T_para_plus_QSS(q2, u, conjugate)) * phi_V(u);
#if FULLNLOQCDF 
    T_amp += N / lambda_B_minus(q2) * T_para_minus_WA(conjugate)* phi_V(u);
#endif    
    return sqrt(q2) * T_amp.real();
}

double MVll::T_para_imag(double q2, double u, bool conjugate)
{
    double N = N_QCDF * (MV / ((MM2pMV2 - q2) / (2. * MM)));

    gslpp::complex T_amp = (N / lambda_B_minus(q2) * (/* + */T_para_minus_O8(q2, u) + T_para_minus_QSS(q2, u, conjugate))
            + N / mySM.getMesons(meson).getLambdaM() * T_para_plus_QSS(q2, u, conjugate)) * phi_V(u);
#if FULLNLOQCDF 
    T_amp += N / lambda_B_minus(q2) * T_para_minus_WA(conjugate) * phi_V(u);
#endif
    return sqrt(q2) * T_amp.imag();
}

double MVll::T_perp_real(double q2, bool conjugate)
{
    FS = convertToGslFunction(boost::bind(&MVll::T_perp_real, &(*this), q2, _1, conjugate));
    gsl_integration_cquad(&FS, 0., 1., 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL);

    return avaSigma;
}

double MVll::T_perp_imag(double q2, bool conjugate)
{
    FS = convertToGslFunction(boost::bind(&MVll::T_perp_imag, &(*this), q2, _1, conjugate));
    gsl_integration_cquad(&FS, 0., 1., 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL);

    return avaSigma;
}

double MVll::T_para_real(double q2, bool conjugate)
{
    FS = convertToGslFunction(boost::bind(&MVll::T_para_real, &(*this), q2, _1, conjugate));
    gsl_integration_cquad(&FS, 0., 1., 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL);

    return avaSigma;
}

double MVll::T_para_imag(double q2, bool conjugate)
{
    FS = convertToGslFunction(boost::bind(&MVll::T_para_imag, &(*this), q2, _1, conjugate));
    gsl_integration_cquad(&FS, 0., 1., 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL);

    return avaSigma;
}

double MVll::QCDF_fit_func(double* x, double* p)
{
    return p[0] + p[1] * x[0] + p[2] * x[0] * x[0] + p[3] * x[0] * x[0] * x[0] + p[4] * x[0] * x[0] * x[0] * x[0] + p[5] * x[0] * x[0] * x[0] * x[0] * x[0] + p[6] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];
}

void MVll::fit_QCDF_func()
{
    int dim = 0;
    for (double i = 0.001; i < 8.6; i += 0.5) {
        myq2.push_back(i);
        Re_T_perp.push_back(T_perp_real(i, false));
        Im_T_perp.push_back(T_perp_imag(i, false));
        Re_T_para.push_back(T_para_real(i, false));
        Im_T_para.push_back(T_para_imag(i, false));

#if COMPUTECP        
        Re_T_perp_conj.push_back(T_perp_real(i, true));
        Im_T_perp_conj.push_back(T_perp_imag(i, true));
        Re_T_para_conj.push_back(T_para_real(i, true));
        Im_T_para_conj.push_back(T_para_imag(i, true));
#endif
        dim++;
    }

    gr1 = TGraph(dim, myq2.data(), Re_T_perp.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Re_T_perp");
    Re_T_perp_res = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Re_T_perp.clear();

    gr1 = TGraph(dim, myq2.data(), Im_T_perp.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Im_T_perp");
    Im_T_perp_res = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Im_T_perp.clear();

    gr1 = TGraph(dim, myq2.data(), Re_T_para.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Re_T_para");
    Re_T_para_res = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Re_T_para.clear();

    gr1 = TGraph(dim, myq2.data(), Im_T_para.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Im_T_para");
    Im_T_para_res = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Im_T_para.clear();

#if COMPUTECP  
    gr1 = TGraph(dim, myq2.data(), Re_T_perp_conj.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Re_T_perp_conj");
    Re_T_perp_res_conj = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Re_T_perp_conj.clear();

    gr1 = TGraph(dim, myq2.data(), Im_T_perp_conj.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Im_T_perp_conj");
    Im_T_perp_res_conj = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Im_T_perp_conj.clear();

    gr1 = TGraph(dim, myq2.data(), Re_T_para_conj.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Re_T_para_conj");
    Re_T_para_res_conj = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Re_T_para_conj.clear();

    gr1 = TGraph(dim, myq2.data(), Im_T_para_conj.data());
    QCDFfit = TF1("QCDFfit", this, &MVll::QCDF_fit_func, 0.001, 8.51, 7, "MVll", "Im_T_para_conj");
    Im_T_para_res_conj = gr1.Fit(&QCDFfit, "SQN0+rob=0.99");
    Im_T_para_conj.clear();
#endif

    myq2.clear();
}

void MVll::spline_QCDF_func()
{
    int dim = GSL_INTERP_DIM;
    int dim_DC = GSL_INTERP_DIM_DC;
    double min = 0.001;
    double interval = (9.9 - min) / ((double) dim);
    double interval_DC = (9.9 - min) / ((double) dim_DC);
    double q2_spline[dim];
    double fq2_Re_T_perp[dim], fq2_Im_T_perp[dim], fq2_Re_T_para[dim], fq2_Im_T_para[dim];
    double q2_spline_DC[dim_DC];
    double fq2_Re_deltaC7_QCDF[dim_DC], fq2_Im_deltaC7_QCDF[dim_DC], fq2_Re_deltaC9_QCDF[dim_DC], fq2_Im_deltaC9_QCDF[dim_DC];
#if COMPUTECP 
    double fq2_Re_T_perp_conj[dim], fq2_Im_T_perp_conj[dim], fq2_Re_T_para_conj[dim], fq2_Im_T_para_conj[dim];
    double fq2_Re_deltaC7_QCDF_conj[dim_DC], fq2_Im_deltaC7_QCDF_conj[dim_DC], fq2_Re_deltaC9_QCDF_conj[dim_DC], fq2_Im_deltaC9_QCDF_conj[dim_DC];
#endif

    for (int i = 0; i < dim; i++) {
        q2_spline[i] = min + (double) i*interval;
        fq2_Re_T_perp[i] = T_perp_real(q2_spline[i], false);
        fq2_Im_T_perp[i] = T_perp_imag(q2_spline[i], false);
        fq2_Re_T_para[i] = T_para_real(q2_spline[i], false);
        fq2_Im_T_para[i] = T_para_imag(q2_spline[i], false);

#if COMPUTECP        
        fq2_Re_T_perp_conj[i] = T_perp_real(q2_spline[i], true);
        fq2_Im_T_perp_conj[i] = T_perp_imag(q2_spline[i], true);
        fq2_Re_T_para_conj[i] = T_para_real(q2_spline[i], true);
        fq2_Im_T_para_conj[i] = T_para_imag(q2_spline[i], true);
#endif
    }
    for (int i = 0; i < dim_DC; i++) {
        q2_spline_DC[i] = min + (double) i*interval_DC;
        fq2_Re_deltaC7_QCDF[i] = deltaC7_QCDF(q2_spline_DC[i], false).real();
        fq2_Im_deltaC7_QCDF[i] = deltaC7_QCDF(q2_spline_DC[i], false).imag();
        fq2_Re_deltaC9_QCDF[i] = deltaC9_QCDF(q2_spline_DC[i], false).real();
        fq2_Im_deltaC9_QCDF[i] = deltaC9_QCDF(q2_spline_DC[i], false).imag();

#if COMPUTECP                
        fq2_Re_deltaC7_QCDF_conj[i] = deltaC7_QCDF(q2_spline_DC[i], true).real();
        fq2_Im_deltaC7_QCDF_conj[i] = deltaC7_QCDF(q2_spline_DC[i], true).imag();
        fq2_Re_deltaC9_QCDF_conj[i] = deltaC9_QCDF(q2_spline_DC[i], true).real();
        fq2_Im_deltaC9_QCDF_conj[i] = deltaC9_QCDF(q2_spline_DC[i], true).imag();
#endif
    }

    gsl_spline_init(spline_Re_T_perp, q2_spline, fq2_Re_T_perp, dim);
    gsl_spline_init(spline_Im_T_perp, q2_spline, fq2_Im_T_perp, dim);
    gsl_spline_init(spline_Re_T_para, q2_spline, fq2_Re_T_para, dim);
    gsl_spline_init(spline_Im_T_para, q2_spline, fq2_Im_T_para, dim);

    gsl_spline_init(spline_Re_deltaC7_QCDF, q2_spline_DC, fq2_Re_deltaC7_QCDF, dim_DC);
    gsl_spline_init(spline_Im_deltaC7_QCDF, q2_spline_DC, fq2_Im_deltaC7_QCDF, dim_DC);
    gsl_spline_init(spline_Re_deltaC9_QCDF, q2_spline_DC, fq2_Re_deltaC9_QCDF, dim_DC);
    gsl_spline_init(spline_Im_deltaC9_QCDF, q2_spline_DC, fq2_Im_deltaC9_QCDF, dim_DC);

#if COMPUTECP
    gsl_spline_init(spline_Re_T_perp_conj, q2_spline, fq2_Re_T_perp_conj, dim);
    gsl_spline_init(spline_Im_T_perp_conj, q2_spline, fq2_Im_T_perp_conj, dim);
    gsl_spline_init(spline_Re_T_para_conj, q2_spline, fq2_Re_T_para_conj, dim);
    gsl_spline_init(spline_Im_T_para_conj, q2_spline, fq2_Im_T_para_conj, dim);

    gsl_spline_init(spline_Re_deltaC7_QCDF_conj, q2_spline_DC, fq2_Re_deltaC7_QCDF_conj, dim_DC);
    gsl_spline_init(spline_Im_deltaC7_QCDF_conj, q2_spline_DC, fq2_Im_deltaC7_QCDF_conj, dim_DC);
    gsl_spline_init(spline_Re_deltaC9_QCDF_conj, q2_spline_DC, fq2_Re_deltaC9_QCDF_conj, dim_DC);
    gsl_spline_init(spline_Im_deltaC9_QCDF_conj, q2_spline_DC, fq2_Im_deltaC9_QCDF_conj, dim_DC);
#endif

}

gslpp::complex MVll::T_minus(double q2, bool conjugate)
{
#if COMPUTECP && SPLINE 
    if (!conjugate) return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (gsl_spline_eval(spline_Re_T_perp, q2, acc_Re_T_perp) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_perp, q2, acc_Im_T_perp));
    else return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (gsl_spline_eval(spline_Re_T_perp_conj, q2, acc_Re_T_perp_conj) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_perp_conj, q2, acc_Im_T_perp_conj));
#elif SPLINE
    return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (gsl_spline_eval(spline_Re_T_perp, q2, acc_Re_T_perp) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_perp, q2, acc_Im_T_perp));
#endif     

#if COMPUTECP  && !SPLINE
    if (!conjugate) return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_perp_res->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_perp_res->GetParams())));
    else return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_perp_res_conj->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_perp_res_conj->GetParams())));
#elif !SPLINE
    return -2. * MM * mb_pole / q2 * (1. - q2 / MM2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_perp_res->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_perp_res->GetParams())));
#endif    
}

gslpp::complex MVll::T_0(double q2, bool conjugate)
{
#if COMPUTECP && SPLINE
    if (!conjugate) return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (gsl_spline_eval(spline_Re_T_para, q2, acc_Re_T_para) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_para, q2, acc_Im_T_para));
    else return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (gsl_spline_eval(spline_Re_T_para_conj, q2, acc_Re_T_para_conj) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_para_conj, q2, acc_Im_T_para_conj));
#elif SPLINE
    return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (gsl_spline_eval(spline_Re_T_para, q2, acc_Re_T_para) + gslpp::complex::i() * gsl_spline_eval(spline_Im_T_para, q2, acc_Im_T_para));
#endif

#if COMPUTECP && !SPLINE
    if (!conjugate) return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_para_res->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_para_res->GetParams())));
    else return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_para_res_conj->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_para_res_conj->GetParams())));
#elif !SPLINE
    return -(1. - q2 / MM2)* (1. - q2 / MM2) * MM * mb_pole / sqrt(q2) * (QCDF_fit_func(&q2, const_cast<double *> (Re_T_para_res->GetParams())) + gslpp::complex::i() * QCDF_fit_func(&q2, const_cast<double *> (Im_T_para_res->GetParams())));
#endif
}

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MVll::H(double q2, double m2, double mu2)
{
    double x = 4. * m2 / q2;
    gslpp::complex par;

    if (x > 1.) par = sqrt(x - 1.) * atan(1. / sqrt(x - 1.));
    else par = sqrt(1. - x) * (log((1. + sqrt(1. - x)) / sqrt(x)) - ihalfMPI);

    return -fournineth * (log(m2 / mu2) - twothird - x) - fournineth * (2. + x) * par;
}

gslpp::complex MVll::H_0(double q2)
{
    return (H_0_pre - fournineth * log(q2 / mu_b2));
}

gslpp::complex MVll::Y(double q2)
{
    return -half * H_0(q2) * H_0_WC + H(q2, mc_pole*mc_pole, mu_b2) * H_c_WC - half * H(q2, mb_pole*mb_pole, mu_b2) * H_b_WC;
}

gslpp::complex MVll::funct_g(double q2)
{
    if (q2 < 4. * Mc * Mc)
        return -8. / 9. * log(Mc / Mb) + 8. / 27. + 16. / 9. * Mc * Mc / q2 - 4. / 9. * (2. + 4. * Mc * Mc / q2) * (sqrt(4. * Mc * Mc / q2 - 1.) * atan(1. / sqrt(4. * Mc * Mc / q2 - 1.)));
    else
        return -8. / 9. * log(Mc / Mb) + 8. / 27. + 16. / 9. * Mc * Mc / q2 - 4. / 9. * (2. + 4. * Mc * Mc / q2) * (sqrt(1. - 4. * Mc * Mc / q2) * (log(1. + sqrt(1. - 4. * Mc * Mc / q2) / sqrt(4. * Mc * Mc / q2)) - gslpp::complex::i() * M_PI_2));
}

gslpp::complex MVll::DeltaC9_KD(double q2, int com)
{
    return ((h_0[com] * (1. - 1. / q2) + h_2[com] / q2) / (1. + h_1[com] * (1. - q2) / mJ2) - (3. * (-0.267) + 1.117) * funct_g(q2))*exp_Phase[com];
    /* C_1 = -0.267 and C_2 = 1.117 in KMPW */
}

gslpp::complex MVll::h_lambda(int hel, double q2)
{
    if (!dispersion) {
        if (h_pole == true) return SU3_breaking * (h_0[hel]+(1. - h_2[hel]) * q2 * (h_1[hel] - h_0[hel]) / (q2 - h_2[hel]));
        else if (hel == 1) return SU3_breaking * (h_0[1] + h_1[1] * q2 + h_2[1] * q2 * q2 + (twoMboMM * h_0[2] * T_p(q2) + h_1[2] * q2 / MM2 * V_p(q2)) / sixteenM_PI2);
        else if (hel == 2) return SU3_breaking * (twoMboMM * h_0[2] * T_m(q2) + h_1[2] * q2 / MM2 * V_m(q2)) / sixteenM_PI2 + h_2[2] * q2 * q2;
        else return SU3_breaking * ((h_0[hel] + h_1[hel] * q2) * sqrt(q2) + (twoMboMM * h_0[2] * T_0t(q2) + h_1[2] * q2 * V_0t(q2) / MM2) / sixteenM_PI2);
    } else {
        if (hel == 0) return SU3_breaking * (-sqrt(q2) / (MM2 * 16. * M_PI * M_PI) * ((MMpMV2 * (MM2mMV2 - q2) * A_1(q2) * DeltaC9_KD(q2, 1) - lambda(q2) * A_2(q2) * DeltaC9_KD(q2, 2)) / (4. * MV * MM * MMpMV)));
        else if (hel == 1) {
            if (q2 == 0.) return SU3_breaking * (-1. / (MM2 * 16. * M_PI * M_PI) * (
                    (MMpMV * A_1(0.)) / (2. * MM) * ((-h_0[1] + h_2[1]) / (1. + h_1[1] / mJ2)) * exp_Phase[1]
                    - sqrt(lambda(0.)) / (2. * MM * MMpMV) * V(0.) * ((-h_0[0] + h_2[0]) / (1. + h_1[0] / mJ2)) * exp_Phase[1]));
            else return SU3_breaking * (-q2 / (MM2 * 16. * M_PI * M_PI) * ((MMpMV * A_1(q2)) / (2. * MM) * DeltaC9_KD(q2, 1) - sqrt(lambda(q2)) / (2. * MM * MMpMV) * V(q2) * DeltaC9_KD(q2, 0)));
        } else {
            if (q2 == 0.) return SU3_breaking * (-1. / (MM2 * 16. * M_PI * M_PI) *
                    ((MMpMV * A_1(0.)) / (2. * MM) * ((-h_0[1] + h_2[1]) / (1. + h_1[1] / mJ2)) * exp_Phase[1]
                    + sqrt(lambda(0.)) / (2. * MM * MMpMV) * V(0.) * ((-h_0[0] + h_2[0]) / (1. + h_1[0] / mJ2)) * exp_Phase[1]));
            else return SU3_breaking * (-q2 / (MM2 * 16. * M_PI * M_PI) * ((MMpMV * A_1(q2)) / (2. * MM) * DeltaC9_KD(q2, 1) + sqrt(lambda(q2)) / (2. * MM * MMpMV) * V(q2) * DeltaC9_KD(q2, 0)));
        }
    }
}

gslpp::complex MVll::H_V_0(double q2, bool bar)
{
    if (!bar) return -gslpp::complex::i() * NN * (((C_9 + deltaC9_QCDF(q2, !bar, SPLINE) /*+ fDeltaC9_0(q2)*/ + Y(q2)) - etaV * pow(-1, angmomV) * C_9p) * V_0t(q2) + T_0(q2, !bar) + MM2 / q2 * (twoMboMM * (C_7 + deltaC7_QCDF(q2, !bar, SPLINE) - etaV * pow(-1, angmomV) * C_7p) * T_0t(q2) - sixteenM_PI2 * h_lambda(0, q2)));
    return -gslpp::complex::i() * NN_conjugate * (((C_9.conjugate() + deltaC9_QCDF(q2, bar, SPLINE) /*+ fDeltaC9_0(q2)*/ + Y(q2)) - etaV * pow(-1, angmomV) * C_9p.conjugate()) * V_0t(q2) + T_0(q2, bar) + MM2 / q2 * (twoMboMM * (C_7.conjugate() + deltaC7_QCDF(q2, bar, SPLINE) - etaV * pow(-1, angmomV) * C_7p.conjugate()) * T_0t(q2) - sixteenM_PI2 * h_lambda(0, q2)));

}

gslpp::complex MVll::H_V_p(double q2, bool bar)
{
    if (!bar) return -gslpp::complex::i() * NN * (((C_9 + deltaC9_QCDF(q2, !bar, SPLINE) /*+ fDeltaC9_p(q2)*/ + Y(q2)) * V_p(q2) - etaV * pow(-1, angmomV) * C_9p * V_m(q2)) + MM2 / q2 * (twoMboMM * ((C_7 + deltaC7_QCDF(q2, !bar, SPLINE)) * T_p(q2) - etaV * pow(-1, angmomV) * C_7p * T_m(q2)) - sixteenM_PI2 * h_lambda(1, q2)));
    return -gslpp::complex::i() * NN_conjugate * (((C_9.conjugate() + deltaC9_QCDF(q2, bar, SPLINE) /*+ fDeltaC9_p(q2)*/ + Y(q2)) * V_p(q2) - etaV * pow(-1, angmomV) * C_9p.conjugate() * V_m(q2)) + MM2 / q2 * (twoMboMM * ((C_7.conjugate() + deltaC7_QCDF(q2, bar, SPLINE)) * T_p(q2) - etaV * pow(-1, angmomV) * C_7p.conjugate() * T_m(q2)) - sixteenM_PI2 * h_lambda(1, q2)));
}

gslpp::complex MVll::H_V_m(double q2, bool bar)
{
    if (!bar) return -gslpp::complex::i() * NN * (((C_9 + deltaC9_QCDF(q2, !bar, SPLINE) /*+ fDeltaC9_m(q2)*/ + Y(q2)) * V_m(q2) - etaV * pow(-1, angmomV) * C_9p * V_p(q2)) + T_minus(q2, !bar) + MM2 / q2 * (twoMboMM * ((C_7 + deltaC7_QCDF(q2, !bar, SPLINE)) * T_m(q2) - etaV * pow(-1, angmomV) * C_7p * T_p(q2)) - sixteenM_PI2 * h_lambda(2, q2)));
    return -gslpp::complex::i() * NN_conjugate * (((C_9.conjugate() + deltaC9_QCDF(q2, bar, SPLINE) /*+ fDeltaC9_m(q2)*/ + Y(q2)) * V_m(q2) - etaV * pow(-1, angmomV) * C_9p.conjugate() * V_p(q2)) + T_minus(q2, bar) + MM2 / q2 * (twoMboMM * ((C_7.conjugate() + deltaC7_QCDF(q2, bar, SPLINE)) * T_m(q2) - etaV * pow(-1, angmomV) * C_7p.conjugate() * T_p(q2)) - sixteenM_PI2 * h_lambda(2, q2)));
}

gslpp::complex MVll::H_A_0(double q2, bool bar)
{
    if (!bar) return gslpp::complex::i() * NN * (-C_10 + etaV * pow(-1, angmomV) * C_10p) * V_0t(q2);
    return gslpp::complex::i() * NN_conjugate * (-C_10.conjugate() + etaV * pow(-1, angmomV) * C_10p.conjugate()) * V_0t(q2);
}

gslpp::complex MVll::H_A_p(double q2, bool bar)
{
    if (!bar) return gslpp::complex::i() * NN * (-C_10 * V_p(q2) + etaV * pow(-1, angmomV) * C_10p * V_m(q2));
    return gslpp::complex::i() * NN_conjugate * (-C_10.conjugate() * V_p(q2) + etaV * pow(-1, angmomV) * C_10p.conjugate() * V_m(q2));
}

gslpp::complex MVll::H_A_m(double q2, bool bar)
{
    if (!bar) return gslpp::complex::i() * NN * (-C_10 * V_m(q2) + etaV * pow(-1, angmomV) * C_10p * V_p(q2));
    return gslpp::complex::i() * NN_conjugate * (-C_10.conjugate() * V_m(q2) + etaV * pow(-1, angmomV) * C_10p.conjugate() * V_p(q2));
}

gslpp::complex MVll::H_S(double q2, bool bar)
{
    if (!bar) return gslpp::complex::i() * NN * MboMW * (C_S - etaV * pow(-1, angmomV) * C_Sp) * S_L(q2);
    return gslpp::complex::i() * NN_conjugate * MboMW * (C_S.conjugate() - etaV * pow(-1, angmomV) * C_Sp.conjugate()) * S_L(q2);
}

gslpp::complex MVll::H_P(double q2, bool bar)
{
    if (!bar) return gslpp::complex::i() * NN * (MboMW * (C_P - etaV * pow(-1, angmomV) * C_Pp) + twoMlepMb / q2 * (C_10 * (1. + etaV * pow(-1, angmomV) * MsoMb) - C_10p * (etaV * pow(-1, angmomV) + MsoMb))) * S_L(q2);
    return gslpp::complex::i() * NN_conjugate * (MboMW * (C_P.conjugate() - etaV * pow(-1, angmomV) * C_Pp.conjugate()) + twoMlepMb / q2 * (C_10.conjugate()*(1. + etaV * pow(-1, angmomV) * MsoMb) - C_10p.conjugate()*(etaV * pow(-1, angmomV) + MsoMb))) * S_L(q2);
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

double MVll::I_1c(double q2, bool bar)
{
    return F(q2, b)*((H_V_0(q2, bar).abs2() + H_A_0(q2, bar).abs2()) / 2. + H_P(q2, bar).abs2() + 2. * Mlep2 / q2 * (H_V_0(q2, bar).abs2()
            - H_A_0(q2, bar).abs2()) + beta2(q2) * H_S(q2, bar).abs2());
}

double MVll::I_1s(double q2, bool bar)
{
    return F(q2, b)*((beta2(q2) + 2.) / 8. * (H_V_p(q2, bar).abs2() + H_V_m(q2, bar).abs2() + H_A_p(q2, bar).abs2() + H_A_m(q2, bar).abs2()) +
            Mlep2 / q2 * (H_V_p(q2, bar).abs2() + H_V_m(q2, bar).abs2() - H_A_p(q2, bar).abs2() - H_A_m(q2, bar).abs2()));
}

double MVll::I_2c(double q2, bool bar)
{
    return -F(q2, b) * beta2(q2) / 2. * (H_V_0(q2, bar).abs2() + H_A_0(q2, bar).abs2());
}

double MVll::I_2s(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 8. * (H_V_p(q2, bar).abs2() + H_V_m(q2, bar).abs2() + H_A_p(q2, bar).abs2() + H_A_m(q2, bar).abs2());
}

double MVll::I_3(double q2, bool bar)
{
    return -F(q2, b) * beta2(q2) / 2. * ((H_V_p(q2, bar) * H_V_m(q2, bar).conjugate()).real() + (H_A_p(q2, bar) * H_A_m(q2, bar).conjugate()).real());
}

double MVll::I_4(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 4. * (((H_V_m(q2, bar) + H_V_p(q2, bar)) * H_V_0(q2, bar).conjugate()).real() + ((H_A_m(q2, bar) + H_A_p(q2, bar)) * H_A_0(q2, bar).conjugate()).real());
}

double MVll::I_5(double q2, bool bar)
{
    return F(q2, b)*(beta(q2) / 2. * (((H_V_m(q2, bar) - H_V_p(q2, bar)) * H_A_0(q2, bar).conjugate()).real() + ((H_A_m(q2, bar) - H_A_p(q2, bar)) * H_V_0(q2, bar).conjugate()).real()) -
            beta(q2) * Mlep / sqrt(q2)*(H_S(q2, bar).conjugate()*(H_V_p(q2, bar) + H_V_m(q2, bar))).real());
}

double MVll::I_6s(double q2, bool bar)
{
    return F(q2, b) * beta(q2)*(H_V_m(q2, bar)*(H_A_m(q2, bar).conjugate()) - H_V_p(q2, bar)*(H_A_p(q2, bar).conjugate())).real();
}

double MVll::I_6c(double q2, bool bar)
{
    return 4. * F(q2, b) * beta(q2) * Mlep / sqrt(q2)*(H_S(q2, bar).conjugate() * H_V_0(q2, bar)).real();
}

double MVll::I_7(double q2, bool bar)
{
    return F(q2, b)*(beta(q2) / 2. * (((H_V_m(q2, bar) + H_V_p(q2, bar)) * H_A_0(q2, bar).conjugate()).imag() + ((H_A_m(q2, bar) + H_A_p(q2, bar)) * H_V_0(q2, bar).conjugate()).imag()) -
            beta(q2) * Mlep / sqrt(q2)*(H_S(q2, bar).conjugate()*(H_V_m(q2, bar) - H_V_p(q2, bar))).imag());
}

double MVll::I_8(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 4. * (((H_V_m(q2, bar) - H_V_p(q2, bar)) * H_V_0(q2, bar).conjugate()).imag() + ((H_A_m(q2, bar) - H_A_p(q2, bar)) * H_A_0(q2, bar).conjugate()).imag());
}

double MVll::I_9(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 2. * ((H_V_p(q2, bar) * H_V_m(q2, bar).conjugate()).imag() + (H_A_p(q2, bar) * H_A_m(q2, bar).conjugate()).imag());
}

double MVll::h_1c(double q2, bool bar)
{
    return F(q2, b)*((H_V_0(q2, bar).abs2() + H_A_0(q2, bar).abs2()) + 2. * H_P(q2, bar).abs2() + 4. * Mlep2 / q2 * (H_V_0(q2, bar).abs2()
            - H_A_0(q2, bar).abs2()) - 2. * beta2(q2) * H_S(q2, bar).abs2());
}

double MVll::h_1s(double q2, bool bar)
{
    return F(q2, b)*((beta2(q2) + 2.) / 2. * ((H_V_p(q2, bar) * H_V_m(q2, bar).conjugate()).real()
            + (H_A_p(q2, bar) * H_A_m(q2, bar).conjugate()).real()) +
            4. * Mlep2 / q2 * ((H_V_p(q2, bar) * H_V_m(q2, bar).conjugate()).real()
            - (H_A_p(q2, bar) * H_A_m(q2, bar).conjugate()).real()));
}

double MVll::h_2c(double q2, bool bar)
{
    return -F(q2, b) * beta2(q2) * (H_V_0(q2, bar).abs2() + H_A_0(q2, bar).abs2());
}

double MVll::h_2s(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 2. * ((H_V_p(q2, bar) * H_V_m(q2, bar).conjugate()).real()
            + (H_A_p(q2, bar) * H_A_m(q2, bar).conjugate()).real());
}

double MVll::h_3(double q2, bool bar)
{
    return -F(q2, b) * beta2(q2) / 2. * (H_V_p(q2, bar).abs2() + H_V_m(q2, bar).abs2()
            + H_A_p(q2, bar).abs2() + H_A_m(q2, bar).abs2());
}

double MVll::h_4(double q2, bool bar)
{
    return F(q2, b) * beta2(q2) / 2. * (((H_V_m(q2, bar) + H_V_p(q2, bar)) * H_V_0(q2, bar).conjugate()).real()
            + ((H_A_m(q2, bar) + H_A_p(q2, bar)) * H_A_0(q2, bar).conjugate()).real());
}

double MVll::h_7(double q2, bool bar)
{
    return F(q2, b)*(beta(q2) * (((H_V_m(q2, bar) + H_V_p(q2, bar)) * H_A_0(q2, bar).conjugate()).imag()
            + ((H_A_m(q2, bar) + H_A_p(q2, bar)) * H_V_0(q2, bar).conjugate()).imag()) -
            beta(q2) * 2. * Mlep / sqrt(q2)*(H_S(q2, bar).conjugate()*(H_V_m(q2, bar) - H_V_p(q2, bar))).imag());
}

double MVll::s_5(double q2, bool bar)
{
    return beta(q2) * (2. * Mlep * ((H_V_m(q2, bar) + H_V_p(q2, bar)) * F(q2, b) * H_S(q2, bar).conjugate()).imag() / sqrt(q2)
            - F(q2, b)*((H_A_m(q2, bar) - H_A_p(q2, bar)) * H_V_0(q2, bar).conjugate()
            + (H_V_m(q2, bar) - H_V_p(q2, bar)) * H_A_0(q2, bar).conjugate()).imag());
}

double MVll::s_6s(double q2, bool bar)
{
    return 2. * beta(q2) * F(q2, b) * (H_A_p(q2, bar) * H_V_m(q2, bar).conjugate() + H_V_p(q2, bar) * H_A_m(q2, bar).conjugate()).imag();
}

double MVll::s_6c(double q2, bool bar)
{
    return -8. * beta(q2) * Mlep * (H_V_0(q2, bar) * F(q2, b) * H_S(q2, bar).conjugate()).imag() / sqrt(q2);
}

double MVll::s_8(double q2, bool bar)
{
    return beta2(q2) * F(q2, b) * ((H_A_m(q2, bar) - H_A_p(q2, bar)) * H_A_0(q2, bar).conjugate()
            + (H_V_m(q2, bar) - H_V_p(q2, bar)) * H_V_0(q2, bar).conjugate()).real() / 2.;
}

double MVll::s_9(double q2, bool bar)
{
    return beta2(q2) * F(q2, b) * (H_A_p(q2, bar).abs2() - H_A_m(q2, bar).abs2()
            + H_V_p(q2, bar).abs2() - H_V_m(q2, bar).abs2()) / 2.;
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
                cacheSigma0[qbin] = avaSigma;
                sigma0Cached[qbin] = 1;
            }
            return cacheSigma0[qbin];
            break;
        case 1:
            if (sigma1Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma1s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma1[qbin] = avaSigma;
                sigma1Cached[qbin] = 1;
            }
            return cacheSigma1[qbin];
            break;
        case 2:
            if (sigma2Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma2c, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma2[qbin] = avaSigma;
                sigma2Cached[qbin] = 1;
            }
            return cacheSigma2[qbin];
            break;
        case 3:
            if (sigma3Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma2s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma3[qbin] = avaSigma;
                sigma3Cached[qbin] = 1;
            }
            return cacheSigma3[qbin];
            break;
        case 4:
            if (sigma4Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma3, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma4[qbin] = avaSigma;
                sigma4Cached[qbin] = 1;
            }
            return cacheSigma4[qbin];
            break;
        case 5:
            if (sigma5Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma4, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma5[qbin] = avaSigma;
                sigma5Cached[qbin] = 1;
            }
            return cacheSigma5[qbin];
            break;
        case 6:
            if (sigma6Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma5, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma6[qbin] = avaSigma;
                sigma6Cached[qbin] = 1;
            }
            return cacheSigma6[qbin];
            break;
        case 7:
            if (sigma7Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma6s, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma7[qbin] = avaSigma;
                sigma7Cached[qbin] = 1;
            }
            return cacheSigma7[qbin];
            break;
        case 9:
            if (sigma9Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma7, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma9[qbin] = avaSigma;
                sigma9Cached[qbin] = 1;
            }
            return cacheSigma9[qbin];
            break;
        case 10:
            if (sigma10Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma8, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma10[qbin] = avaSigma;
                sigma10Cached[qbin] = 1;
            }
            return cacheSigma10[qbin];
            break;
        case 11:
            if (sigma11Cached[qbin] == 0) {
                FS = convertToGslFunction(boost::bind(&MVll::getSigma9, &(*this), _1));
                if (gsl_integration_cquad(&FS, q_min, q_max, 1.e-2, 1.e-1, w_sigma, &avaSigma, &errSigma, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheSigma11[qbin] = avaSigma;
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
            throw std::runtime_error("MVll::getSigma: index " + out.str() + " not implemented");
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
                FD = convertToGslFunction(boost::bind(&MVll::getDelta1c, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta0[qbin] = avaDelta;
                delta0Cached[qbin] = 1;
            }
            return cacheDelta0[qbin];
            break;
        case 1:
            if (delta1Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta1s, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta1[qbin] = avaDelta;
                delta1Cached[qbin] = 1;
            }
            return cacheDelta1[qbin];
            break;
        case 2:
            if (delta2Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta2c, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta2[qbin] = avaDelta;
                delta2Cached[qbin] = 1;
            }
            return cacheDelta2[qbin];
            break;
        case 3:
            if (delta3Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta2s, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta3[qbin] = avaDelta;
                delta3Cached[qbin] = 1;
            }
            return cacheDelta3[qbin];
            break;
        case 6:
            if (delta6Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta5, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta6[qbin] = avaDelta;
                delta6Cached[qbin] = 1;
            }
            return cacheDelta6[qbin];
            break;
        case 7:
            if (delta7Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta6s, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta7[qbin] = avaDelta;
                delta7Cached[qbin] = 1;
            }
            return cacheDelta7[qbin];
            break;
        case 8:
            if (delta8Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta6c, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta8[qbin] = avaDelta;
                delta8Cached[qbin] = 1;
            }
            return cacheDelta8[qbin];
            break;
        case 10:
            if (delta10Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta8, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta10[qbin] = avaDelta;
                delta10Cached[qbin] = 1;
            }
            return cacheDelta10[qbin];
            break;
        case 11:
            if (delta11Cached[qbin] == 0) {
                FD = convertToGslFunction(boost::bind(&MVll::getDelta9, &(*this), _1));
                if (gsl_integration_cquad(&FD, q_min, q_max, 1.e-2, 1.e-1, w_delta, &avaDelta, &errDelta, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
                cacheDelta11[qbin] = avaDelta;
                delta11Cached[qbin] = 1;
            }
            return cacheDelta11[qbin];
            break;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVll::integrateDelta: index " + out.str() + " not implemented");
    }

    gsl_set_error_handler(old_handler);

}
