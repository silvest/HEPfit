/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVgamma.h"
#include "StandardModel.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include "F_1.h"
#include "F_2.h"
#include "AmpDB2.h"
#include <boost/bind.hpp>
#include <limits>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_gegenbauer.h>

MVgamma::MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: SM(SM_i), myF_1(new F_1()), myF_2(new F_2())
{
    meson = meson_i;
    vectorM = vector_i;
    dispersion = false;
    FixedWCbtos = false;
    mJ2 = 3.096*3.096;
    
    w_GSL = gsl_integration_cquad_workspace_alloc (100);
}

MVgamma::~MVgamma()
{}

std::vector<std::string> MVgamma::initializeMVgammaParameters()
{
    dispersion = SM.getFlavour().getFlagUseDispersionRelation();
    FixedWCbtos = SM.getFlavour().getFlagFixedWCbtos();
    
#if NFPOLARBASIS_MVGAMMA
    if (vectorM == StandardModel::PHI) mVgammaParameters = make_vector<std::string>() << 
            "a_0T1phi" << "a_0A1phi" << "a_0Vphi" << 
            "absh_p" << "absh_m" << "argh_p" << "argh_m" << "SU3_breaking_abs" << "SU3_breaking_arg";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1" << "a_0A1" << "a_0V" << 
            "absh_p" << "absh_m" << "argh_p" << "argh_m";
    else if (vectorM == StandardModel::RHO || vectorM == StandardModel::RHO_P) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1rho" << "a_0A1rho" << "a_0Vrho" << 
            "absh_p" << "absh_m" << "argh_p" << "argh_m";
    else if (vectorM == StandardModel::OMEGA) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1omega" << "a_0A1omega" << "a_0Vomega" << 
            "absh_p" << "absh_m" << "argh_p" << "argh_m";
#else
    if (vectorM == StandardModel::PHI) mVgammaParameters = make_vector<std::string>() << 
            "a_0T1phi" << "a_0A1phi" << "a_0Vphi" << 
            "reh_p" << "reh_m" << "imh_p" << "imh_m" << "SU3_breaking_abs" << "SU3_breaking_arg";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1" << "a_0A1" << "a_0V" << 
            "reh_p" << "reh_m" << "imh_p" << "imh_m";
    else if (vectorM == StandardModel::RHO || vectorM == StandardModel::RHO_P) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1rho" << "a_0A1rho" << "a_0Vrho" << 
            "reh_p" << "reh_m" << "imh_p" << "imh_m";
    else if (vectorM == StandardModel::OMEGA) 
        mVgammaParameters = make_vector<std::string>() << "a_0T1omega" << "a_0A1omega" << "a_0Vomega" << 
            "reh_p" << "reh_m" << "imh_p" << "imh_m";
#endif
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    if (dispersion) {
        mVgammaParameters.clear();
        if (vectorM == StandardModel::PHI) mVgammaParameters = make_vector<std::string>() << 
                "a_0T1phi" << "a_0A1phi" << "a_0Vphi" << 
                "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1" << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2" << "SU3_breaking_abs" << "SU3_breaking_arg";
        else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) 
            mVgammaParameters = make_vector<std::string>() << "a_0T1" << "a_0A1" << "a_0V" <<
                "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1" << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2";
        else if (vectorM == StandardModel::RHO || vectorM == StandardModel::RHO_P) 
            mVgammaParameters = make_vector<std::string>() << "a_0T1rho" << "a_0A1rho" << "a_0Vrho" << 
                "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1" << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2";
        else if (vectorM == StandardModel::OMEGA) 
            mVgammaParameters = make_vector<std::string>() << "a_0T1omega" << "a_0A1omega" << "a_0Vomega" << 
                "r1_1" << "r2_1" << "deltaC9_1" << "phDC9_1" << "r1_2" << "r2_2" << "deltaC9_2" << "phDC9_2";
        else {
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
        }
    }
    
    if (FixedWCbtos) mVgammaParameters.push_back("C7_SM" );
    
    SM.initializeMeson(meson);
    SM.initializeMeson(vectorM);
    return mVgammaParameters;
}

void MVgamma::updateParameters()
{
    if (!SM.getFlavour().getUpdateFlag(meson, vectorM, QCD::NOLEPTON)) return;
    
    GF = SM.getGF();
    ale = SM.getAle();
    MM = SM.getMesons(meson).getMass();
    MM2 = MM * MM;
    MV = SM.getMesons(vectorM).getMass();
    Mb = SM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    mb_pole = SM.Mbar2Mp(Mb); /* Conversion to pole mass*/
    mc_pole = SM.Mbar2Mp(SM.getQuarks(QCD::CHARM).getMass()); /* Conversion to pole mass*/
    Ms = SM.getQuarks(QCD::STRANGE).getMass();
    MW = SM.Mw();
    mu_b = SM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    fB = SM.getMesons(meson).getDecayconst();
    width = SM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV, 2.);
    alpha_s_mub = SM.Als(mu_b, FULLNLO); /* Used for QCDF @ NLO */
    
    switch (vectorM) {
        case StandardModel::K_star:
            a_0T1 = SM.getOptionalParameter("a_0T1");
            a_0A1 = SM.getOptionalParameter("a_0A1");
            a_0V = SM.getOptionalParameter("a_0V");
            lambda_t = SM.getCKM().computelamt_s();
            lambda_u = SM.getCKM().computelamu_s();
            spectator_charge = SM.getQuarks(QCD::DOWN).getCharge();
            SU3_breaking = 1.;
            break;
        case StandardModel::K_star_P:
            a_0T1 = SM.getOptionalParameter("a_0T1");
            a_0A1 = SM.getOptionalParameter("a_0A1");
            a_0V = SM.getOptionalParameter("a_0V");
            lambda_t = SM.getCKM().computelamt_s();
            lambda_u = SM.getCKM().computelamu_s();
            spectator_charge = SM.getQuarks(QCD::UP).getCharge();
            SU3_breaking = 1.;
            break;
        case StandardModel::PHI:
            a_0T1 = SM.getOptionalParameter("a_0T1phi");
            a_0A1 = SM.getOptionalParameter("a_0A1phi");
            a_0V = SM.getOptionalParameter("a_0Vphi");
            lambda_t = SM.getCKM().computelamt_s();
            lambda_u = SM.getCKM().computelamu_s();
            spectator_charge = SM.getQuarks(QCD::STRANGE).getCharge();
            SU3_breaking = gslpp::complex(1. + SM.getOptionalParameter("SU3_breaking_abs"),
                    SM.getOptionalParameter("SU3_breaking_arg"), true);
            break;
        case StandardModel::RHO:
            a_0T1 = SM.getOptionalParameter("a_0T1rho");
            a_0A1 = SM.getOptionalParameter("a_0A1rho");
            a_0V = SM.getOptionalParameter("a_0Vrho");
            lambda_t = SM.getCKM().computelamt_d();
            lambda_u = SM.getCKM().computelamu_d();
            spectator_charge = SM.getQuarks(QCD::DOWN).getCharge();
            SU3_breaking = 1.;
            break;
        case StandardModel::RHO_P:
            a_0T1 = SM.getOptionalParameter("a_0T1rho");
            a_0A1 = SM.getOptionalParameter("a_0A1rho");
            a_0V = SM.getOptionalParameter("a_0Vrho");
            lambda_t = SM.getCKM().computelamt_d();
            lambda_u = SM.getCKM().computelamu_d();
            spectator_charge = SM.getQuarks(QCD::UP).getCharge();
            SU3_breaking = 1.;
            break;
        case StandardModel::OMEGA:
            a_0T1 = SM.getOptionalParameter("a_0T1omega");
            a_0A1 = SM.getOptionalParameter("a_0A1omega");
            a_0V = SM.getOptionalParameter("a_0Vomega");
            lambda_t = SM.getCKM().computelamt_d();
            lambda_u = SM.getCKM().computelamu_d();
            spectator_charge = SM.getQuarks(QCD::DOWN).getCharge();
            SU3_breaking = 1.;
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    fpara = SM.getMesons(vectorM).getDecayconst();
    fperp = SM.getMesons(vectorM).getDecayconst_p();
    
    double ms_over_mb = SM.Mrun(mu_b, SM.getQuarks(QCD::STRANGE).getMass_scale(), 
                        SM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO)
                       /SM.Mrun(mu_b, SM.getQuarks(QCD::BOTTOM).getMass_scale(), 
                        SM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    
    if (!dispersion) {
#if NFPOLARBASIS_MVGAMMA
        h[0] = gslpp::complex(SM.getOptionalParameter("absh_p"), SM.getOptionalParameter("argh_p"), true); //h_plus
        h[1] = gslpp::complex(SM.getOptionalParameter("absh_m"), SM.getOptionalParameter("argh_m"), true); //h_minus
        h[1] *= 2. * (Mb / MM) / (16. * M_PI * M_PI) * (T_1() * lambda / MM2) ;
        h[0] += ms_over_mb * h[1] ;
        
        r1_1 = 0.;
        r1_2 = 0.;
        r2_1 = 0.;
        r2_2 = 0.;
        deltaC9_1 = 0.;
        deltaC9_2 = 0.;
        exp_Phase_1 = 0.;
        exp_Phase_2 = 0.;
#else
        h[0] = gslpp::complex(SM.getOptionalParameter("reh_p"), SM.getOptionalParameter("imh_p"), false); //h_plus
        h[1] = gslpp::complex(SM.getOptionalParameter("reh_m"), SM.getOptionalParameter("imh_m"), false); //h_minus
        h[1] *= 2. * (Mb / MM) / (16. * M_PI * M_PI) * (T_1() * lambda / MM2) ;
        h[0] += ms_over_mb * h[1] ;
        
        r1_1 = 0.;
        r1_2 = 0.;
        r2_1 = 0.;
        r2_2 = 0.;
        deltaC9_1 = 0.;
        deltaC9_2 = 0.;
        exp_Phase_1 = 0.;
        exp_Phase_2 = 0.;
#endif
    } else {
        //gslpp::complex DC7_1 = SM.getOptionalParameter("deltaC7_1")*exp(gslpp::complex::i()*SM.getOptionalParameter("phDC7_1"));
        //gslpp::complex DC7_2 = SM.getOptionalParameter("deltaC7_2")*exp(gslpp::complex::i()*SM.getOptionalParameter("phDC7_2"));
        //h[0] = (-(2.*Mb)/(MM*16.*M_PI*M_PI) * lambda/(2.*MM2) * T_1()*(DC7_2 - DC7_1)).abs();
        //h[1] = (-(2.*Mb)/(MM*16.*M_PI*M_PI) * lambda/(2.*MM2) * T_1()*(DC7_2 + DC7_1)).abs();
        r1_1 = SM.getOptionalParameter("r1_1");
        r1_2 = SM.getOptionalParameter("r1_2");
        r2_1 = SM.getOptionalParameter("r2_1");
        r2_2 = SM.getOptionalParameter("r2_2");
        deltaC9_1 = SM.getOptionalParameter("deltaC9_1");
        deltaC9_2 = SM.getOptionalParameter("deltaC9_2");
        exp_Phase_1 = exp(gslpp::complex::i()*SM.getOptionalParameter("phDC9_1"));
        exp_Phase_2 = exp(gslpp::complex::i()*SM.getOptionalParameter("phDC9_2"));
        
        h[0] = h_lambda(0);
        h[1] = h_lambda(1);
    }

#if UNIFIEDBTOS       
    allcoeff = SM.getFlavour().ComputeCoeffBMll(mu_b, QCD::MU); //check the mass scale, scheme fixed to NDR. QCD::MU does not make any difference to the WC necessary here.
    allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu_b, QCD::MU); //check the mass scale, scheme fixed to NDR. QCD::MU does not make any difference to the WC necessary here.

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_2 = (*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1);
    C_3 = (*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2);
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 = (*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4);
    C_6 = (*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5);
    C_8 = (*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7);
    
    if (FixedWCbtos) {/** NOTE: ComputeCoeff with different argumetns cannot be mixed. They have to be called sequentially. **/
        allcoeff_noSM = SM.getFlavour().ComputeCoeffBMll(mu_b, StandardModel::NOLEPTON, true); //check the mass scale, scheme fixed to NDR
        C_7 = SM.getOptionalParameter("C7_SM") + ((*(allcoeff_noSM[LO]))(6) + (*(allcoeff_noSM[NLO]))(6));
    }
    else C_7 = ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6));
    C_7p = ms_over_mb * ((*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6));
//    C_7p -= ms_over_mb * C_7;
    /* Done in the dirty way to remove from the effective basis since the effective C7p does not involve the non-primed C_1 to C_6.*/
    C_7p += ms_over_mb * (-C_7 - 1. / 3. * C_3 - 4 / 9 * C_4 - 20. / 3. * C_5 - 80. / 9. * C_6);
#else   
    allcoeff = SM.getFlavour().ComputeCoeffsgamma(mu_b);
    allcoeffprime = SM.getFlavour().ComputeCoeffprimesgamma(mu_b);

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_2 = (*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1);
    C_3 = (*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2);
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 = (*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4);
    C_6 = (*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5);
    C_8 = (*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7);

    if (FixedWCbtos) {/** NOTE: ComputeCoeff with different argumetns cannot be mixed. They have to be called sequentially. **/
        allcoeff_noSM = SM.getFlavour().ComputeCoeffsgamma(mu_b, true); //check the mass scale, scheme fixed to NDR
        C_7 = SM.getOptionalParameter("C7_SM") + ((*(allcoeff_noSM[LO]))(6) + (*(allcoeff_noSM[NLO]))(6));
    }
    else C_7 = ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6)););
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
    /* Done in the dirty way to remove from the effective basis since the effective C7p does not involve the non-primed C_1 to C_6.*/
    C_7p += -ms_over_mb * C_7 - 1. / 3. * C_3 - 4 / 9 * C_4 - 20. / 3. * C_5 - 80. / 9. * C_6;
#endif    
    DC7_QCDF = deltaC7_QCDF(false);
    DC7_QCDF_bar = deltaC7_QCDF(true);
    
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    
    f_GSL = convertToGslFunction(boost::bind(&MVgamma::getT_perp_integrand_real, &(*this), _1));
    if (gsl_integration_cquad(&f_GSL, 0., 1., 1.e-2, 1.e-1, w_GSL, &average, &error, NULL) != 0) T_perp_real = std::numeric_limits<double>::quiet_NaN();
    T_perp_real = average;
    
    f_GSL = convertToGslFunction(boost::bind(&MVgamma::getT_perp_integrand_imag, &(*this), _1));
    if (gsl_integration_cquad(&f_GSL, 0., 1., 1.e-2, 1.e-1, w_GSL, &average, &error, NULL) != 0) T_perp_imag = std::numeric_limits<double>::quiet_NaN();
    T_perp_imag = average;
    
    f_GSL = convertToGslFunction(boost::bind(&MVgamma::getT_perp_bar_integrand_real, &(*this), _1));
    if (gsl_integration_cquad(&f_GSL, 0., 1., 1.e-2, 1.e-1, w_GSL, &average, &error, NULL) != 0) T_perp_bar_real = std::numeric_limits<double>::quiet_NaN();
    T_perp_bar_real = average;
    
    f_GSL = convertToGslFunction(boost::bind(&MVgamma::getT_perp_bar_integrand_imag, &(*this), _1));
    if (gsl_integration_cquad(&f_GSL, 0., 1., 1.e-2, 1.e-1, w_GSL, &average, &error, NULL) != 0) T_perp_bar_imag = std::numeric_limits<double>::quiet_NaN();
    T_perp_bar_imag = average;
    
    gsl_set_error_handler(old_handler);
    
    SM.getFlavour().setUpdateFlag(meson, vectorM, QCD::NOLEPTON, false);
    
}

/*******************************************************************************
 * Form Factor                                                     *
 * ****************************************************************************/
double MVgamma::T_1()
{
    return ( a_0T1);
}

gslpp::complex MVgamma::deltaC7_QCDF(bool conjugate)
{
    double muh = mu_b/mb_pole;
    double z = mc_pole*mc_pole/mb_pole/mb_pole;

#if FULLNLOQCDF_MVGAMMA    
    gslpp::complex A_Seidel = 1./729. * (833. + 120.*gslpp::complex::i()*M_PI - 312. * log(mb_pole*mb_pole/mu_b/mu_b)); /* hep-ph/0403185v2.*/
    gslpp::complex Fu_17 = -A_Seidel; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
    gslpp::complex Fu_27 = 6. * A_Seidel; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
#endif    
    gslpp::complex F_17 = myF_1->F_17re(muh, z, 0.00001, 20) + gslpp::complex::i() * myF_1->F_17im(muh, z, 0.00001, 20); /*q^2 = 0 gives nan. Independent of how small q^2 is. arXiv:0810.4077*/
    gslpp::complex F_27 = myF_2->F_27re(muh, z, 0.00001, 20) + gslpp::complex::i() * myF_2->F_27im(muh, z, 0.00001, 20); /*q^2 = 0 gives nan. Independent of how small q^2 is. arXiv:0810.4077*/
    gslpp::complex F_87 = (-4.*(33. + 24.*log(muh) + 6.*gslpp::complex::i()*M_PI - 2.*M_PI*M_PI))/27.; 
    
    if (!conjugate) {
        gslpp::complex delta = C_1 * F_17 + C_2 * F_27;
        gslpp::complex delta_t = C_8 * F_87 + delta;
#if FULLNLOQCDF_MVGAMMA        
        gslpp::complex delta_u = delta + C_1 * Fu_17 + C_2 * Fu_27;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - lambda_u / lambda_t * delta_u);
#else
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    } else {
        gslpp::complex delta = C_1.conjugate() * F_17 + C_2.conjugate() * F_27;
        gslpp::complex delta_t = C_8.conjugate() * F_87 + delta;
#if FULLNLOQCDF_MVGAMMA        
        gslpp::complex delta_u = delta + C_1.conjugate() * Fu_17 + C_2.conjugate() * Fu_27;
        return -alpha_s_mub / (4. * M_PI) * (delta_t - (lambda_u / lambda_t).conjugate() * delta_u);
#else        
        return -alpha_s_mub / (4. * M_PI) * delta_t;
#endif        
    }
}

gslpp::complex MVgamma::Cq34(bool conjugate)
{
    gslpp::complex T_t = -C_3 + 4./3.*(C_4 + 12.*C_5 + 16.*C_6);
    gslpp::complex T_u = 0.; /* 0 for K*0, phi*/
    if (meson == QCD::B_P) T_u = -3.*C_2;
    else if (vectorM == QCD::PHI) T_t = T_t + 6.*(C_3 + 10.*C_5);
    if (!conjugate) return T_t + lambda_u / lambda_t * T_u;
    else return T_t + (lambda_u / lambda_t).conjugate() * T_u;
}

gslpp::complex MVgamma::T_perp_WA_1()
{
    return -spectator_charge * 4./mb_pole * (C_3 + 4./3.*(C_4 + 3.*C_5 + 4.*C_6));
}

gslpp::complex MVgamma::T_perp_WA_2(bool conjugate)
{
    return spectator_charge * 2./mb_pole * Cq34(conjugate);
}

gslpp::complex MVgamma::L1(gslpp::complex x)
{
    if (x == 0.) return -(M_PI*M_PI/6.);
    if (x == 1.) return 0.;
    else return log((x-1.)/x)*log(1.-x) - (M_PI*M_PI/6.) + dilog(x/(x-1.));
}

double MVgamma::phi_V(double u)
{
    return 6.* u * (1. - u) * (1. + SM.getMesons(vectorM).getGegenalpha(0) * gsl_sf_gegenpoly_1(3./2., (2.*u - 1.)) + SM.getMesons(vectorM).getGegenalpha(1) * gsl_sf_gegenpoly_2(3./2., (2.*u - 1.)));
}

gslpp::complex MVgamma::t_perp(double u, double m)
{
    double ubar = 1. - u;
    gslpp::complex x0 = sqrt(0.25 - (m*m - gslpp::complex::i()*1.e-10)/(ubar * MM2));
    gslpp::complex xp = 0.5 + x0;
    gslpp::complex xm = 0.5 - x0;

    return 4./ubar * (1. + 2.*(m*m - gslpp::complex::i()*1.e-10)/(ubar*MM2) * (L1(xp) + L1(xm)));
}

gslpp::complex MVgamma::T_perp_plus_QSS(double u, bool conjugate)
{
    gslpp::complex t_perp_mc = t_perp(u, mc_pole);
    double eu = 2./3.;
#if FULLNLOQCDF_MVGAMMA     
    gslpp::complex t_perp_0 = t_perp(u, 0.);    
    double ed = -1./3.;
    gslpp::complex T_t = (alpha_s_mub/(3.*M_PI))*MM/(2.*mb_pole)*(eu * t_perp_mc * (C_1/6. + C_2 + 6.*C_6)
        + ed * t_perp(u, mb_pole) * (C_3 - C_4/6. + 16.*C_5 + 10.*C_6/3. + mb_pole/MM*(C_3 + C_4/6. - 4.*C_5 + 2.*C_6/3.))
        + ed * t_perp_0  * (-C_3 + C_4/6. - 16.*C_5 + 8.*C_6/3.));
    
    gslpp::complex T_u = ((alpha_s_mub/(3.*M_PI))*eu*MM/(2.*mb_pole)*(t_perp_mc - t_perp_0)*(C_2 - C_1/6.));
    if (!conjugate) return T_t + lambda_u / lambda_t * T_u;
    else return T_t + (lambda_u / lambda_t).conjugate() * T_u;
#else        
    return (alpha_s_mub/(3.*M_PI))*MM/(2.*mb_pole)*(eu * t_perp_mc * (C_1/6. + C_2 + 6.*C_6));
#endif    
}

gslpp::complex MVgamma::T_perp_plus_O8(double u) 
{   
    return -(alpha_s_mub/(3.*M_PI))*4.*(-1./3.)*C_8/u;
}

gslpp::complex MVgamma::T_perp(double u, bool conjugate) 
{
    double N = M_PI*M_PI/3.*fB*fperp/MM;
    gslpp::complex T_amp = N/SM.getMesons(meson).getLambdaM() * phi_V(u) * (T_perp_plus_O8(u) + T_perp_plus_QSS(u, conjugate));
#if FULLNLOQCDF_MVGAMMA    
    double ubar = 1. - u;
    T_amp += N * phi_V(u)/ubar * T_perp_WA_1() + N/SM.getMesons(meson).getLambdaM() * fpara/fperp * MV * T_perp_WA_2(conjugate); 
            /*last term proportional to T_perp_WA_2 is a constant but is included in the integral because u is integrated over the range [0,1]*/
#endif            
    return T_amp;
}

gslpp::complex MVgamma::T_QCDF_minus(bool conjugate)
{
    if (vectorM == StandardModel::RHO || vectorM == StandardModel::RHO_P || vectorM == StandardModel::OMEGA) return 0.; // Temporary
    if (!conjugate) return (T_perp_real + gslpp::complex::i() * T_perp_imag);
    else return (T_perp_bar_real + gslpp::complex::i() * T_perp_bar_imag);
}

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/

gslpp::complex MVgamma::h_lambda(int hel) 
{
    if (hel == 0) {
        return SU3_breaking * ( -1./(MM2*16.*M_PI*M_PI) * (
                ((MM+MV)*a_0A1) / (2.*MM) * ((- r1_2 + deltaC9_2) / (1. + r2_2 / mJ2) )*exp_Phase_2
                - lambda / (2.*MM*(MM+MV))*a_0V * ((- r1_1 + deltaC9_1) / (1. + r2_1 / mJ2) )*exp_Phase_1 ) );
    }
    else if (hel == 1) {
        return SU3_breaking * (-1./(MM2*16.*M_PI*M_PI) *
                (((MM+MV)*a_0A1) / (2.*MM) * ((- r1_2 + deltaC9_2) / (1. + r2_2 / mJ2) )*exp_Phase_2
                + lambda / (2.*MM*(MM+MV))*a_0V * ((- r1_1 + deltaC9_1) / (1. + r2_1 / mJ2) )*exp_Phase_1 ) );
    }
    else {
        std::stringstream out;
        out << hel;
        throw std::runtime_error("MVgamma: hel " + out.str() + " not implemented, can only be 1 (+) or 2 (-)");
    }
}

gslpp::complex MVgamma::H_V_m()
{
    return lambda_t * (((C_7 + DC7_QCDF) * T_1() + MM2/(MM2 - MV*MV) * T_QCDF_minus(false)) * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p()
{
    return lambda_t * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

gslpp::complex MVgamma::H_V_m_bar()
{
    return lambda_t.conjugate() * (((C_7 + DC7_QCDF_bar) * T_1() + MM2/(MM2 - MV*MV) * T_QCDF_minus(true)) * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p_bar()
{
    return lambda_t.conjugate() * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MVgamma::BR_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i), myAmpDB2(*(new AmpDB2(SM_i)))
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double BR_MVgamma::computeBR_MVgamma(QCD::meson meson, QCD::meson vector)
{
    QCD::meson meson_i = meson;
    QCD::meson vector_i = vector;
    
    SM.getFlavour().getMVgamma(meson_i, vector_i).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_p_bar();
    
    switch (vector_i) {
        case StandardModel::K_star:
        case StandardModel::K_star_P:
        case StandardModel::RHO:
        case StandardModel::RHO_P:
        case StandardModel::OMEGA:
            arg = myAmpDB2.getAmpBd(FULLNLO).arg();
            t_int = 1.;
            break;
        case StandardModel::PHI:
            arg = myAmpDB2.getAmpBs(FULLNLO).arg();
            /* For correctly defined polarization the numerator should be H_V_p().conjugate()*H_V_p_bar() + H_V_m().conjugate()*H_V_m_bar(). Switched to keep consistency with K*ll.*/
            /* See discussion around eq.53 in hep-ph/0510104*/
            ADG = 2.*(exp(gslpp::complex::i()*arg)*(HVp.conjugate()*HVm_bar + HVm.conjugate()*HVp_bar)).real() / (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2());
            ys = SM.getMesons(QCD::B_S).getDgamma_gamma()/2.;
            t_int = (1. - ADG * ys)/(1. - ys*ys);
            break;
        default:
            std::stringstream out;
            out << vector_i;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
    
    double GF = SM.getGF();
    double ale = SM.getAle();
    double MM = SM.getMesons(meson_i).getMass();
    double MM2 = MM * MM;
    double Mb = SM.getQuarks(QCD::BOTTOM).getMass();
    double MV = SM.getMesons(vector_i).getMass();
    double width = SM.getMesons(meson_i).computeWidth();
    double lambda = MM2 - pow(MV, 2.);
    
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda / (4. * width) * (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2()) * t_int;
}

double BR_MVgamma::computeThValue()
{
    return computeBR_MVgamma(meson, vectorM);
}

R_MVgamma::R_MVgamma(const StandardModel& SM_i, QCD::meson meson_1, QCD::meson vector_1, QCD::meson meson_2, QCD::meson vector_2)
: BR_MVgamma(SM_i, meson_1, vector_1)
{
    meson1 = meson_1;
    meson2 = meson_2;
    vector1 = vector_1;
    vector2 = vector_2;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson1, vector1).initializeMVgammaParameters());
    setParametersForObservable(SM.getFlavour().getMVgamma(meson2, vector2).initializeMVgammaParameters());
}

double R_MVgamma::computeThValue()
{
    return computeBR_MVgamma(meson1, vector1)/computeBR_MVgamma(meson2, vector2);
}

D0p_MVgamma::D0p_MVgamma(const StandardModel& SM_i, QCD::meson meson_1, QCD::meson vector_1, QCD::meson meson_2, QCD::meson vector_2)
: BR_MVgamma(SM_i, meson_1, vector_1)
{
    meson1 = meson_1;
    meson2 = meson_2;
    vector1 = vector_1;
    vector2 = vector_2;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson1, vector1).initializeMVgammaParameters());
    setParametersForObservable(SM.getFlavour().getMVgamma(meson2, vector2).initializeMVgammaParameters());
}

double D0p_MVgamma::computeThValue()
{
    return (computeBR_MVgamma(meson1, vector1) - computeBR_MVgamma(meson2, vector2)) / 
            (computeBR_MVgamma(meson1, vector1) + computeBR_MVgamma(meson2, vector2));
}

ACP_MVgamma::ACP_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) 
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ACP_MVgamma::computeACP_MVgamma(QCD::meson meson, QCD::meson vector)
{
    QCD::meson meson_i = meson;
    QCD::meson vector_i = vector;
    
    SM.getFlavour().getMVgamma(meson_i, vector_i).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson_i, vector_i).H_V_p_bar();
    double CC = ((HVp.abs2() + HVm.abs2() - HVp_bar.abs2() - HVm_bar.abs2())) / (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2());
    return -CC;          
}

double ACP_MVgamma::computeThValue()
{
    return computeACP_MVgamma(meson, vectorM);
}

DACP_MVgamma::DACP_MVgamma(const StandardModel& SM_i, QCD::meson meson_1, QCD::meson vector_1, QCD::meson meson_2, QCD::meson vector_2)
: ACP_MVgamma(SM_i, meson_1, vector_1)
{
    meson1 = meson_1;
    meson2 = meson_2;
    vector1 = vector_1;
    vector2 = vector_2;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson1, vector1).initializeMVgammaParameters());
    setParametersForObservable(SM.getFlavour().getMVgamma(meson2, vector2).initializeMVgammaParameters());
}

double DACP_MVgamma::computeThValue()
{
    return (computeACP_MVgamma(meson1, vector1) - computeACP_MVgamma(meson2, vector2));
}

C_MVgamma::C_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) 
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double C_MVgamma::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p_bar();
    /* REMEMBER: ACP = -C by definition in neutral B mesons.*/
    double CC = ((HVp.abs2() + HVm.abs2() - HVp_bar.abs2() - HVm_bar.abs2())) / (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2());
    if (meson == QCD::B_P) return -CC;
    else return CC;
            
}

S_MVgamma::S_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) 
: ThObservable(SM_i), myAmpDB2(*(new AmpDB2(SM_i)))
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double S_MVgamma::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p_bar();
    
    switch (vectorM) {
        case StandardModel::K_star:
        case StandardModel::RHO:
            arg = myAmpDB2.getAmpBd(FULLNLO).arg();
            break;
        case StandardModel::PHI:
            arg = myAmpDB2.getAmpBs(FULLNLO).arg();
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    /* For correctly defined polarization the numerator should be H_V_p().conjugate()*H_V_p_bar() + H_V_m().conjugate()*H_V_m_bar(). Switched to keep consistency with K*ll.*/
    /* See discussion around eq.53 in hep-ph/0510104*/
    return 2.*(exp(gslpp::complex::i()*arg)*(HVp.conjugate()*HVm_bar + HVm.conjugate()*HVp_bar)).imag() / (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2());
}

ADG_MVgamma::ADG_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) 
: ThObservable(SM_i), myAmpDB2(*(new AmpDB2(SM_i)))
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ADG_MVgamma::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p_bar();
    
    switch (vectorM) {
        case StandardModel::K_star:
            arg = myAmpDB2.getAmpBd(FULLNLO).arg();
            break;
        case StandardModel::PHI:
            arg = myAmpDB2.getAmpBs(FULLNLO).arg();
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    /* For correctly defined polarization the numerator should be H_V_p().conjugate()*H_V_p_bar() + H_V_m().conjugate()*H_V_m_bar(). Switched to keep consistency with K*ll.*/
    /* See discussion around eq.53 in hep-ph/0510104*/
    return 2.*(exp(gslpp::complex::i()*arg)*(HVp.conjugate()*HVm_bar + HVm.conjugate()*HVp_bar)).real() / (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2());
}

DC7_1::DC7_1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double DC7_1::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[1] - SM.getFlavour().getMVgamma(meson, vectorM).h[0])).abs();
}

DC7_2::DC7_2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double DC7_2::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[1] + SM.getFlavour().getMVgamma(meson, vectorM).h[0])).abs();
}

hp0_hm0::hp0_hm0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
   
   setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double hp0_hm0::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return SM.getFlavour().getMVgamma(meson, vectorM).h[0].abs()/SM.getFlavour().getMVgamma(meson, vectorM).h[1].abs();
}

AbsDC7_L::AbsDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double AbsDC7_L::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[1])).abs();
}

AbsDC7_R::AbsDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double AbsDC7_R::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[0])).abs();
}

ReDC7_L::ReDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ReDC7_L::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[1])).real();
}

ReDC7_R::ReDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ReDC7_R::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[0])).real();
}

ImDC7_L::ImDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ImDC7_L::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[1])).imag();
}

ImDC7_R::ImDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ImDC7_R::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    return ( (8. * M_PI * M_PI * SM.getFlavour().getMVgamma(meson, vectorM).MM2 * SM.getFlavour().getMVgamma(meson, vectorM).MM) / (SM.getFlavour().getMVgamma(meson, vectorM).lambda * SM.getFlavour().getMVgamma(meson, vectorM).Mb * SM.getFlavour().getMVgamma(meson, vectorM).T_1())*(SM.getFlavour().getMVgamma(meson, vectorM).h[0])).imag();
}

AbsDC7_QCDF::AbsDC7_QCDF(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double AbsDC7_QCDF::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(false)/T1 ).abs();
}

AbsDC7_QCDF_bar::AbsDC7_QCDF_bar(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double AbsDC7_QCDF_bar::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF_bar + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(true)/T1 ).abs();
}

ReDC7_QCDF::ReDC7_QCDF(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ReDC7_QCDF::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(false)/T1 ).real();
}

ReDC7_QCDF_bar::ReDC7_QCDF_bar(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ReDC7_QCDF_bar::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF_bar + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(true)/T1 ).real();
}

ImDC7_QCDF::ImDC7_QCDF(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ImDC7_QCDF::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(false)/T1 ).imag();
}

ImDC7_QCDF_bar::ImDC7_QCDF_bar(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).initializeMVgammaParameters());
}

double ImDC7_QCDF_bar::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double MV = SM.getMesons(vectorM).getMass();
    double T1 = SM.getFlavour().getMVgamma(meson, vectorM).T_1();
    
    return ( SM.getFlavour().getMVgamma(meson, vectorM).DC7_QCDF_bar + MM2/(MM2 - MV*MV) * SM.getFlavour().getMVgamma(meson, vectorM).T_QCDF_minus(true)/T1 ).imag();
}
