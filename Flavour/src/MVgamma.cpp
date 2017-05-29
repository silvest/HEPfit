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
#include <boost/bind.hpp>
#include <limits>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_gegenbauer.h>

MVgamma::MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: SM(SM_i), myBXqll(SM_i, QCD::BOTTOM, QCD::MU)
{
    meson = meson_i;
    vectorM = vector_i;
#if NFPOLARBASIS_MVGAMMA
    if (vectorM == StandardModel::PHI) mVgammaParameters = make_vector<std::string>() << "a_0T1phi" << "absh_p" << "absh_m" << "argh_p" << "argh_m";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) mVgammaParameters = make_vector<std::string>() << "a_0T1" << "absh_p" << "absh_m" << "argh_p" << "argh_m";
#else
    if (vectorM == StandardModel::PHI) mVgammaParameters = make_vector<std::string>() << "a_0T1phi" << "reh_p" << "reh_m" << "imh_p" << "imh_m";
    else if (vectorM == StandardModel::K_star || vectorM == StandardModel::K_star_P) mVgammaParameters = make_vector<std::string>() << "a_0T1" << "reh_p" << "reh_m" << "imh_p" << "imh_m";
#endif
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
        
    
    w_GSL = gsl_integration_cquad_workspace_alloc (100);
}

MVgamma::~MVgamma()
{
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
    Mc = SM.getQuarks(QCD::CHARM).getMass();
    Ms = SM.getQuarks(QCD::STRANGE).getMass();
    MW = SM.Mw();
    lambda_t = SM.computelamt_s();
    lambda_u = SM.computelamu_s();
    mu_b = SM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    fB = SM.getMesons(meson).getDecayconst();
    width = SM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV, 2.);
    alpha_s_mub = SM.Als(mu_b);

    switch (vectorM) {
        case StandardModel::K_star:
            a_0T1 = SM.getOptionalParameter("a_0T1");
            fperp = SM.getFKstarp();
            spectator_charge = SM.getQuarks(QCD::DOWN).getCharge();
            break;
        case StandardModel::K_star_P:
            a_0T1 = SM.getOptionalParameter("a_0T1");
            fperp = SM.getFKstarPp();
            spectator_charge = SM.getQuarks(QCD::UP).getCharge();
            break;
        case StandardModel::PHI:
            a_0T1 = SM.getOptionalParameter("a_0T1phi");
            fperp = SM.getFphip();
            spectator_charge = SM.getQuarks(QCD::STRANGE).getCharge();
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    fpara = SM.getMesons(vectorM).getDecayconst();
    
#if NFPOLARBASIS_MVGAMMA
        h[0] = gslpp::complex(SM.getOptionalParameter("absh_p"), SM.getOptionalParameter("argh_p"), true); //h_plus
        h[1] = gslpp::complex(SM.getOptionalParameter("absh_m"), SM.getOptionalParameter("argh_m"), true); //h_minus
#else
        h[0] = gslpp::complex(SM.getOptionalParameter("reh_p"), SM.getOptionalParameter("imh_p"), false); //h_plus
        h[1] = gslpp::complex(SM.getOptionalParameter("reh_m"), SM.getOptionalParameter("imh_m"), false); //h_minus
#endif
    
    allcoeff = SM.getFlavour().ComputeCoeffBMll(mu_b, QCD::MU); //check the mass scale, scheme fixed to NDR. QCD::MU does not make any difference to the WC necessary here.
    allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu_b, QCD::MU); //check the mass scale, scheme fixed to NDR. QCD::MU does not make any difference to the WC necessary here.
    
    double ms_over_mb = SM.Mrun(mu_b, SM.getQuarks(QCD::STRANGE).getMass_scale(), 
                        SM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO)
                       /SM.Mrun(mu_b, SM.getQuarks(QCD::BOTTOM).getMass_scale(), 
                        SM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);

    C_1 = (*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0);
    C_2 = (*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1);
    C_3 = (*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2);
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 = (*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4);
    C_6 = (*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5);
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    /* Defined with a -ve sign since Jager et. al. 2013 define C7prime with a -ve sign while others define C7 with a +ve sign in the amplitude. See Altmannshofer et. al. 2008.*/
    /* Done in the dirty way to remove from the effective basis since the effective C7p does not involve the non-primed C_1 to C_6.*/
    C_7p = ms_over_mb * (((*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6)) - C_7 - 1./3. * C_3 - 4/9 * C_4 - 20./3. * C_5 - 80./9. * C_6);
    C_1_bar = C_1/2.;
    C_2_bar = C_2 - C_1/6.;
    C_8 = (*(allcoeff[LO]))(7) + (*(allcoeff[NLO]))(7);

    allcoeffh = SM.getFlavour().ComputeCoeffBMll(mu_h, QCD::MU); //check the mass scale, scheme fixed to NDR

    C_2h_bar = (*(allcoeffh[LO]))(1) - (*(allcoeffh[LO]))(0)/6.;
    C_8h = (*(allcoeffh[LO]))(7);
    
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

/*
 * QCDF alpha_s corrections
 * The following commented out part is the old LO implementation. The NLO follows it.
 * Retaining code for cleaning up later.
 */
/* gslpp::complex MVgamma::G1(double s)
{
    double logs = log(s);
    double M2 = M_PI*M_PI;
    double M3 = M_PI * M_PI*M_PI;

    return -104. / 27. * log(mu_b / Mb) - 833. / 162. - 20. * gslpp::complex::i() * M_PI / 27. +
            8. * M2 / 9. * pow(s, 1.5) + 2. / 9. * (48. + 30. * gslpp::complex::i() * M_PI -
            5. * M2 - 2. * gslpp::complex::i() * M3 - 36. * gsl_sf_zeta(3.)+
            (36. + 6. * gslpp::complex::i() * M_PI - 9. * M2) * logs + (3. + 6. *
            gslpp::complex::i() * M_PI) * logs * logs + logs * logs * logs) * s + 2. / 9. * (
            18. + 2. * M2 - 2. * gslpp::complex::i() * M3 + (12. - 6. * M2) * logs + 6. *
            gslpp::complex::i() * M_PI * logs * logs + logs * logs * logs) * s * s + 1. / 27. * (-9. + 112. *
            gslpp::complex::i() * M_PI - 14. * M2 + (182. - 48. * gslpp::complex::i() * M_PI) * logs -
            126. * logs * logs) * s * s*s;
}

gslpp::complex MVgamma::G8()
{
    return 8. / 3. * log(mu_b / Mb) + 11. / 3. - 2. / 9. * M_PI * M_PI + 2. / 3. * gslpp::complex::i() * M_PI;
}

gslpp::complex MVgamma::H1(double s)
{
    gslpp::complex c0, c1, c2;

    c0 = 1.05171 + 1.02281 + gslpp::complex::i()*2.75305;
    c1 = 1.41919 + 0.413974 - gslpp::complex::i()*1.85404;
    c2 = 0.269769 - 1.73577 - gslpp::complex::i()*1.50017;
    
    return -2. * M_PI * M_PI / 9. * SM.getMesons(meson).getDecayconst() *
            fperp / T_1() / MM / SM.getMesons(meson).getLambdaM()*
            (c0 + c1 * SM.getMesons(vectorM).getGegenalpha(0) + c2 * SM.getMesons(vectorM).getGegenalpha(1));
}

gslpp::complex MVgamma::H8()
{
    return 4. * M_PI * M_PI / 9. * SM.getMesons(meson).getDecayconst() *
            fperp / T_1() / MM / SM.getMesons(meson).getLambdaM()*
            (3. - 3. * SM.getMesons(vectorM).getGegenalpha(0) + 3. * SM.getMesons(vectorM).getGegenalpha(1));
}*/

/*
 * QCDF NLO
 */

gslpp::complex MVgamma::deltaC7_QCDF(bool conjugate)
{
    double muh = mu_b/Mb;
    double z = Mc*Mc/Mb/Mb;
    
    gslpp::complex A_Seidel = 1./729. * (833. + 120.*gslpp::complex::i()*M_PI - 312. * log(Mb*Mb/mu_b/mu_b)); /* hep-ph/0403185v2.*/
    gslpp::complex Fu_17 = -A_Seidel; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
    gslpp::complex Fu_27 = 6. * A_Seidel; /* sign different from hep-ph/0403185v2 but consistent with hep-ph/0412400 */
    gslpp::complex F_17 = myBXqll.F_17re(muh, z, 0.00001, 20) + gslpp::complex::i() * myBXqll.F_17im(muh, z, 0.00001, 20); /*q^2 = 0 gives nan. Independent of how small q^2 is. arXiv:0810.4077*/
    gslpp::complex F_27 = myBXqll.F_27re(muh, z, 0.00001, 20) + gslpp::complex::i() * myBXqll.F_27im(muh, z, 0.00001, 20); /*q^2 = 0 gives nan. Independent of how small q^2 is. arXiv:0810.4077*/
    gslpp::complex F_87 = (-4.*(33. + 24.*log(muh) + 6.*gslpp::complex::i()*M_PI - 2.*M_PI*M_PI))/27.; 
    
    if (!conjugate) {
        gslpp::complex delta = C_1 * F_17 + C_2 * F_27;
        gslpp::complex delta_t = C_8 * F_87 + delta;
        gslpp::complex delta_u = delta + C_1 * Fu_17 + C_2 * Fu_27;

        return -SM.Als(mu_b) / (4. * M_PI) * (delta_t - lambda_u / lambda_t * delta_u);
    } else {
        gslpp::complex delta = C_1.conjugate() * F_17 + C_2.conjugate() * F_27;
        gslpp::complex delta_t = C_8.conjugate() * F_87 + delta;
        gslpp::complex delta_u = delta + C_1.conjugate() * Fu_17 + C_2.conjugate() * Fu_27;

        return -SM.Als(mu_b) / (4. * M_PI) * (delta_t - (lambda_u / lambda_t).conjugate() * delta_u);
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
    return -spectator_charge * 4./Mb * (C_3 + 4./3.*(C_4 + 3.*C_5 + 4.*C_6));
}

gslpp::complex MVgamma::T_perp_WA_2(bool conjugate)
{
    return spectator_charge * 2./Mb * Cq34(conjugate);
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
    gslpp::complex t_perp_mc = t_perp(u, Mc);
    gslpp::complex t_perp_0 = t_perp(u, 0.);
    double eu = 2./3.;
    double ed = -1./3.;
    
    gslpp::complex T_t = (alpha_s_mub/(3.*M_PI))*MM/(2.*Mb)*(eu * t_perp_mc * (C_1/6. + C_2 + 6.*C_6)
        + ed * t_perp(u, Mb) * (C_3 - C_4/6. + 16.*C_5 + 10.*C_6/3. + Mb/MM*(C_3 + C_4/6. - 4.*C_5 + 2.*C_6/3.))
        + ed * t_perp_0  * (C_3 + C_4/6. - 16.*C_5 + 8.*C_6/3.));
    
    gslpp::complex T_u = ((alpha_s_mub/(3.*M_PI))*eu*MM/(2.*Mb)*(t_perp_mc - t_perp_0)*(C_2 - C_1/6.));
    
    if (!conjugate) return T_t + lambda_u / lambda_t * T_u;
    else return T_t + (lambda_u / lambda_t).conjugate() * T_u;
    
}

gslpp::complex MVgamma::T_perp_plus_O8(double u) 
{   
    return -(alpha_s_mub/(3.*M_PI))*4.*(-1./3.)*C_8/u;
}

gslpp::complex MVgamma::T_perp(double u, bool conjugate) 
{
    double N = M_PI*M_PI/3.*fB*fperp/MM;
    double ubar = 1. - u;
    gslpp::complex T_amp = N/SM.getMesons(meson).getLambdaM() * phi_V(u) * (T_perp_plus_O8(u) + T_perp_plus_QSS(u, conjugate)) + N * phi_V(u)/ubar * T_perp_WA_1() + N/SM.getMesons(meson).getLambdaM() * fpara/fperp * MV * T_perp_WA_2(conjugate); /*last term proportional to T_perp_WA_2 is a constant but is included in the integral because u is integrated over the range [0,1]*/
    return T_amp;
}

gslpp::complex MVgamma::T_QCDF_minus(bool conjugate)
{
    if (!conjugate) return (T_perp_real + gslpp::complex::i() * T_perp_imag);
    else return (T_perp_bar_real + gslpp::complex::i() * T_perp_bar_imag);
}

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MVgamma::H_V_m()
{
//    double s = Mc * Mc / Mb / Mb;
    
//    return lambda_t * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2_bar * G1(s) + C_8 * G8())
//            + SM.Als(mu_h) / 3. / M_PI * (C_2h_bar * H1(s) + C_8h * H8())) * T_1()
//            * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
    return lambda_t * (((C_7 + DC7_QCDF) * T_1() + MM2/(MM2 - MV*MV) * T_QCDF_minus(false)) * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p()
{
    return lambda_t * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

gslpp::complex MVgamma::H_V_m_bar()
{
//    double s = Mc * Mc / Mb / Mb;

//    return lambda_t.conjugate() * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2_bar * G1(s) + C_8 * G8())
//            + SM.Als(mu_h) / 3. / M_PI * (C_2h_bar * H1(s) + C_8h * H8())) * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
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
: ThObservable(SM_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
}

double BR_MVgamma::computeThValue()
{
    SM.getFlavour().getMVgamma(meson, vectorM).updateParameters();
    
    gslpp::complex HVm = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m();
    gslpp::complex HVm_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_m_bar();
    gslpp::complex HVp = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p();
    gslpp::complex HVp_bar = SM.getFlavour().getMVgamma(meson, vectorM).H_V_p_bar();
    
    switch (vectorM) {
        case StandardModel::K_star:
        case StandardModel::K_star_P:
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
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }
    
    double GF = SM.getGF();
    double ale = SM.getAle();
    double MM = SM.getMesons(meson).getMass();
    double MM2 = MM * MM;
    double Mb = SM.getQuarks(QCD::BOTTOM).getMass();
    double MV = SM.getMesons(vectorM).getMass();
    double width = SM.getMesons(meson).computeWidth();
    double lambda = MM2 - pow(MV, 2.);
    
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda / (4. * width) * (HVp.abs2() + HVm.abs2() + HVp_bar.abs2() + HVm_bar.abs2()) * t_int;
}

C_MVgamma::C_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) 
: ThObservable(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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

S_MVgamma::S_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) :  ThObservable(SM_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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

ADG_MVgamma::ADG_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) :  ThObservable(SM_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
   
   setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
    
    setParametersForObservable(SM.getFlavour().getMVgamma(meson, vectorM).getMVgammaParameters());
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
