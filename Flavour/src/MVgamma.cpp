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
#include <gsl/gsl_sf_zeta.h>

MVgamma::MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) : ThObservable(SM_i), myBXqll(SM_i, QCD::BOTTOM, QCD::MU)
{
    meson = meson_i;
    vectorM = vector_i;
    if (vectorM == StandardModel::PHI) setParametersForObservable(make_vector<std::string>() << "a_0T1phi" << "absh_p" << "absh_m" << "argh_p" << "argh_m");
    if (vectorM == StandardModel::K_star) setParametersForObservable(make_vector<std::string>() << "a_0T1" << "absh_p" << "absh_m" << "argh_p" << "argh_m");
}

MVgamma::~MVgamma()
{
}

void MVgamma::updateParameters()
{
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
    width = SM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV, 2.);

    switch (vectorM) {
        case StandardModel::K_star:
            a_0T1 = SM.getOptionalParameter("a_0T1");
            
            fperp = SM.getFKstarp();

            break;
        case StandardModel::PHI:
            a_0T1 = SM.getOptionalParameter("a_0T1phi");
            
            fperp = SM.getFphip();

            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }


    h[0] = gslpp::complex(SM.getOptionalParameter("absh_p"), SM.getOptionalParameter("argh_p"), true); //h_plus
    h[1] = gslpp::complex(SM.getOptionalParameter("absh_m"), SM.getOptionalParameter("argh_m"), true); //h_minus
    
    allcoeff = SM.getFlavour().ComputeCoeffBMll(mu_b); //check the mass scale, scheme fixed to NDR
    allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu_b); //check the mass scale, scheme fixed to NDR
    
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
    /* Done in the dirty way to remove from the effective basis since C7p is not in the effective basis according to EOS.*/
    C_7p = ms_over_mb * (((*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6)) - C_7 - 1./3. * C_3 - 4/9 * C_4 - 20./3. * C_5 - 80./9. * C_6);
    C_1_bar = C_1/2.;
    C_2_bar = C_2 - C_1/6.;
    C_8 = (*(allcoeff[LO]))(7);

    allcoeffh = SM.getFlavour().ComputeCoeffBMll(mu_h); //check the mass scale, scheme fixed to NDR

    C_2h_bar = (*(allcoeffh[LO]))(1) - (*(allcoeffh[LO]))(0)/6.;
    C_8h = (*(allcoeffh[LO]))(7);
    
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
 */
gslpp::complex MVgamma::G1(double s)
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
}

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
/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MVgamma::H_V_m()
{
//    double s = Mc * Mc / Mb / Mb;
    
//    return lambda_t * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2_bar * G1(s) + C_8 * G8())
//            + SM.Als(mu_h) / 3. / M_PI * (C_2h_bar * H1(s) + C_8h * H8())) * T_1()
//            * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
    return lambda_t * ((C_7 + deltaC7_QCDF(false)) * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
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
    return lambda_t.conjugate() * ((C_7 + deltaC7_QCDF(true)) * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p_bar()
{
    return lambda_t.conjugate() * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MVgamma::BR_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double BR_MVgamma::computeThValue()
{
    updateParameters();
    
    switch (vectorM) {
        case StandardModel::K_star:
            arg = myAmpDB2.getAmpBd(FULLNLO).arg();
            t_int = 1.;
            break;
        case StandardModel::PHI:
            arg = myAmpDB2.getAmpBs(FULLNLO).arg();
            /* For correctly defined polarization the numerator should be H_V_p().conjugate()*H_V_p_bar() + H_V_m().conjugate()*H_V_m_bar(). Switched to keep consistency with K*ll.*/
            /* See discussion around eq.53 in hep-ph/0510104*/
            ADG = 2.*(exp(gslpp::complex::i()*arg)*(H_V_p().conjugate()*H_V_m_bar() + H_V_m().conjugate()*H_V_p_bar())).real() / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
            ys = SM.getMesons(QCD::B_S).getDgamma_gamma()/2.;
            t_int = (1. - ADG * ys)/(1. - ys*ys);
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }

    
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda / (4. * width) * (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2()) * t_int;
}

C_MVgamma::C_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) : MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double C_MVgamma::computeThValue()
{
    updateParameters();

    return ((H_V_p().abs2() + H_V_m().abs2() - H_V_p_bar().abs2() - H_V_m_bar().abs2())) / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

S_MVgamma::S_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) : MVgamma(SM_i, meson_i, vector_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double S_MVgamma::computeThValue()
{
    updateParameters();
    
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
    return 2.*(exp(gslpp::complex::i()*arg)*(H_V_p().conjugate()*H_V_m_bar() + H_V_m().conjugate()*H_V_p_bar())).imag() / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

ADG_MVgamma::ADG_MVgamma(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i) : MVgamma(SM_i, meson_i, vector_i), myAmpDB2(SM_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ADG_MVgamma::computeThValue()
{
    updateParameters();
    
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
    return 2.*(exp(gslpp::complex::i()*arg)*(H_V_p().conjugate()*H_V_m_bar() + H_V_m().conjugate()*H_V_p_bar())).real() / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

DC7_1::DC7_1(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double DC7_1::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[1] - h[0])).abs();
}

DC7_2::DC7_2(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double DC7_2::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[1] + h[0])).abs();
}

hp0_hm0::hp0_hm0(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
   meson = meson_i;
   vectorM = vector_i;
}

double hp0_hm0::computeThValue()
{
    updateParameters();
    return h[0].abs()/h[1].abs();
}

AbsDC7_L::AbsDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double AbsDC7_L::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[1])).abs();
}

AbsDC7_R::AbsDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double AbsDC7_R::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[0])).abs();
}

ReDC7_L::ReDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ReDC7_L::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[1])).real();
}

ReDC7_R::ReDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ReDC7_R::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[0])).real();
}

ImDC7_L::ImDC7_L(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ImDC7_L::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[1])).imag();
}

ImDC7_R::ImDC7_R(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ImDC7_R::computeThValue()
{
    updateParameters();
    return ( (8. * M_PI * M_PI * MM2 * MM) / (lambda * Mb * T_1())*(h[0])).imag();
}
