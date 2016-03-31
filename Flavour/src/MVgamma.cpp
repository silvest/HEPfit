/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVgamma.h"
#include <gsl/gsl_sf_zeta.h>

MVgamma::MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : ThObservable(SM_i)
{
    if (SM.getModelName().compare("StandardModel") != 0 && SM.getModelName().compare("FlavourWilsonCoefficient") != 0) std::cout << "\nWARNING: B to V gamma not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
    meson = meson_i;
    vectorM = vector_i;
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
    mu_b = SM.getMub();
    mu_h = sqrt(mu_b * .5); // From Beneke Neubert
    width = SM.getMesons(meson).computeWidth();
    lambda = MM2 - pow(MV, 2.);

    switch (vectorM) {
        case StandardModel::K_star:
            a_0T1 = SM.geta_0T1();
            
            fperp = SM.getFKstarp();

            break;
        case StandardModel::PHI:
            a_0T1 = SM.geta_0T1phi();
            
            fperp = SM.getFphip();

            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVgamma: vector " + out.str() + " not implemented");
    }


    h[0] = SM.geth_p(); //h_plus
    h[1] = SM.geth_m(); //h_minus

    allcoeff = SM.getMyFlavour()->ComputeCoeffBMll(mu_b); //check the mass scale, scheme fixed to NDR
    allcoeffprime = SM.getMyFlavour()->ComputeCoeffprimeBMll(mu_b); //check the mass scale, scheme fixed to NDR
    
    double ms_over_mb = SM.Mrun(mu_b, SM.getQuarks(QCD::STRANGE).getMass_scale(), 
                        SM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO)
                       /SM.Mrun(mu_b, SM.getQuarks(QCD::BOTTOM).getMass_scale(), 
                        SM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);

    C_3 = (*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2);
    C_4 = (*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3);
    C_5 = (*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4);
    C_6 = (*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5);
    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    /* Defined with a -ve sign since Jager et. al. 2013 define C7prime with a -ve sign while others define C7 with a _ve sign in the amplitude. See Altmannshofer et. al. 2008.*/
    /* Done in the dirty way to remove from the effective basis since C7p is not in the effective basis according to EOS.*/
    C_7p = ms_over_mb * (((*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6)) - C_7 - 1./3. * C_3 - 4/9 * C_4 - 20./3. * C_5 - 80./9. * C_6);
    C_2 = (*(allcoeff[LO]))(1);
    C_8 = (*(allcoeff[LO]))(7);

    allcoeffh = SM.getMyFlavour()->ComputeCoeffBMll(mu_h); //check the mass scale, scheme fixed to NDR

    C_2h = (*(allcoeffh[LO]))(1);
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

/*******************************************************************************
 * Helicity amplitudes                                                         *
 * ****************************************************************************/
gslpp::complex MVgamma::H_V_m()
{
    double s = Mc * Mc / Mb / Mb;
    return lambda_t * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2 * G1(s) + C_8 * G8())
            + SM.Als(mu_h) / 3. / M_PI * (C_2h * H1(s) + C_8h * H8())) * T_1()
            * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p()
{
    return lambda_t * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

gslpp::complex MVgamma::H_V_m_bar()
{
    double s = Mc * Mc / Mb / Mb;
    return lambda_t.conjugate() * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2 * G1(s) + C_8 * G8())
            + SM.Als(mu_h) / 3. / M_PI * (C_2h * H1(s) + C_8h * H8())) * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p_bar()
{
    return lambda_t.conjugate() * (-C_7p * T_1() * lambda / MM2 - MM / (2. * Mb)*16. * M_PI * M_PI * h[0]);
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


BR_MVgamma::BR_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i)
: MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double BR_MVgamma::computeThValue()
{
    updateParameters();
    
    return ale * pow(GF * Mb / (4 * M_PI * M_PI), 2.) * MM * lambda / (4. * width) * (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

C_MVgamma::C_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double C_MVgamma::computeThValue()
{
    updateParameters();

    return ((H_V_p().abs2() + H_V_m().abs2() - H_V_p_bar().abs2() - H_V_m_bar().abs2())) / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
}

S_MVgamma::S_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : MVgamma(SM_i, meson_i, vector_i), myAmpDB2(SM_i)
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

DC7_1::DC7_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i)
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

DC7_2::DC7_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i)
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