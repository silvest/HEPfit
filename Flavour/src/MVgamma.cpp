/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "MVgamma.h"
#include <gsl/gsl_sf_zeta.h>

MVgamma::MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : ThObservable(SM_i)
{
    if (SM.ModelName().compare("StandardModel") != 0 && SM.ModelName().compare("FlavourWilsonCoefficient") != 0) std::cout << "\nWARNING: B to V gamma not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
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

            break;
        case StandardModel::PHI:
            a_0T1 = SM.geta_0T1phi();

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

    C_7 = (*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6);
    C_7p = (*(allcoeffprime[LO]))(6) + (*(allcoeffprime[NLO]))(6);
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
            SM.getMesons(vectorM).getDecayconst() / T_1() / MM / SM.getMesons(meson).getLambdaM()*
            (c0 + c1 * SM.getMesons(vectorM).getGegenalpha(0) + c2 * SM.getMesons(vectorM).getGegenalpha(1));
}

gslpp::complex MVgamma::H8()
{
    return 4. * M_PI * M_PI / 9. * SM.getMesons(meson).getDecayconst() *
            SM.getMesons(vectorM).getDecayconst() / T_1() / MM / SM.getMesons(meson).getLambdaM()*
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
            * lambda / MM2 - MM / (2 * Mb)*16 * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p()
{
    return lambda_t * (-C_7p * T_1() * lambda / MM2 - MM / (2 * Mb)*16 * M_PI * M_PI * h[0]);
}

gslpp::complex MVgamma::H_V_m_bar()
{
    double s = Mc * Mc / Mb / Mb;
    return lambda_t.conjugate() * ((C_7 + SM.Als(mu_b) / 3. / M_PI * (C_2 * G1(s) + C_8 * G8())
            + SM.Als(mu_h) / 3. / M_PI * (C_2h * H1(s) + C_8h * H8())) * T_1() * lambda / MM2 - MM / (2 * Mb)*16 * M_PI * M_PI * h[1]);
}

gslpp::complex MVgamma::H_V_p_bar()
{
    return lambda_t.conjugate() * (-C_7p * T_1() * lambda / MM2 - MM / (2 * Mb)*16 * M_PI * M_PI * h[0]);
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

ACP_MVgamma::ACP_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i) : MVgamma(SM_i, meson_i, vector_i)
{
    meson = meson_i;
    vectorM = vector_i;
}

double ACP_MVgamma::computeThValue()
{
    updateParameters();

    return ((H_V_p().abs2() + H_V_m().abs2() - H_V_p_bar().abs2() - H_V_m_bar().abs2())) / (H_V_p().abs2() + H_V_m().abs2() + H_V_p_bar().abs2() + H_V_m_bar().abs2());
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