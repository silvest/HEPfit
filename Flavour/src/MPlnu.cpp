/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "MPlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <gsl/gsl_sf_zeta.h>
#include <boost/bind.hpp>
#include <limits>
#include <TFitResult.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_expint.h>

MPlnu::MPlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson pseudoscalar_i, QCD::lepton lep_i)
: mySM(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    pseudoscalarM = pseudoscalar_i;
    CLNflag = false;
    btocNPpmflag = false;

    w_J = gsl_integration_cquad_workspace_alloc(100);

    checkcache_int_tau = 0;
    checkcache_int_mu = 0;
    checkcache_int_el = 0;
}

MPlnu::~MPlnu()
{}

std::vector<std::string> MPlnu::initializeMPlnuParameters()
{
    CLNflag = mySM.getFlavour().getFlagCLN();
    btocNPpmflag = mySM.getFlavour().getbtocNPpmflag();

    if (btocNPpmflag) {
        if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
            << "af0_1" << "af0_2" << "afplus_0" << "afplus_1" << "afplus_2"
            << "mBc1m_1" << "mBc1m_2" << "mBc1m_3" << "mBc1m_4"
            << "mBc0p_1" << "mBc0p_2" << "chitildeT" << "chiL" << "nI"
            << "CS_NP" << "CP_NP" << "CV_NP" << "CA_NP" << "CT_NP";
        else {
            std::stringstream out;
            out << pseudoscalarM;
            throw std::runtime_error("MPlnu: vector " + out.str() + " not implemented");
        }

        if (CLNflag) {
            mplnuParameters.clear();
            if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
                << "fplusz0" << "rho1to2"
                << "CS_NP" << "CP_NP" << "CV_NP" << "CA_NP" << "CT_NP";
        }
    } else {
        if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
            << "af0_1" << "af0_2" << "afplus_0" << "afplus_1" << "afplus_2"
            << "mBc1m_1" << "mBc1m_2" << "mBc1m_3" << "mBc1m_4"
            << "mBc0p_1" << "mBc0p_2" << "chitildeT" << "chiL" << "nI"
            << "CSL_NP" << "CSR_NP" << "CVL_NP" << "CVR_NP" << "CT_NP";
        else {
            std::stringstream out;
            out << pseudoscalarM;
            throw std::runtime_error("MPlnu: vector " + out.str() + " not implemented");
        }

        if (CLNflag) {
            mplnuParameters.clear();
            if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
                << "fplusz0" << "rho1to2" 
                << "CSL_NP" << "CSR_NP" << "CVL_NP" << "CVR_NP" << "CT_NP";
        }
    }

    mySM.initializeMeson(meson);
    mySM.initializeMeson(pseudoscalarM);
    return mplnuParameters;
}

void MPlnu::updateParameters() 
{
    if (!mySM.getFlavour().getUpdateFlag(meson, pseudoscalarM, lep)) return;
    
    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
    Mnu = 0.; // neutrinos assumed to be massless
    MM = mySM.getMesons(meson).getMass();
    MP = mySM.getMesons(pseudoscalarM).getMass();
    width = mySM.getMesons(meson).computeWidth();
    w0 = (MM*MM+MP*MP)/(2.*MM*MP);
    RV = 2.*sqrt(MM*MP)/(MM+MP);
    mu_b = mySM.getMub();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass(); // add the PS b mass
    Vcb = mySM.getCKM().getV_cb(); // mySM.getOptionalParameter("AbsVcb");
    ale_mub = mySM.Ale(mu_b,FULLNLO);
    /* Amplitude propto 4*GF*Vij/sqrt(2) & kinematics requires 1/(2^9 pi^3 MB^3) */
    amplsq_factor = GF*GF*Vcb.abs2()/(64.*M_PI*M_PI*M_PI*MM*MM*MM);
    q2min = Mlep*Mlep;
    q2max = (MM-MP)*(MM-MP);
    
    /* SM Wilson coefficients */
    eta_EW = 1.0066;
    CV_SM = 1./2.;
    CV = CV_SM*eta_EW;
    CA = -CV_SM*eta_EW;
    CVp = 0.;
    CAp = 0.;
    CS = 0.;
    CSp = 0.;
    CP = 0.;
    CPp = 0.;
    C7 = 0.;
    C7p = 0.;
    CT = 0.;
    CTp = 0.;

    /* SM + NP Wilson coefficients */
    if (lep == StandardModel::TAU) {
        if (btocNPpmflag) {
            CV += (mySM.getOptionalParameter("CV_NP") - mySM.getOptionalParameter("CA_NP")) / 2. / M_SQRT2;
            CVp = (mySM.getOptionalParameter("CV_NP") + mySM.getOptionalParameter("CA_NP")) / 2. / M_SQRT2;
            CA -= (mySM.getOptionalParameter("CV_NP") - mySM.getOptionalParameter("CA_NP")) / 2. / M_SQRT2;
            CAp = -(mySM.getOptionalParameter("CV_NP") + mySM.getOptionalParameter("CA_NP")) / 2. / M_SQRT2;
            CS = (mySM.getOptionalParameter("CS_NP") - mySM.getOptionalParameter("CP_NP")) / 2. / M_SQRT2;
            CSp = (mySM.getOptionalParameter("CS_NP") + mySM.getOptionalParameter("CP_NP")) / 2. / M_SQRT2;
            CP = -(mySM.getOptionalParameter("CS_NP") - mySM.getOptionalParameter("CP_NP")) / 2. / M_SQRT2;
            CPp = -(mySM.getOptionalParameter("CS_NP") + mySM.getOptionalParameter("CP_NP")) / 2. / M_SQRT2;
            CTp = mySM.getOptionalParameter("CT_NP");
        } else {
            CV += mySM.getOptionalParameter("CVL_NP") / 2.;
            CVp = mySM.getOptionalParameter("CVR_NP") / 2.;
            CA -= mySM.getOptionalParameter("CVL_NP") / 2.;
            CAp = -mySM.getOptionalParameter("CVR_NP") / 2.;
            CS = mySM.getOptionalParameter("CSL_NP") / 2.;
            CSp = mySM.getOptionalParameter("CSR_NP") / 2.;
            CP = -mySM.getOptionalParameter("CSL_NP") / 2.;
            CPp = -mySM.getOptionalParameter("CSR_NP") / 2.;
            CTp = mySM.getOptionalParameter("CT_NP");
        }
    }
    
    switch (pseudoscalarM) {
        case StandardModel::D_P:
            if (CLNflag) {
                fplusz0 = mySM.getOptionalParameter("fplusz0");
                rho1to2 = mySM.getOptionalParameter("rho1to2");
                af0_1 = 0.;
                af0_2 = 0.;
                afplus_1 = 0.;
                afplus_2 = 0.;
                mBc1m_1 = 0.;
                mBc1m_2 = 0.;
                mBc1m_3 = 0.;
                mBc1m_4 = 0.;
                mBc0p_1 = 0.;
                mBc0p_2 = 0.;
                chitildeT = 0.;
                chiL = 0.;
                nI = 0.;
            } else {
                fplusz0 = 0.;
                rho1to2 = 0.;
                af0_1 = mySM.getOptionalParameter("af0_1");
                af0_2 = mySM.getOptionalParameter("af0_2");
                afplus_0 = mySM.getOptionalParameter("afplus_0");
                afplus_1 = mySM.getOptionalParameter("afplus_1");
                afplus_2 = mySM.getOptionalParameter("afplus_2");
                mBc1m_1 = mySM.getOptionalParameter("mBc1m_1");
                mBc1m_2 = mySM.getOptionalParameter("mBc1m_2");
                mBc1m_3 = mySM.getOptionalParameter("mBc1m_3");
                mBc1m_4 = mySM.getOptionalParameter("mBc1m_4");
                mBc0p_1 = mySM.getOptionalParameter("mBc0p_1");
                mBc0p_2 = mySM.getOptionalParameter("mBc0p_2");
                chitildeT = mySM.getOptionalParameter("chitildeT");
                chiL = mySM.getOptionalParameter("chiL");
                nI = mySM.getOptionalParameter("nI");
            }
            break;
        default:
            std::stringstream out;
            out << pseudoscalarM;
            throw std::runtime_error("MPlnu: vector " + out.str() + " not implemented");
    }
    
    z1m_1 = sqrt((MM+MP)*(MM+MP)-mBc1m_1*mBc1m_1)-sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP));
    z1m_1 /= (sqrt((MM+MP)*(MM+MP)-mBc1m_1*mBc1m_1)+sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP)));
    z1m_2 = sqrt((MM+MP)*(MM+MP)-mBc1m_2*mBc1m_2)-sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP));
    z1m_2 /= (sqrt((MM+MP)*(MM+MP)-mBc1m_2*mBc1m_2)+sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP)));
    z1m_3 = sqrt((MM+MP)*(MM+MP)-mBc1m_3*mBc1m_3)-sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP));
    z1m_3 /= (sqrt((MM+MP)*(MM+MP)-mBc1m_3*mBc1m_3)+sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP)));

    z0p_1 = sqrt((MM+MP)*(MM+MP)-mBc0p_1*mBc0p_1)-sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP));
    z0p_1 /= (sqrt((MM+MP)*(MM+MP)-mBc0p_1*mBc0p_1)+sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP)));
    z0p_2 = sqrt((MM+MP)*(MM+MP)-mBc0p_2*mBc0p_2)-sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP));
    z0p_2 /= (sqrt((MM+MP)*(MM+MP)-mBc0p_2*mBc0p_2)+sqrt((MM+MP)*(MM+MP)-(MM-MP)*(MM-MP)));

    if ((fplusz0 != fplusz0_cache) || (rho1to2 != rho1to2_cache)
            || (af0_1 != af0_1_cache) || (af0_2 != af0_2_cache)
            || (afplus_0 != afplus_0_cache) || (afplus_1 != afplus_1_cache)
            || (afplus_2 != afplus_2_cache)
            || (CS != CS_cache) || (CSp != CSp_cache)
            || (CP != CP_cache) || (CPp != CPp_cache)
            || (CV != CV_cache) || (CVp != CVp_cache)
            || (CA != CA_cache) || (CAp != CAp_cache)
            || (CT != CT_cache) || (CTp != CTp_cache)) {
        checkcache_int_tau = 0;
        checkcache_int_mu = 0;
        checkcache_int_el = 0;
    }

    if ((checkcache_int_tau == 0) || (checkcache_int_mu == 0) || (checkcache_int_el == 0)) {
        if (lep == StandardModel::TAU) {
            cached_intJ1_tau = integrateJ(1, q2min, q2max);
            cached_intJ3_tau = integrateJ(3, q2min, q2max);
            cached_intJ2_tau = 0.;
            // cached_intJ2_tau = integrateJ(2,q2min,q2max);
            checkcache_int_tau = 1;
        }
        if (lep == StandardModel::MU) {
            cached_intJ1_mu = integrateJ(1, q2min, q2max);
            cached_intJ3_mu = integrateJ(3, q2min, q2max);
            cached_intJ2_mu = 0.;
            // cached_intJ2_mu = integrateJ(2,q2min,q2max);
            checkcache_int_mu = 1;
        }
        if (lep == StandardModel::ELECTRON) {
            cached_intJ1_el = integrateJ(1, q2min, q2max);
            cached_intJ3_el = integrateJ(3, q2min, q2max);
            cached_intJ2_el = 0.;
            // cached_intJ2_el = integrateJ(2,q2min,q2max);
            checkcache_int_el = 1;
        }
    }

    fplusz0_cache = fplusz0;
    rho1to2_cache = rho1to2;

    af0_1_cache = af0_1;
    af0_2_cache = af0_2;
    afplus_0_cache = afplus_0;
    afplus_1_cache = afplus_1;
    afplus_2_cache = afplus_2;
    /* f+(q2=0) = f0(q2=0) */
    z0 = (sqrt(w0+1.)-sqrt(2.))/(sqrt(w0+1.)+sqrt(2.));
    if (CLNflag) af0_0 = 0.;
    else {
        af0_0 = fplus(0.);
        af0_0 -= (af0_1 * z0 + af0_2 * z0 * z0) / phi_f0(z0) / ((z0 - z0p_1) / (1. - z0 * z0p_1)*(z0 - z0p_2) / (1. - z0 * z0p_2));
        af0_0 *= phi_f0(z0)*((z0 - z0p_1) / (1. - z0 * z0p_1)*(z0 - z0p_2) / (1. - z0 * z0p_2));
    }
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

    mySM.getFlavour().setUpdateFlag(meson, pseudoscalarM, lep, false);

    return;
    
}

/*******************************************************************************
 * Kinematic functions                                                          *
 * ****************************************************************************/

double MPlnu::lambda_half(double a, double b, double c) 
{   
    return sqrt(a*a+b*b+c*c-2.*(a*b+a*c+b*c));
}
/*******************************************************************************
 * Form factors                                                                *
 * ****************************************************************************/

double MPlnu::phi_fplus(double z)
{
    // chitildeT in GeV-2
    double prefac = 8. * (MP / MM)*(MP / MM) / MM * sqrt(8. * nI / (3. * M_PI * chitildeT));
    double num = (1. + z)*(1. + z) * sqrt(1. - z);
    double den = (1. + MP / MM)*(1. - z) + 2. * sqrt(MP / MM)*(1. + z);
    double den5 = den * den * den * den*den;
    return prefac * num / den5;
}

double MPlnu::phi_f0(double z)
{
    // chiL dimensionless
    double prefac = (MP / MM)*(1. - (MP / MM)*(MP / MM)) * sqrt(8. * nI / (M_PI * chiL));
    double num = (1. - z * z) * sqrt((1. - z));
    double den = (1. + MP / MM)*(1. - z) + 2. * sqrt(MP / MM)*(1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MPlnu::fplus(double q2)
{
    double w = w0 - q2 / (2. * MM * MP);
    double z = (sqrt(w + 1.) - sqrt(2.)) / (sqrt(w + 1.) + sqrt(2.));
    if (CLNflag) {
        return fplusz0 * (1. - 8. * rho1to2 * z + (51. * rho1to2 - 10.) * z * z - (252. * rho1to2 - 84.) * z * z * z);
    } else {
        double P_fplus = (z - z1m_1) / (1. - z * z1m_1)*(z - z1m_2) / (1. - z * z1m_2)*(z - z1m_3) / (1. - z * z1m_3);
        return (afplus_0 + afplus_1 * z + afplus_2 * z * z) / phi_fplus(z) / P_fplus;
    }
}

double MPlnu::f0(double q2)
{
    double w = w0 - q2 / (2. * MM * MP);
    double z = (sqrt(w + 1.) - sqrt(2.)) / (sqrt(w + 1.) + sqrt(2.));
    if (CLNflag) {
        double prefac = fplus(q2)*2. * sqrt(MP / MM) / (1. + MP / MM)*2. * sqrt(MP / MM) / (1. + MP / MM)*(1. + w) / 2.;
        return prefac * 1.0036 * (1. - 0.0068 * (w - 1.) + 0.0017 * (w - 1.)*(w - 1.) - 0.0013 * (w - 1.)*(w - 1.)*(w - 1.));
    } else {
        double P_f0 = (z - z0p_1) / (1. - z * z0p_1)*(z - z0p_2) / (1. - z * z0p_2);
        return (af0_0 + af0_1 * z + af0_2 * z * z) / phi_f0(z) / P_f0;
    }
}

double MPlnu::fT(double q2)
{
    return 0.;
}
/********************************************************************************
 * Helicity amplitudes  (normalization such that all H \propto (mass scale)^-1) *
 * *****************************************************************************/

gslpp::complex MPlnu::HV(double q2)
{
    return lambda_half(MM*MM, MP*MP, q2) / 2. / sqrt(q2)*((CV + CVp) * fplus(q2) + 2. * Mb / (MM + MP)*(C7 + C7p) * fT(q2));
}

gslpp::complex MPlnu::HA(double q2)
{
    return lambda_half(MM*MM, MP*MP, q2) / 2. / sqrt(q2)*(CA + CAp) * fplus(q2);
}

gslpp::complex MPlnu::HP(double q2)
{
    return (MM * MM - MP * MP) / 2. * ((CP + CPp) / (Mb - Mc)+(Mlep + Mnu) / q2 * (CA + CAp)) * f0(q2);
}

gslpp::complex MPlnu::HS(double q2)
{
    return (MM * MM - MP * MP) / 2. * ((CS + CSp) / (Mb - Mc)+(Mlep - Mnu) / q2 * (CV + CVp)) * f0(q2);
}

gslpp::complex MPlnu::HT(double q2)
{
    return -gslpp::complex::i() * lambda_half(MM*MM, MP*MP, q2) / sqrt(2.) / (MM + MP)*(CT - CTp) * fT(q2);
}

gslpp::complex MPlnu::HTt(double q2)
{
    return -gslpp::complex::i() * lambda_half(MM*MM, MP*MP, q2) / 2. / (MM + MP)*(CT + CTp) * fT(q2);
}
/*******************************************************************************
 * Generalized angular coefficients  (see 1506.03970)                          *
 * ****************************************************************************/

gslpp::complex MPlnu::G0(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MP*MP, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep * Mlep + lambda_lep2 / (4. * q2));
    double Enu = sqrt(Mnu * Mnu + lambda_lep2 / (4. * q2));
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * ((4. * (Elep * Enu + Mlep * Mnu) + lambda_lep2 / (3. * q2)) * HV(q2).abs2()+
            (4. * (Elep * Enu - Mlep * Mnu) + lambda_lep2 / (3. * q2)) * HA(q2).abs2()+
            (4. * (Elep * Enu - Mlep * Mnu) + lambda_lep2 / q2) * HS(q2).abs2()+
            (4. * (Elep * Enu + Mlep * Mnu) + lambda_lep2 / q2) * HP(q2).abs2()+
            (8. * (Elep * Enu - Mlep * Mnu) - lambda_lep2 / (12. * q2)) * HT(q2).abs2()+
            (16. * (Elep * Enu + Mlep * Mnu) - lambda_lep2 / (12. * q2)) * HTt(q2).abs2() +
            8. * sqrt(2.)*(Enu * Mlep - Elep * Mnu)*(HA(q2) * HT(q2).conjugate()).imag() +
            16. * (Enu * Mlep + Elep * Mnu)*(HV(q2) * HTt(q2).conjugate()).imag());
}

gslpp::complex MPlnu::G1(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MP*MP, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 4. * lambda_lep * ((HA(q2)*(Mlep - Mnu) * HP(q2).conjugate()
            + HV(q2)*(Mlep + Mnu) * HS(q2).conjugate()).real() / sqrt(q2)
            -(sqrt(2.) * HT(q2) * HP(q2).conjugate() + 2. * HTt(q2) * HS(q2).conjugate()).imag());
}

gslpp::complex MPlnu::G2(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MP*MP, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * (4. * lambda_lep2 / (3. * q2))*(HA(q2).abs2() + HV(q2).abs2() - 2. * HT(q2).abs2() - 4. * HTt(q2).abs2());
}
/***************************************************************************
 * 12 independent J angular coefficients  (see again 1506.03970)           *
 * ************************************************************************/

double MPlnu::J1(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MP)*(MM - MP))) return 0.;
    return amplsq_factor * (G0(q2) - G2(q2) / 2).real();
}

double MPlnu::J2(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MP)*(MM - MP))) return 0.;
    return amplsq_factor * G1(q2).real();
}

double MPlnu::J3(double q2) 
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MP)*(MM - MP))) return 0.;
    return amplsq_factor * 3. / 2. * G2(q2).real();
}

double MPlnu::dGammadq2(double q2)
{
    if ((q2 < q2min) or (q2 > (MM - MP)*(MM - MP))) return 0.;
    double sqlambdaB = lambda_half(q2,MM*MM,MP*MP);
    double prefac = (CV-CA)*(CV-CA)*GF*GF*Vcb.abs2()*MM*sqlambdaB/192./M_PI/M_PI/M_PI;
    double coeff_fp = (1.+Mlep*Mlep/(2.*q2))*sqlambdaB*sqlambdaB/MM/MM/MM/MM;
    double coeff_f0 = (1.-MP*MP/MM/MM)*(1.-MP*MP/MM/MM)*3*Mlep*Mlep/(2.*q2);
    double TotAmp2 = coeff_fp*fplus(q2)*fplus(q2)+coeff_f0*f0(q2)*f0(q2);
    return prefac*(1.-Mlep*Mlep/q2)*(1.-Mlep*Mlep/q2)*TotAmp2;
 }

/***************************************************************************
 * Integration of angular coefficients Js                                  *
 * ************************************************************************/

double MPlnu::integrateJ(int i, double q2_min, double q2_max)
{
    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 1:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1_el;
            FJ = convertToGslFunction(boost::bind(&MPlnu::J1, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 2:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2_el;
            FJ = convertToGslFunction(boost::bind(&MPlnu::J2, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 3:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_el;
            FJ = convertToGslFunction(boost::bind(&MPlnu::J3, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
        case 4:
            FJ = convertToGslFunction(boost::bind(&MPlnu::dGammadq2, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        default:
            gsl_set_error_handler(old_handler);
            std::stringstream out;
            out << i;
            throw std::runtime_error("MPlnu::integrateJ: index " + out.str() + " not implemented");
    }
}

/* misleading name, needs to be changed: this is like dGammadq2 ... */
double MPlnu::dGammadw(double q2)
{
    updateParameters();

    return 2. * (J1(q2) + J3(q2) / 3.);
}

/* misleading name, needs to be changed: this is like DeltaGamma only ... */
double MPlnu::getDeltaGammaDeltaw(double w_min, double w_max)
{
    updateParameters();
    
    double q2_min = (2. * MM * MP)*(w0 - w_max); // min is Mlep*Mlep;
    double q2_max = (2. * MM * MP)*(w0 - w_min); // max is (MM-MP)*(MM-MP);
    /*
    double intJ1 = integrateJ(1, q2_min, q2_max);
    double intJ3 = integrateJ(3, q2_min, q2_max);

    return 2. * (intJ1 + intJ3 / 3.);
    */
    return integrateJ(4, q2_min, q2_max);
}

double MPlnu::get_unitarity_1min_BGL()
{
    updateParameters();

    return afplus_0 * afplus_0 + afplus_1 * afplus_1 + afplus_2 * afplus_2;

}

double MPlnu::get_unitarity_0plus_BGL()
{
    updateParameters();

    return af0_0 * af0_0 + af0_1 * af0_1 + af0_2*af0_2;
}

double MPlnu::get_fplus(double q2)
{
    updateParameters();

    return fplus(q2);
}

double MPlnu::get_f0(double q2)
{
    updateParameters();

    return f0(q2);
}

double MPlnu::get_fT(double q2)
{
    updateParameters();

    return fT(q2);
}
