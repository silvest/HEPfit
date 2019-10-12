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
#include <limits>

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
    
    double max_double = std::numeric_limits<double>::max();
    
    fplusz0_cache = max_double;
    rho1to2_cache = max_double;
    N_0_cache = max_double;
    alpha_0_cache = max_double;
    alpha_p_cache = max_double;
    beta_0_cache = max_double;
    beta_p_cache = max_double;
    gamma_0_cache = max_double;
    gamma_p_cache = max_double;
    af0_1_cache = max_double;
    af0_2_cache = max_double;
    afplus_0_cache = max_double;
    afplus_1_cache = max_double;
    afplus_2_cache = max_double;
    
#if NBGL == 3
    af0_3_cache = max_double;
    afplus_3_cache = max_double;
#endif 
    
    CS_cache = max_double;
    CSp_cache = max_double;
    CP_cache = max_double;
    CPp_cache = max_double;
    CV_cache = max_double;
    CVp_cache = max_double;
    CA_cache = max_double;
    CAp_cache = max_double;
    CT_cache = max_double;
    CTp_cache = max_double;
}

MPlnu::~MPlnu()
{}

std::vector<std::string> MPlnu::initializeMPlnuParameters()
{
    CLNflag = mySM.getFlavour().getFlagCLN();
    btocNPpmflag = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0);
    NPanalysis = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0 || mySM.getModelName().compare("RealWeakEFTCC") == 0);

    if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
        << "af0_1" << "af0_2" << "afplus_0" << "afplus_1" << "afplus_2"
#if NBGL == 3
        << "af0_3" << "afplus_3"
#endif                
        << "mBc1m_1" << "mBc1m_2" << "mBc1m_3" << "mBc1m_4"
        << "mBc0p_1" << "mBc0p_2" << "chitildeT" << "chiL" << "nI";
    else {
        std::stringstream out;
        out << pseudoscalarM;
        throw std::runtime_error("MPlnu: vector " + out.str() + " not implemented");
    }

    if (CLNflag) {
        mplnuParameters.clear();
        if (pseudoscalarM == StandardModel::D_P) mplnuParameters = make_vector<std::string>()
            << "fplusz0" << "rho1to2"
            << "N_0" << "alpha_0" << "alpha_p" << "beta_0" << "beta_p" << "gamma_0" << "gamma_p";
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
    mu_b = MM; // mySM.getMub();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass(); // add the PS b mass
    Vcb = mySM.getCKM().getV_cb(); // mySM.getOptionalParameter("AbsVcb");
    ale_mub = mySM.Ale(mu_b,FULLNLO);
    /* Amplitude propto 4*GF*Vij/sqrt(2) & kinematics requires 1/(2^9 pi^3 MB^3) */
    amplsq_factor = GF*GF*Vcb.abs2()/(64.*M_PI*M_PI*M_PI*MM*MM*MM);
    q2min = Mlep*Mlep;
    q2max = (MM-MP)*(MM-MP);
    
    /* SM Wilson coefficients */
    eta_EW = (1.+ale_mub/M_PI*log(mySM.getMz()/mu_b));
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
    if (NPanalysis) {
        if (lep == StandardModel::TAU) {
            if (btocNPpmflag) {
                CV += (mySM.getCCC3() - mySM.getCCC4()) / 2. / M_SQRT2;
                CVp = (mySM.getCCC3() + mySM.getCCC4()) / 2. / M_SQRT2;
                CA -= (mySM.getCCC3() - mySM.getCCC4()) / 2. / M_SQRT2;
                CAp = -(mySM.getCCC3() + mySM.getCCC4()) / 2. / M_SQRT2;
                CS = (mySM.getCCC1() - mySM.getCCC2()) / 2. / M_SQRT2;
                CSp = (mySM.getCCC1() + mySM.getCCC2()) / 2. / M_SQRT2;
                CP = -(mySM.getCCC1() - mySM.getCCC2()) / 2. / M_SQRT2;
                CPp = -(mySM.getCCC1() + mySM.getCCC2()) / 2. / M_SQRT2;
                CTp = mySM.getCCC5();
            } else {
                CV += mySM.getCCC3() / 2.;
                CVp = mySM.getCCC4() / 2.;
                CA -= mySM.getCCC3() / 2.;
                CAp = -mySM.getCCC4() / 2.;
                CS = mySM.getCCC1() / 2.;
                CSp = mySM.getCCC2() / 2.;
                CP = -mySM.getCCC1() / 2.;
                CPp = -mySM.getCCC2() / 2.;
                CTp = mySM.getCCC5();
            }
        }
    }
    
    switch (pseudoscalarM) {
        case StandardModel::D_P:
            if (CLNflag) {
                fplusz0 = mySM.getOptionalParameter("fplusz0");
                rho1to2 = mySM.getOptionalParameter("rho1to2");
                N_0 = mySM.getOptionalParameter("N_0");
                alpha_0 = mySM.getOptionalParameter("alpha_0"); 
                alpha_p = mySM.getOptionalParameter("alpha_p");
                beta_0 = mySM.getOptionalParameter("beta_0");
                beta_p = mySM.getOptionalParameter("beta_p");
                gamma_0 = mySM.getOptionalParameter("gamma_0");
                gamma_p = mySM.getOptionalParameter("gamma_p");
                af0_1 = 0.;
                af0_2 = 0.;
                afplus_1 = 0.;
                afplus_2 = 0.;
#if NBGL == 3
                af0_3 = 0.;
                afplus_3 = 0.;
#endif                
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
                N_0 = 0.;
                alpha_0 = 0.;
                alpha_p = 0.;
                beta_0 = 0.;
                beta_p = 0.;
                gamma_0 = 0.;
                gamma_p = 0.;
                af0_1 = mySM.getOptionalParameter("af0_1");
                af0_2 = mySM.getOptionalParameter("af0_2");
                afplus_0 = mySM.getOptionalParameter("afplus_0");
                afplus_1 = mySM.getOptionalParameter("afplus_1");
                afplus_2 = mySM.getOptionalParameter("afplus_2");
#if NBGL == 3
                af0_3 = mySM.getOptionalParameter("af0_3");
                afplus_3 = mySM.getOptionalParameter("afplus_3");
#endif                
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
            || (N_0 != N_0_cache) 
            || (alpha_0 != alpha_0_cache) || (alpha_p != alpha_p_cache)
            || (beta_0 != beta_0_cache) || (beta_p != beta_p_cache)
            || (gamma_0 != gamma_0_cache) || (gamma_p != gamma_p_cache)
            || (af0_1 != af0_1_cache) || (af0_2 != af0_2_cache)
            || (afplus_0 != afplus_0_cache) || (afplus_1 != afplus_1_cache)
            || (afplus_2 != afplus_2_cache)
#if NBGL == 3
            || (af0_3 != af0_3_cache) || (afplus_3 != afplus_3_cache)
#endif            
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
    N_0_cache = N_0;
    alpha_0_cache = alpha_0;
    alpha_p_cache = alpha_p;
    beta_0_cache = beta_0;
    beta_p_cache = beta_p;
    gamma_0_cache = gamma_0;
    gamma_p_cache = gamma_p;

    af0_1_cache = af0_1;
    af0_2_cache = af0_2;
    afplus_0_cache = afplus_0;
    afplus_1_cache = afplus_1;
    afplus_2_cache = afplus_2;
#if NBGL == 3
    af0_3_cache = af0_3;
    afplus_3_cache = afplus_3;
#endif    
    /* f+(q2=0) = f0(q2=0) */
    z0 = (sqrt(w0+1.)-sqrt(2.))/(sqrt(w0+1.)+sqrt(2.));
    if (CLNflag) af0_0 = 0.;
    else {
        af0_0 = fplus(0.);
#if NBGL == 3
        af0_0 -= (af0_1 * z0 + af0_2 * z0 * z0 + af0_3 * z0 * z0 *z0) / phi_f0(z0) / ((z0 - z0p_1) / (1. - z0 * z0p_1)*(z0 - z0p_2) / (1. - z0 * z0p_2));
#else        
        af0_0 -= (af0_1 * z0 + af0_2 * z0 * z0) / phi_f0(z0) / ((z0 - z0p_1) / (1. - z0 * z0p_1)*(z0 - z0p_2) / (1. - z0 * z0p_2));
#endif        
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
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    if (CLNflag) {
        return fplusz0 * N_0 * (1. - alpha_p*8. * rho1to2 * z + beta_p*(51. * rho1to2 - 10.) * z * z - gamma_p*(252. * rho1to2 - 84.) * z * z * z);
    } else {
        double P_fplus = (z - z1m_1) / (1. - z * z1m_1)*(z - z1m_2) / (1. - z * z1m_2)*(z - z1m_3) / (1. - z * z1m_3);
#if NBGL == 3        
        return (afplus_0 + afplus_1 * z + afplus_2 * z * z +  afplus_3 * z * z * z) / phi_fplus(z) / P_fplus;
#else
        return (afplus_0 + afplus_1 * z + afplus_2 * z * z) / phi_fplus(z) / P_fplus;
#endif        
    }
}

double MPlnu::f0(double q2)
{
    double w = w0 - q2 / (2. * MM * MP);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    if (CLNflag) {
        double prefac0 = 2. * (MP / MM) / (1. + MP / MM)/ (1. + MP / MM)*(1. + w0);
        double norm = prefac0 * (1. - alpha_0*0.0068 * (w0 - 1.) + beta_0*0.0017 * (w0 - 1.)*(w0 - 1.) - gamma_0*0.0013 * (w0 - 1.)*(w0 - 1.)*(w0 - 1.));
        double prefac = fplus(q2)* prefac0/(1. + w0)*(1. + w);
        // norm introduced to respect f+(q2=0)=f0(q2=0) exactly
        return prefac/norm * (1. - alpha_0*0.0068 * (w - 1.) + beta_0*0.0017 * (w - 1.)*(w - 1.) - gamma_0*0.0013 * (w - 1.)*(w - 1.)*(w - 1.));
    } else {
        double P_f0 = (z - z0p_1) / (1. - z * z0p_1)*(z - z0p_2) / (1. - z * z0p_2);
#if NBGL == 3
        return (af0_0 + af0_1 * z + af0_2 * z * z + af0_3 * z * z * z) / phi_f0(z) / P_f0;
#else        
        return (af0_0 + af0_1 * z + af0_2 * z * z) / phi_f0(z) / P_f0;
#endif        
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
            8. * M_SQRT2 *(Enu * Mlep - Elep * Mnu)*(HA(q2) * HT(q2).conjugate()).imag() +
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
    double coeff_f0 = (1.-MP*MP/MM/MM)*(1.-MP*MP/MM/MM)*3.*Mlep*Mlep/(2.*q2);
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
    
    double intJ1 = integrateJ(1, q2_min, q2_max);
    double intJ3 = integrateJ(3, q2_min, q2_max);

    return 2. * (intJ1 + intJ3 / 3.);
    
    // x-check of the SM computation
    //return integrateJ(4, q2_min, q2_max);
}

double MPlnu::get_unitarity_1min_BGL()
{
    updateParameters();

#if NBGL == 3
    return afplus_0 * afplus_0 + afplus_1 * afplus_1 + afplus_2 * afplus_2 + afplus_3 * afplus_3;
#else    
    return afplus_0 * afplus_0 + afplus_1 * afplus_1 + afplus_2 * afplus_2;
#endif    
}

double MPlnu::get_unitarity_0plus_BGL()
{
    updateParameters();

#if NBGL == 3
    return af0_0 * af0_0 + af0_1 * af0_1 + af0_2 * af0_2 + af0_3 * af0_3;
#else    
    return af0_0 * af0_0 + af0_1 * af0_1 + af0_2*af0_2;
#endif    
}

double MPlnu::get_strong_unitarity_BGL()
{
    updateParameters();

#if NBGL == 3
    return 1707.54 * afplus_0 * afplus_0 + 1299.57 * afplus_0 * afplus_1 + 442.82 * afplus_1 * afplus_1 - 356.01 * afplus_0 * afplus_2 
               - 101.62 * afplus_1 * afplus_2 + 34.947 * afplus_2 * afplus_2 - 206.767 * afplus_0 * afplus_3 - 127.668 * afplus_1 * afplus_3 
               + 33.234 * afplus_2 * afplus_3 + 16.475 * afplus_3 * afplus_3;
#else    
    return 442.82 * afplus_0 * afplus_0 - 101.619 * afplus_0 * afplus_1 + 34.947* afplus_1 * afplus_1 - 127.668 * afplus_0 * afplus_2 + 33.234 * afplus_1 * afplus_2 + 16.4754 * afplus_2 * afplus_2;
#endif    
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
