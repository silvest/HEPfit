/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "MVlnu.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <boost/bind.hpp>
#include <limits>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_expint.h>
#include <limits>

MVlnu::MVlnu(const StandardModel& SM_i, QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i)
: mySM(SM_i)
{
    lep = lep_i;
    meson = meson_i;
    vectorM = vector_i;
    CLNflag = false;
    btocNPpmflag = false;
    
    w_J = gsl_integration_cquad_workspace_alloc (100); 
    
    checkcache_int_tau = 0;
    checkcache_int_mu = 0;
    checkcache_int_el = 0;
    
    double max_double = std::numeric_limits<double>::max();
    
    hA1w1_cache = max_double;
    rho2_cache = max_double;
    R1w1_cache = max_double;
    R2w1_cache = max_double;
    N_A_cache = max_double;
    N_1_cache =max_double;
    N_2_cache = max_double;
    j_A_cache= max_double;
    j_0_cache = max_double;
    j_1_cache = max_double;
    j_2_cache = max_double;
    k_A_cache = max_double;
    k_0_cache = max_double;
    k_1_cache = max_double;
    k_2_cache = max_double;
    l_A_cache = max_double;

    af0_cache = max_double;
    af1_cache = max_double;
    af2_cache = max_double;
    ag0_cache = max_double;
    ag1_cache = max_double;
    ag2_cache = max_double;
    aF11_cache = max_double;
    aF12_cache = max_double;
    aF21_cache = max_double;
    aF22_cache = max_double;

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

MVlnu::~MVlnu() {
}

std::vector<std::string> MVlnu::initializeMVlnuParameters()
{
    CLNflag = mySM.getFlavour().getFlagCLN();
    btocNPpmflag = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0);
    NPanalysis = (mySM.getModelName().compare("RealWeakEFTCCPM") == 0 || mySM.getModelName().compare("RealWeakEFTCC") == 0);

    if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>()
        << "af0" << "af1" << "af2" << "ag0" << "ag1" << "ag2"
        << "aF11" << "aF12" << "aF21" << "aF22"
        << "mBcstV1" << "mBcstV2" << "mBcstV3" << "mBcstV4"
        << "mBcstA1" << "mBcstA2" << "mBcstA3" << "mBcstA4"
        << "mBcstP1" << "mBcstP2" << "mBcstP3"
        << "chiTV" << "chiTA" << "chiTP" << "nI";
    else {
        std::stringstream out;
        out << vectorM;
        throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }

    if (CLNflag) {
        mvlnuParameters.clear();
        if (vectorM == StandardModel::D_star_P) mvlnuParameters = make_vector<std::string>()
            << "hA1w1" << "rho2" << "R1w1" << "R2w1"
            << "N_A" << "N_1" << "N_2" << "j_A" << "j_0" << "j_1" << "j_2"
            << "k_A" << "k_0" << "k_1" << "k_2" << "l_A";
    }    

    mySM.initializeMeson(meson);
    mySM.initializeMeson(vectorM);
    return mvlnuParameters;
}

void MVlnu::updateParameters() 
{
    if (!mySM.getFlavour().getUpdateFlag(meson, vectorM, lep)) return;

    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
    Mnu = 0.; // neutrinos assumed to be massless
    MM = mySM.getMesons(meson).getMass();
    MV = mySM.getMesons(vectorM).getMass();
    width = mySM.getMesons(meson).computeWidth();
    w0 = (MM*MM+MV*MV)/(2.*MM*MV);
    RV = 2.*sqrt(MM*MV)/(MM+MV);
    mu_b = MM; // mySM.getMub();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass(); // add the PS b mass
    Vcb = mySM.getCKM().getV_cb(); // mySM.getOptionalParameter("AbsVcb");
    ale_mub = mySM.Ale(mu_b,FULLNLO);
    /* Amplitude propto 4*GF*Vij/sqrt(2) & kinematics requires 1/(2^9 pi^3 MB^3) */
    amplsq_factor = GF*GF*Vcb.abs2()/(64.*M_PI*M_PI*M_PI*MM*MM*MM);
    q2min = Mlep*Mlep;
    q2max = (MM-MV)*(MM-MV);
    
    MV_o_MM = MV / MM;
    sqrtMV_o_MM = sqrt(MV_o_MM);
    
    /* SM Wilson coefficients */
    CV_SM = 1./2.*(1.+ale_mub/M_PI*log(mySM.getMz()/mu_b));
    CV = CV_SM;
    CA = -CV_SM;
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
                CTp = mySM.getOptionalParameter("CT_NP");
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

    switch (vectorM) {
        case StandardModel::D_star_P:
            if (CLNflag) {
                hA1w1 = mySM.getOptionalParameter("hA1w1");
                rho2 = mySM.getOptionalParameter("rho2");
                R1w1 = mySM.getOptionalParameter("R1w1");
                R2w1 = mySM.getOptionalParameter("R2w1");
                N_A = mySM.getOptionalParameter("N_A");
                N_1 = mySM.getOptionalParameter("N_1");
                N_2 = mySM.getOptionalParameter("N_2"); 
                j_A = mySM.getOptionalParameter("j_A");
                j_0 = mySM.getOptionalParameter("j_0");
                j_1 = mySM.getOptionalParameter("j_1"); 
                j_2 = mySM.getOptionalParameter("j_2");
                k_A = mySM.getOptionalParameter("k_A");
                k_0 = mySM.getOptionalParameter("k_0");
                k_1 = mySM.getOptionalParameter("k_1");
                k_2 = mySM.getOptionalParameter("k_2");
                l_A = mySM.getOptionalParameter("l_A");
                af0 = 0.;
                af1 = 0.;
                af2 = 0.;
                ag0 = 0.;
                ag1 = 0.;
                ag2 = 0.;
                aF11 = 0.;
                aF12 = 0.;
                aF21 = 0.;
                aF22 = 0.;
                mBcstV1 = 0.;
                mBcstV2 = 0.;
                mBcstV3 = 0.;
                mBcstV4 = 0.;
                mBcstA1 = 0.;
                mBcstA2 = 0.;
                mBcstA3 = 0.;
                mBcstA4 = 0.;
                mBcstP1 = 0.;
                mBcstP2 = 0.;
                mBcstP3 = 0.;
                chiTV = 0.;
                chiTA = 0.;
                chiTP = 0.;
                nI = 0.;
            } else {
                hA1w1 = 0.;
                rho2 = 0.;
                R1w1 = 0.;
                R2w1 = 0.;
                N_A = 0.;
                N_1 = 0.;
                N_2 = 0.;
                j_A = 0.;
                j_0 = 0.;
                j_1 = 0.;
                j_2 = 0.;
                k_A = 0.;
                k_0 = 0.;
                k_1 = 0.;
                k_2 = 0.;
                l_A = 0.;
                af0 = mySM.getOptionalParameter("af0");
                af1 = mySM.getOptionalParameter("af1");
                af2 = mySM.getOptionalParameter("af2");
                ag0 = mySM.getOptionalParameter("ag0");
                ag1 = mySM.getOptionalParameter("ag1");
                ag2 = mySM.getOptionalParameter("ag2");
                aF11 = mySM.getOptionalParameter("aF11");
                aF12 = mySM.getOptionalParameter("aF12");
                aF21 = mySM.getOptionalParameter("aF21");
                aF22 = mySM.getOptionalParameter("aF22");
                mBcstV1 = mySM.getOptionalParameter("mBcstV1");
                mBcstV2 = mySM.getOptionalParameter("mBcstV2");
                mBcstV3 = mySM.getOptionalParameter("mBcstV3");
                mBcstV4 = mySM.getOptionalParameter("mBcstV4");
                mBcstA1 = mySM.getOptionalParameter("mBcstA1");
                mBcstA2 = mySM.getOptionalParameter("mBcstA2");
                mBcstA3 = mySM.getOptionalParameter("mBcstA3");
                mBcstA4 = mySM.getOptionalParameter("mBcstA4");
                mBcstP1 = mySM.getOptionalParameter("mBcstP1");
                mBcstP2 = mySM.getOptionalParameter("mBcstP2");
                mBcstP3 = mySM.getOptionalParameter("mBcstP3");
                chiTV = mySM.getOptionalParameter("chiTV");
                chiTA = mySM.getOptionalParameter("chiTA");
                chiTP = mySM.getOptionalParameter("chiTP");
                nI = mySM.getOptionalParameter("nI");
            }
            break;
        default:
            std::stringstream out;
            out << vectorM;
            throw std::runtime_error("MVlnu: vector " + out.str() + " not implemented");
    }
    
    zV1 = sqrt((MM+MV)*(MM+MV)-mBcstV1*mBcstV1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zV1 /= (sqrt((MM+MV)*(MM+MV)-mBcstV1*mBcstV1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zV2 = sqrt((MM+MV)*(MM+MV)-mBcstV2*mBcstV2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zV2 /= (sqrt((MM+MV)*(MM+MV)-mBcstV2*mBcstV2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zV3 = sqrt((MM+MV)*(MM+MV)-mBcstV3*mBcstV3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zV3 /= (sqrt((MM+MV)*(MM+MV)-mBcstV3*mBcstV3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zV4 = sqrt((MM+MV)*(MM+MV)-mBcstV4*mBcstV4)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zV4 /= (sqrt((MM+MV)*(MM+MV)-mBcstV4*mBcstV4)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));

    zA1 = sqrt((MM+MV)*(MM+MV)-mBcstA1*mBcstA1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zA1 /= (sqrt((MM+MV)*(MM+MV)-mBcstA1*mBcstA1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zA2 = sqrt((MM+MV)*(MM+MV)-mBcstA2*mBcstA2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zA2 /= (sqrt((MM+MV)*(MM+MV)-mBcstA2*mBcstA2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zA3 = sqrt((MM+MV)*(MM+MV)-mBcstA3*mBcstA3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zA3 /= (sqrt((MM+MV)*(MM+MV)-mBcstA3*mBcstA3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zA4 = sqrt((MM+MV)*(MM+MV)-mBcstA4*mBcstA4)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zA4 /= (sqrt((MM+MV)*(MM+MV)-mBcstA4*mBcstA4)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    
    zP1 = sqrt((MM+MV)*(MM+MV)-mBcstP1*mBcstP1)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zP1 /= (sqrt((MM+MV)*(MM+MV)-mBcstP1*mBcstP1)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zP2 = sqrt((MM+MV)*(MM+MV)-mBcstP2*mBcstP2)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zP2 /= (sqrt((MM+MV)*(MM+MV)-mBcstP2*mBcstP2)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));
    zP3 = sqrt((MM+MV)*(MM+MV)-mBcstP3*mBcstP3)-sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV));
    zP3 /= (sqrt((MM+MV)*(MM+MV)-mBcstP3*mBcstP3)+sqrt((MM+MV)*(MM+MV)-(MM-MV)*(MM-MV)));

    if ((hA1w1 != hA1w1_cache) || (rho2 != rho2_cache) || (R1w1 != R1w1_cache) || (R2w1 != R2w1_cache)
            || (N_A != N_A_cache) || (N_1 != N_1_cache) || (N_2 != N_2_cache)
            || (j_A != j_A_cache) || (j_0 != j_0_cache) || (j_1 != j_1_cache) || (j_2 != j_2_cache)
            || (k_A != k_A_cache) || (k_0 != k_0_cache) || (k_1 != k_1_cache) || (k_2 != k_2_cache)
            || (l_A != l_A_cache)
            || (af0 != af0_cache) || (af1 != af1_cache) || (af2 != af2_cache)
            || (ag0 != ag0_cache) || (ag1 != af1_cache) || (ag2 != af2_cache)
            || (aF11 != aF11_cache) || (aF12 != aF12_cache)
            || (aF21 != aF21_cache) || (aF22 != aF22_cache)
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
            cached_intJ1s_tau = integrateJ(1, q2min, q2max);
            cached_intJ1c_tau = integrateJ(2, q2min, q2max);
            cached_intJ2s_tau = integrateJ(3, q2min, q2max);
            cached_intJ2c_tau = integrateJ(4, q2min, q2max);
            cached_intJ3_tau = integrateJ(5, q2min, q2max);
            cached_intJ6s_tau = integrateJ(8, q2min, q2max);
            cached_intJ6c_tau = integrateJ(9, q2min, q2max);
            cached_intJ9_tau = integrateJ(12, q2min, q2max);
            cached_intJ4_tau = 0.;
            cached_intJ5_tau = 0.;
            cached_intJ7_tau = 0.;
            cached_intJ8_tau = 0.;
            // not needed at present
            /*
            cached_intJ4_tau = integrateJ(6,q2min,q2max);
            cached_intJ5_tau = integrateJ(7,q2min,q2max);
            cached_intJ7_tau = integrateJ(10,q2min,q2max);
            cached_intJ8_tau = integrateJ(11,q2min,q2max);
             */
            checkcache_int_tau = 1;
        }
        if (lep == StandardModel::MU) {
            cached_intJ1s_mu = integrateJ(1, q2min, q2max);
            cached_intJ1c_mu = integrateJ(2, q2min, q2max);
            cached_intJ2s_mu = integrateJ(3, q2min, q2max);
            cached_intJ2c_mu = integrateJ(4, q2min, q2max);
            cached_intJ3_mu = integrateJ(5, q2min, q2max);
            cached_intJ6s_mu = integrateJ(8, q2min, q2max);
            cached_intJ6c_mu = integrateJ(9, q2min, q2max);
            cached_intJ9_mu = integrateJ(12, q2min, q2max);
            cached_intJ4_mu = 0.;
            cached_intJ5_mu = 0.;
            cached_intJ7_mu = 0.;
            cached_intJ8_mu = 0.;
            // not needed at present
            /*
            cached_intJ4_mu = integrateJ(6,q2min,q2max);
            cached_intJ5_mu = integrateJ(7,q2min,q2max);
            cached_intJ7_mu = integrateJ(10,q2min,q2max);
            cached_intJ8_mu = integrateJ(11,q2min,q2max);
             */
            checkcache_int_mu = 1;
        }
        if (lep == StandardModel::ELECTRON) {
            cached_intJ1s_el = integrateJ(1, q2min, q2max);
            cached_intJ1c_el = integrateJ(2, q2min, q2max);
            cached_intJ2s_el = integrateJ(3, q2min, q2max);
            cached_intJ2c_el = integrateJ(4, q2min, q2max);
            cached_intJ3_el = integrateJ(5, q2min, q2max);
            cached_intJ6s_el = integrateJ(8, q2min, q2max);
            cached_intJ6c_el = integrateJ(9, q2min, q2max);
            cached_intJ9_el = integrateJ(12, q2min, q2max);
            cached_intJ4_el = 0.;
            cached_intJ5_el = 0.;
            cached_intJ7_el = 0.;
            cached_intJ8_el = 0.;
            // not needed at present
            /*
            cached_intJ4_el = integrateJ(6,q2min,q2max);
            cached_intJ5_el = integrateJ(7,q2min,q2max);
            cached_intJ7_el = integrateJ(10,q2min,q2max);
            cached_intJ8_el = integrateJ(11,q2min,q2max);
             */
            checkcache_int_el = 1;
        }
    }

    hA1w1_cache = hA1w1;
    rho2_cache = rho2;
    R1w1_cache = R1w1;
    R2w1_cache = R2w1;
    N_A_cache = N_A;
    N_1_cache =N_1;
    N_2_cache = N_2;
    j_A_cache= j_A;
    j_0_cache = j_0;
    j_1_cache = j_1;
    j_2_cache = j_2;
    k_A_cache = k_A;
    k_0_cache = k_0;
    k_1_cache = k_1;
    k_2_cache = k_2;
    l_A_cache = l_A;

    af0_cache = af0;
    af1_cache = af1;
    af2_cache = af2;
    ag0_cache = ag0;
    ag1_cache = ag1;
    ag2_cache = ag2;
    aF11_cache = aF11;
    aF12_cache = aF12;
    aF21_cache = aF21;
    aF22_cache = aF22;

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

    mySM.getFlavour().setUpdateFlag(meson, vectorM, lep, false);

    return;
    
}

/*******************************************************************************
 * Kinematic functions                                                          *
 * ****************************************************************************/

double MVlnu::lambda_half(double a, double b, double c) 
{   
    return sqrt(a*a+b*b+c*c-2.*(a*b+a*c+b*c));
}
/*******************************************************************************
 * Form factors                                                                *
 * ****************************************************************************/

double MVlnu::phi_f(double z)
{
    double prefac = 4. * (MV_o_MM) / MM / MM * sqrt(nI / (3. * M_PI * chiTA));
    double num = (1. + z) * sqrt((1. - z)*(1. - z)*(1. - z));
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::f_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double Pfacf = (z - zA1) / (1. - z * zA1)*(z - zA2) / (1. - z * zA2)*(z - zA3) / (1. - z * zA3)*(z - zA4) / (1. - z * zA4);
    double phif = phi_f(z);
    return (af0 + af1 * z + af2 * z * z) / phif / Pfacf;
}

double MVlnu::phi_g(double z)
{
    double prefac = sqrt(nI / (3. * M_PI * chiTV));
    double num = 16. * (MV_o_MM)*(MV_o_MM)*(1. + z)*(1. + z) / sqrt(1. - z);
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::g_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double Pfacg = (z - zV1) / (1. - z * zV1)*(z - zV2) / (1. - z * zV2)*(z - zV3) / (1. - z * zV3)*(z - zV4) / (1. - z * zV4);
    double phig = phi_g(z);
    return (ag0 + ag1 * z + ag2 * z * z) / phig / Pfacg;
}

double MVlnu::phi_F1(double z)
{
    double prefac = 4. * (MV_o_MM) / MM / MM / MM * sqrt(nI / (6. * M_PI * chiTA));
    double num = (1. + z) * sqrt((1. - z)*(1. - z)*(1. - z)*(1. - z)*(1. - z));
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den5 = den * den * den * den*den;
    return prefac * num / den5;
}

double MVlnu::F1_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double PfacF1 = (z - zA1) / (1. - z * zA1)*(z - zA2) / (1. - z * zA2)*(z - zA3) / (1. - z * zA3)*(z - zA4) / (1. - z * zA4);
    double phiF1 = phi_F1(z);
    double aF10 = (MM - MV)*(phi_F1(0.) / phi_f(0.)) * af0; // F1(z=0) = (MM-MV)*f(z=0)
    return (aF10 + aF11 * z + aF12 * z * z) / phiF1 / PfacF1;
}

double MVlnu::phi_F2(double z)
{
    double prefac = 8. * sqrt(2.)*(MV_o_MM)*(MV_o_MM) * sqrt(nI / (M_PI * chiTP));
    double num = (1. + z)*(1. + z) / sqrt(1. - z);
    double den = (1. + MV_o_MM)*(1. - z) + 2. * sqrtMV_o_MM * (1. + z);
    double den4 = den * den * den*den;
    return prefac * num / den4;
}

double MVlnu::F2_BGL(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    double z0 = (sqrt(w0 + 1.) - M_SQRT2) / (sqrt(w0 + 1.) + M_SQRT2);
    double PfacF2 = (z - zP1) / (1. - z * zP1)*(z - zP2) / (1. - z * zP2)*(z - zP3) / (1. - z * zP3);
    double PfacF2z0 = (z0 - zP1) / (1. - z0 * zP1)*(z0 - zP2) / (1. - z0 * zP2)*(z0 - zP3) / (1. - z0 * zP3);
    double phiF2 = phi_F2(z);
    double phiF2z0 = phi_F2(z0);
    double aF20 = PfacF2z0 * phiF2z0 * 2. * F1_BGL(0.) / (MM * MM - MV * MV) - aF21 * z0 - aF22 * z0*z0; // F2(q2=0) = 2.*F1(q2=0)/(MM*MM-MV*MV)
    return (aF20 + aF21 * z + aF22 * z * z) / phiF2 / PfacF2;
}

double MVlnu::hA1(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    double z = (sqrt(w + 1.) - M_SQRT2) / (sqrt(w + 1.) + M_SQRT2);
    if (CLNflag) return hA1w1 * N_A * (1. - j_A * 8. * rho2 * z + k_A * (53. * rho2 - 15.) * z * z - l_A * (231. * rho2 - 91.) * z * z * z);
    else return f_BGL(q2) / sqrt(MM * MV) / (1. + w);
}

double MVlnu::R1(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    return N_1 * R1w1 - j_1 * 0.12 * (w - 1.) + k_1 * 0.05 * (w - 1.)*(w - 1.);
}

double MVlnu::R2(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    return N_2 * R2w1 + j_2 * 0.11 * (w - 1.) - k_2 * 0.06 * (w - 1.)*(w - 1.);
}

double MVlnu::R0(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    /* form factor relation among A0, A1 and A2 at q2=0 */
    double R2q2at0 = R2(0.);
    double R0q2at0 = (MM + MV - (MM - MV) * R2q2at0) / (2. * MV);
    // caveat: HQET rel at the kinematic endpoint, q2 = 0 ...
    double R0w1 = R0q2at0 + j_0 * 0.11 * (w0 - 1.) - k_0 * 0.01 * (w0 - 1.)*(w0 - 1.);
    // one may consider "lattice" R0w1 = 1.14 +- O(10%) + consistency rel at q2 = 0 ...
    return R0w1 - j_0 * 0.11 * (w - 1.) + k_0 * 0.01 * (w - 1.)*(w - 1.);
}

double MVlnu::V(double q2)
{
    if (CLNflag) return R1(q2) / RV * hA1(q2);
    else return (MM + MV) * g_BGL(q2) / 2.;
}

double MVlnu::A0(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    if (CLNflag) return R0(q2) / RV * hA1(q2);
    else return F2_BGL(q2) / RV / (1. + w);
}

double MVlnu::A1(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    return (w + 1.)*RV / 2. * hA1(q2);
}

double MVlnu::A2(double q2)
{
    double w = w0 - q2 / (2. * MM * MV);
    if (CLNflag) return R2(q2) / RV * hA1(q2);
    else return (MM + MV) / 2. / (w * w - 1.) / MM / MV * ((w - MV_o_MM) * f_BGL(q2) - F1_BGL(q2) / MM);
}

double MVlnu::T1(double q2)
{
    double delta_T1 = 0.;
    return (Mb + Mc) / (MM + MV) * V(q2)*(1. + delta_T1);
}

double MVlnu::T2(double q2)
{
    double delta_T2 = 0.;
    return (Mb - Mc) / (MM - MV) * A1(q2)*(1. + delta_T2);
}

double MVlnu::A12(double q2)
{
    return (A1(q2)*(MM + MV)*(MM + MV)*(MM * MM - MV * MV - q2) -
            A2(q2)*(MM * MM * MM * MM + (MV * MV - q2)*(MV * MV - q2) -
            2. * MM * MM * (MV * MV + q2))) / (16. * MM * MV * MV * (MM + MV));
}

double MVlnu::T23(double q2)
{
    double delta_T23 = 0.;
    return ((Mb - Mc)*((MM - MV)*(MM - MV) - q2)*((MM + MV)*(MM + MV) - q2) * A0(q2) +
            8 * MM * MV * (MV * MV - MM * MM) * A12(q2)) / (4. * MM * (MV - MM) * MV * q2)*(1. + delta_T23);
}
/********************************************************************************
 * Helicity amplitudes  (normalization such that all H \propto (mass scale)^-1) *
 * *****************************************************************************/

gslpp::complex MVlnu::HV0(double q2)
{
    return 4. * gslpp::complex::i() * MM * MV / (sqrt(q2)*(MM + MV))*((CV - CVp)*(MM + MV) * A12(q2) + Mb * (C7 - C7p) * T23(q2));
}

gslpp::complex MVlnu::HVp(double q2)
{
    return gslpp::complex::i()*((((CV + CVp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CV - CVp) * A1(q2))) / (2. * (MM + MV))
            + (Mb / q2)*((C7 + C7p) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(C7 - C7p)*(MM * MM - MV * MV) * T2(q2)));
}

gslpp::complex MVlnu::HVm(double q2)
{
    return gslpp::complex::i()*(((-(CV + CVp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CV - CVp) * A1(q2))) / (2. * (MM + MV))
            + (Mb / q2)*(-(C7 + C7p) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(C7 - C7p)*(MM * MM - MV * MV) * T2(q2)));
}

gslpp::complex MVlnu::HAp(double q2)
{
    return gslpp::complex::i()*((CA + CAp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CA - CAp) * A1(q2)) / (2. * (MM + MV));
}

gslpp::complex MVlnu::HAm(double q2)
{
    return gslpp::complex::i()*(-(CA + CAp) * lambda_half(MM*MM, MV*MV, q2) * V(q2)-(MM + MV)*(MM + MV)*(CA - CAp) * A1(q2)) / (2. * (MM + MV));
}

gslpp::complex MVlnu::HA0(double q2)
{
    return 4. * gslpp::complex::i() * MV * MM / (sqrt(q2))*(CA - CAp) * A12(q2);
}

gslpp::complex MVlnu::HP(double q2)
{
    return gslpp::complex::i() * lambda_half(MM*MM, MV*MV, q2) / 2. * ((CP - CPp) / (Mb + Mc)+(Mlep + Mnu) / q2 * (CA - CAp)) * A0(q2);
}

gslpp::complex MVlnu::HS(double q2)
{
    return gslpp::complex::i() * lambda_half(MM*MM, MV*MV, q2) / 2. * ((CS - CSp) / (Mb + Mc)+(Mlep - Mnu) / q2 * (CV - CVp)) * A0(q2);
}

gslpp::complex MVlnu::HT0(double q2)
{
    return 2. * M_SQRT2 * (MM * MV) / (MM + MV)*(CT + CTp) * T23(q2);
}

gslpp::complex MVlnu::HT0t(double q2)
{
    return 2. * (MM * MV) / (MM + MV)*(CT - CTp) * T23(q2);
}

gslpp::complex MVlnu::HTp(double q2)
{
    return ((CT - CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT + CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}

gslpp::complex MVlnu::HTpt(double q2)
{
    return ((CT + CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT - CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}

gslpp::complex MVlnu::HTm(double q2)
{
    return (-(CT - CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT + CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}

gslpp::complex MVlnu::HTmt(double q2)
{
    return (-(CT + CTp) * lambda_half(MM*MM, MV*MV, q2) * T1(q2)-(CT - CTp)*(MM * MM - MV * MV) * T2(q2)) / (sqrt(2. * q2));
}
/*******************************************************************************
 * Generalized angular coefficients  (see 1506.03970)                          *
 * ****************************************************************************/

gslpp::complex MVlnu::G000(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep * Mlep + lambda_lep2 / (4. * q2));
    double Enu = sqrt(Mnu * Mnu + lambda_lep2 / (4. * q2));
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * (4. / 9. * (3. * Elep * Enu + lambda_lep2 / (4. * q2))*(HVp(q2).abs2() + HVm(q2).abs2() + HV0(q2).abs2() + HAp(q2).abs2() + HAm(q2).abs2() + HA0(q2).abs2()) +
            4. * Mlep * Mnu / 3. * (HVp(q2).abs2() + HVm(q2).abs2() + HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() - HA0(q2).abs2()) +
            4. / 3. * ((Elep * Enu - Mlep * Mnu + lambda_lep2 / (4. * q2)) * HS(q2).abs2()+(Elep * Enu + Mlep * Mnu + lambda_lep2 / (4. * q2)) * HP(q2).abs2()) +
            16. / 9. * (3. * (Elep * Enu + Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTpt(q2).abs2() + HTmt(q2).abs2() + HT0t(q2).abs2()) +
            8. / 9. * (3. * (Elep * Enu - Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTpt(q2).abs2() + HTmt(q2).abs2() + HT0t(q2).abs2()) +
            16. / 3. * (Mlep * Enu + Mnu * Elep)*(HVp(q2) * HTpt(q2).conjugate() + HVm(q2) * HTmt(q2).conjugate() + HV0(q2) * HT0t(q2).conjugate()).imag() +
            8. * M_SQRT2 / 3. * (Mlep * Enu - Mnu * Elep)*(HAp(q2) * HTp(q2).conjugate() + HAm(q2) * HTm(q2).conjugate() + HA0(q2) * HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G010(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / 3. * lambda_lep * ((HVp(q2) * HAp(q2).conjugate() - HVm(q2) * HAm(q2).conjugate()).real()+
            (2. * M_SQRT2) / q2 * (Mlep * Mlep - Mnu * Mnu)*(HTp(q2) * HTpt(q2).conjugate() - HTm(q2) * HTmt(q2).conjugate()).real() +
            2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HTpt(q2).conjugate() - HAm(q2) * HTmt(q2).conjugate()).imag() +
            M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HTp(q2).conjugate() - HVm(q2) * HTm(q2).conjugate()).imag()-
            (Mlep - Mnu) / sqrt(q2)*(HA0(q2) * HP(q2).conjugate()).real()-(Mlep + Mnu) / sqrt(q2)*(HV0(q2) * HS(q2).conjugate()).real()+
            (M_SQRT2 * HT0(q2) * HP(q2).conjugate() + 2. * HT0t(q2) * HS(q2).conjugate()).imag());

}

gslpp::complex MVlnu::G020(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 2. / 9. * lambda_lep2 / q2 * ((-HVp(q2).abs2() - HVm(q2).abs2() + 2. * HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() + 2. * HA0(q2).abs2()) -
            2. * (2. * HT0(q2).abs2() - HTp(q2).abs2() - HTm(q2).abs2()) - 4. * (2. * HT0t(q2).abs2() - HTpt(q2).abs2() - HTmt(q2).abs2()));
}

gslpp::complex MVlnu::G200(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Elep = sqrt(Mlep * Mlep + lambda_lep2 / (4. * q2));
    double Enu = sqrt(Mnu * Mnu + lambda_lep2 / (4. * q2));
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * (-4. / 9. * (3. * Elep * Enu + lambda_lep2 / (4. * q2))*(HVp(q2).abs2() + HVm(q2).abs2() - 2. * HV0(q2).abs2() + HAp(q2).abs2() + HAm(q2).abs2() - 2. * HA0(q2).abs2()) -
            4. / 3. * Mlep * Mnu * (HVp(q2).abs2() + HVm(q2).abs2() - 2. * HV0(q2).abs2() - HAp(q2).abs2() - HAm(q2).abs2() + 2. * HA0(q2).abs2()) +
            8. / 3. * (Elep * Enu - Mlep * Mnu + lambda_lep2 / (4. * q2)) * HS(q2).abs2() + 8. / 3. * (Elep * Enu + Mlep * Mnu + lambda_lep2 / (4. * q2)) * HP(q2).abs2() -
            16. / 9. * (3. * (Elep * Enu + Mlep * Mnu) - lambda_lep2 / (4. * q2))*(HTpt(q2).abs2() + HTmt(q2).abs2() - 2. * HT0t(q2).abs2()) - 8. / 9. * (3. * (Elep * Enu - Mlep * Mnu) - lambda_lep2 / (4. * q2))*
            (HTp(q2).abs2() + HTm(q2).abs2() - 2 * HT0(q2).abs2()) - 16. / 3. * (Mlep * Enu + Mnu * Elep)*(HVp(q2) * HTpt(q2).conjugate() + HVm(q2) * HTmt(q2).conjugate() - 2. * HV0(q2) * HT0t(q2).conjugate()).imag() -
            8. * M_SQRT2 / 3. * (Mlep * Enu - Mnu * Elep)*(HAp(q2) * HTp(q2).conjugate() + HAm(q2) * HTm(q2).conjugate() - 2. * HA0(q2) * HT0(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G210(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 4. * lambda_lep / 3. * ((HVp(q2) * HAp(q2).conjugate() - HVm(q2) * HAm(q2).conjugate()).real() +
            2. * M_SQRT2 * (Mlep * Mlep - Mnu * Mnu) / q2 * (HTp(q2) * HTpt(q2).conjugate() - HTm(q2) * HTmt(q2).conjugate()).real() +
            2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HTpt(q2).conjugate() - HAm(q2) * HTmt(q2).conjugate()).imag() +
            M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HTp(q2).conjugate() - HVm(q2) * HTm(q2).conjugate()).imag() +
            2. * (Mlep - Mnu) / sqrt(q2)*(HA0(q2) * HP(q2).conjugate()).real() + 2. * (Mlep + Mnu) / sqrt(q2)*(HV0(q2) * HS(q2).conjugate()).real() -
            2. * (M_SQRT2 * HT0(q2) * HP(q2).conjugate() + 2. * HT0t(q2) * HS(q2).conjugate()).imag());
}

gslpp::complex MVlnu::G220(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 2. / 9. * lambda_lep2 / q2 * (HVp(q2).abs2() + HVm(q2).abs2() + 4. * HV0(q2).abs2() + HAp(q2).abs2() +
            HAm(q2).abs2() + 4. * HA0(q2).abs2() - 2. * (HTp(q2).abs2() + HTm(q2).abs2() + 4. * HT0(q2).abs2()) -
            4. * (HTpt(q2).abs2() + HTmt(q2).abs2() + 4. * HT0t(q2).abs2()));
}

gslpp::complex MVlnu::G211(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / sqrt(3.) * lambda_lep * (HVp(q2) * HA0(q2).conjugate() + HAp(q2) * HV0(q2).conjugate() -
            HV0(q2) * HAm(q2).conjugate() - HA0(q2) * HVm(q2).conjugate()+(Mlep + Mnu) / sqrt(q2)*(HVp(q2) * HS(q2).conjugate() + HS(q2) * HVm(q2).conjugate()) -
            gslpp::complex::i() * M_SQRT2 * (HP(q2) * HTm(q2).conjugate() - HTp(q2) * HP(q2).conjugate() + M_SQRT2 * (HS(q2) * HTmt(q2).conjugate() - HTpt(q2) * HS(q2).conjugate()))+
            (Mlep - Mnu) / sqrt(q2)*(HAp(q2) * HP(q2).conjugate() + HP(q2) * HAm(q2).conjugate()) -
            gslpp::complex::i()*2. * (Mlep + Mnu) / sqrt(q2)*(HAp(q2) * HT0t(q2).conjugate() + HT0t(q2) * HAm(q2).conjugate() - HTpt(q2) * HA0(q2).conjugate() - HA0(q2) * HTmt(q2).conjugate()) -
            gslpp::complex::i() * M_SQRT2 * (Mlep - Mnu) / sqrt(q2)*(HVp(q2) * HT0(q2).conjugate() + HT0(q2) * HVm(q2).conjugate() - HTp(q2) * HV0(q2).conjugate() - HV0(q2) * HTm(q2).conjugate()) +
            2. * M_SQRT2 * (Mlep * Mlep - Mnu * Mnu) / q2 * (HTp(q2) * HT0t(q2).conjugate() + HTpt(q2) * HT0(q2).conjugate() - HT0(q2) * HTmt(q2).conjugate() - HT0t(q2) * HTm(q2).conjugate()));
}

gslpp::complex MVlnu::G221(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return Gprefactor * 4. / 3. * lambda_lep2 / q2 * (HVp(q2) * HV0(q2).conjugate() + HV0(q2) * HVm(q2).conjugate() +
            HAp(q2) * HA0(q2).conjugate() + HA0(q2) * HAm(q2).conjugate() - 2. * (HTp(q2) * HT0(q2).conjugate() +
            HT0(q2) * HTm(q2).conjugate() + 2. * (HTpt(q2) * HT0t(q2).conjugate() + HT0t(q2) * HTmt(q2).conjugate())));
}

gslpp::complex MVlnu::G222(double q2)
{
    double lambda_MM = lambda_half(MM*MM, MV*MV, q2);
    double lambda_lep = lambda_half(Mlep*Mlep, Mnu*Mnu, q2);
    double lambda_lep2 = lambda_lep*lambda_lep;
    double Gprefactor = lambda_MM * lambda_lep / q2;

    return -Gprefactor * 8. / 3. * lambda_lep2 / q2 * (HVp(q2) * HVm(q2).conjugate() + HAp(q2) * HAm(q2).conjugate() -
            2. * (HTp(q2) * HTm(q2).conjugate() + 2. * HTpt(q2) * HTmt(q2).conjugate()));
}
/***************************************************************************
 * 12 independent J angular coefficients  (see again 1506.03970)           *
 * ************************************************************************/

double MVlnu::J1s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (8. * G000(q2) + 2. * G020(q2) - 4. * G200(q2) - G220(q2)).real() / 3.;
}

double MVlnu::J1c(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (8. * G000(q2) + 2. * G020(q2) + 8. * G200(q2) + 2. * G220(q2)).real() / 3.;
}

double MVlnu::J2s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * G020(q2) - G220(q2)).real();
}

double MVlnu::J2c(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * (G020(q2) + G220(q2))).real();
}

double MVlnu::J3(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (G222(q2).real());
}

double MVlnu::J4(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (G221(q2).real());
}

double MVlnu::J5(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (2. * G211(q2).real() / sqrt(3.));
}

double MVlnu::J6s(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (4. * (2. * G010(q2) - G210(q2)).real() / 3.);
}

double MVlnu::J6c(double q2)
{
    if (q2 < Mlep * Mlep) return 0.;
    if (q2 > (MM - MV)*(MM - MV)) return 0.;
    return -amplsq_factor * (8. * (G010(q2) + G210(q2)).real() / 3.);
}

double MVlnu::J7(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (2. * sqrt(3.)*(G211(q2).imag()) / 3.);
}

double MVlnu::J8(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return amplsq_factor * (G221(q2).imag());
}

double MVlnu::J9(double q2)
{
    if ((q2 < Mlep * Mlep) or (q2 > (MM - MV)*(MM - MV))) return 0.;
    return -amplsq_factor * (G222(q2).imag());
}
/***************************************************************************
 * Integration of angular coefficients Js                                  *
 * ************************************************************************/

double MVlnu::integrateJ(int i, double q2_min, double q2_max)
{
    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 1:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1s_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J1s, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 2:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ1c_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J1c, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 3:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2s_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J2s, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 4:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ2c_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J2c, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 5:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ3_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J3, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 6:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ4_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J4, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 7:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ5_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J5, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 8:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6s_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J6s, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 9:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_mu;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ6c_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J6c, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            return J_res;
            break;
        case 10:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ7_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J7, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 11:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ8_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J8, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 12:
            if (lep == StandardModel::TAU) if ((checkcache_int_tau == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_tau;
            if (lep == StandardModel::MU) if ((checkcache_int_mu == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_mu;
            if (lep == StandardModel::ELECTRON) if ((checkcache_int_el == 1) && (q2_min == q2min) && (q2_max == q2max)) return cached_intJ9_el;
            FJ = convertToGslFunction(boost::bind(&MVlnu::J9, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        default:
            gsl_set_error_handler(old_handler);
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVlnu::integrateJ: index " + out.str() + " not implemented");
    }
}

double MVlnu::dGammadw(double q2)
{
    updateParameters();

    return 3. / 4. * (2. * J1s(q2) + J1c(q2)) - 1. / 4. * (2. * J2s(q2) + J2c(q2));
}

double MVlnu::getDeltaGammaDeltaw(double w_min, double w_max)
{
    updateParameters();

    double q2_min = (2. * MM * MV)*(w0 - w_max); // min is Mlep*Mlep;
    double q2_max = (2. * MM * MV)*(w0 - w_min); // max is (MM-MV)*(MM-MV);

    double intJ1s = integrateJ(1, q2_min, q2_max);
    double intJ1c = integrateJ(2, q2_min, q2_max);
    double intJ2s = integrateJ(3, q2_min, q2_max);
    double intJ2c = integrateJ(4, q2_min, q2_max);

    return 3. / 4. * (2. * intJ1s + intJ1c) - 1. / 4. * (2. * intJ2s + intJ2c);
}

double MVlnu::dGammadcldq2(double q2, double cl)
{
    updateParameters();

    return 3. / 8. * ((J1s(q2) + 2. * J1c(q2)) + cl * (J6s(q2) + 2. * J6c(q2))+(2. * cl * cl - 1.)*(J2s(q2) + 2. * J2c(q2)));
}

double MVlnu::dGammadcl(double cl)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ6s = integrateJ(8, q2min, q2max);
    double intJ6c = integrateJ(9, q2min, q2max);

    return 3. / 8. * ((intJ1s + 2. * intJ1c) + cl * (intJ6s + 2. * intJ6c)+(2. * cl * cl - 1.)*(intJ2s + 2. * intJ2c));
}

double MVlnu::getDeltaGammaDeltacl(double cl_min, double cl_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ6s = integrateJ(8, q2min, q2max);
    double intJ6c = integrateJ(9, q2min, q2max);

    return 3. / 8. * ((cl_max - cl_min)*(intJ1c + 2. * intJ1s)+
            (cl_max * cl_max - cl_min * cl_min) / 2. * (intJ6c + 2. * intJ6s)+
            (2. / 3. * (cl_max * cl_max * cl_max - cl_min * cl_min * cl_min)-(cl_max - cl_min))*(intJ2c + 2. * intJ2s));
}

double MVlnu::dGammadcVdq2(double q2, double cl)
{
    updateParameters();

    return 3. / 8. * ((J1s(q2) + 2. * J1c(q2)) + cl * (J6s(q2) + 2. * J6c(q2))+(2. * cl * cl - 1.)*(J2s(q2) + 2. * J2c(q2)));
}

double MVlnu::dGammadcV(double cV)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);

    return 3. / 8. * (cV * cV * (3. * intJ1c - intJ2c)+(1. - cV * cV)*(3. * intJ1s - intJ2s));
}

double MVlnu::getDeltaGammaDeltacV(double cV_min, double cV_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);

    return 3. / 8. * ((cV_max * cV_max * cV_max - cV_min * cV_min * cV_min) / 3. * (3. * intJ1c - intJ2c)+
            ((cV_max - cV_min)-(cV_max * cV_max * cV_max - cV_min * cV_min * cV_min) / 3.)*(3. * intJ1s - intJ2s));
}

double MVlnu::dGammadchidq2(double q2, double chi)
{
    updateParameters();

    return (3. * J1c(q2) + 6. * J1s(q2) - J2c(q2) - 2. * J2s(q2)) / 8. / M_PI +
            cos(2. * chi) / 2. / M_PI * J3(q2) + sin(2. * chi) / 2. / M_PI * J9(q2);
}

double MVlnu::dGammadchi(double chi)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ3 = integrateJ(5, q2min, q2max);
    double intJ9 = integrateJ(12, q2min, q2max);

    return ((3. * intJ1c + 6. * intJ1s - intJ2c - 2. * intJ2s) / 4. +
            cos(2. * chi) * intJ3 + sin(2. * chi) * intJ9) / 2. / M_PI;
}

double MVlnu::getDeltaGammaDeltachi(double chi_min, double chi_max)
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    double intJ3 = integrateJ(5, q2min, q2max);
    double intJ9 = integrateJ(12, q2min, q2max);

    return ((chi_max - chi_min)*(3. * intJ1c + 6. * intJ1s - intJ2c - 2. * intJ2s) / 4. +
            (sin(2. * chi_max) - sin(2. * chi_min)) / 2. * intJ3 -
            (cos(2. * chi_max) - cos(2. * chi_min)) / 2. * intJ9) / (2. * M_PI);
}

double MVlnu::getFL()
{
    updateParameters();

    double intJ1s = integrateJ(1, q2min, q2max);
    double intJ1c = integrateJ(2, q2min, q2max);
    double intJ2s = integrateJ(3, q2min, q2max);
    double intJ2c = integrateJ(4, q2min, q2max);
    
    double DeltaJL = (3. * intJ1c - intJ2c) / 4.;
    double DeltaJ = 3. / 4. * (2. * intJ1s + intJ1c) - 1. / 4. * (2. * intJ2s + intJ2c);
    return DeltaJL/DeltaJ;
    
}

double MVlnu::get_unitarity_V_BGL()
{
    updateParameters();

    return ag0 * ag0 + ag1 * ag1 + ag2*ag2;

}

double MVlnu::get_unitarity_A_BGL()
{
    updateParameters();

    double aF10 = (MM - MV)*(phi_F1(0.) / phi_f(0.)) * af0;
    return af0 * af0 + af1 * af1 + af2 * af2 + aF10 * aF10 + aF11 * aF11 + aF12*aF12;
}

double MVlnu::get_unitarity_P_BGL()
{
    updateParameters();

    double z0 = (sqrt(w0 + 1.) - M_SQRT2) / (sqrt(w0 + 1.) + M_SQRT2);
    double PfacF2z0 = (z0 - zP1) / (1. - z0 * zP1)*(z0 - zP2) / (1. - z0 * zP2)*(z0 - zP3) / (1. - z0 * zP3);
    double phiF2z0 = phi_F2(z0);
    double aF20 = PfacF2z0 * phiF2z0 * 2. * F1_BGL(0.) / (MM * MM - MV * MV) - aF21 * z0 - aF22 * z0*z0;
    return aF20 * aF20 + aF21 * aF21 + aF22*aF22;
}

double MVlnu::get_hA1w1()
{
    updateParameters();

    return hA1(q2max);
}

double MVlnu::get_hA1(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return hA1(q2);
}

double MVlnu::get_hA2(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return (hA1(q2) * (1. + w) - RV * (A0(q2) * (1. + MV_o_MM) - A2(q2) * (MV_o_MM - w))) /  (1. + MV_o_MM*MV_o_MM - 2.* MV_o_MM * w);
    // return (hA1(q2) * (1. + w) - A0(q2) * RV * (1. + MV_o_MM) - A2(q2) * RV * (w - MV_o_MM)) / (1. + MV_o_MM*MV_o_MM - 2.*MV_o_MM*w) ;
}

double MVlnu::get_hA3(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return A2(q2) * RV - MV_o_MM * get_hA2(w) ;
}

double MVlnu::get_hV(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    return V(q2) * RV;
}

double MVlnu::get_R1(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R1(q2);
    else return V(q2) * RV / hA1(q2);
}

double MVlnu::get_R2(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R2(q2);
    else return A2(q2) * RV / hA1(q2);
}

double MVlnu::get_R0(double w)
{
    updateParameters();
    double q2 = (2. * MM * MV)*(w0 - w);

    if (CLNflag) return R0(q2);
    else return A0(q2) * RV / hA1(q2);
}


/***************************************************************************
 * SM computation  ... lep hel asymmetry, see 1203.2654, 1707.09509        *
 * ************************************************************************/

double MVlnu::Hplus(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return (MM+MV)*A1(q2)-2.*MM/(MM+MV)*abs_p*V(q2);
}

double MVlnu::Hminus(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return (MM+MV)*A1(q2)+2.*MM/(MM+MV)*abs_p*V(q2);
}

double MVlnu::H0(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return ((MM*MM-MV*MV-q2)*(MM+MV)*A1(q2)-4.*MM*MM*abs_p*abs_p/(MM+MV)*A2(q2))/(2.*MV*sqrt(q2));
}

double MVlnu::H0t(double q2){
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    return 2.*MM*abs_p*A0(q2)/sqrt(q2);
}

double MVlnu::dGpdq2(double q2){
        
    updateParameters();
        
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    double lep_factor = (1.-Mlep*Mlep/q2)*(1.-Mlep*Mlep/q2)*Mlep*Mlep/(2.*q2);
    return 2./3.*amplsq_factor*abs_p*q2*lep_factor*(Hplus(q2)*Hplus(q2)+Hminus(q2)*Hminus(q2)+H0(q2)*H0(q2)
            +3.*H0t(q2)*H0t(q2));
}

double MVlnu::dGmdq2(double q2){
    
    updateParameters(); 
    
    double abs_p = lambda_half(MM*MM,MV*MV,q2);
    double lep_factor = (1.-Mlep*Mlep/q2)*(1.-Mlep*Mlep/q2);
    return 2./3.*amplsq_factor*abs_p*q2*lep_factor*(Hplus(q2)*Hplus(q2)+Hminus(q2)*Hminus(q2)+H0(q2)*H0(q2));
}

double MVlnu::integrateGpm(int i, double q2_min, double q2_max)
{
    old_handler = gsl_set_error_handler_off();

    switch (i) {
        case 1:
            FJ = convertToGslFunction(boost::bind(&MVlnu::dGpdq2, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        case 2:
            FJ = convertToGslFunction(boost::bind(&MVlnu::dGmdq2, &(*this), _1));
            if (gsl_integration_cquad(&FJ, q2_min, q2_max, 1.e-2, 1.e-1, w_J, &J_res, &J_err, NULL) != 0) std::numeric_limits<double>::quiet_NaN();
            gsl_set_error_handler(old_handler);
            return J_res;
            break;
        default:
            gsl_set_error_handler(old_handler);
            std::stringstream out;
            out << i;
            throw std::runtime_error("MVlnu::integrateGpm: index " + out.str() + " not implemented");
    }
}

double MVlnu::getPlep()
{
    updateParameters();
     
    double DeltaGammaPlus = integrateGpm(1,q2min,q2max);
    double DeltaGammaMinus = integrateGpm(2,q2min,q2max);
    
    return (DeltaGammaPlus-DeltaGammaMinus)/(DeltaGammaPlus+DeltaGammaMinus);
    
}