/*
 * Copyright (C) 2026 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6flavchain.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

const std::string NPSMEFTd6flavchain::NPSMEFTd6flavchainVars[NNPSMEFTd6flavchainVars] = {
    "MVLQ", "lamu1", "lamu2", "lamd1", "lamd2", "lamd3", "eps", 
    "x1u", "x2u", "x1d", "x2d", "x3d", "x4d", "x5d", "x6d", "phi", "Ytop"
};

// BEGIN-CHAIN-MATCHING (self-contained: only gslpp and <cmath> are used)
namespace {

    typedef gslpp::complex cplx;

    /**
     * Exact tree-level matching of the explicit chain model (7 U + 9 D VLQs), all orders in eps.
     * Output in the chain (= UP) basis, doublet index first:
     *   Yu, Yd : dimension-four Yukawas (Yd with Cholesky right-handed normalization)
     *   Cu, Cd : M^2 Chat_U, M^2 Chat_D (dimensionless, hermitian)
     * Notation and equations as in Sections 5-6 of the chain_flavor_matching notes.
     */
    void chainMatching(double e, double lu1, double lu2, double ld1, double ld2, double ld3,
            double xu1, double xu2, double x1, double x2, double x3, double x4, double x5, double x6,
            double phi, double yt,
            gslpp::matrix<cplx>& Yu, gslpp::matrix<cplx>& Yd, gslpp::matrix<cplx>& Cu, gslpp::matrix<cplx>& Cd)
    {
        const double e2 = e * e, e3 = e2 * e, e4 = e2 * e2, e6 = e4 * e2;
        const cplx P(cos(phi), sin(phi)); // e^{i phi}

        // ---- Step 1: chain constants ----
        auto S = [e2](int n) { double s = 0., t = 1.; for (int k = 0; k < n; k++) { s += t; t *= e2; } return s; };
        auto T = [&S](int n) { double s = 0.; for (int c = 1; c <= n; c++) s += S(c); return s; };
        auto R = [&S, e2](int n) { double s = 0.; for (int c = 1; c <= n; c++) s += pow(e2, n - c) * S(c) * S(c); return s; };
        const double S2 = S(2), S3 = S(3), S4 = S(4), S5 = S(5);

        // ---- Up sector: everything diagonal, exact closed forms (Eqs. Yuexpl, ChatU) ----
        auto cn = [&](double lam, double x, int n) {
            const double z = x * x * e2 * S(n), e2n = pow(e2, n);
            return lam * lam * (S(n) - 2. * x * x * e2n * T(n) / (1. + z)
                    + x * x * x * x * e2n * e2 * R(n) / ((1. + z) * (1. + z)));
        };
        Yu = gslpp::matrix<cplx>(3, 3, 0.);
        Cu = gslpp::matrix<cplx>(3, 3, 0.);
        Yu.assign(0, 0, -xu1 * lu1 * e4 * e / sqrt(1. + xu1 * xu1 * e2 * S5));
        Yu.assign(1, 1, xu2 * lu2 * e2 / sqrt(1. + xu2 * xu2 * e2 * S2));
        Yu.assign(2, 2, yt);
        Cu.assign(0, 0, cn(lu1, xu1, 5));
        Cu.assign(1, 1, cn(lu2, xu2, 2));

        // ---- Down sector ----
        const double T4 = 1. + S2 + S3 + S4, k3 = 1. + S2 + (1. + e2) * S3;
        const double R4 = S4 * S4 + e2 * S3 * S3 + e4 * S2 * S2 + e6;
        const double h33 = (1. + e2) * S3 * S3 + e2 * S2 * S2 + e4;
        const double h34 = -e * (S3 * S4 + S3 * S3 + e2 * S2 * S2 + e4);

        // Step 2: UV blocks (M^2 units for A, B, C)
        cplx Y[3][3], A[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { Y[i][j] = 0.; A[i][j] = 0.; }
        Y[0][1] = ld1 * x1 * e4;
        Y[0][2] = P * (ld1 * x2 * e4);
        Y[1][0] = ld2 * x5 * e4;
        Y[1][1] = -ld2 * x3 * e3;
        Y[1][2] = -ld2 * x4 * e3;
        Y[2][2] = -ld3 * x6 * e;
        A[0][0] = e2 * x5 * x5 * R4;
        A[0][1] = e2 * x3 * x5 * h34;
        A[0][2] = e2 * x4 * x5 * h34;
        A[1][1] = e2 * (x1 * x1 * R4 + x3 * x3 * h33);
        A[1][2] = e2 * (P * (x1 * x2 * R4) + x3 * x4 * h33);
        A[2][2] = e2 * (x2 * x2 * R4 + x4 * x4 * h33 + x6 * x6);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < i; j++) A[i][j] = A[j][i].conjugate();
        const double C11 = ld1 * ld1 * S4, C22 = ld2 * ld2 * S4, C33 = ld3 * ld3;

        // Step 3: Cholesky factor R of 1+Z (R^dagger R = 1+Z) and its inverse r
        const double D = 1. + e2 * x5 * x5 * S4;
        const double R11 = sqrt(D);
        const double R12 = -e3 * x3 * x5 * S3 / R11;
        const double R13 = -e3 * x4 * x5 * S3 / R11;
        const double R22 = sqrt(1. + e2 * (x1 * x1 * S4 + x3 * x3 * S3) - e6 * x3 * x3 * x5 * x5 * S3 * S3 / D);
        const cplx R23 = (P * (e2 * x1 * x2 * S4) + e2 * x3 * x4 * S3 - e6 * x3 * x4 * x5 * x5 * S3 * S3 / D) / R22;
        const double R33 = sqrt(1. + e2 * (x2 * x2 * S4 + x4 * x4 * S3 + x6 * x6)
                - e6 * x4 * x4 * x5 * x5 * S3 * S3 / D - R23.abs2());
        cplx r[3][3];
        r[0][0] = 1. / R11;
        r[0][1] = -R12 / (R11 * R22);
        r[0][2] = (R12 * R23 - R13 * R22) / (R11 * R22 * R33);
        r[1][0] = 0.;
        r[1][1] = 1. / R22;
        r[1][2] = -R23 / (R22 * R33);
        r[2][0] = 0.;
        r[2][1] = 0.;
        r[2][2] = 1. / R33;

        // Step 4: canonical blocks Yhat = Y r, Bhat = B r, Ahat = r^dagger A r
        // B has the zero pattern of Y, weighted by the exit site: T4 (site 4), k3 (site 3), 1 (one-site chain)
        cplx B[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) B[i][j] = 0.;
        B[0][1] = -T4 * Y[0][1];
        B[0][2] = -T4 * Y[0][2];
        B[1][0] = -T4 * Y[1][0];
        B[1][1] = -k3 * Y[1][1];
        B[1][2] = -k3 * Y[1][2];
        B[2][2] = -Y[2][2];
        cplx yh[3][3], bh[3][3], ah[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                yh[i][j] = 0.;
                bh[i][j] = 0.;
                for (int k = 0; k <= j; k++) {
                    yh[i][j] += Y[i][k] * r[k][j];
                    bh[i][j] += B[i][k] * r[k][j];
                }
            }
        for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) {
                ah[k][l] = 0.;
                for (int m = 0; m <= k; m++)
                    for (int n = 0; n <= l; n++)
                        ah[k][l] += r[m][k].conjugate() * A[m][n] * r[n][l];
            }

        // Step 5: Chat_D = C + Bhat Yhat^dagger + Yhat Bhat^dagger + Yhat Ahat Yhat^dagger
        cplx Ch[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Ch[i][j] = 0.;
                for (int k = 0; k < 3; k++) {
                    Ch[i][j] += bh[i][k] * yh[j][k].conjugate() + yh[i][k] * bh[j][k].conjugate();
                    for (int l = 0; l < 3; l++)
                        Ch[i][j] += yh[i][k] * ah[k][l] * yh[j][l].conjugate();
                }
            }
        Ch[0][0] += C11;
        Ch[1][1] += C22;
        Ch[2][2] += C33;

        Yd = gslpp::matrix<cplx>(3, 3, 0.);
        Cd = gslpp::matrix<cplx>(3, 3, 0.);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Yd.assign(i, j, yh[i][j]);
                // enforce exact hermiticity
                Cd.assign(i, j, i == j ? cplx(Ch[i][i].real(), 0.) : (i < j ? Ch[i][j] : Ch[j][i].conjugate()));
            }
    }

} // namespace
// END-CHAIN-MATCHING

NPSMEFTd6flavchain::NPSMEFTd6flavchain()
: NPSMEFTd6General(),
{
    setModelName("NPSMEFTd6flavchain");

    // All the Wilson coefficients follow from the UV parameters: no *_LNP input is accepted
    restrictWilsonCoefficientInputs(NPSMEFTd6flavchainVars);

    ModelParamMap.insert(std::make_pair("MVLQ", std::cref(MVLQ)));
    ModelParamMap.insert(std::make_pair("lamu1", std::cref(lamu1)));
    ModelParamMap.insert(std::make_pair("lamu2", std::cref(lamu2)));
    ModelParamMap.insert(std::make_pair("lamd1", std::cref(lamd1)));
    ModelParamMap.insert(std::make_pair("lamd2", std::cref(lamd2)));
    ModelParamMap.insert(std::make_pair("lamd3", std::cref(lamd3)));
    ModelParamMap.insert(std::make_pair("eps", std::cref(eps)));
    ModelParamMap.insert(std::make_pair("x1u", std::cref(xu1)));
    ModelParamMap.insert(std::make_pair("x2u", std::cref(xu2)));
    ModelParamMap.insert(std::make_pair("x1d", std::cref(xd1)));
    ModelParamMap.insert(std::make_pair("x2d", std::cref(xd2)));
    ModelParamMap.insert(std::make_pair("x3d", std::cref(xd3)));
    ModelParamMap.insert(std::make_pair("x4d", std::cref(xd4)));
    ModelParamMap.insert(std::make_pair("x5d", std::cref(xd5)));
    ModelParamMap.insert(std::make_pair("x6d", std::cref(xd6)));
    ModelParamMap.insert(std::make_pair("phi", std::cref(phi)));
    ModelParamMap.insert(std::make_pair("Ytop", std::cref(Ytop)));
}

void NPSMEFTd6flavchain::setParameter(const std::string name, const double& value)
{
    if (name.compare("MVLQ") == 0) {
        MVLQ = value;
    } else if (name.compare("lamu1") == 0)
        lamu1 = value;
    else if (name.compare("lamu2") == 0)
        lamu2 = value;
    else if (name.compare("lamd1") == 0)
        lamd1 = value;
    else if (name.compare("lamd2") == 0)
        lamd2 = value;
    else if (name.compare("lamd3") == 0)
        lamd3 = value;
    else if (name.compare("eps") == 0)
        eps = value;
    else if (name.compare("x1u") == 0)
        xu1 = value;
    else if (name.compare("x2u") == 0)
        xu2 = value;
    else if (name.compare("x1d") == 0)
        xd1 = value;
    else if (name.compare("x2d") == 0)
        xd2 = value;
    else if (name.compare("x3d") == 0)
        xd3 = value;
    else if (name.compare("x4d") == 0)
        xd4 = value;
    else if (name.compare("x5d") == 0)
        xd5 = value;
    else if (name.compare("x6d") == 0)
        xd6 = value;
    else if (name.compare("phi") == 0)
        phi = value;
    else if (name.compare("Ytop") == 0)
        Ytop = value;
    else
        NPSMEFTd6General::setParameter(name, value);
}

bool NPSMEFTd6flavchain::setFlag(const std::string name, const bool value)
{
    if (name.compare("U2") == 0) {
        flag_U2 = value;
        return true;
    } else
        return NPSMEFTd6General::setFlag(name, value);
}

bool NPSMEFTd6flavchain::setFlagStr(const std::string name, const std::string value)
{
    if (name.compare("SMEFTBasisFlag") == 0 && value.compare("UP") != 0)
        throw std::runtime_error("NPSMEFTd6flavchain: the model is defined in the UP basis, SMEFTBasisFlag = "
            + value + " is not allowed");
    return NPSMEFTd6General::setFlagStr(name, value);
}

bool NPSMEFTd6flavchain::CheckParameters(const std::map<std::string, double>& DPars)
{
    std::vector<std::string> mandatory = {"MVLQ", "lamu1", "lamu2", "lamd1", "lamd2", "lamd3", "eps", 
        "x1u", "x2u", "x1d", "x2d", "x3d", "x4d", "x5d", "x6d", "phi", "Ytop"};


    for (const std::string& par : mandatory) {
        if (DPars.find(par) == DPars.end()) {
            std::cout << "ERROR: missing mandatory NPSMEFTd6flavchain parameter " << par << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(par);
        }
    }
    return (NPSMEFTd6General::CheckParameters(DPars));
}

void NPSMEFTd6flavchain::computeMatrices()
{
    chainMatching(eps, lamu1, lu2, lamd1, ld2, lamd3, xu1, xu2, xd1, xd2, xd3, xd4, xd5, xd6, phi, Ytop,
            Yuhat, Ydhat, Cuhat, Cdhat);
}



bool NPSMEFTd6flavchain::PostUpdate()
{
    if (SMEFTBasisFlag.compare("UP") != 0)
        throw std::runtime_error("NPSMEFTd6flavchain::PostUpdate(): SMEFTBasisFlag must be UP");
    if (!(MVLQ > 0.) || !(eps > 0.)) {
        std::cout << "ERROR: NPSMEFTd6flavchain requires MVLQ > 0 and eps > 0" << std::endl;
        return (false);
    }

    // SM initial conditions (Yukawas in the UP basis) at the matching scale (gauge and Higgs couplings). The quark Yukawas
    // generated here from the SM inputs are replaced below by those of the chain model.
    NPSMEFTd6General::GenerateSMInitialConditions();

    computeMatrices();

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            SMEFTEvolEW.setCoefficient("YuR", i, j, Yuhat(j, i).real());
            SMEFTEvolEW.setCoefficient("YuI", i, j,-Yuhat(j, i).imag());
            SMEFTEvolEW.setCoefficient("YdR", i, j, Ydhat(j, i).real());
            SMEFTEvolEW.setCoefficient("YdI", i, j,-Ydhat(j, i).imag());
        }

   // Coefficients in units of 1/Lambda_NP^2 with Lambda_NP = MVLQ: Cuhat, Cdhat = M^2 Chat
    const gslpp::matrix<gslpp::complex> CuH = Cuhat * Yuhat * 0.5;
    const gslpp::matrix<gslpp::complex> CdH = Cdhat * Ydhat * 0.5;

    // C_Hq^(1) = (Chat_U - Chat_D)/4, C_Hq^(3) = -(Chat_U + Chat_D)/4 (hermitian)
    CHq1_11r_LNP = 0.25 * (Cuhat(0, 0) - Cdhat(0, 0)).real();
    CHq1_22r_LNP = 0.25 * (Cuhat(1, 1) - Cdhat(1, 1)).real();
    CHq1_33r_LNP = 0.25 * (Cuhat(2, 2) - Cdhat(2, 2)).real();
    CHq1_12r_LNP = 0.25 * (Cuhat(0, 1) - Cdhat(0, 1)).real();
    CHq1_12i_LNP = 0.25 * (Cuhat(0, 1) - Cdhat(0, 1)).imag();
    CHq1_13r_LNP = 0.25 * (Cuhat(0, 2) - Cdhat(0, 2)).real();
    CHq1_13i_LNP = 0.25 * (Cuhat(0, 2) - Cdhat(0, 2)).imag();
    CHq1_23r_LNP = 0.25 * (Cuhat(1, 2) - Cdhat(1, 2)).real();
    CHq1_23i_LNP = 0.25 * (Cuhat(1, 2) - Cdhat(1, 2)).imag();

    CHq3_11r_LNP = -0.25 * (Cuhat(0, 0) + Cdhat(0, 0)).real();
    CHq3_22r_LNP = -0.25 * (Cuhat(1, 1) + Cdhat(1, 1)).real();
    CHq3_33r_LNP = -0.25 * (Cuhat(2, 2) + Cdhat(2, 2)).real();
    CHq3_12r_LNP = -0.25 * (Cuhat(0, 1) + Cdhat(0, 1)).real();
    CHq3_12i_LNP = -0.25 * (Cuhat(0, 1) + Cdhat(0, 1)).imag();
    CHq3_13r_LNP = -0.25 * (Cuhat(0, 2) + Cdhat(0, 2)).real();
    CHq3_13i_LNP = -0.25 * (Cuhat(0, 2) + Cdhat(0, 2)).imag();
    CHq3_23r_LNP = -0.25 * (Cuhat(1, 2) + Cdhat(1, 2)).real();
    CHq3_23i_LNP = -0.25 * (Cuhat(1, 2) + Cdhat(1, 2)).imag();

    // C_uH = Chat_U Yu^dagger / 2: (C_uH)_ij = CU_i (Yu_ji)^* / 2, diagonal in the UP basis
    CuH_11r_LNP = CuH(0, 0).real();
    CuH_11i_LNP = CuH(0, 0).imag();
    CuH_12r_LNP = CuH(0, 1).real();
    CuH_12i_LNP = CuH(0, 1).imag();
    CuH_13r_LNP = CuH(0, 2).real();
    CuH_13i_LNP = CuH(0, 2).imag();
    CuH_21r_LNP = CuH(1, 0).real();
    CuH_21i_LNP = CuH(1, 0).imag();
    CuH_22r_LNP = CuH(1, 1).real();
    CuH_22i_LNP = CuH(1, 1).imag();
    CuH_23r_LNP = CuH(1, 2).real();
    CuH_23i_LNP = CuH(1, 2).imag();
    CuH_31r_LNP = CuH(2, 0).real();
    CuH_31i_LNP = CuH(2, 0).imag();
    CuH_32r_LNP = CuH(2, 1).real();
    CuH_32i_LNP = CuH(2, 1).imag();
    CuH_33r_LNP = CuH(2, 2).real();
    CuH_33i_LNP = CuH(2, 2).imag();

    CdH_11r_LNP = CdH(0, 0).real();
    CdH_11i_LNP = CdH(0, 0).imag();
    CdH_12r_LNP = CdH(0, 1).real();
    CdH_12i_LNP = CdH(0, 1).imag();
    CdH_13r_LNP = CdH(0, 2).real();
    CdH_13i_LNP = CdH(0, 2).imag();
    CdH_21r_LNP = CdH(1, 0).real();
    CdH_21i_LNP = CdH(1, 0).imag();
    CdH_22r_LNP = CdH(1, 1).real();
    CdH_22i_LNP = CdH(1, 1).imag();
    CdH_23r_LNP = CdH(1, 2).real();
    CdH_23i_LNP = CdH(1, 2).imag();
    CdH_31r_LNP = CdH(2, 0).real();
    CdH_31i_LNP = CdH(2, 0).imag();
    CdH_32r_LNP = CdH(2, 1).real();
    CdH_32i_LNP = CdH(2, 1).imag();
    CdH_33r_LNP = CdH(2, 2).real();
    CdH_33i_LNP = CdH(2, 2).imag();

    // NPSMEFTd6General::PostUpdate() does not regenerate the SM initial conditions (already computed),
    // runs to the EW scale and sets quark masses and CKM from the diagonalized mass matrices.
 
    if (!NPSMEFTd6General::PostUpdate()) return (false);

    return (true);
}
