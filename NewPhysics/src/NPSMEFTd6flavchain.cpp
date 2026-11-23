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
    "MVLQ", "lamu1", "lamu2", "lamd1", "lamd2", "lamd3", "eps", "deltaU2u", "deltaU2d"
};

NPSMEFTd6flavchain::NPSMEFTd6flavchain()
: NPSMEFTd6General(),
  MVLQ(0.), lamu1(0.), lamu2(0.), lamd1(0.), lamd2(0.), lamd3(0.), eps(0.), deltaU2u(0.), deltaU2d(0.),
  lamu1_eff(0.), lamu2_eff(0.), lamd1_eff(0.), lamd2_eff(0.), lamd3_eff(0.),
  ChatU11(0.), ChatU22(0.), ChatD11(0.), ChatD22(0.), ChatD33(0.),
  flag_U2(false), flag_minimal(false)
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
    ModelParamMap.insert(std::make_pair("deltaU2u", std::cref(deltaU2u)));
    ModelParamMap.insert(std::make_pair("deltaU2d", std::cref(deltaU2d)));

    // Derived quantities, only to be used as observables
    ModelParamMap.insert(std::make_pair("lamu1_eff", std::cref(lamu1_eff)));
    ModelParamMap.insert(std::make_pair("lamu2_eff", std::cref(lamu2_eff)));
    ModelParamMap.insert(std::make_pair("lamd1_eff", std::cref(lamd1_eff)));
    ModelParamMap.insert(std::make_pair("lamd2_eff", std::cref(lamd2_eff)));
    ModelParamMap.insert(std::make_pair("lamd3_eff", std::cref(lamd3_eff)));
    ModelParamMap.insert(std::make_pair("ChatU11", std::cref(ChatU11)));
    ModelParamMap.insert(std::make_pair("ChatU22", std::cref(ChatU22)));
    ModelParamMap.insert(std::make_pair("ChatD11", std::cref(ChatD11)));
    ModelParamMap.insert(std::make_pair("ChatD22", std::cref(ChatD22)));
    ModelParamMap.insert(std::make_pair("ChatD33", std::cref(ChatD33)));
}

void NPSMEFTd6flavchain::setParameter(const std::string name, const double& value)
{
    if (name.compare("MVLQ") == 0) {
        MVLQ = value;
        Lambda_NP = value; // matching scale
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
    else if (name.compare("deltaU2u") == 0)
        deltaU2u = value;
    else if (name.compare("deltaU2d") == 0)
        deltaU2d = value;
    else if (name.compare("Lambda_NP") == 0) {
        // Lambda_NP is fixed to MVLQ, any input value is ignored
    } else
        NPSMEFTd6General::setParameter(name, value);
}

bool NPSMEFTd6flavchain::setFlag(const std::string name, const bool value)
{
    if (name.compare("U2") == 0) {
        flag_U2 = value;
        return true;
    } else if (name.compare("MinimalCouplings") == 0) {
        flag_minimal = value;
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
    std::vector<std::string> mandatory = {"MVLQ"};
    std::vector<std::string> ignored = {"Lambda_NP"};

    if (flag_minimal) {
        mandatory.push_back("eps");
        ignored.insert(ignored.end(), {"lamu1", "lamu2", "lamd1", "lamd2", "lamd3"});
    } else {
        mandatory.insert(mandatory.end(), {"lamu1", "lamd1", "lamd3"});
        if (flag_U2)
            ignored.insert(ignored.end(), {"lamu2", "lamd2"});
        else
            mandatory.insert(mandatory.end(), {"lamu2", "lamd2"});
        ignored.push_back("eps");
    }
    if (flag_U2)
        mandatory.insert(mandatory.end(), {"deltaU2u", "deltaU2d"});
    else
        ignored.insert(ignored.end(), {"deltaU2u", "deltaU2d"});

    for (const std::string& par : mandatory) {
        if (DPars.find(par) == DPars.end()) {
            std::cout << "ERROR: missing mandatory NPSMEFTd6flavchain parameter " << par << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(par);
        }
    }
    for (const std::string& par : ignored) {
        if (DPars.find(par) != DPars.end())
            std::cout << "WARNING: NPSMEFTd6flavchain parameter " << par
                << " is not used with the current flags (Lambda_NP is always set to MVLQ)" << std::endl;
    }

    return (NPSMEFTd6General::CheckParameters(DPars));
}

void NPSMEFTd6flavchain::computeCouplings(const gslpp::matrix<gslpp::complex>& Yu,
        const gslpp::matrix<gslpp::complex>& Yd)
{
    double lu1 = lamu1, lu2 = lamu2, ld1 = lamd1, ld2 = lamd2, ld3 = lamd3;

    if (flag_minimal) {
        // UP basis: Yu = diag(yu), Yd = diag(yd) V^dagger, hence Yd(j,i) = yd_j V_ij^*
        double yu[3], yd[3];
        for (int i = 0; i < 3; i++) {
            yu[i] = Yu(i, i).abs();
            double s = 0.;
            for (int j = 0; j < 3; j++)
                s += Yd(i, j).abs2();
            yd[i] = sqrt(s);
        }
        const double aVus = Yd(1, 0).abs() / yd[1];
        const double aVub = Yd(2, 0).abs() / yd[2];
        const double aVcb = Yd(2, 1).abs() / yd[2];

        // Leading-order minimal couplings (largest hopping coefficient of each chain set to one)
        const double e2 = eps * eps, e3 = e2 * eps, e4 = e3 * eps, e5 = e4 * eps;
        lu1 = yu[0] / e5;
        lu2 = yu[1] / e2;
        ld1 = std::max(aVus * yd[1], aVub * yd[2]) / e4;
        ld2 = std::max({yd[0] / (aVus * e4), yd[1] / e3, aVcb * yd[2] / e3});
        ld3 = yd[2] / eps;

        if (flag_U2) {
            lu1 = lu2 = std::max(lu1, lu2);
            ld1 = ld2 = std::max(ld1, ld2);
        }
    } else if (flag_U2) {
        lu2 = lu1;
        ld2 = ld1;
    }

    lamu1_eff = lu1;
    lamu2_eff = lu2;
    lamd1_eff = ld1;
    lamd2_eff = ld2;
    lamd3_eff = ld3;

    // M^2 Chat in the UP basis (diagonal), with the optional U(2) breaking of the (2,2) entries
    ChatU11 = lu1 * lu1;
    ChatU22 = lu2 * lu2 * (flag_U2 ? 1. - deltaU2u : 1.);
    ChatD11 = ld1 * ld1;
    ChatD22 = ld2 * ld2 * (flag_U2 ? 1. - deltaU2d : 1.);
    ChatD33 = ld3 * ld3;
}

void NPSMEFTd6flavchain::setChainWilsonCoefficients(const gslpp::matrix<gslpp::complex>& Yu,
        const gslpp::matrix<gslpp::complex>& Yd)
{
    // Coefficients in units of 1/Lambda_NP^2 with Lambda_NP = MVLQ, i.e. directly in terms of M^2 Chat
    const double CU[3] = {ChatU11, ChatU22, 0.};
    const double CD[3] = {ChatD11, ChatD22, ChatD33};

    // C_Hq^(1) = (Chat_U - Chat_D)/4, C_Hq^(3) = -(Chat_U + Chat_D)/4 (diagonal and real)
    CHq1_11r_LNP = 0.25 * (CU[0] - CD[0]);
    CHq1_22r_LNP = 0.25 * (CU[1] - CD[1]);
    CHq1_33r_LNP = 0.25 * (CU[2] - CD[2]);
    CHq3_11r_LNP = -0.25 * (CU[0] + CD[0]);
    CHq3_22r_LNP = -0.25 * (CU[1] + CD[1]);
    CHq3_33r_LNP = -0.25 * (CU[2] + CD[2]);

    // C_uH = Chat_U Yu^dagger / 2: (C_uH)_ij = CU_i (Yu_ji)^* / 2, diagonal in the UP basis
    CuH_11r_LNP = 0.5 * CU[0] * Yu(0, 0).real();
    CuH_22r_LNP = 0.5 * CU[1] * Yu(1, 1).real();
    // CuH_33 vanishes since the top quark is not generated by a chain

    // C_dH = Chat_D Yd^dagger / 2: (C_dH)_ij = CD_i (Yd_ji)^* / 2, a full complex matrix
    double* CdHr[3][3] = {
        {&CdH_11r_LNP, &CdH_12r_LNP, &CdH_13r_LNP},
        {&CdH_21r_LNP, &CdH_22r_LNP, &CdH_23r_LNP},
        {&CdH_31r_LNP, &CdH_32r_LNP, &CdH_33r_LNP}
    };
    double* CdHi[3][3] = {
        {&CdH_11i_LNP, &CdH_12i_LNP, &CdH_13i_LNP},
        {&CdH_21i_LNP, &CdH_22i_LNP, &CdH_23i_LNP},
        {&CdH_31i_LNP, &CdH_32i_LNP, &CdH_33i_LNP}
    };
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            *CdHr[i][j] = 0.5 * CD[i] * Yd(j, i).real();
            *CdHi[i][j] = -0.5 * CD[i] * Yd(j, i).imag();
        }
}

bool NPSMEFTd6flavchain::PostUpdate()
{
    if (SMEFTBasisFlag.compare("UP") != 0)
        throw std::runtime_error("NPSMEFTd6flavchain::PostUpdate(): SMEFTBasisFlag must be UP");

    if (flag_minimal && !(eps > 0.)) {
        std::cout << "ERROR: NPSMEFTd6flavchain requires eps > 0 with MinimalCouplings = true" << std::endl;
        return (false);
    }

    Lambda_NP = MVLQ;

    // SM initial conditions (Yukawas in the UP basis) at the matching scale
    NPSMEFTd6General::GenerateSMInitialConditions();

    gslpp::matrix<gslpp::complex> Yu(3, 3, 0.), Yd(3, 3, 0.);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            Yu.assign(i, j, gslpp::complex(getRawSMEFTCoeffEW("YuR", i, j), getRawSMEFTCoeffEW("YuI", i, j)));
            Yd.assign(i, j, gslpp::complex(getRawSMEFTCoeffEW("YdR", i, j), getRawSMEFTCoeffEW("YdI", i, j)));
        }

    computeCouplings(Yu, Yd);
    setChainWilsonCoefficients(Yu, Yd);

    if (!NPSMEFTd6General::PostUpdate()) return (false);

    return (true);
}
