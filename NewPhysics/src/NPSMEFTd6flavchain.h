/*
 * Copyright (C) 2026 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6FLAVCHAIN_H
#define NPSMEFTD6FLAVCHAIN_H

#include "NPSMEFTd6General.h"

/**
 * @class NPSMEFTd6flavchain
 * @ingroup NewPhysics
 * @brief Chain model of quark flavour (SU(2)-singlet vector-like quarks, arXiv:2602.17754),
 * matched at tree level onto the SMEFT, exactly in @f$\epsilon=\mu/M@f$ (dimension-six truncation only).
 * @details The UV parameters are sampled directly. The VLQ mass is @f$M = M_{\rm VLQ}\,\Lambda_{\rm NP}@f$,
 * with MVLQ a dimensionless input and Lambda_NP a fixed input that sets the starting scale of the RG
 * evolution. At @f$\Lambda_{\rm NP}@f$ the model
 * - replaces the SMEFT Yukawa matrices of RGESolver with those of the chain model, so that quark
 *   masses and CKM matrix are predictions (NPSMEFTd6General::PostUpdate() rediagonalizes the
 *   quark mass matrices at the electroweak scale and sets masses and CKM from them);
 * - sets the non-vanishing dimension-six coefficients
 * @f[
 * C_{Hq}^{(1)} = \tfrac14(\hat C_U-\hat C_D),\quad
 * C_{Hq}^{(3)} = -\tfrac14(\hat C_U+\hat C_D),\quad
 * C_{uH} = \tfrac12\hat C_U\hat Y_u,\quad
 * C_{dH} = \tfrac12\hat C_D\hat Y_d .
 * @f]
 * Conventions (chain_flavor_matching notes): @f$\mathcal L\supset-\bar q\tilde H\hat Y_uu-\bar qH\hat Y_dd+{\rm h.c.}@f$,
 * i.e. @f$\hat Y@f$ has the doublet index first. RGESolver stores the hermitian conjugate,
 * @f$\mathcal Y=\hat Y^\dagger@f$ (right-handed row, left-handed column). The flavour basis is
 * the chain basis, which is the UP basis of the SMEFT (@f$\hat Y_u@f$ diagonal); the right-handed
 * down quarks are canonically normalized with a Cholesky factor. Cuhat, Cdhat are stored as
 * @f$M^2\hat C@f$ (dimensionless), so the *_LNP inputs of %NPSMEFTd6General are
 * @f$\Lambda_{\rm NP}^2\hat C_{U,D}@f$ = Cuhat/MVLQ^2, Cdhat/MVLQ^2.
 * The coefficients and Yukawas are set at @f$\Lambda_{\rm NP}@f$: the running between @f$M@f$ and
 * @f$\Lambda_{\rm NP}@f$ is neglected, so @f$\Lambda_{\rm NP}@f$ should be chosen close to the
 * scanned range of @f$M@f$.
 * Contributions quadratic in @f$\hat C@f$ (one-loop VLQ boxes, double Z exchange) are not included.
 *
 * @anchor NPSMEFTd6flavchainParameters
 * <h3>Model parameters</h3>
 * <table class="model">
 * <tr><th>Label</th><th>Description</th></tr>
 * <tr><td>Lambda_NP</td><td>Fixed scale in GeV at which the matching conditions are imposed and the RG evolution starts.</td></tr>
 * <tr><td>MVLQ</td><td>Dimensionless ratio @f$M/\Lambda_{\rm NP}@f$ of the common VLQ mass to Lambda_NP.</td></tr>
 * <tr><td>eps</td><td>Hopping parameter @f$\epsilon=\mu/M@f$.</td></tr>
 * <tr><td>lamu1, lamu2</td><td>Couplings of q1, q2 to the first site of the up chains (B5, B2).</td></tr>
 * <tr><td>lamd1, lamd2, lamd3</td><td>Couplings of q1, q2, q3 to the first site of the down chains (B4, B4, B1).</td></tr>
 * <tr><td>x1u, x2u</td><td>Mixings of the last up sites with u1, u2.</td></tr>
 * <tr><td>x1d ... x6d</td><td>Mixings of the down sites with the light singlets: chain 1 site 4 to d2 (x1d) and d3 (x2d e^{i phi});
 * chain 2 site 3 to d2 (x3d) and d3 (x4d); chain 2 site 4 to d1 (x5d); chain 3 site 1 to d3 (x6d).</td></tr>
 * <tr><td>phi</td><td>The physical phase of the down-sector chain loop.</td></tr>
 * <tr><td>Ytop</td><td>Elementary top Yukawa coupling at Lambda_NP.</td></tr>
 * </table>
 *
 * <h3>Model flags</h3>
 * <table class="modelflags">
 * <tr><th>Label</th><th>Values</th><th>Description</th></tr>
 * <tr><td>U2</td><td>true/false</td><td>U(2)_q limit: lamu2 and lamd2 are ignored and set equal to lamu1 and lamd1.</td></tr>
 * </table>
 */
class NPSMEFTd6flavchain : public NPSMEFTd6General {
public:

    static const int NNPSMEFTd6flavchainVars = 18;

    static const std::string NPSMEFTd6flavchainVars[NNPSMEFTd6flavchainVars];

    NPSMEFTd6flavchain();

    /**
     * @brief Computes the chain-model matrices, replaces the SMEFT Yukawas at MVLQ and sets the
     * non-vanishing Wilson coefficients, then calls NPSMEFTd6General::PostUpdate().
     */
    virtual bool PostUpdate() override;

    virtual bool CheckParameters(const std::map<std::string, double>& DPars) override;

    virtual bool setFlag(const std::string name, const bool value) override;

    /**
     * @brief Only SMEFTBasisFlag = UP is accepted.
     */
    virtual bool setFlagStr(const std::string name, const std::string value) override;

    /** @brief Matrices at Lambda_NP in the chain (= UP) basis, notes conventions (Cuhat, Cdhat = M^2 Chat). */
    const gslpp::matrix<gslpp::complex>& getYuhat() const { return Yuhat; }
    const gslpp::matrix<gslpp::complex>& getYdhat() const { return Ydhat; }
    const gslpp::matrix<gslpp::complex>& getCuhat() const { return Cuhat; }
    const gslpp::matrix<gslpp::complex>& getCdhat() const { return Cdhat; }

protected:

    virtual void setParameter(const std::string name, const double& value) override;

    // Input parameters, set from the configuration file
    double MVLQ = 0.; ///< Ratio M/Lambda_NP of the common vector-like quark mass to Lambda_NP (dimensionless).
    double lamu1 = 0., lamu2 = 0.; ///< Up-type light-heavy couplings.
    double lamd1 = 0., lamd2 = 0., lamd3 = 0.; ///< Down-type light-heavy couplings.
    double eps = 0.; ///< mu/M.
    double xu1 = 0., xu2 = 0.; ///< Hopping coefficients of the up-type chains.
    double xd1 = 0., xd2 = 0., xd3 = 0., xd4 = 0., xd5 = 0., xd6 = 0.; ///< Hopping coefficients of the down-type chains.
    double phi = 0.; ///< Phase of the down-type chain
    double Ytop = 0.; ///< Top Yukawa coupling at Lambda_NP.

    gslpp::matrix<gslpp::complex> Yuhat, Ydhat; ///< Yukawa matrices at Lambda_NP (doublet index first).
    gslpp::matrix<gslpp::complex> Cuhat, Cdhat; ///< M^2 Chat_U, M^2 Chat_D (dimensionless, hermitian).

    /**
     * @brief Computes Yuhat, Ydhat, Cuhat, Cdhat exactly in eps (chain_flavor_matching notes, Sects. 5-6).
     */
    void computeMatrices();

private:

    bool flag_U2 = false; ///< U(2)_q limit.

};

#endif /* NPSMEFTD6FLAVCHAIN_H */