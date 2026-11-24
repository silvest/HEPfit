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
 * matched at tree level onto the SMEFT at the scale @f$\mu = M_{\rm VLQ}@f$, exactly in
 * @f$\epsilon=\mu/M@f$ (dimension-six truncation only).
 * @details The UV parameters are sampled directly. At @f$\Lambda_{\rm NP}=M_{\rm VLQ}@f$ the model
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
 * down quarks are canonically normalized with a Cholesky factor. Yuhat, Ydhat, Cuhat, Cdhat are
 * stored in units such that Cuhat, Cdhat @f$=M^2\hat C@f$ are dimensionless, which with
 * @f$\Lambda_{\rm NP}=M_{\rm VLQ}@f$ are directly the *_LNP inputs of %NPSMEFTd6General.
 * Contributions quadratic in @f$\hat C@f$ (one-loop VLQ boxes, double Z exchange) are not included.
 *
 * @anchor NPSMEFTd6flavchainParameters
 * <h3>Model parameters</h3>
 * <table class="model">
 * <tr><th>Label</th><th>Description</th></tr>
 * <tr><td>MVLQ</td><td>Common VLQ mass in GeV; it is also the matching scale (Lambda_NP = MVLQ).</td></tr>
 * <tr><td>eps</td><td>Hopping parameter @f$\epsilon=\mu/M@f$.</td></tr>
 * <tr><td>lamu1, lamu2</td><td>Couplings of q1, q2 to the first site of the up chains (B5, B2).</td></tr>
 * <tr><td>lamd1, lamd2, lamd3</td><td>Couplings of q1, q2, q3 to the first site of the down chains (B4, B4, B1).</td></tr>
 * <tr><td>x1u, x2u</td><td>Mixings of the last up sites with u1, u2.</td></tr>
 * <tr><td>x1d ... x6d</td><td>Mixings of the down sites with the light singlets: chain 1 site 4 to d2 (x1d) and d3 (x2d e^{i phi});
 * chain 2 site 3 to d2 (x3d) and d3 (x4d); chain 2 site 4 to d1 (x5d); chain 3 site 1 to d3 (x6d).</td></tr>
 * <tr><td>phi</td><td>The physical phase of the down-sector chain loop.</td></tr>
 * <tr><td>Ytop</td><td>Elementary top Yukawa coupling at MVLQ.</td></tr>
 * </table>
 *
 * <h3>Model flags</h3>
 * <table class="modelflags">
 * <tr><th>Label</th><th>Values</th><th>Description</th></tr>
 * <tr><td>U2</td><td>true/false</td><td>lamu2 = lamu1 and lamd2 = lamd1 (U(2)_q limit), up to deltaU2u, deltaU2d.</td></tr>
 * <tr><td>MinimalCouplings</td><td>true/false</td><td>The lambda's are computed from eps and from the SMEFT
 * Yukawas at MVLQ with the leading-order minimal-coupling formulas (largest hopping coefficient of each chain = 1).</td></tr>
 * </table>
 */
class NPSMEFTd6flavchain : public NPSMEFTd6General {
public:

    static const int NNPSMEFTd6flavchainVars = 17;

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

    /** @brief Matrices at MVLQ in the chain (= UP) basis, notes conventions. */
    const gslpp::matrix<gslpp::complex>& getYuhat() const { return Yuhat; }
    const gslpp::matrix<gslpp::complex>& getYdhat() const { return Ydhat; }
    const gslpp::matrix<gslpp::complex>& getCuhat() const { return Cuhat; }
    const gslpp::matrix<gslpp::complex>& getCdhat() const { return Cdhat; }

protected:

    virtual void setParameter(const std::string name, const double& value) override;

    const char* matchingScaleInputAlias() const override { return "MVLQ"; }

    // Input parameters, set from the configuration file
    double MVLQ = 0.; ///< Common vector-like quark mass and matching scale [GeV].
    double lamu1 = 0., lamu2 = 0.; ///< Up-type light-heavy couplings.
    double lamd1 = 0., lamd2 = 0., lamd3 = 0.; ///< Down-type light-heavy couplings.
    double eps = 0.; ///< mu/M, used only with MinimalCouplings.
    double xu1 = 0., xu2 = 0.; ///< Hopping coefficients of the up-type chains, used only with MinimalCouplings.
    double xd1 = 0., xd2 = 0., xd3 = 0., xd4 = 0., xd5 = 0., xd6 = 0.; ///< Hopping coefficients of the down-type chains, used only with MinimalCouplings.
    double phi = 0.; ///< Phase of the down-type chain
    double Ytop = 0.; ///< Top Yukawa coupling.

    gslpp::matrix<gslpp::complex> Yuhat, Ydhat; ///< Yukawa matrices at MVLQ (doublet index first).
    gslpp::matrix<gslpp::complex> Cuhat, Cdhat; ///< M^2 Chat_U, M^2 Chat_D at MVLQ (dimensionless, hermitian).

    /**
     * @brief Computes Yuhat, Ydhat, Cuhat, Cdhat exactly in eps (chain_flavor_matching notes, Sects. 5-6).
     */
    void computeMatrices();

private:

    bool flag_U2; ///< U(2)_q limit.

};

#endif /* NPSMEFTD6FLAVCHAIN_H */
