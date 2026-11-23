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
 * @brief Chain models of quark flavour (SU(2)-singlet vector-like quarks, arXiv:2602.17754)
 * matched onto the SMEFT at the scale @f$\mu = M_{\rm VLQ}@f$.
 * @details The model sets directly the non-vanishing Wilson coefficients of %NPSMEFTd6General,
 * in the UP flavour basis (the only basis allowed for this model: SMEFTBasisFlag = UP).
 * In this basis the Yukawa matrices of RGESolver are
 * @f$\mathcal Y_u = {\rm diag}(y_u,y_c,y_t)@f$, @f$\mathcal Y_d = {\rm diag}(y_d,y_s,y_b)V^\dagger@f$,
 * with @f$\mathcal L_Y = -(\bar d\,\mathcal Y_d\, q\, H^\dagger + \bar u\,\mathcal Y_u\, q\,\tilde H^\dagger) + {\rm h.c.}@f$,
 * and the chain model is flavour diagonal:
 * @f$M^2\hat C_U = {\rm diag}(\lambda^{u\,2}_1,\lambda^{u\,2}_2,0)@f$,
 * @f$M^2\hat C_D = {\rm diag}(\lambda^{d\,2}_1,\lambda^{d\,2}_2,\lambda^{d\,2}_3)@f$.
 * With @f$\Lambda_{NP} = M_{\rm VLQ}@f$ the coefficients @f$C_i^{LNP}=\Lambda_{NP}^2 C_i@f$ are
 * @f[
 * C_{Hq}^{(1)} = \frac{M^2}{4}(\hat C_U-\hat C_D),\quad
 * C_{Hq}^{(3)} = -\frac{M^2}{4}(\hat C_U+\hat C_D),\quad
 * C_{uH} = \frac{M^2}{2}\hat C_U\mathcal Y_u^\dagger,\quad
 * C_{dH} = \frac{M^2}{2}\hat C_D\mathcal Y_d^\dagger ,
 * @f]
 * All other Wilson coefficients vanish at the matching scale. No *_LNP input is accepted.
 * The top-Yukawa enhanced running into @f$\Delta F=2@f$ operators is obtained with RGEci = true
 * (flag of %NPSMEFTd6General). Contributions quadratic in @f$\hat C@f$ (one-loop VLQ boxes,
 * double Z exchange) are not included.
 *
 * @anchor NPSMEFTd6flavchainParameters
 * <h3>Model parameters</h3>
 * <table class="model">
 * <tr><th>Label</th><th>Description</th></tr>
 * <tr><td>MVLQ</td><td>Common mass of the vector-like quarks in GeV; it is also the matching scale (Lambda_NP).</td></tr>
 * <tr><td>lamu1, lamu2</td><td>Light-heavy couplings of q1, q2 to the first sites of the up-type chains.</td></tr>
 * <tr><td>lamd1, lamd2, lamd3</td><td>Light-heavy couplings of q1, q2, q3 to the first sites of the down-type chains.</td></tr>
 * <tr><td>eps</td><td>Ratio mu/M of the hopping masses, used only with MinimalCouplings = true.</td></tr>
 * <tr><td>deltaU2u, deltaU2d</td><td>Relative U(2) breaking, \f$\hat C_{22}=\hat C_{11}(1-\delta)\f$, used only with U2 = true.</td></tr>
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

    static const int NNPSMEFTd6flavchainVars = 9;

    static const std::string NPSMEFTd6flavchainVars[NNPSMEFTd6flavchainVars];

    NPSMEFTd6flavchain();

    /**
     * @brief Computes the couplings and sets the non-vanishing Wilson coefficients of
     * %NPSMEFTd6General before calling NPSMEFTd6General::PostUpdate().
     */
    virtual bool PostUpdate() override;

    virtual bool CheckParameters(const std::map<std::string, double>& DPars) override;

    virtual bool setFlag(const std::string name, const bool value) override;

    /**
     * @brief Only SMEFTBasisFlag = UP is accepted.
     */
    virtual bool setFlagStr(const std::string name, const std::string value) override;

    /** @brief The couplings actually used in the matching (after the U2/MinimalCouplings options). */
    double getlamu1_eff() const { return lamu1_eff; }
    double getlamu2_eff() const { return lamu2_eff; }
    double getlamd1_eff() const { return lamd1_eff; }
    double getlamd2_eff() const { return lamd2_eff; }
    double getlamd3_eff() const { return lamd3_eff; }

protected:

    virtual void setParameter(const std::string name, const double& value) override;

    const char* matchingScaleInputAlias() const override { return "MVLQ"; }

    // Input parameters, set from the configuration file
    double MVLQ; ///< Common vector-like quark mass and matching scale [GeV].
    double lamu1, lamu2; ///< Up-type light-heavy couplings.
    double lamd1, lamd2, lamd3; ///< Down-type light-heavy couplings.
    double eps; ///< mu/M, used only with MinimalCouplings.
    double deltaU2u, deltaU2d; ///< Relative U(2) breaking of Chat_22, used only with U2.

    // Derived quantities (can be used as observables through ModelParameter)
    double lamu1_eff, lamu2_eff, lamd1_eff, lamd2_eff, lamd3_eff;
    double ChatU11, ChatU22; ///< Diagonal of M^2 Chat_U (UP basis).
    double ChatD11, ChatD22, ChatD33; ///< Diagonal of M^2 Chat_D (UP basis).

    /**
     * @brief Sets the effective couplings according to the flags.
     * @param[in] Yu, Yd the RGESolver Yukawa matrices at MVLQ in the UP basis
     */
    void computeCouplings(const gslpp::matrix<gslpp::complex>& Yu, const gslpp::matrix<gslpp::complex>& Yd);

    /**
     * @brief Sets the non-vanishing *_LNP coefficients of %NPSMEFTd6General.
     */
    void setChainWilsonCoefficients(const gslpp::matrix<gslpp::complex>& Yu, const gslpp::matrix<gslpp::complex>& Yd);

private:

    bool flag_U2; ///< U(2)_q limit.
    bool flag_minimal; ///< Minimal couplings from eps.

};

#endif /* NPSMEFTD6FLAVCHAIN_H */
