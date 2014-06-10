/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE_H
#define	NPEFFECTIVE_H

#include <string.h>
#include <stdexcept>
#include "NPbase.h"

/**
 * @class NPEffective
 * @brief The auxiliary base model class for NPEffective1 and NPEffective2 classes.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is an auxiliary Model class containing parameters and functions 
 * associated with the dimension-six effective Lagrangian in the basis of
 * \cite Barbieri:1999tm,
 * @f[
 * \mathcal{L}_\mathrm{eff}
 * = \mathcal{L}_\mathrm{SM}
 *   + \sum_i \frac{C_i}{\Lambda^2} \mathcal{O}_i
 * @f]
 * where only the SM gauge-invariant and flavour-diagonal dimension-six operators
 * that are relevant to the electroweak precision observables are considered:
 * @f{align}{
 * \mathcal{O}_{WB} &= (H^\dagger\tau^a H) W^a_{\mu\nu}B^{\mu\nu},
 * &
 * \mathcal{O}_{H} &= | H^\dagger D_\mu H|^2\,,
 * \nonumber\\
 * \mathcal{O}_{LL} &= \frac{1}{2}(\overline{L}\gamma_\mu\tau^a L)^2\,,
 * &
 * \mathcal{O}^\prime_{HL} &= i(H^\dagger D_\mu\tau^a H)(\overline{L}\gamma^\mu\tau^a L)\,,
 * \nonumber\\
 * \mathcal{O}^\prime_{HQ} &= i(H^\dagger D_\mu\tau^a H)(\overline{Q}\gamma^\mu\tau^a Q)\,,
 * &
 * \mathcal{O}_{HL} &= i(H^\dagger D_\mu H)(\overline{L}\gamma^\mu L)\,,
 * \nonumber\\
 * \mathcal{O}_{HQ} &= i(H^\dagger D_\mu H)(\overline{Q}\gamma^\mu Q)\,,
 * &
 * \mathcal{O}_{HE} &= i(H^\dagger D_\mu H)(\overline{E}\gamma^\mu E)\,,
 * \nonumber\\
 * \mathcal{O}_{HU} &= i(H^\dagger D_\mu H)(\overline{U}\gamma^\mu U)\,,
 * &
 * \mathcal{O}_{HD} &= i(H^\dagger D_\mu H)(\overline{D}\gamma^\mu D)\,,
 * @f}
 * as well as the Hermitian conjugate for operators @f$\mathcal{O}^\prime_{HL}@f$
 * to @f$\mathcal{O}_{HD}@f$. 
 *
 * 
 * @anchor NPEffectiveInitialization
 * <h3>Initialization</h3>
 *
 * This class is intended to be used with an inherited model class.
 *
 *
 * @anchor NPEffectiveParameters
 * <h3>%Model parameters</h3>
 *
 * There is no model parameter in the current class.
 *
 *
 * @anchor NPEffectiveFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPEffectiveFunctions
 * <h3>Important member functions</h3>
 *
 * Compared to the base class NPbase, the functions for the following quantities
 * are reimplemented in the current class:
 *
 * @li @f$v@f$&nbsp;&nbsp;(with v()),
 * @li @f$M_{W}^{\mathrm{tree}}@f$&nbsp;&nbsp;(with Mw_tree()),
 * @li @f$\Delta G@f$&nbsp;&nbsp;(with DeltaGF()),
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW()),
 * @li @f$\delta g_V^f@f$&nbsp;&nbsp;(with deltaGVl() and deltaGVq()),
 * @li @f$\delta g_A^f@f$&nbsp;&nbsp;(with deltaGAl() and deltaGAq()),
 *
 * and also the get methods getCoeff() to retrieve the value of each coefficient.
 * 
 */
class NPEffective : public NPbase {
public:

    /**
     * @brief The default constructor.
     */
    NPEffective();

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A get method to retrieve the value of the coefficient for a
     * dimension-six operator.
     * @param[in] name name of the coefficient to be retrieved:
     * cWB, cH, cL1L1, cL1L2, etc.
     * @return the coefficient of the dimension-six operator 
     */
    double getCoeff(const std::string name) const;

    /**
     * @brief A get method to retrieve the value of the new physics scale \f$\Lambda\f$.
     * @return the value of the new physics scale \f$\Lambda\f$ in GeV
     */
    double getLambdaNP() const
    {
        return LambdaNP;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The SM Higgs vacuum expectation value \f$v\f$.
     * @details
     * @f[
     *    v = \left(\frac{1}{\sqrt{2} G_F}\right)^{1/2}.
     * @f]
     * @return the SM value of the Higgs VEV in GeV
     * 
     */
    virtual double v() const;

    /**
     * @brief @copybrief StandardModel::Mw_tree()
     * @details The tree-level mass is obtained by subtracting the effect
     * of @f$\Delta G@f$. 
     * @return @f$M_W^{\mathrm{tree}}@f$ in GeV.
     */
    virtual double Mw_tree() const;

    /**
     * @brief @copybrief NPbase::DeltaGF()
     * @details The new physics contribution @f$\Delta G@f$ is defined as
     * @f[
     * G_\mu = G_{\mu,\mathrm{SM}}(1+\Delta G)\,,
     * @f]
     * where @f$G_\mu@f$ is the experimentl value measured through muon decays,
     * and @f$G_{\mu,\mathrm{SM}}@f$ is the Fermi constant in the SM.
     * The dimension-six operators yield the following corrections:
     * @f[
     *  \Delta G_F
     *  = - C_{LL}\left(\frac{v}{\Lambda}\right)^2
     *    + C^\prime_{HL_1}\left(\frac{v}{\Lambda}\right)^2
     *    + C^\prime_{HL_2}\left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return @f$\Delta G@f$
     */
    virtual double DeltaGF() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @brief @copybrief NPbase::obliqueS()
     * @details The operator @f$\mathcal{O}_{WB}@f$ contributes to the @f$S@f$
     * parameter:
     * @f[
     * S = \frac{4s_Wc_W\, C_{WB} }{\alpha(M_Z^2)} \left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @return the value of @f$S@f$
     */
    virtual double obliqueS() const;

    /**
     * @brief @copybrief NPbase::obliqueT()
     * @details The operator @f$\mathcal{O}_{H}@f$ contributes to the @f$T@f$
     * parameter:
     * @f[
     * T = - \frac{C_H}{2\alpha(M_Z^2)}\left(\frac{v}{\Lambda}\right)^2.
     * @f]
     * 
     * See @cite Ciuchini:2013pca and references therein.
     * @return the value of @f$T@f$
     */
    virtual double obliqueT() const;

    /**
     * @brief @copybrief NPbase::obliqueU()
     * @return @f$U=0@f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief @copybrief NPbase::GammaW()
     * @details The width of the @f$W@f$ boson receives the corrections from the
     * operators @f$\mathcal{O}^\prime_{HL}@f$ and @f$\mathcal{O}^\prime_{HQ}@f$,
     * in addition to the corrections via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$:
     * @f[
     * \Gamma_W
     * =
     * \Gamma_{W,\mathrm{SM}}
     * \left[
     * 1 - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *     \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     *   - \frac{1+c_W^2}{2(c_W^2-s_W^2)}\, \Delta G
     *   + \left( C_{HL_1}^\prime + C_{HL_2}^\prime  + C_{HL_3}^\prime
     *          + C_{HQ_1}^\prime + C_{HQ_2}^\prime\right)
     * \left(\frac{v}{\Lambda}\right)^2
     * \right].
     * @f]
     * 
     * See @cite Ciuchini:2013pca and references therein.
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief @copybrief NPbase::deltaGVl()
     * @details New physics contribution to the neutral-current vector
     * coupling @f$g_V^l@f$ is given by
     * @f[
     * \delta g_V^l = \delta g_L^l + \delta g_R^l.
     * @f]
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_V^l@f$
     */
    virtual double deltaGV_f(const Particle p) const;

    /**
     * @brief @copybrief NPbase::deltaGAl()
     * @details New physics contribution to the neutral-current axial-vector
     * coupling @f$g_A^l@f$ is given by
     * @f[
     * \delta g_A^l = \delta g_L^l - \delta g_R^l.
     * @f]
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_A^l@f$
     */    
    virtual double deltaGA_f(const Particle p) const;
    
    ////////////////////////////////////////////////////////////////////////    
protected:

    double cWB;///< The dimension-6 operator coefficient \f$C_{WB}\f$.
    double cH;///< The dimension-6 operator coefficient \f$C_{H}\f$.
    double cL1L1;///< The dimension-6 operator coefficient \f$C_{L_1L_1}\f$.
    double cL1L2;///< The dimension-6 operator coefficient \f$C_{L_1L_2}\f$.
    double cL1L3;///< The dimension-6 operator coefficient \f$C_{L_1L_3}\f$.
    double cL2L2;///< The dimension-6 operator coefficient \f$C_{L_2L_2}\f$.
    double cL2L3;///< The dimension-6 operator coefficient \f$C_{L_2L_3}\f$.
    double cL3L3;///< The dimension-6 operator coefficient \f$C_{L_3L_3}\f$.
    double cHL1p;///< The dimension-6 operator coefficient \f$C_{HL_1}^\prime\f$.
    double cHL2p;///< The dimension-6 operator coefficient \f$C_{HL_2}^\prime\f$.
    double cHL3p;///< The dimension-6 operator coefficient \f$C_{HL_3}^\prime\f$.
    double cHQ1p;///< The dimension-6 operator coefficient \f$C_{HQ_1}^\prime\f$.
    double cHQ2p;///< The dimension-6 operator coefficient \f$C_{HQ_2}^\prime\f$.
    double cHQ3p;///< The dimension-6 operator coefficient \f$C_{HQ_3}^\prime\f$.
    double cHL1;///< The dimension-6 operator coefficient \f$C_{HL_1}\f$.
    double cHL2;///< The dimension-6 operator coefficient \f$C_{HL_2}\f$.
    double cHL3;///< The dimension-6 operator coefficient \f$C_{HL_3}\f$.
    double cHQ1;///< The dimension-6 operator coefficient \f$C_{HQ_1}\f$.
    double cHQ2;///< The dimension-6 operator coefficient \f$C_{HQ_2}\f$.
    double cHQ3;///< The dimension-6 operator coefficient \f$C_{HQ_3}\f$.
    double cHE1;///< The dimension-6 operator coefficient \f$C_{HE_1}\f$.
    double cHE2;///< The dimension-6 operator coefficient \f$C_{HE_2}\f$.
    double cHE3;///< The dimension-6 operator coefficient \f$C_{HE_3}\f$.
    double cHU1;///< The dimension-6 operator coefficient \f$C_{HU_1}\f$.
    double cHU2;///< The dimension-6 operator coefficient \f$C_{HU_2}\f$.
    double cHU3;///< The dimension-6 operator coefficient \f$C_{HU_3}\f$.
    double cHD1;///< The dimension-6 operator coefficient \f$C_{HD_1}\f$.
    double cHD2;///< The dimension-6 operator coefficient \f$C_{HD_2}\f$.
    double cHD3;///< The dimension-6 operator coefficient \f$C_{HD_3}\f$.
    double LambdaNP;///< The new physics scale \f$\Lambda\f$.
    
    ////////////////////////////////////////////////////////////////////////
private:

    /**
     * @brief New physics contribution to @f$g_L^l@f$.
     * @details New physics contributions to the neutral-current left-handed
     * coupling @f$g_L^l@f$ from the operators @f$\mathcal{O}_{HL_i}^\prime@f$
     * and @f$\mathcal{O}_{HL_i}@f$ are given by
     * @f[
     * \delta g_L^{\nu_i}
     * = \frac{C_{HL_i}^\prime-C_{HL_i}}{2} \left(\frac{v}{\Lambda}\right)^2,
     * \qquad
     * \delta g_L^{e_i}
     * = -\frac{C_{HL_i}^\prime+C_{HL_i}}{2} \left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_L^l@f$
     *
     * @attention The new physics contribution via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$ are not included in this function.
     */
    double deltaGL_f_tmp(const Particle p) const;

    /**
     * @brief New physics contribution to @f$g_R^l@f$.
     * @details New physics contributions to the neutral-current right-handed
     * coupling @f$g_R^l@f$ for @f$l=e_i@f$ from the operators
     * @f$\mathcal{O}_{HE_i}@f$ are given by
     * @f[
     * \delta g_R^{e_i} = -\frac{C_{HE_i}}{2} \left(\frac{v}{\Lambda}\right)^2.
     * @f]
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\delta g_R^l@f$
     *
     * @attention The new physics contribution via @f$S@f$, @f$T@f$, @f$U@f$ and
     * @f$\Delta G@f$ are not included in this function.
     */
    double deltaGR_f_tmp(const Particle p) const;
    
};

#endif	/* NPEFFECTIVE_H */

