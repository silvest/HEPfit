/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSUSY_H
#define	EWSUSY_H

#include <gslpp.h>
#include <PVfunctions.h>
#include "SUSY.h"

using namespace gslpp;

/**
 * @class EWSUSY
 * @ingroup SUSY
 * @brief A class for SUSY contributions to the EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used for the calculations of SUSY contributions to
 * the EW precision observables, where Rosiek's notation is adopted internally.
 * The conversions from Rosiek's notation to SLHA one are implemented in 
 * EWSUSY::SetRosiekParameters(), which is called from SUSY::PostUpdate().
 * @par References
 * <A HREF="http://inspirehep.net/record/363948?ln=en">
 * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>;
 * @n
 * <A HREF="http://inspirehep.net/record/344599?ln=en">
 * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>;
 * @n
 * <A HREF="http://inspirehep.net/record/401937?ln=en">Rosiek, hep-ph/9511250</A>,
 * where an updated version is available at author's webpage. 
 */
class EWSUSY {
public:

    /**
     * @brief An EWSUSY constructor.
     * @param[in] SUSY_in A reference to a SUSY object. 
     */
    EWSUSY(const SUSY& SUSY_in);

    /**
     * @brief Sets parameters in Rosiek's notation. 
     *  <table border=0>
     *  <tr>
     *  <td>
     *    <table border>
     *    <tr>
     *      <td align="center"><b> Rosiek </b></td>
     *      <td align="center"><b> SLHA </b></td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Y_u@f$ </td>
     *      <td align="center"> @f$Y_U@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Y_d@f$ </td>
     *      <td align="center"> @f$-Y_D@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Y_\ell@f$ </td>
     *      <td align="center"> @f$-Y_E@f$ </td>
     *    </tr>
     *    </table>
     *  </td>
     *  <td>
     *    <table border>
     *    <tr>
     *      <td align="center"><b> Rosiek </b></td>
     *      <td align="center"><b> SLHA </b></td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$A_u@f$ </td>
     *      <td align="center"> @f$-T_U@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$A_d@f$ </td>
     *      <td align="center"> @f$T_D@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$A_\ell@f$ </td>
     *      <td align="center"> @f$T_E@f$ </td>
     *    </tr>
     *    </table>
     *  </td>
     *  <td>
     *    <table border>
     *    <tr>
     *      <td align="center"><b> Rosiek </b></td>
     *      <td align="center"><b> SLHA </b></td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_-@f$ </td>
     *      <td align="center"> @f$U^\dagger@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_+@f$ </td>
     *      <td align="center"> @f$V^\dagger@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_N@f$ </td>
     *      <td align="center"> @f$N^\dagger@f$ </td>
     *    </tr>
     *    </table>
     *  </td>
     *  <td>
     *    <table border>
     *    <tr>
     *      <td align="center"><b> Rosiek </b></td>
     *      <td align="center"><b> SLHA </b></td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_U@f$ </td>
     *      <td align="center"> @f$R_u^\dagger@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_D@f$ </td>
     *      <td align="center"> @f$R_d^T@f$ </td>
     *    </tr>
     *    <tr>
     *      <td align="center"> @f$Z_L@f$ </td>
     *      <td align="center"> @f$R_e^T@f$ </td>
     *    </tr>
     *    </table>
     *  </td>
     *  </tr>
     *  </table>
     * @f$Z_R=\left(\begin{array}{cc} \cos\alpha & -\sin\alpha \\ \sin\alpha & \cos\alpha \end{array}\right)@f$,
     * @f$Z_H=\left(\begin{array}{cc} \sin\beta & -\cos\beta \\ \cos\beta & \sin\beta \end{array}\right)@f$,
     */
    void SetRosiekParameters();

    /**
     * @brief Fermionic contribuiton to the transverse part of a gauge-boson 
     * self-energy, @f$F_A^{ab}(p^2,m_i,m_j,c_V^{aij},c_V^{bji},c_A^{aij},c_A^{bji})@f$.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] p2 The momentum squared @f$p^2@f$.
     * @param[in] mi The mass of a fermion @f$(i)@f$ running in the loop.
     * @param[in] mj The mass of a fermion @f$(j)@f$ running in the loop.
     * @param[in] cV_aij The vector coupling @f$c_V^{aij}@f$ for a vertex with
     * an incoming vector meson @f$(a)@f$, an incoming fermion @f$(i)@f$ and
     * an outgoing fermion @f$(j)@f$.
     * @param[in] cV_bji The vector coupling @f$c_V^{bji}@f$ for a vertex with
     * an incoming vector meson @f$(b)@f$, an incoming fermion @f$(j)@f$ and
     * an outgoing fermion @f$(i)@f$.
     * @param[in] cA_aij The axial-vector coupling @f$c_A^{aij}@f$ for a vertex with
     * an incoming vector meson @f$(a)@f$, an incoming fermion @f$(i)@f$ and
     * an outgoing fermion @f$(j)@f$.
     * @param[in] cA_bji The axial-vector coupling @f$c_A^{bji}@f$ for a vertex with
     * an incoming vector meson @f$(b)@f$, an incoming fermion @f$(j)@f$ and
     * an outgoing fermion @f$(i)@f$.
     * @return @f$F_A^{ab}(p^2,m_i,m_j,c_V^{aij},c_V^{bji},c_A^{aij},c_A^{bji})@f$
     * renormalized at the scale @f$\mu@f$.
     * @par References
     * Eq. (A.7) in [<A HREF="http://inspirehep.net/record/344599?ln=en">
     * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>].
     */
    complex FA(const double mu, const double p2, const double mi, const double mj,
               const complex cV_aij, const complex cV_bji,
               const complex cA_aij, const complex cA_bji) const;

    /**
     * @brief The transverse part of the Z-boson self-energy, @f$\Pi_Z^T(p^2)@f$,
     * in the 't Hooft-Feynman gauge.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] p2 The momentum squared @f$p^2@f$.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\Pi_Z^T(p^2)@f$ renormalized at the scale @f$\mu@f$ in the
     * 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.15) in [<A HREF="http://inspirehep.net/record/344599?ln=en">
     * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>].
     */
    complex PiT_Z(const double mu, const double p2, const double Mw_i) const;

    /**
     * @brief The transverse part of the self-energy, @f$\Pi_{\gamma Z}^T(p^2)@f$,
     * for the mixing between photon and Z boson in the 't Hooft-Feynman gauge.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] p2 The momentum squared @f$p^2@f$.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\Pi_{\gamma Z}^T(p^2)@f$ renormalized at the scale @f$\mu@f$
     * in the 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.18) in [<A HREF="http://inspirehep.net/record/344599?ln=en">
     * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>].
     */
    complex PiT_AZ(const double mu, const double p2, const double Mw_i) const;

    /**
     * @brief The transverse part of the W-boson self-energy, @f$\Pi_W^T(p^2)@f$,
     * in the 't Hooft-Feynman gauge.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] p2 The momentum squared @f$p^2@f$.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\Pi_W^T(p^2)@f$ renormalized at the scale @f$\mu@f$ in the
     * 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.20) in [<A HREF="http://inspirehep.net/record/344599?ln=en">
     * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>].
     */
    complex PiT_W(const double mu, const double p2, const double Mw_i) const;

    /**
     * @brief The derivative of the transverse part of the photon self-energy
     * with respect to @f$p^2@f$, @f$\Pi_{\gamma}^{T\prime}(p^2)@f$,
     * in the 't Hooft-Feynman gauge.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] p2 The momentum squared @f$p^2@f$.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\Pi_{\gamma}^{T\prime}(p^2)@f$ renormalized at the scale @f$\mu@f$
     * in the 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.17) in [<A HREF="http://inspirehep.net/record/344599?ln=en">
     * Chankowski, Pokorski and Rosiek, NPB 423 (1994) 437</A>].
     */
    complex PiTp_A(const double mu, const double p2, const double Mw_i) const;

    /**
     * @brief The renormalized transverse W-boson self-energy at zero momentum 
     * transefer in the 't Hooft-Feynman gauge.
     * @param [in] mu The renormalization scale @f$\mu@f$.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\hat{\Pi}_W^T(0)@f$ in the 't Hooft-Feynman gauge.
     */
    double PiThat_W_0(const double Mw_i) const;

    /**
     * @brief The SM one-loop renormalized vertex and box corrections
     * to @f$\Delta r@f$ in the 't Hooft-Feynman gauge. 
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return The SM one-loop renormalized vertex and box corrections
     * to @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @par References
     * Eq. (4) in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>],
     * in which only the finite contribution is presented. 
     */
    double DeltaR_rem_SM(const double Mw_i) const;

    /**
     * @brief The LL SUSY box corrections to @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return The LL SUSY box corrections to @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.17) in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>].
     */
    double DeltaR_boxLL_SUSY(const double Mw_i) const;

    /**
     * @brief The LR SUSY box corrections to @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return The LR SUSY box corrections to @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     */
    double DeltaR_boxLR_SUSY(const double Mw_i) const;
    
    /**
     * @brief
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] M A charged lepton.
     * @param[in] J A neutrino. 
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$v(M,J)@f$.
     * @par References
     * Eq. (A.19) in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>],
     * in which only the finite contribution is presented.
     */
    complex v(const double mu, const StandardModel::lepton M,
              const StandardModel::lepton J, const double Mw_i) const;

    /**
     * @brief
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] M A charged lepton.
     * @param[in] J A neutrino.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\delta v(M,J)@f$.
     * @par References
     * Eq. (A.21) in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>],
     * in which only the finite contribution is presented.
     */
    complex delta_v(const double mu, const StandardModel::lepton M,
                    const StandardModel::lepton J, const double Mw_i) const;

    /**
     * @brief The renormalized SUSY vertex corrections to @f$\Delta r@f$
     * in the 't Hooft-Feynman gauge.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return The renormalized SUSY vertex corrections to @f$\Delta r@f$
     * in the 't Hooft-Feynman gauge.
     * @par References
     * Eqs. (A.19), (A.21) and (A.23)
     * in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>].
     */
    double DeltaR_vertex_SUSY(const double Mw_i) const;

    /**
     * @brief The SUSY neutrino self-energy at zero momentum transfer
     * in the 't Hooft-Feynman gauge.
     * @param[in] mu The renormalization scale @f$\mu@f$.
     * @param[in] I A neutrino.
     * @param[in] J A neutrino.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return @f$\Sigma_\nu(0,I,J)@f$ in the 't Hooft-Feynman gauge.
     * @par References
     * Eq. (A.24) in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>].
     */
    complex Sigma_nu_0(const double mu, const StandardModel::lepton I,
                       const StandardModel::lepton J, const double Mw_i) const;

    /**
     * @brief The renormalized SUSY neutrino wave-function contribution to
     * @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @param[in] Mw_i The W-boson mass @f$M_W@f$.
     * @return The renormalized SUSY neutrino wave-function contribution to
     * @f$\Delta r@f$ in the 't Hooft-Feynman gauge.
     * @par References
     * Eqs. (A.21), (A.24) and (A.25)
     * in [<A HREF="http://inspirehep.net/record/363948?ln=en">
     * Chankowski, Dabelstein, Hollik, Mosle, Pokorski and Rosiek, NPB 417 (1994) 101</A>].
     */
    double DeltaR_neutrino_SUSY(const double Mw_i) const;

    /**
     * @brief The one-loop contribution to @f$\Delta r@f$ in the MSSM. 
     * @return @f$\Delta r_{\rm MSSM}^{\alpha} = \Delta r_{\rm SM}^{\alpha} + \Delta r_{\rm SUSY}^{\alpha}@f$.
     */
    double DeltaR_MSSM_EW1(const double Mw_i) const;

    /**
     * @brief The one-loop SUSY contribution to @f$\Delta r@f$.
     * @return @f$\Delta r_{\rm SUSY}^{\alpha}@f$.
     */
    double DeltaR_SUSY_EW1(const double Mw_i) const;


private:

    const PVfunctions PV; ///< An object of PVfunctions class.
    const SUSY& mySUSY; ///< A reference to the SUSY object passed to the constructor.

    matrix<complex> Yu; ///< The Yukawa coupling matrix for up-type quarks in Rosiek's notation.
    matrix<complex> Yd; ///< The Yukawa coupling matrix for down-type quarks in Rosiek's notation.
    matrix<complex> Yl; ///< The Yukawa coupling matrix for charged leptons in Rosiek's notation.

    matrix<complex> Au; ///< The trilinear coupling matrix for up-type squarks in Rosiek's notation.
    matrix<complex> Ad; ///< The trilinear coupling matrix for down-type squarks in Rosiek's notation.
    matrix<complex> Al; ///< The trilinear coupling matrix for charged sleptons in Rosiek's notation.

    matrix<complex> Zm; ///< The rotation matrix for negative charginos in Rosiek's notation.
    matrix<complex> Zp; ///< The rotation matrix for positive charginos in Rosiek's notation.
    matrix<complex> ZN; ///< The rotation matrix for neutralinos in Rosiek's notation.
    matrix<complex> ZU; ///< The rotation matrix for up-type squarks in Rosiek's notation.
    matrix<complex> ZD; ///< The rotation matrix for down-type squarks in Rosiek's notation.
    matrix<complex> ZL; ///< The rotation matrix for charged sleptons in Rosiek's notation.
    matrix<complex> Zne; ///< The rotation matrix for sneutrinos in Rosiek's notation.
    matrix<double> ZR; ///< The rotation matrix for CP-even neutral Higgses in Rosiek's notation.
    matrix<double> ZH; ///< The rotation matrix for charged (CP-odd neutral) Higgs and charged (neutral) Goldstone boson in Rosiek's notation.

};

#endif	/* EWSUSY_H */

