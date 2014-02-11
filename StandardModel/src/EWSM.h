/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSM_H
#define	EWSM_H

#include <cstring>
#include <stdexcept>
#include <gslpp.h>
#include "StandardModel.h"
#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"
#include "EWSMTwoLoopQCD.h"
#include "EWSMThreeLoopQCD.h"
#include "EWSMTwoLoopEW.h"
#include "EWSMThreeLoopEW2QCD.h"
#include "EWSMThreeLoopEW.h"
#include "EWSMApproximateFormulae.h"
#include "EWSMOneLoopEW_HV.h"
#include "EWSMTwoFermionsLEP2.h"

using namespace gslpp;

class EWSM_Output; // forward reference to EWSM_Output class

/**
 * @class EWSM
 * @ingroup StandardModel
 * @brief A class for radiative corrections to the %EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains ....
 *
 *
 * the on-shell scheme:
 * @cite Bardin:1980fe
 * @cite Bardin:1981sv
 *
 *
 * @f{align}{
 * {\cal L}
&= 
%\left( \sqrt{2}G_\mu M_Z^2 \right)^{1/2}
\frac{e}{2 s_W c_W}\,
Z_\mu \sum_f \bar{f}
\left( g_{V}^f\gamma_\mu - g_{A}^f \gamma_\mu\gamma_5 \right) f\,,
\nonumber\\
&=
%\left( \sqrt{2}G_\mu M_Z^2 \rho_Z^f \right)^{1/2}
\frac{e}{2 s_W c_W}
\sum_f \sqrt{\rho_Z^f}\,
Z_\mu \bar{f}
\left[( I_3^f - 2Q_f\kappa_Z^f\sthw)\gamma^\mu
  - I_3^f\gamma^\mu\gamma_5\right]f\,,
\label{eq:L_NCint}
\end{align}
where the couplings satisfy the relations
\begin{align}
g_V^f &= \sqrt{\rho_Z^f} I_3^f (1 - 4|Q_f|\kappa_Z^f\sthw)
= \sqrt{\rho_Z^f} (I_3^f - 2Q_f\kappa_Z^f\sthw)\,,
\\
g_A^f &= \sqrt{\rho_Z^f} I_3^f\,,
\label{eq:effectivecouplings_gV_gA}
\end{align}
and
\begin{align}
\rho_Z^f &= \left( \frac{g_A^f}{I_3^f} \right)^2,
\\
\kappa_Z^f &= \frac{1}{4|Q_f|\sthw}
\left( 1 - \frac{g_V^{f}}{g_A^{f}}\right).
 * @f}
 *
 * cashing method
 *
 *
 * experimental/running-width
 *
 * complex-pole/fixed-width
 *
 */
class EWSM {
public:

    friend class EWSM_Output;///< Friend class of the current class.

    /**
     * @brief The target accuracy of the iterative calculation of the
     * @f$W@f$-boson mass in units of GeV.
     */
    static const double Mw_error;
    
    /**
     * @brief The order of radiative corrections to %EW precision observables.
     */
    enum orders_EW
    {
        EW1=0,
        EW1QCD1,
        EW1QCD2,
        EW2,
        EW2QCD1,
        EW3,
        orders_EW_size ///< The size of this enum.
    };

    /**
     * @brief The number of the SM parameters that are relevant to the %EW
     * precision observalbles.
     * @details This constant is used for the cashing method.
     * See also checkSMparams().
     */
    static const int NumSMParams = 27;

    /**
     * @brief A method to check if a given scheme name in string form is valid.
     * @param[in] scheme scheme name for @f$M_W@f$, @f$\rho_Z^f@f$ or @f$\kappa_Z^f@f$
     * @return a boolean that is true if the scheme name is valid
     */
    bool checkEWPOscheme(const std::string scheme) const
    {
        if (scheme.compare("NORESUM") == 0
                || scheme.compare("OMSI") == 0
                || scheme.compare("INTERMEDIATE") == 0
                || scheme.compare("OMSII") == 0
                || scheme.compare("APPROXIMATEFORMULA") == 0)
            return true;
        else
            return false;
    }

    /**
     * @brief A method to convert a given scheme name in string form into a 
     * floating-point number with double precision.
     * @details This method is used in EWSM::checkSMparams() for caching the
     * schemes used in computing @f$M_W@f$, @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$.
     * @param[in] scheme scheme name that is used in computing @f$M_W@f$,
     * @f$\rho_Z^f@f$ or @f$\kappa_Z^f@f$
     * @return a floating-point number with double precision corresponding to
     * the given scheme name 
     */
    double SchemeToDouble(const std::string scheme) const
    {
        if (scheme.compare("NORESUM") == 0) 
            return 0.0;
        else if (scheme.compare("OMSI") == 0)
            return 1.0;
        else if (scheme.compare("INTERMEDIATE") == 0)
            return 2.0;
        else if (scheme.compare("OMSII") == 0)
            return 3.0;
        else if (scheme.compare("APPROXIMATEFORMULA") == 0)
            return 4.0;
        else
            throw std::runtime_error("EWSM::SchemeToDouble: bad scheme");
    }
        
    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief Constructor. 
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSM(const StandardModel& SM_i);

    /**
     * @brief The default destructor.
     */
    virtual ~EWSM();

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A set method to change the model flag CacheInEWSM in StandardModel.
     * @details Setting CacheInEWSM to false, the caching methods defined in the
     * current class are not employed in numerical computations. The flag is set
     * to true in the constructor EWSM() by default.
     *
     * See also @ref StandardModelFlags "the description of the StandardModel flags".
     * @param[in] FlagCacheInEWSM a boolean flag for caching
     */
    void setFlagCacheInEWSM(bool FlagCacheInEWSM)
    {
        this->FlagCacheInEWSM = FlagCacheInEWSM;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMcache.
     * @return myCache
     */
    EWSMcache* getMyCache() const
    {
        return myCache;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMOneLoopEW,
     * @return myOneLoopEW
     */
    EWSMOneLoopEW* getMyOneLoopEW() const
    {
        return myOneLoopEW;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMApproximateFormulae.
     * @return myApproximateFormulae
     */
    EWSMApproximateFormulae* getMyApproximateFormulae() const
    {
        return myApproximateFormulae;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMTwoFermionsLEP2.
     * @return myTwoFermionsLEP2
     */
    EWSMTwoFermionsLEP2* getMyTwoFermionsLEP2() const 
    {
        return myTwoFermionsLEP2;
    }

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief
     * @details
     * @param[in,out] Params_cache
     * @param[in] bUpdate
     * @return
     */
    bool checkSMparams(double Params_cache[], const bool bUpdate=true) const;

    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @param[in] s invariant mass squared 
     * @return the leptonic corrections to alpha
     */
    double DeltaAlphaLepton(const double s) const;    

    /**
     * @return the sum of the leptonic and hadronic corrections to alpha at Mz
     */
    double DeltaAlphaL5q() const;

    /**
     * @param[in] s invariant mass squared 
     * @return the top-quark corrections to alpha
     */
    double DeltaAlphaTop(const double s) const;    

    /**
     * @return the total (leptonic+hadronic+top) corrections to alpha at Mz
     */
    double DeltaAlpha() const;
    
    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double alphaMz() const;

   
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The square of the sine of the weak angle @f$s_0^2@f$ defined
     * without weak radiative corrections.
     * @details The quantity @f$s_0^2@f$ is defined through
     * @f[
     * s_0^2 c_0^2 = \frac{\pi\,\alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2}
     * \ \ \rightarrow\ \
     * s_0^2 = \frac{1}{2}
     * \left(1 - \sqrt{1
     *  - \frac{4\pi \alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2}}\  \right)\,.
     * @f]
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$s_0^2@f$
     */
    double s02() const;

    /**
     * @brief The square of the cosine of the weak angle @f$c_0^2@f$ defined
     * without weak radiative corrections.
     * @details The quantity @f$c_0^2@f$ is given by
     * @f[
     * c_0^2 = 1 - s_0^2\,,
     * @f]
     * where @f$s_0^2@f$ is defined in s02(). 
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$s_0^2@f$
     */
    double c02() const;


    ////////////////////////////////////////////////////////////////////////     
        
    /**
     * @brief The SM prediction for the @f$W@f$-boson mass in the on-shell scheme,
     * @f$M_{W,\mathrm{SM}}@f$.
     * @details If the model flag @ref StandardModelFlags "CacheInEWSM"
     * in StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     *
     * When the model flag @ref StandardModelFlags "Mw" in StandardModel is set
     * to APPROXIMATEFORMULA, the current function uses
     * the two-loop approximate formula in EWSMApproximateFormulae::Mw(),
     * which includes the full two-loop %EW contribution of
     * @f${\cal O}(\alpha^2)@f$ as well as the leading
     * @f${\cal O}(G_\mu^2\alpha_s m_t^4)@f$ and @f${\cal O}(G_\mu^3m_t^6)@f$
     * contributions.
     *
     * When @ref StandardModelFlags "Mw" is not set to APPROXIMATEFORMULA,
     * the @f$W@f$-boson mass is computed from @f$\Delta r(M_W)@f$ with an
     * iterative procedure. The target accuracy of the iterative calculation
     * is specified with the constant #Mw_error.
     * See also resumMw(), in which @f$M_W@f$ is computed with a given
     * @f$\Delta r@f$, equivalently with @f$\Delta\rho@f$ and
     * @f$\Delta r_{\mathrm{rem}}@f$
     *
     * @return @f$M_{W,\mathrm{SM}}@f$ in GeV
     */
    double Mw_SM() const;
    
    /** 
     * @brief The SM prediction for @f$\Delta r@f$ derived from @f$M_{W,\mathrm{SM}}@f$. 
     * @details 
     * If the model flag @ref StandardModelFlags "Mw" in StandardModel is set
     * to NORESUM, the quantity @f$\Delta r@f$ is computed by using the following
     * relation:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta r)\,.
     * @f]
     * Otherwise, the following relation is employed instead:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu(1-\Delta r)}\,.
     * @f]
     * @sa The corresponding quantity in the so-called complex-pole scheme
     * (instead of the experimental/running-widthr scheme)
     * is defined in DeltaRbar_SM().
     * @return @f$\Delta r_{\mathrm{SM}}@f$
     */
    double DeltaR_SM() const;

    /**
     * @brief The SM prediction for the square of the cosine of the weak mixing
     * angle in the on-shell scheme, @f$c_{W,\mathrm{SM}}^2@f$.
     * @details
     * @f[
     *   c_{W,\mathrm{SM}}^2 = \frac{M_{W,\mathrm{SM}}^2}{M_Z^2}\,.
     * @f]
     * @return @f$c_{W,\mathrm{SM}}^2@f$
     */
    double cW2_SM() const;
    
    /**
     * @brief The SM prediction for the square of the sine of the weak mixing
     * angle in the on-shell scheme, @f$s_{W,\mathrm{SM}}^2@f$.
     * @details
     * @f[
     *   s_{W,\mathrm{SM}}^2 = 1 - c_{W,\mathrm{SM}}^2\,,
     * @f]
     * where @f$c_{W,\mathrm{SM}}^2@f$ is given by cW2_SM().
     * @return @f$s_{W,\mathrm{SM}}^2@f$
     */
    double sW2_SM() const;

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The @f$Z@f$-boson mass @f$\overline{M}_Z@f$
     * in the complex-pole/fixed-width scheme.
     * @details The mass parameter @f$\overline{M}_Z@f$ in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt is given by 
     * @f[
     *   \overline{M}_{Z} = M_{Z} - \frac{\Gamma_{Z}^2}{2M_{Z}}\,,
     * @f]
     * where @f$M_Z@f$ and @f$\Gamma_{Z}@f$ are the mass and width of the
     * @f$Z@f$ boson in the experimental/running-width scheme:
     * @f{align}{
     * \Gamma(Z\to f\bar{f})
     * = \frac{G_\mu M_Z^3}{24\sqrt{2}\pi}
     * \left[
     * \left( \frac{v_f}{a_f} \right)^2 + 1
     * \right]
     * \times
     * \left\{
     * \begin{array}{ll}
     * 1 & \mathrm{for}\quad f=\ell\,,
     * \\[2mm]
     * \displaystyle
     * N_c \left( 1 + \frac{\alpha_s(M_Z^2)}{\pi} \right)
     * & \mathrm{for}\quad f=q
     * \end{array}
     * \right.
     * @f}
     * with @f$v_f/a_f=1-4|Q_f|s_{W,\mathrm{tree}}^2@f$.
     * @return @f$\overline{M}_Z@f$ in GeV
     */
    double Mzbar() const;

    /**
     * @brief A method to convert the @f$W@f$-boson mass
     * from the experimental/running-width scheme 
     * to the complex-pole/fixed-width scheme.
     * @details The mass parameter @f$\overline{M}_W@f$ in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt is given by
     * @f[
     *   \overline{M}_{W} = M_{W} - \frac{\Gamma_{W}^2}{2M_{W}}\,,
     * @f]
     * where @f$M_W@f$ and @f$\Gamma_{W}@f$ are the mass and width of the
     * @f$W@f$ boson in the experimental/running-width scheme:
     * @f[
     * \Gamma_W
     * =
     *   \frac{3G_\mu M_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(M_W^2)}{3\pi} \right)\,.
     * @f]
     * @param[in] Mw the @f$W@f$-boson mass in the experimental/running-width scheme
     * @return @f$\overline{M}_W@f$ in GeV
     */
    double MwbarFromMw(const double Mw) const;

    /**
     * @brief A method to convert the @f$W@f$-boson mass
     * from the complex-pole/fixed-width scheme 
     * to the experimental/running-width scheme.
     * @details The experimental mass @f$M_W@f$ is derived
     * @f[
     *   M_W = \overline{M}_W + \frac{\Gamma_{W}^2}{2\overline{M}_{W}}\,, 
     * @f]
     * where @f$\overline{M}_W@f$ is the mass parameter in the
     * complex-pole/fixed-width scheme, and @f$\Gamma_{W}@f$ is the @f$W@f$-boson
     * width in the experimental/running-width scheme:
     * @f[
     * \Gamma_W
     * =
     *   \frac{3G_\mu M_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(M_W^2)}{3\pi} \right)
     * \approx
     *   \frac{3G_\mu \overline{M}_W^3}{2\sqrt{2}\pi}
     *   \left( 1 + \frac{2\alpha_s(\overline{M}_W^2)}{3\pi} \right)\,.
     * @f]
     * @param[in] Mwbar the @f$W@f$-boson mass in the complex-pole/fixed-width scheme 
     * @return @f$M_W@f$ in GeV
     */
    double MwFromMwbar(const double Mwbar) const;

    /**
     * @brief The SM prediction for @f$\Delta \overline{r}@f$ derived from
     * @f$\overline{M}_{W,\mathrm{SM}}@f$
     * in the complex-pole/fixed-width scheme.
     * @details The quantity @f$\Delta \overline{r}@f$ is computed by using
     * the following relation:
     * @f[
     * \overline{s}_W^2 \overline{M}_W^2
     * = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta \overline{r})\,,
     * @f]
     * where @f$\overline{M}_W@f$ and @f$\overline{s}_W@f$ are the @f$W@f$-boson
     * mass and the sine of the weak mixing angle in the complex-mass scheme. 
     * @sa The corresponding quantity in the experimental/running-width 
     * scheme is defined in DeltaR_SM().
     * @return @f$\Delta \overline{r}_{\mathrm{SM}}@f$
     */
    double DeltaRbar_SM() const;


    ////////////////////////////////////////////////////////////////////////
    // SM contribution to the effective couplings

    /**
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$ in the SM.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\rho_Z^l@f$ in the SM
     */
    complex rhoZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$ in the SM.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\rho_Z^q@f$ in the SM
     */
    complex rhoZ_q_SM(const QCD::quark q) const;    
    
    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$ in the SM.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\kappa_Z^l@f$ in the SM
     */
    complex kappaZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$ in the SM.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\kappa_Z^q@f$ in the SM
     */
    complex kappaZ_q_SM(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$ in the SM.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_V^l@f$ in the SM
     */
    complex gVl_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current vector coupling @f$g_V^q@f$ in the SM.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_V^q@f$ in the SM
     */
    complex gVq_SM(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$ in the SM.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_A^l@f$ in the SM
     */
    complex gAl_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current axial-vector coupling @f$g_A^q@f$ in the SM.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_A^q@f$ in the SM
     */
    complex gAq_SM(const QCD::quark q) const;
    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\rho_Z^l@f$ 
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\rho_Z^q@f$ 
     */
    virtual complex rhoZ_q(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\kappa_Z^l@f$ 
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\kappa_Z^q@f$ 
     */
    virtual complex kappaZ_q(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_V^l@f$
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current vector coupling @f$g_V^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_V^q@f$ 
     */
    virtual complex gVq(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_A^l@f$ 
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current axial-vector coupling @f$g_A^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_A^q@f$ 
     */
    virtual complex gAq(const QCD::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return the top-quark corrections to the Z-b-bbar vertex
     */
    double taub() const;    
    

    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] l lepton
     * @return flavor-dependent correction to rho_Z^l with respect to that for the charged leptons
     */
    complex rhoZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark
     * @return flavor-dependent correction to rho_Z^q with respect to that for the charged leptons
     */
    complex rhoZ_q_SM_FlavorDep(QCD::quark q) const;

    /**
     * @param[in] l lepton
     * @return flavor-dependent correction to kappa_Z^l with respect to that for the charged leptons
     */
    complex kappaZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark
     * @return flavor-dependent correction to kappa_Z^q with respect to that for the charged leptons
     */
    complex kappaZ_q_SM_FlavorDep(QCD::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_1@f$.
     * @details The parameters @f$\varepsilon_1@f$ is defined as
     * @f[
     * \varepsilon_1 = \Delta\rho'\,,
     * @f]
     * where @f$\Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right)@f$.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{1,\mathrm{SM}}@f$
     */
    double epsilon1_SM() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_2@f$.
     * @details The parameters @f$\varepsilon_2@f$ is computed via the formula: 
     * @f[
     * \varepsilon_2 = c_0^2  \Delta\rho'
     * + \frac{s_0^2}{c_0^2 - s_0^2} \Delta r_W
     * - 2 s_0^2 \Delta\kappa'\,,
     * @f]
     * where @f$\Delta\rho'@f$, @f$\Delta r_W@f$ and @f$\Delta\kappa'@f$ are 
     * defined as 
     * @f{align}{
     * \Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right),\qquad
     * \Delta r_W = 1 - \frac{\pi\,\alpha(M_Z^2)}{\sqrt{2}\,G_\mu M_Z^2 s_W^2 c_W^2},\qquad
     * \Delta\kappa' = \frac{\sin^2\theta_{\mathrm{eff}}^e}{s_0^2} - 1\,,
     * @f}
     * and @f$s_0^2@f$ and @f$c_0^2@f$ are given in s02() and c02(), respectively.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{2,\mathrm{SM}}@f$
     */
    double epsilon2_SM() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_3@f$.
     * @details The parameters @f$\varepsilon_3@f$ is computed via the formula:
     * @f[
     * \varepsilon_3 = c_0^2\Delta\rho' + (c_0^2-s_0^2)\Delta\kappa'\,,
     * @f]
     * where @f$\Delta\rho'@f$ and @f$\Delta\kappa'@f$ are
     * defined as
     * @f{align}{
     * \Delta\rho'=2\left(\sqrt{{\rm Re}(\rho_Z^e)}-1\right),\qquad
     * \Delta\kappa' = \frac{\sin^2\theta_{\mathrm{eff}}^e}{s_0^2} - 1\,,
     * @f}
     * and @f$s_0^2@f$ and @f$c_0^2@f$ are given in s02() and c02(), respectively.
     *
     * See @cite Altarelli:1990zd and @cite Altarelli:1991fk.
     * @return @f$\varepsilon_{3,\mathrm{SM}}@f$
     */
    double epsilon3_SM() const;

    /**
     * @brief The SM contribution to the epsilon parameter @f$\varepsilon_b@f$.
     * @details The parameters @f$\varepsilon_b@f$ is computed via the formula:
     * @f[
     * \epsilon_b = 
     * \frac{ {\rm Re}\left[ \kappa_Z^e + \Delta\kappa_Z^b \right]}
     * {{\rm Re}(\kappa_Z^b)} - 1\,,
     * @f]
     * where @f$\Delta\kappa_Z^b@f$, representing flavour non-universal vertex
     * corrections to the @f$Zb\bar{b}@f$ vertex, is neglected when the
     * model flag WithoutNonUniversalVC in StandardModel is set to true.
     * 
     * See @cite Altarelli:1990zd, @cite Altarelli:1991fk and @cite Altarelli:1993sz
     * for the @f$\varepsilon@f$ parameterization
     * and @cite Ciuchini:2013pca for the flavour non-universal vertex corrections. 
     * @return @f$\varepsilon_{b,\mathrm{SM}}@f$
     */
    double epsilonb_SM() const;
    

    ////////////////////////////////////////////////////////////////////////     
    // The W-boson decay width

    /**
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return rho_ij^W for Gamma_W in the SM
     * @attention Fermion masses are neglected. 
     */    
    double rho_GammaW_l_SM(const StandardModel::lepton li, 
                           const StandardModel::lepton lj) const;
    
    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return rho_ij^W for Gamma_W in the SM
     * @attention Fermion masses are neglected. 
     */    
    double rho_GammaW_q_SM(const QCD::quark qi, 
                           const QCD::quark qj) const;
    
    /**
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return the partial width of W^+ decay into an @f$l_i\bar{l_j}@f$ pair in the SM
     * @attention Mixings in the lepton sector are neglected. 
     * @attention Fermion masses are neglected. 
     */
    double GammaW_l_SM(const StandardModel::lepton li, 
                       const StandardModel::lepton lj) const;

    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return the partial width of W^+ decay into an @f$q_i\bar{q_j}@f$ pair in the SM
     * @attention Fermion masses are neglected. 
     */
    double GammaW_q_SM(const QCD::quark qi, 
                       const QCD::quark qj) const;    
    
    /**
     * @return the total width of the W boson in the SM
     */
    double GammaW_SM() const;   


    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;///< A reference to an object of type StandardModel.

    /**
     * @brief computes Delta rho
     * @param[in] Mw_i the W boson mass
     * @param[out] DeltaRho[]
     */
    void ComputeDeltaRho(const double Mw_i, double DeltaRho[orders_EW_size]) const;

    /**
     * @brief computes Delta r_rem
     * @param[in] Mw_i the W boson mass
     * @param[out] DeltaR_rem
     */
    void ComputeDeltaR_rem(const double Mw_i, double DeltaR_rem[orders_EW_size]) const;


    ////////////////////////////////////////////////////////////////////////         
private:
    bool flag_order[orders_EW_size]; 
    
    EWSMcache* myCache;
    EWSMOneLoopEW* myOneLoopEW;
    EWSMTwoLoopQCD* myTwoLoopQCD;
    EWSMThreeLoopQCD* myThreeLoopQCD;
    EWSMTwoLoopEW* myTwoLoopEW;
    EWSMThreeLoopEW2QCD* myThreeLoopEW2QCD;
    EWSMThreeLoopEW* myThreeLoopEW; 
    EWSMApproximateFormulae* myApproximateFormulae;
    
    EWSMTwoFermionsLEP2* myTwoFermionsLEP2;
    
    
    ////////////////////////////////////////////////////////////////////////     
    //caches
    
    bool FlagCacheInEWSM; // true for caching
    
    mutable double DeltaAlphaLepton_params_cache[NumSMParams];
    mutable double DeltaAlphaLepton_cache;
    
    mutable double DeltaAlpha_params_cache[NumSMParams];
    mutable double DeltaAlpha_cache;
    
    mutable double Mw_params_cache[NumSMParams];
    mutable double Mw_cache;

    mutable double rhoZ_l_params_cache[6][NumSMParams];
    mutable complex rhoZ_l_cache[6];

    mutable double rhoZ_q_params_cache[6][NumSMParams];
    mutable complex rhoZ_q_cache[6];

    mutable double kappaZ_l_params_cache[6][NumSMParams];
    mutable complex kappaZ_l_cache[6];

    mutable double kappaZ_q_params_cache[6][NumSMParams];
    mutable complex kappaZ_q_cache[6];
    
    mutable double GammaW_params_cache[NumSMParams];
    mutable double GammaW_cache;
    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] Mw_i the W boson mass
     * @param[in] DeltaRho
     * @param[in] DeltaR_rem
     * @return resummed Mw
     */
    double resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
                   const double DeltaR_rem[orders_EW_size]) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaRho_rem 
     * @param[in] DeltaRbar_rem
     * @param[in] bool_Zbb true for f=b
     * @return resummed Re[rho_Z^f]
     */
    double resumRhoZ(const double DeltaRho[orders_EW_size],
                     const double deltaRho_rem[orders_EW_size], 
                     const double DeltaRbar_rem, const bool bool_Zbb) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaKappa_rem 
     * @param[in] DeltaRbar_rem
     * @param[in] bool_Zbb true for f=b
     * @return resummed Re[kappa_Z^f]
     */
    double resumKappaZ(const double DeltaRho[orders_EW_size],
                       const double deltaKappa_rem[orders_EW_size],
                       const double DeltaRbar_rem, const bool bool_Zbb) const;
    
};

#endif	/* EWSM_H */

