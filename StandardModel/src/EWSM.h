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
 * @details This class contains functions to compute the radiative corrections
 * to the electromagnetic couplings, to the mass and the width of the @f$W@f$
 * boson and to the effective @f$Zf\bar{f}@f$ couplings in the Standard %Model
 * (SM).
 *
 *
 * <h3>Notation</h3>
 *
 * The on-mass-shell renormalization scheme @cite Sirlin:1980nh,
 * @cite Marciano:1980pb, @cite Bardin:1980fe, @cite Bardin:1981sv is adopted
 * for UV divergences, and the weak mixing angle is defined in terms of the
 * physical masses of the gauge bosons:
 * @f[
 * s_W^2 \equiv \sin^2\theta_W = 1 - \frac{M_W^2}{M_Z^2}\,,
 * @f]
 * and @f$c_W^2=1-s_W^2@f$.

 * The Fermi constant @f$G_\mu@f$ in @f$\mu@f$ decay is taken as an input
 * quantity instead of the @f$W@f$-boson mass, since the latter has not been
 * measured very precisely compared to the former. The relation between
 * @f$G_\mu@f$ and @f$M_W@f$ is written as
 * @f[
 * G_\mu = \frac{\pi\,\alpha}{\sqrt{2} s_W^2 M_W^2} (1+\Delta r)\,,
 * @f]
 * where @f$\Delta r@f$ represents radiative corrections. From this relation,
 * the @f$W@f$-boson mass is calculated as
 * @f[
 * M_W^2
 * = \frac{M_Z^2}{2}
 * \left( 1+\sqrt{1-\frac{4\pi\alpha}{\sqrt{2}G_\mu M_Z^2}\,(1+\Delta r)}\
 * \right).
 * @f]
 *
 * The interaction between the @f$Z@f$ boson and the neutral current can be
 * written in terms of the effective @f$Zf\bar{f}@f$ couplings @f$g_{V}^f@f$
 * and @f$g_{A}^f@f$, of @f$g_{R}^f@f$ and @f$g_{L}^f@f$, or of @f$\rho_Z^f@f$
 * and @f$\kappa_Z^f@f$:
 * @f{eqnarray}{
 * \mathcal{L}
 * &=&
 * \frac{e}{2 s_W c_W}\,
 * Z_\mu \sum_f \bar{f}
 * \left( g_{V}^f\gamma_\mu - g_{A}^f \gamma_\mu\gamma_5 \right)\, f\,,
 * \\
 * &=&
 * \frac{e}{2s_W c_W}\,
 * Z_\mu \sum_f \bar{f}
 * \left[ g_{R}^f \gamma_\mu (1 + \gamma_5)
 * + g_{L}^f \gamma_\mu (1 - \gamma_5) \right]\, f\,,
 * \\
 * &=&
 * \frac{e}{2 s_W c_W}\sqrt{\rho_Z^f}\,
 * Z_\mu \sum_f \bar{f}
 * \left[( I_3^f - 2Q_f\kappa_Z^f s_W^2)\gamma^\mu
 *   - I_3^f\gamma^\mu\gamma_5\right]\,f\,,
 * @f}
 * where @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$ are related to
 * @f$g_{V}^f@f$ and @f$g_{A}^f@f$ as the relations:
 * @f{eqnarray}{
 * g_V^f
 * &=&
 * \sqrt{\rho_Z^f} I_3^f (1 - 4|Q_f|\kappa_Z^fs_W^2)
 * = \sqrt{\rho_Z^f} (I_3^f - 2Q_f\kappa_Z^fs_W^2)\,,
 * \qquad
 * g_A^f
 * &=&
 * \sqrt{\rho_Z^f} I_3^f\,,
 * @f}
 * and 
 * @f{eqnarray}{
 * \rho_Z^f &=& \left( \frac{g_A^f}{I_3^f} \right)^2,
 * \qquad
 * \kappa_Z^f &=& \frac{1}{4|Q_f|s_W^2}
 * \left( 1 - \frac{g_V^{f}}{g_A^{f}}\right).
 * @f}
 *
 *
 * <h3>Important member functions</h3>
 *
 * The current class handles the following quantities:
 *
 * @li @f$\Delta\alpha_{\mathrm{lept}}(s)@f$&nbsp;&nbsp; (with DeltaAlphaLepton()),
 * @li @f$\Delta\alpha^{\ell+5q}(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlphaL5q()),
 * @li @f$\Delta\alpha_{\mathrm{top}}(s)@f$&nbsp;&nbsp; (with DeltaAlphaTop()),
 * @li @f$\Delta\alpha(M_Z^2)@f$&nbsp;&nbsp; (with DeltaAlpha()),
 * @li @f$\alpha(M_Z^2)@f$&nbsp;&nbsp; (with alphaMz()),
 *
 * @li @f$M_W@f$&nbsp;&nbsp; (with Mw_SM()),
 * @li @f$\Delta r@f$&nbsp;&nbsp; (with DeltaR_SM()),
 * @li @f$c_W^2@f$ and @f$s_W^2@f$&nbsp;&nbsp; (with cW2_SM() and sW2_SM()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp; (with GammaW_SM()),
 *
 * @li @f$\rho_Z^f@f$&nbsp;&nbsp; (with rhoZ_l() and rhoZ_q()),
 * @li @f$\kappa_Z^f@f$&nbsp;&nbsp; (with kappaZ_l() and kappaZ_q()),
 * @li @f$g_V^f@f$&nbsp;&nbsp; (with gVl() and gVq()),
 * @li @f$g_A^f@f$&nbsp;&nbsp; (with gAl() and gAq()),
 *
 * @li @f$\varepsilon_{1,2,3,b}@f$&nbsp;&nbsp; (with epsilon1_SM(), epsilon2_SM(),
 * epsilon3_SM() and epsilonb_SM()).
 *
 * Moreover, the functions Mzbar(), MwbarFromMw(), MwFromMwbar() and DeltaRbar_SM()
 * can be used for the quantities in the complex-pole/fixed-width scheme.
 *
 * 
 * <h3>Schemes</h3>
 *
 * The formulae used for the @f$W@f$-boson mass @f$M_W@f$ and the effective
 * couplings @f$\rho_Z^f@f$ and @f$\kappa_Z^f@f$ are controlled with the model
 * flags @ref StandardModelFlags "Mw", @ref StandardModelFlags "RhoZ" and
 * @ref StandardModelFlags "KappaZ" of StandardModel. For each flag, the
 * available schemes are as follows:
 *
 * @li NORESUM:&nbsp;&nbsp; No resummation is considered;
 * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
 * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
 * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
 * @li APPROXIMATEFORMULA:&nbsp;&nbsp; the approximate two-loop formula given
 * in EWSMApproximateFormulae class is employed. 
 *
 * The scheme APPROXIMATEFORMULA provides the most accurate SM predictions for
 * @f$M_W@f$ and @f$\kappa_Z^f@f$, while the approximate two-loop formula is
 * not available for @f$\rho_Z^f@f$. 
 *
 * See resumMw(), resumRhoZ() and resumKappaZ() for details on the other schemes.
 *
 *
 * <h3>Cashes</h3>
 *
 * This class contains caching methods for the following functions:
 * DeltaAlphaLepton(), DeltaAlpha(), Mw_SM(), GammaW_SM(), 
 * rhoZ_l_SM(), rhoZ_q_SM(), kappaZ_l_SM() and kappaZ_q_SM(), to improve the
 * performance of the Monte Carlo run. The caching methods are implemented
 * with the function checkSMparams(). 
 *
 * The use of the caching methods can be controlled with the model flag
 * @ref StandardModelFlags "CacheInEWSM" of StandardModel.
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
     * @brief An enumerated type representing perturbative orders of radiative
     * corrections to %EW precision observables.
     */
    enum orders_EW
    {
        EW1=0,///< One-loop of @f$\mathcal{O}(\alpha)@f$.
        EW1QCD1,///< Two-loop of @f$\mathcal{O}(\alpha\alpha_s)@f$.
        EW1QCD2,///< Three-loop of @f$\mathcal{O}(\alpha\alpha_s^2)@f$.
        EW2,///< Two-loop of @f$\mathcal{O}(\alpha^2)@f$.
        EW2QCD1,///< Three-loop of @f$\mathcal{O}(\alpha^2\alpha_s)@f$.
        EW3,///< Three-loop of @f$\mathcal{O}(\alpha^3)@f$.
        orders_EW_size ///< The size of this enum.
    };

    /**
     * @brief The number of the SM parameters that are relevant to the %EW
     * precision observalbles.
     * @details This constant is used for the cashing method.
     *
     * @sa checkSMparams()
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
     * @brief A set method to change 
     * the model flag @ref StandardModelFlags "CacheInEWSM" of StandardModel.
     * @details Setting CacheInEWSM to false, the caching methods defined in the
     * current class are not employed in numerical computations. The flag is set
     * to true in the constructor EWSM() by default.
     * @param[in] FlagCacheInEWSM true (false) if the caching methods are turned
     * on (off);
     *
     * @sa @ref StandardModelFlags "the description of the StandardModel flags"
     */
    void setFlagCacheInEWSM(bool FlagCacheInEWSM)
    {
        this->FlagCacheInEWSM = FlagCacheInEWSM;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMcache.
     * @return the pointer #myCache
     */
    EWSMcache* getMyCache() const
    {
        return myCache;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMOneLoopEW,
     * @return the pointer #myOneLoopEW
     */
    EWSMOneLoopEW* getMyOneLoopEW() const
    {
        return myOneLoopEW;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMApproximateFormulae.
     * @return the pointer #myApproximateFormulae
     */
    EWSMApproximateFormulae* getMyApproximateFormulae() const
    {
        return myApproximateFormulae;
    }

    /**
     * @brief A get method to retrieve the member pointer of type EWSMTwoFermionsLEP2.
     * @return the pointer #myTwoFermionsLEP2
     */
    EWSMTwoFermionsLEP2* getMyTwoFermionsLEP2() const 
    {
        return myTwoFermionsLEP2;
    }

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief A method to check whether the values of the parameters stored
     * in the given cache are all identical to those of the corresponding
     * model parameters in StandardModel.
     * @details This function is used for the cashing methods implemented in
     * the current class:
     * DeltaAlphaLepton(), DeltaAlpha(), Mw_SM(), rhoZ_l_SM(), rhoZ_q_SM(),
     * kappaZ_l_SM(), kappaZ_q_SM() and GammaW_SM().
     * When the values of the StandardModel parameters are updated in the Monte
     * Carlo run and differ from those stored in the given cache, Params_cache,
     * this function updates the cache, and returns false.
     * @param[in,out] Params_cache the cache of the parameters to be checked
     * @return a boolean that is true if the values of the parameters stored in
     * the given cache differ from those of the corresponding model parameters
     * in StandardModel
     *
     * @sa NumSMParams
     */
    bool checkSMparams(double Params_cache[]) const;

    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief @copybrief StandardModel::DeltaAlphaLepton()
     * @copydetails StandardModel::DeltaAlphaLepton()
     * 
     * @sa checkSMparams()
     * @attention This function uses the cache when the model flag
     * @ref StandardModelFlags "CacheInEWSM" of StandardModel is set to true.
     */
    double DeltaAlphaLepton(const double s) const;    

    /**
     * @brief @copybrief StandardModel::DeltaAlphaL5q()
     * @copydetails StandardModel::DeltaAlphaL5q()
     */
    double DeltaAlphaL5q() const;

    /**
     * @brief Top-quark contribution to the electromagnetic coupling @f$\alpha@f$,
     * denoted as @f$\Delta\alpha_{\mathrm{top}}(s)@f$.
     * @param[in] s invariant mass squared
     * @return @f$\Delta\alpha_{\mathrm{top}}(s)@f$
     */
    double DeltaAlphaTop(const double s) const;    

    /**
     * @brief @copybrief StandardModel::DeltaAlpha()
     * @copydetails StandardModel::DeltaAlpha()
     */
    double DeltaAlpha() const;

    /**
     * @brief @copybrief StandardModel::alphaMz()
     * @copydetails StandardModel::alphaMz()
     */
    double alphaMz() const;

   
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The square of the sine of the weak mixing angle @f$s_0^2@f$ defined
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
     * @brief The square of the cosine of the weak mixing angle @f$c_0^2@f$ defined
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
     * @details 
     * When the model flag @ref StandardModelFlags "Mw" of StandardModel is set
     * to APPROXIMATEFORMULA, the current function uses
     * the two-loop approximate formula in EWSMApproximateFormulae::Mw(),
     * which includes the full two-loop %EW contribution of
     * @f${\cal O}(\alpha^2)@f$ as well as the leading
     * @f${\cal O}(G_\mu^2\alpha_s m_t^4)@f$ and @f${\cal O}(G_\mu^3m_t^6)@f$
     * contributions.
     *
     * When the model flag @ref StandardModelFlags "Mw" is not set to APPROXIMATEFORMULA,
     * the @f$W@f$-boson mass is computed from @f$\Delta r(M_W)@f$ with an
     * iterative procedure. The target accuracy of the iterative calculation
     * is specified with the constant #Mw_error.
     * This function calls resumMw(), in which @f$M_W@f$ is computed with a given
     * @f$\Delta r@f$, equivalently with @f$\Delta\rho@f$ and
     * @f$\Delta r_{\mathrm{rem}}@f$
     * @return @f$M_{W,\mathrm{SM}}@f$ in GeV
     *
     * @sa resumMw()
     * @attention If the model flag @ref StandardModelFlags "CacheInEWSM"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    double Mw_SM() const;
    
    /** 
     * @brief The SM prediction for @f$\Delta r@f$ derived from that for the
     * @f$W@f$ boson mass.
     * @details 
     * If the model flag @ref StandardModelFlags "Mw" of StandardModel is set
     * to NORESUM or APPROXIMATEFORMULA, the quantity @f$\Delta r@f$ is computed
     * by using the following relation:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta r)\,.
     * @f]
     * Otherwise, the following relation is employed instead:
     * @f[
     * s_W^2 M_W^2 = \frac{\pi\,\alpha}{\sqrt{2}G_\mu(1-\Delta r)}\,,
     * @f]
     * where the resummation for @f$\Delta r@f$ is considered.
     * @return @f$\Delta r_{\mathrm{SM}}@f$
     * 
     * @sa The corresponding quantity in the complex-pole/fixed-width scheme
     * (instead of the experimental/running-widthr scheme)
     * is defined in DeltaRbar_SM().
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
     * in the experimental/running-width scheme
     * to that in the complex-pole/fixed-width scheme.
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
     * in the complex-pole/fixed-width scheme
     * to that in the experimental/running-width scheme.
     * @details The experimental mass @f$M_W@f$ is derived
     * @f[
     *   M_W = \overline{M}_W + \frac{\Gamma_{W}^2}{2\overline{M}_{W}}\,, 
     * @f]
     * where @f$\overline{M}_W@f$ is the mass parameter in the
     * complex-pole/fixed-width scheme @cite Bardin:1988xt, and @f$\Gamma_{W}@f$
     * is the @f$W@f$-boson width in the experimental/running-width scheme:
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
     * that for the @f$W@f$-boson mass.
     * @details The quantity @f$\Delta \overline{r}@f$ is computed by using
     * the following relation:
     * @f[
     * \overline{s}_W^2 \overline{M}_W^2
     * = \frac{\pi\,\alpha}{\sqrt{2}G_\mu}(1+\Delta \overline{r})\,,
     * @f]
     * where @f$\overline{M}_W@f$ and @f$\overline{s}_W@f$ are the @f$W@f$-boson
     * mass and the sine of the weak mixing angle in the complex-pole/fixed-width 
     * scheme @cite Bardin:1988xt.
     * @return @f$\Delta \overline{r}_{\mathrm{SM}}@f$
     *
     * @sa DeltaR_SM(), defining the corresponding quantity in the
     * experimental/running-width scheme.
     */
    double DeltaRbar_SM() const;


    ////////////////////////////////////////////////////////////////////////
    // SM contribution to the effective couplings

    /**
     * @brief The effective leptonic neutral-current coupling @f$\rho_Z^l@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\rho_Z^l@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumRhoZ(), while only the one-loop contribution is kept
     * in the imaginary part. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\rho_{Z,\,\mathrm{SM}}^l@f$
     * 
     * @sa resumRhoZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInEWSM"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    complex rhoZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\rho_Z^q@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\rho_Z^q@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumRhoZ(), while only the one-loop contribution is kept
     * in the imaginary part.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\rho_{Z,\,\mathrm{SM}}^q@f$
     *
     * @sa resumRhoZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInEWSM"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    complex rhoZ_q_SM(const QCD::quark q) const;    
    
    /**
     * @brief The effective leptonic neutral-current coupling @f$\kappa_Z^l@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\kappa_Z^l@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumKappaZ(), while only the one-loop contribution is kept
     * in the imaginary part.
     *
     * As a part of the two-loop %EW contribution, a correction associated with
     * the product of the imaginary part of @f$\Delta\alpha@f$ and that of
     * @f$\Pi_{Z\gamma}@f$ is included @cite Bardin:1999ak, @cite Bardin:1999yd :
     * @f{eqnarray}{
     * \Delta \kappa_Z^l
     * = - \frac{1}{s_W^2}\left( \frac{\alpha(M_Z^2)}{4\pi} \right)^2
     * {\rm Im}\,\overline{\Pi}_{\gamma\gamma}^{\rm fer}(M_Z^2)\,\,
     * {\rm Im}\,\overline{\Pi}_{Z\gamma}^{\rm fer}(M_Z^2)
     * = \frac{35\alpha^2(M_Z^2)}{18 s_W^2}\,
     * \left( 1 - \frac{8}{3}\, {\rm Re}(\kappa_Z^l) s_W^2 \right).
     * @f}
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\kappa_{Z,\,\mathrm{SM}}^l@f$
     *
     * @sa resumKappaZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInEWSM"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    complex kappaZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current coupling @f$\kappa_Z^q@f$ in the SM.
     * @details This function collects the radiative corrections to @f$\kappa_Z^q@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. The real part is computed
     * with the function resumKappaZ(), while only the one-loop contribution is kept
     * in the imaginary part.
     *
     * As a part of the two-loop %EW contribution, a correction associated with 
     * the product of the imaginary part of @f$\Delta\alpha@f$ and that of
     * @f$\Pi_{Z\gamma}@f$ is included @cite Bardin:1999ak, @cite Bardin:1999yd :
     * @f{eqnarray}{
     * \Delta \kappa_Z^q
     * = - \frac{1}{s_W^2}\left( \frac{\alpha(M_Z^2)}{4\pi} \right)^2
     * {\rm Im}\,\overline{\Pi}_{\gamma\gamma}^{\rm fer}(M_Z^2)\,\,
     * {\rm Im}\,\overline{\Pi}_{Z\gamma}^{\rm fer}(M_Z^2)
     * = \frac{35\alpha^2(M_Z^2)}{18 s_W^2}\,
     * \left( 1 - \frac{8}{3}\, {\rm Re}(\kappa_Z^q) s_W^2 \right).
     * @f}
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\kappa_{Z,\,\mathrm{SM}}^q@f$
     *
     * @sa resumKappaZ()
     * @attention If the model flag @ref StandardModelFlags "CacheInEWSM"
     * of StandardModel is set to true, the caching method implemented in the
     * current class is employed.
     */
    complex kappaZ_q_SM(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current vector coupling @f$g_V^l@f$ in the SM.
     * @details
     * @f[
     * g_V^l = g_A^l (1 - 4|Q_l|\kappa_Z^l s_W^2)\,.
     * @f]
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_{V,\,\mathrm{SM}}^l@f$
     */
    complex gVl_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current vector coupling @f$g_V^q@f$ in the SM.
     * @details
     * @f[
     * g_V^q = g_A^q (1 - 4|Q_q|\kappa_Z^q s_W^2)\,.
     * @f]
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_{V,\,\mathrm{SM}}^q@f$
     */
    complex gVq_SM(const QCD::quark q) const;

    /**
     * @brief The effective leptonic neutral-current axial-vector coupling @f$g_A^l@f$ in the SM.
     * @details
     * @f[
     * g_A^l = \sqrt{\rho_Z^l}\, I_3^l\,.
     * @f]
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$g_{A,\,\mathrm{SM}}^l@f$
     */
    complex gAl_SM(const StandardModel::lepton l) const;

    /**
     * @brief The effective quark neutral-current axial-vector coupling @f$g_A^q@f$ in the SM.
     * @details
     * @f[
     * g_A^q = \sqrt{\rho_Z^q}\, I_3^q\,.
     * @f]
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$g_{A,\,\mathrm{SM}}^q@f$
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
     * @brief Top-mass corrections to the @f$Zb\bar{b}@f$ vertex, denoted by
     * @f$\tau_b@f$. 
     * @details The large top-quark mass gives important corrections to the
     * %EW observables through the gauge-boson self-energies, i.e.,
     * @f$\Delta\rho@f$, and through the @f$Zb\bar{b}@f$ vertex. The latter
     * contribution is parameterised by the quantity @f$\tau_b@f$:
     * @f[
     * \tau_{b} =
     * -2\, X_t^{G_\mu}
     * \left[ 1 - \frac{\pi}{3}\alpha_s(M^2_t)
     *   + X_t^{G_\mu} \tau^{(2)}
     *     \left( \frac{M_t^2}{m_h^2} \right)
     * \right],
     * @f]
     * where the @f$O(G_\mu\alpha_s m_t^2)@f$ term was calculated in
     * @cite Fleischer:1992fq, @cite Buchalla:1992zm, @cite Degrassi:1993ij,
     * @cite Chetyrkin:1993jp, and the @f$O(G_\mu^2 m_t^4)@f$ term can be found
     * in @cite Barbieri:1992nz, @cite Barbieri:1992dq, @cite Fleischer:1993ub,
     * @cite Fleischer:1994cb.
     * @return @f$\tau_b@f$
     */
    double taub() const;    
    

    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Flavour non-universal vertex corrections to @f$\rho_Z^l@f$,
     * denoted by @f$\Delta\rho_Z^l@f$.
     * @details The non-universal contribution @f$\Delta\rho_Z^l@f$ is given by
     * @f[
     * \Delta \rho_Z^l = \rho_Z^l - \rho_Z^e
     * = \frac{\alpha}{2\pi s_W^2}\left(u_l - u_e\right),
     * @f]
     * where @f$u_l@f$ is defined as
     * @f[
     * u_l = \frac{3v_l^2+a_l^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^l(M_Z^2)
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_l = I_3^l - 2Q_l s_W^2@f$ and @f$a_l = I_3^l@f$ and the form factors,
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^l@f$.
     *
     * See @cite Ciuchini:2013pca and references therein. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\Delta\rho_Z^l@f$
     */
    complex rhoZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @brief Flavour non-universal vertex corrections to @f$\rho_Z^q@f$,
     * denoted by @f$\Delta\rho_Z^q@f$.
     * @details The non-universal contribution @f$\Delta\rho_Z^q@f$ is given by
     * @f[
     * \Delta \rho_Z^q = \rho_Z^q - \rho_Z^e
     * = \frac{\alpha}{2\pi s_W^2}\left(u_q - u_e\right),
     * @f]
     * where @f$u_q@f$ is defined as
     * @f[
     * u_q = \frac{3v_q^2+a_q^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^q(M_Z^2)
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_q = I_3^q - 2Q_q s_W^2@f$ and @f$a_q = I_3^q@f$ and the form factors,
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^q@f$.
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Delta\rho_Z^q@f$
     */
    complex rhoZ_q_SM_FlavorDep(QCD::quark q) const;

    /**
     * @brief Flavour non-universal vertex corrections to @f$\kappa_Z^l@f$,
     * denoted by @f$\Delta\kappa_Z^l@f$.
     * @details The non-universal contribution @f$\Delta\kappa_Z^l@f$ is given by
     * @f[
     * \Delta \kappa_Z^l = \kappa_Z^l - \kappa_Z^e
     * = \frac{\alpha}{4\pi s_W^2}
     *  \left( \frac{\delta_l^2-\delta_e^2}{4c_W^2}\,\mathcal{F}_Z(M_Z^2)
     *  -u_l+u_e\right),
     * @f]
     * where @f$u_l@f$ and @f$\delta_l@f$ are defined as
     * @f[
     * u_l = \frac{3v_l^2+a_l^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^l(M_Z^2)\,,
     * \qquad
     * \delta_l = v_l - a_l
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_l = I_3^l - 2Q_l s_W^2@f$ and @f$a_l = I_3^l@f$, and the form factors
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^l@f$.
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$\Delta\kappa_Z^l@f$
     */
    complex kappaZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @brief Flavour non-universal vertex corrections to @f$\kappa_Z^q@f$,
     * denoted by @f$\Delta\kappa_Z^q@f$.
     * @details The non-universal contribution @f$\Delta\kappa_Z^q@f$ is given by
     * @f[
     * \Delta \kappa_Z^q = \kappa_Z^q - \kappa_Z^e
     * = \frac{\alpha}{4\pi s_W^2}
     *  \left( \frac{\delta_q^2-\delta_e^2}{4c_W^2}\,\mathcal{F}_Z(M_Z^2)
     *  -u_q+u_e\right),
     * @f]
     * where @f$u_q@f$ and @f$\delta_q@f$ are defined as
     * @f[
     * u_q = \frac{3v_q^2+a_q^2}{4c_W^2}\mathcal{F}_Z(M_Z^2) + \mathcal{F}_W^q(M_Z^2)\,,
     * \qquad
     * \delta_q = v_q - a_q
     * @f]
     * with the tree-level vector and axial-vector couplings
     * @f$v_q = I_3^q - 2Q_q s_W^2@f$ and @f$a_q = I_3^q@f$, and the form factors
     * @f$\mathcal{F}_Z@f$ and @f$\mathcal{F}_W^q@f$.
     *
     * See @cite Ciuchini:2013pca and references therein.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$\Delta\kappa_Z^q@f$
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
     * model flag WithoutNonUniversalVC of StandardModel is set to true.
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
     * @brief @copybrief EWSMOneLoopEW::rho_GammaW_l()
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return @f$\rho^W_{ij}@f$
     *
     * @sa EWSMOneLoopEW::rho_GammaW_l()
     */    
    double rho_GammaW_l_SM(const StandardModel::lepton li, 
                           const StandardModel::lepton lj) const;
    
    /**
     * @brief @copybrief EWSMOneLoopEW::rho_GammaW_q()
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return @f$\rho^W_{ij}@f$
     *
     * @sa EWSMOneLoopEW::rho_GammaW_q()
     */    
    double rho_GammaW_q_SM(const QCD::quark qi, 
                           const QCD::quark qj) const;
    
    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a lepton pair.
     * @details
     * @f[
     * \Gamma^W_{ij} 
     * =
     * |U_{ij}|^2\,\frac{G_\mu M_W^3}{6\sqrt{2}\,\pi}\,\rho^W_{ij}
     * @f]
     * where @f$U@f$ denotes the %MNS matrix, and @f$\rho^W_{ij}@f$ represents
     * %EW radiative corrections.
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return @f$\Gamma^W_{ij}@f$
     *
     * @sa rho_GammaW_l_SM()
     * @attention Fermion masses are neglected. 
     */
    double GammaW_l_SM(const StandardModel::lepton li, 
                       const StandardModel::lepton lj) const;

    /**
     * @brief A partial decay width of the @f$W@f$ boson decay into a quark pair.
     * @details
     * @f[
     * \Gamma^W_{ij}
     * =
     * 3 |V_{ij}|^2\,\frac{G_\mu M_W^3}{6\sqrt{2}\,\pi}\,\rho^W_{ij}
     * \left( 1 + \frac{\alpha_s(M_W^2)}{\pi} \right).
     * @f]
     * where @f$V@f$ denotes the %CKM matrix, and @f$\rho^W_{ij}@f$ represents
     * %EW radiative corrections.
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return @f$\Gamma^W_{ij}@f$
     *
     * @sa rho_GammaW_q_SM()
     * @attention Fermion masses are neglected. 
     */
    double GammaW_q_SM(const QCD::quark qi, 
                       const QCD::quark qj) const;    
    
    /**
     * @brief The total widht of the @f$W@f$ boson in the SM.
     * @details
     * @f[
     * \Gamma_W
     * =
     * \Gamma^W_{\nu_e e} + \Gamma^W_{\nu_\mu \mu} + \Gamma^W_{\nu_\tau \tau}
     * + \Gamma^W_{ud} + \Gamma^W_{us} + \Gamma^W_{ub}
     * + \Gamma^W_{cd} + \Gamma^W_{cs} + \Gamma^W_{cb}\,.
     * @f]
     * @return @f$\Gamma_W@f$
     *
     * @sa GammaW_l_SM() and GammaW_q_SM()
     */
    double GammaW_SM() const;   


    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;///< A reference to an object of type StandardModel.

    /**
     * @brief A method to collect @f$\Delta\rho@f$ computed via subclasses.
     * @details This function collects @f$\Delta\rho@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[out] DeltaRho Array of @f$\Delta\rho@f$
     */
    void ComputeDeltaRho(const double Mw_i, double DeltaRho[orders_EW_size]) const;

    /**
     * @brief A method to collect @f$\Delta r_{\mathrm{rem}}@f$ computed via subclasses.
     * @details This function collects @f$\Delta r_{\mathrm{rem}}@f$
     * computed via EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
     * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[out] DeltaR_rem Array of @f$\Delta r_{\mathrm{rem}}@f$
     */
    void ComputeDeltaR_rem(const double Mw_i, double DeltaR_rem[orders_EW_size]) const;


    ////////////////////////////////////////////////////////////////////////         
private:

    /**
     * @brief An array of internal flags controlling the inclusions of higher-order
     * corrections.
     * @details These flags are prepared for debugging.
     * The flags are initialized in the constructor EWSM().
     */
    bool flag_order[orders_EW_size];
    
    EWSMcache* myCache;///< A pointer to an object of type EWSMcache.
    EWSMOneLoopEW* myOneLoopEW;///< A pointer to an object of type EWSMOneLoopEW.
    EWSMTwoLoopQCD* myTwoLoopQCD;///< A pointer to an object of type EWSMTwoLoopQCD.
    EWSMThreeLoopQCD* myThreeLoopQCD;///< A pointer to an object of type EWSMThreeLoopQCD.
    EWSMTwoLoopEW* myTwoLoopEW;///< A pointer to an object of type EWSMTwoLoopEW.
    EWSMThreeLoopEW2QCD* myThreeLoopEW2QCD;///< A pointer to an object of type EWSMThreeLoopEW2QCD.
    EWSMThreeLoopEW* myThreeLoopEW; ///< A pointer to an object of type EWSMThreeLoopEW.
    EWSMApproximateFormulae* myApproximateFormulae;///< A pointer to an object of type EWSMApproximateFormulae.
    
    EWSMTwoFermionsLEP2* myTwoFermionsLEP2;///< A pointer to an object of type EWSMTwoFermionsLEP2.
    
    
    ////////////////////////////////////////////////////////////////////////     
    //caches

    bool FlagCacheInEWSM;///< A flag for caching (true by default).

    /**
     * @brief A cache array of a set of SM parameters, used together with #DeltaAlphaLepton_cache.
     */
    mutable double DeltaAlphaLepton_params_cache[NumSMParams];
    mutable double DeltaAlphaLepton_cache;///< A cache of the value of @f$\Delta\alpha_{\mathrm{lept}}(M_Z^2)@f$.
    
    /**
     * @brief A cache array of a set of SM parameters, used together with #DeltaAlpha_cache.
     */
    mutable double DeltaAlpha_params_cache[NumSMParams];
    mutable double DeltaAlpha_cache;///< A cache of the value of @f$\Delta\alpha(M_Z^2)@f$.
    
    /**
     * @brief A cache array of a set of SM parameters, used together with #Mw_cache.
     */
    mutable double Mw_params_cache[NumSMParams];
    mutable double Mw_cache;///< A cache of the value of @f$M_W@f$.

    /**
     * @brief A cache array of a set of SM parameters, used together with #rhoZ_l_cache.
     */
    mutable double rhoZ_l_params_cache[6][NumSMParams];
    mutable complex rhoZ_l_cache[6];///< A cache of the value of @f$\rho_Z^l@f$.

    /**
     * @brief A cache array of a set of SM parameters, used together with #rhoZ_q_cache.
     */
    mutable double rhoZ_q_params_cache[6][NumSMParams];
    mutable complex rhoZ_q_cache[6];///< A cache of the value of @f$\rho_Z^q@f$.

    /**
     * @brief A cache array of a set of SM parameters, used together with #kappaZ_l_cache.
     */
    mutable double kappaZ_l_params_cache[6][NumSMParams];
    mutable complex kappaZ_l_cache[6];///< A cache of the value of @f$\kappa_Z^l@f$.

    /**
     * @brief A cache array of a set of SM parameters, used together with #kappaZ_q_cache.
     */
    mutable double kappaZ_q_params_cache[6][NumSMParams];
    mutable complex kappaZ_q_cache[6];///< A cache of the value of @f$\kappa_Z^q@f$.
    
    /**
     * @brief A cache array of a set of SM parameters, used together with #GammaW_cache.
     */
    mutable double GammaW_params_cache[NumSMParams];
    mutable double GammaW_cache;///< A cache of the value of @f$\Gamma_W@f$.
    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A method to compute the @f$W@f$-boson mass from 
     * @f$\Delta\rho@f$ and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes the @f$W@f$-boson mass without or with
     * resummation of @f$\Delta r@f$, 
     * depending on the model flag @ref StandardModelFlags "Mw" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER. 
     * The full two-loop %EW contribution is included in the case of "NORESUM",
     * while the large-@f$m_t@f$ expansion for the two-loop contribution is
     * adopted in the other cases.
     *
     * In the case of "NORESUM", the two-loop %EW contribution to @f$\Delta r@f$
     * is calculated via the function EWSMApproximateFormulae::DeltaR_TwoLoopEW_rem(),
     * given in the complex-pole/fixed-width scheme. The @f$W@f$-boson mass in
     * the complex-pole/fixed-width scheme, obtained from @f$\Delta r@f$, is
     * converted into the one in the experimental/running-width scheme with the
     * function MwFromMwbar().
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] DeltaR_rem Array of @f$\Delta r_{\mathrm{rem}}@f$
     * @return @f$M_W@f$
     */
    double resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
                   const double DeltaR_rem[orders_EW_size]) const;
    
    /**
     * @brief A method to compute the real part of the effective coupling
     * @f$\rho_Z^f@f$ from @f$\Delta\rho@f$, @f$\delta\rho_{\rm rem}^{f}@f$
     * and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes @f$\rho_Z^f@f$ without or with
     * resummation of @f$\Delta\rho@f$, depending on
     * the model flag @ref StandardModelFlags "RhoZ" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER. 
     * In all the cases, the two-loop %EW corrections are calculated in the
     * large-@f$m_t@f$ expansion. 
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] deltaRho_rem Array of @f$\delta\rho_{\rm rem}^{f}@f$
     * @param[in] DeltaRbar_rem Array of @f$\Delta \bar{r}_{\rm rem}@f$
     * @param[in] bool_Zbb true for @f$Zb\bar{b}@f$
     * @return @f$\mathrm{Re}(\rho_Z^f)@f$
     */
    double resumRhoZ(const double DeltaRho[orders_EW_size],
                     const double deltaRho_rem[orders_EW_size], 
                     const double DeltaRbar_rem, const bool bool_Zbb) const;
    
    /**
     * @brief A method to compute the real part of the effetvive coupling 
     * @f$\kappa_Z^f@f$ from @f$\Delta\rho@f$, @f$\delta\rho_{\rm rem}^{f}@f$
     * and @f$\Delta r_{\mathrm{rem}}@f$.
     * @details This function computes @f$\kappa_Z^f@f$ without or with
     * resummation of @f$\Delta\rho@f$, depending on
     * the model flag @ref StandardModelFlags "KappaZ" of StandardModel:
     *
     * @li NORESUM (recommended):&nbsp;&nbsp; no resummation is considered;
     * @li OMSI:&nbsp;&nbsp; the so-called OMS-I scheme is adopted;
     * @li INTERMEDIATE:&nbsp;&nbsp; an intermediate scheme between OMS-I and OMS-II is adopted;
     * @li OMSII:&nbsp;&nbsp; the so-called OMS-II scheme is adopted;
     * @li APPROXIMATEFORMULA:&nbsp;&nbsp; this is not applicable to the current function.
     *
     * where the OMS-I, INTERMEDIATE and OMS-II schemes are adopted in ZFITTER
     * @cite Bardin:1999yd (see also @cite Degrassi:1996mg, @cite Degrassi:1996ps,
     * @cite Degrassi:1999jd, @cite Bardin:1999ak),
     * and used for making comparisons to the outputs of ZFITTER. 
     * In all the cases, the two-loop %EW corrections are calculated in the
     * large-@f$m_t@f$ expansion.
     * @param[in] DeltaRho Array of @f$\Delta\rho@f$
     * @param[in] deltaKappa_rem Array of @f$\delta\kappa_{\rm rem}^{f}@f$
     * @param[in] DeltaRbar_rem Array of @f$\Delta \bar{r}_{\rm rem}@f$
     * @param[in] bool_Zbb true for @f$Zb\bar{b}@f$
     * @return @f$\mathrm{Re}(\kappa_Z^f)@f$
     */
    double resumKappaZ(const double DeltaRho[orders_EW_size],
                       const double deltaKappa_rem[orders_EW_size],
                       const double DeltaRbar_rem, const bool bool_Zbb) const;
    
};

#endif	/* EWSM_H */

