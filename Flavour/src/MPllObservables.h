/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MPLLOBSERVABLES_H
#define	MPLLOBSERVABLES_H

#include "MPll.h"
#include <StandardModel.h>
#include <ThObservable.h>




/**
 * @class BR_MPll
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<BR>@f$ in @f$M \to P l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<BR>@f$ in 
 * @f$M \to P l^+l^-@f$ in terms of the binned CP averaged helicity coefficients 
 * @f$<\Sigma_i>@f$, computed in the MPll class, and the meson width @f$W_M@f$:
 * @f[
 * <BR>= \frac {<\Gamma'>}{W_M} \,.
 * @f]
 */
class BR_MPll : public ThObservable{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    BR_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
     * @brief A method to compute the binned observable @f$<BR>@f$ in @f$M \to P l^+l^-@f$ in a generic bin.
     * @param[in] qmin minimal value of the bin
     * @param[in] qmax maximal value of the bin
     * @param[in] lep final leptons of the decay
     * @return @f$<BR>_{[qmin,qmax]}@f$
     */
    double computeBR_MPll(double qmin, double qmax, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<BR>@f$ in @f$M \to P l^+l^-@f$.
    * @return @f$<BR>@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson pseudoscalar; /**< Final pseudoscalar meson type. */
};


/**
 * @class R_MPll
 * @ingroup Flavour
 * @brief A class for the binned observable ratio @f$\frac {<BR>_{M \to P l_1^+l_1^-}}{<BR>_{M \to P l_2^+l_2^-}}@f$ in @f$M \to P l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable ratio 
 * @f$\frac {<BR>_{M \to P l_1^+l_1^-}}{<BR>_{M \to P l_2^+l_2^-}}@f$ 
 * in @f$M \to P l^+l^-@f$ in terms of the binned observables @f$<\Gamma'>_{M \to P l_i^+l_i^-}@f$:
 * @f[
 * <R>_{M \to P l^+ l^-}= \frac {<\Gamma'>_{M \to P l_1^+l_1^-}}{<\Gamma'>_{M \to P l_2^+l_2^-}} \,.
 * @f]
 */
class R_MPll : public BR_MPll{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_1 first final leptons of the decay
     * @param[in] lep_1 second final leptons of the decay
     */
    R_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_1, StandardModel::lepton lep_2);
    
    /**
    * @brief The binned observable ratio @f$\frac {<BR>_{M \to P l_1^+l_1^-}}{<BR>_{M \to P l_2^+l_2^-}}@f$ in @f$M \to P l^+l^-@f$.
    * @return @f$<R>_{M \to P l^+ l^-}@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep1; /**< First final leptons type. */
    StandardModel::lepton lep2; /**< Second final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson pseudoscalar; /**< Final pseudoscalar meson type. */
};


/**
 * @class ACP_MPll
 * @ingroup Flavour
 * @brief A class for the binned observable @f$<A_{CP}>@f$ in @f$M \to P l^+l^-@f$. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the binned observable @f$<A_{CP}>@f$ in 
 * @f$M \to P l^+l^-@f$ in terms of the binned CP asymmetry helicity coefficients 
 * @f$<\Delta_i>@f$, computed in the MVll class, and the @f$\Gamma'@f$:
 * @f[
 * <A_{CP}>= -\frac {3<\Delta_{1c}> - <\Delta_{2c}> + 2(3<\Delta_{1s}> - <\Delta_{2s}>)}{4<\Gamma'>} \,.
 * @f]
 */
class ACP_MPll : public BR_MPll{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] pseudoscalar_i final pseudoscalar meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    ACP_MPll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i);
    
    /**
    * @brief The binned observable @f$<A_{CP}>@f$ in @f$M \to P l^+l^-@f$.
    * @return @f$<A_{CP}>@f$
    */
    double computeThValue ();
    
private:
    StandardModel::lepton lep; /**< Final leptons type. */
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson pseudoscalar; /**< Final pseudoscalar meson type. */
};


#endif	/* MPLLOBSERVABLES_H */

    