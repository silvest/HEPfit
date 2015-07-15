/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AFBLEPTON_H
#define	AFBLEPTON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class AFBlepton
 * @ingroup EW 
 * @brief An observable class for the forward-backward asymmetry in
 * @f$e^+ e^-\to Z\to \ell^+ \ell^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the forward-backward asymmetry for 
 * an @f$e^+ e^-\to Z\to \ell^+ \ell^-@f$ process at the @f$Z@f$ pole in terms
 * of the asymmetry parameters @f$\mathcal{A}_f@f$, computed in the EW class:
 * @f[
 * A_{\mathrm{FB}}^{0,\ell}=\frac 34 \mathcal{A}_e \mathcal{A}_\ell\,,
 * @f]
 * where @f$\ell@f$ denotes a charged lepton, and lepton-flavour universality
 * is assumed.
 *
 * @sa EW_NPZff::AFBlepton() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class AFBlepton : public ThObservable {
public:
 
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AFBlepton(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The forward-backward asymmetry for @f$e^+ e^-\to Z \to \ell^+ \ell^-@f$
     * at the @f$Z@f$ pole, @f$A^{0,\ell}_{\mathrm{FB}}@f$
     * @return @f$A^{0,\ell}_{\mathrm{FB}}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* AFBLEPTON_H */

