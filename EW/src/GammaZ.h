/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAZ_H
#define	GAMMAZ_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class GammaZ
 * @ingroup EW 
 * @brief A class for the total decay width of the \f$Z\f$ boson 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the total decay width of the \f$Z\f$ boson 
 * (in GeV),
 * \f[\Gamma_Z=\Gamma_h+\Gamma_e+\Gamma_\mu+\Gamma_\tau+\Gamma_\mathrm{inv},\f] 
 * with @f$\Gamma_h=\sum_{q\not =t}\Gamma_q@f$ the total hadronic width, and 
 * @f$\Gamma_\mathrm{inv}=3\Gamma_\nu@f$ the invisible width.
 */
class GammaZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    GammaZ(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The theory prediction for \f$\Gamma_Z\f$.
     * @return the total width of the \f$Z\f$ boson in GeV
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;///< A reference to an object of the EW class.
};

#endif	/* GAMMAZ_H */

