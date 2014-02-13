/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class PtauPol
 * @ingroup EW 
 * @brief An observable class for the @f$\tau@f$ polarization in 
 * @f$e^+e^-\to Z\to \tau^+\tau^-@f$ at the @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\tau@f$ polarization in the 
 * @f$e^+e^-\to Z\to \tau^+\tau^-@f$ process at the @f$Z@f$ pole,
 * @f[
 * P_\tau^{\mathrm{pol}}=\mathrm{A}_\tau\,.
 * @f]
 */
class PtauPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    PtauPol(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The @f$\tau@f$ polarization in @f$e^+e^-\to Z\to \tau^+\tau^-@f$
     * at the @f$Z@f$ pole, @f$P_\tau^{\mathrm{pol}}@f$.
     * @return @f$P_\tau^{\mathrm{pol}}@f$
     */
    double computeThValue();

    
private:
    const EW& myEW;///< A reference to an object of type EW.

};

#endif	/* PTAUPOL_H */

