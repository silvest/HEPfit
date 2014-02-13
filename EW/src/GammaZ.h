/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
 * @brief An observable class for the total decay width of the @f$Z@f$ boson. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the total decay width of the @f$Z@f$
 * boson,
 * @f[
 * \Gamma_Z
 * = 3\,\Gamma_\nu + \Gamma_{e} + \Gamma_{\mu} + \Gamma_{\tau} + \Gamma_h\,,
 * @f]
 * where @f$\Gamma_h=\sum_{q\neq t}\Gamma_q@f$ is the total hadronic width.
 */
class GammaZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    GammaZ(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The total decay width of the @f$Z@f$ boson, @f$\Gamma_Z@f$,
     * in units of GeV.
     * @return @f$\Gamma_Z@f$ in units of GeV
     */
    double computeThValue();

    
private:
    const EW& myEW;///< A reference to an object of type EW.

};

#endif	/* GAMMAZ_H */

