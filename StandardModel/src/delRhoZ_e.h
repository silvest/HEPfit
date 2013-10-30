/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELRHOZ_E_H
#define	DELRHOZ_E_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"
#include "EWSM.h"

/**
 * @class delRhoZ_nu
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\delta\rho_Z^\ell@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class delRhoZ_e : public ThObservable {
public:

    delRhoZ_e(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double computeThValue()
    {
        return SM.getDelRhoZ_e();
    };

private:

};

#endif	/* DELRHOZ_E_H */