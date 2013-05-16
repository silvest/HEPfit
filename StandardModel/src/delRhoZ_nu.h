/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELRHOZ_NU_H
#define	DELRHOZ_NU_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"
#include "EWSM.h"

/**
 * @class delRhoZ_nu
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\delta\rho_Z^\nu@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class delRhoZ_nu : public ThObservable {
public:

    delRhoZ_nu(const ThObsType& ObsType)
    : ThObservable(ObsType)
    {
    };

    double getThValue()
    {
        return SM.getDelRhoZ_nu();
    };

private:

};

#endif	/* DELRHOZ_NU_H */
