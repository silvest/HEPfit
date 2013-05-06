/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELRHOZ_U_H
#define	DELRHOZ_U_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"
#include "EWSM.h"

/**
 * @class delRhoZ_u
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\delta\rho_Z^u@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class delRhoZ_u : public ThObservable {
public:

    delRhoZ_u(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getEWSM()->delRhoZ_q(SM.UP);
    };

private:

};

#endif	/* DELRHOZ_U_H */

