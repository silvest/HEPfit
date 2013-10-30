/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELRHOZ_B_H
#define	DELRHOZ_B_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"
#include "EWSM.h"

/**
 * @class delRhoZ_b
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\delta\rho_Z^b@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class delRhoZ_b : public ThObservable {
public:

    delRhoZ_b(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double computeThValue()
    {
        return SM.getEWSM()->delRhoZ_q(SM.BOTTOM);
    };

private:

};

#endif	/* DELRHOZ_B_H */

