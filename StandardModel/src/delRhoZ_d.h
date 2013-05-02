/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELRHOZ_D_H
#define	DELRHOZ_D_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"
#include "EWSM.h"

/**
 * @class delRhoZ_d
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\delta\rho_Z^d@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class delRhoZ_d : public ThObservable {
public:

    delRhoZ_d(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( SM.StandardModel::rhoZ_q(SM.DOWN).real()
                 - SM.getEWSM()->rhoZ_q_SM(SM.DOWN).real() );
    };

private:

};

#endif	/* DELRHOZ_D_H */

