/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTARHOZB_H
#define	DELTARHOZB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

/**
 * @class deltaRhoZb
 * @brief A class for @f$\delta\rho_Z^b@f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaRhoZb : public ThObservable {
public:

    deltaRhoZb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( SM.rhoZ_q(SM.BOTTOM).real() 
                 - SM.StandardModel::rhoZ_q(SM.BOTTOM).real() );
    };
    
private:

};

#endif	/* DELTARHOZB_H */

