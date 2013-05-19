/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAKAPPAZB_H
#define	DELTAKAPPAZB_H

#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class deltaKappaZb
 * @brief A class for @f$\delta\kappa_Z^b@f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaKappaZb : public ThObservable {
public:

    deltaKappaZb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        if (SM.IsFlagNotLinearizedNP())
            return ( SM.kappaZ_q(SM.BOTTOM).real() 
                     - SM.StandardModel::kappaZ_q(SM.BOTTOM).real() );
        else {
            complex gVb = SM.StandardModel::gVq(SM.BOTTOM) + SM.deltaGVq(SM.BOTTOM);
            complex gAb = SM.StandardModel::gAq(SM.BOTTOM) + SM.deltaGAq(SM.BOTTOM);
            double Qb = SM.getQuarks(SM.BOTTOM).getCharge();
            double kappaZb_full = (1.0 - (gVb/gAb).real())/(4.0*fabs(Qb)*SM.sW2());
            return ( kappaZb_full
                     - SM.StandardModel::kappaZ_q(SM.BOTTOM).real() );
        }
    };
    
private:

};

#endif	/* DELTAKAPPAZB_H */

