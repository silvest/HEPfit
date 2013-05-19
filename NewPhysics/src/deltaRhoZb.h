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
        if (SM.IsFlagApproximateGqOverGb() 
                && !SM.IsFlagRhoZbFromGuOverGb() 
                && !SM.IsFlagRhoZbFromGdOverGb()
                && !SM.IsFlagTestSubleadingTwoLoopEW())
            // SM prediction for rho_Z^b is missing!
            throw std::runtime_error("deltaRhoZb::getThValue() cannot be used!");
        else
        if (SM.IsFlagNotLinearizedNP())
            return ( SM.rhoZ_q(SM.BOTTOM).real() 
                     - SM.StandardModel::rhoZ_q(SM.BOTTOM).real() );
        else {
            complex gAb = SM.StandardModel::gAq(SM.BOTTOM) + SM.deltaGAq(SM.BOTTOM);
            double I3b = SM.getQuarks(SM.BOTTOM).getIsospin();
            double rhoZb_full = (gAb*gAb/I3b/I3b).real();
            return ( rhoZb_full
                     - SM.StandardModel::rhoZ_q(SM.BOTTOM).real() );
        }
    };
    
private:

};

#endif	/* DELTARHOZB_H */

