/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAKAPPAZB_H
#define	DELTAKAPPAZB_H

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
        else
            throw std::runtime_error("ERROR: deltaKappaZb::getThValue() cannot be used with flag NotLinearizedNP=1.");
    };
    
private:

};

#endif	/* DELTAKAPPAZB_H */

