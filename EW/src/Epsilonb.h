/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONB_H
#define	EPSILONB_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Epsilonb 
 * @ingroup EW
 * @brief An observable class for the @f$\epsilon_b@f$ parameter.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details To be added
 */
class Epsilonb : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Epsilonb(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief To be added
     * @return @f$\epsilon_b@f$
     */
    double computeThValue() 
    {
        return SM.epsilonb();
    }

    
private:


};

#endif	/* EPSILONB_H */
