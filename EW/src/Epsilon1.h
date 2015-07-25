/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON1_H
#define	EPSILON1_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Epsilon1 
 * @ingroup EW
 * @brief An observable class for the @f$\epsilon_1@f$ parameter.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details To be added
 */
class Epsilon1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Epsilon1(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief To be added
     * @return @f$\epsilon_1@f$
     */
    double computeThValue() 
    {
        return SM.epsilon1();
    }

    
private:


};

#endif	/* EPSILON1_H */
