/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON3_H
#define	EPSILON3_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Epsilon3 
 * @ingroup EW
 * @brief An observable class for the @f$\epsilon_3@f$ parameter.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details To be added
 */
class Epsilon3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Epsilon3(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief To be added
     * @return @f$\epsilon_3@f$
     */
    double computeThValue() 
    {
        return SM.epsilon3();
    }

    
private:


};

#endif	/* EPSILON3_H */
