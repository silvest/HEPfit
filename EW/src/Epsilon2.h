/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILON2_H
#define	EPSILON2_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Epsilon2 
 * @ingroup EW
 * @brief An observable class for the @f$\epsilon_2@f$ parameter.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details To be added
 */
class Epsilon2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i a reference to an object of type EW
     */
    Epsilon2(const StandardModel& SM_i) 
    : ThObservable(SM_i)  
    {
    };

    /**
     * @brief To be added
     * @return @f$\epsilon_2@f$
     */
    double computeThValue() 
    {
        return SM.epsilon2();
    }

    
private:


};

#endif	/* EPSILON2_H */
