/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALPHAMZ_H
#define	ALPHAMZ_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class AlphaEmMz 
 * @brief An observable class for the em constant evaluated at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the em constant evaluated at the @f$Z@f$ pole.
 *
 */
class AlphaEmMz : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    AlphaEmMz(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The em constant evaluated at the @f$Z@f$ pole, @f$\alpha(M_Z)@f$.
     * @return @f$\alpha(M_Z)@f$
     */
    double computeThValue();

    
private:


};

#endif	/* ALPHAMZ_H */

