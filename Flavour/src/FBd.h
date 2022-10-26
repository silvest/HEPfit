/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FBD_H
#define FBD_H

#include "ThObservable.h"

/**
 * @class FBd
 * @ingroup Flavour
 * @brief A class for @f$F_{B_d}@f$, the @f$B_d@f$ decay constant
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to retrieve the value of
 * @f$F_{B_d}@f$.
 */
class FBd : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    FBd(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$F_{B_d}@f$ 
     */
    virtual double computeThValue();
};

#endif /* FBD_H */
