/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BBD_H
#define BBD_H

#include "ThObservable.h"

/**
 * @class BBd
 * @ingroup Flavour
 * @brief A class for @f$B_{B_d}@f$, the B-parameter for @f$B_d@f$ mesons
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to retrieve the value of
 * @f$B_{B_d}@f$.
 */
class BBd : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    BBd(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$B_{B_d}@f$ 
     */
    virtual double computeThValue();
};

#endif /* BBD_H */
