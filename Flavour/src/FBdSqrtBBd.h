/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FBDSQRTBBD_H
#define FBDSQRTBBD_H

#include "ThObservable.h"

/**
 * @class FBdSqrtBBd
 * @ingroup Flavour
 * @brief A class for @f$F_{B_d}\sqrt{B_{B_d}}@f$, relevant for @f$B_d@f$ mesons mixing in the Standard Model
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to retrieve the value of
 * @f$F_{B_d}\sqrt{B_{B_d}}@f$.
 */
class FBdSqrtBBd : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    FBdSqrtBBd(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$F_{B_d}\sqrt{B_{B_d}}@f$ 
     */
    virtual double computeThValue();
};

#endif /* FBDSQRTBBD_H */

