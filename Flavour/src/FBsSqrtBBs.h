/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FBSSQRTBBS_H
#define FBSSQRTBBS_H

#include "ThObservable.h"

/**
 * @class FBsSqrtBBs
 * @ingroup Flavour
 * @brief A class for @f$F_{B_s}\sqrt{B_{B_s}}@f$, relevant for @f$B_s@f$ mesons mixing in the Standard Model
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to retrieve the value of
 * @f$F_{B_s}\sqrt{B_{B_s}}@f$.
 */
class FBsSqrtBBs : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    FBsSqrtBBs(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$F_{B_s}\sqrt{B_{B_s}}@f$ 
     */
    virtual double computeThValue();
};

#endif /* FBSSQRTBBS_H */

