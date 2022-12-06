/* 
 * Copyright (C) 2022 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef XI_H
#define XI_H

#include "ThObservable.h"

/**
 * @class xi
 * @ingroup Flavour
 * @brief A class for @f$\xi \equiv \frac{F_{B_s}\sqrt{B_{B_s}}}{F_{B_d}\sqrt{B_{B_d}}}@f$, relevant for @f$B_{s,d}@f$ mesons mixing in the Standard Model
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to retrieve the value of
 * @f$\xi@f$.
 */
class xi : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    xi(const StandardModel& SM_i);

    /**
     * 
     * @return value of @f$\xi@f$ 
     */
    virtual double computeThValue();
};

#endif /* XI_H */

