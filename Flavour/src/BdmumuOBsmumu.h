/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BDMUMUOBSMUMU_H
#define BDMUMUOBSMUMU_H

class StandardModel;
#include "ThObservable.h"

/**
 * @class BdmumuOBsmumu
 * @ingroup Flavour
 * @brief A class for the ratio of the branching ratios of \f$B_d\to \mu\mu\f$ over \f$B_s\to \mu\mu\f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * the  ratio of the branching ratios of \f$B_d\to \mu\mu\f$ over \f$B_s\to \mu\mu\f$.
 */
class BdmumuOBsmumu : public ThObservable {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    BdmumuOBsmumu(const StandardModel& SM_i);

    /**
     * 
     * @return theoretical value of |\f$ BR(B_d\to \mu\mu) / BR(B_s\to \mu\mu) \f$|
     */
    double computeThValue();
};

#endif /* BDMUMUOBSMUMU_H */

