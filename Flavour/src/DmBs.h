/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBS_H
#define	DMBS_H

#include "ThObservable.h"

/**
 * @class DMBs
 * @ingroup Flavour
 * @brief A class for @f$\Delta M_{B_s}@f$, the difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * @f$B_d@f$ system.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\Delta M_{B_s}@f$. The short distance contributions to this
 * parameter is calculable both in the SM and in NP models. In the SM
 * the short distance contribution is dominated by the top loop in the
 * box diagram.
 */
class DmBs : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    DmBs(const StandardModel& SM_i) : ThObservable(SM_i) {};

    /**
     *
     * @return theoretical value of @f$\Delta M_{B_s}@f$
     */

    double computeThValue();
};


/**
 * @class RmBs
 * @ingroup Flavour
 * @brief A class for @f$R(\Delta M_{B_s})@f$, the ratio of the total difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * @f$B_d@f$ system, over the SM difference.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$R(\Delta M_{B_s})@f$. The short distance contributions to this
 * parameter is calculable both in the SM and in NP models. In the SM
 * the short distance contribution is dominated by the top loop in the
 * box diagram.
 */
class RmBs : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    RmBs(const StandardModel& SM_i) : ThObservable(SM_i){};

    /**
     *
     * @return theoretical value of @f$R (\Delta M_{B_s})@f$
     */

    double computeThValue();
};

#endif	/* DMBS_H */
