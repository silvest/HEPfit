/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMBD_H
#define	DMBD_H

#include "ThObservable.h"

/**
 * @class DmBd
 * @ingroup Flavour
 * @brief A class for @f$\Delta M_{B_d}@f$, the difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * @f$B_d@f$ system.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\Delta M_{B_d}@f$. The short distance contributions to this
 * parameter is calculable both in the SM and in NP models. In the SM
 * the short distance contribution is dominated by the top loop in the
 * box diagram.
 */
class DmBd : public ThObservable {
public:
    
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    DmBd(const StandardModel& SM_i);
    
    /**
     *
     * @return theoretical value of @f$\Delta M_{B_d}@f$
     */
    
    double computeThValue();
};


/**
 * @class RmBd
 * @ingroup Flavour
 * @brief A class for @f$R(\Delta M_{B_d})@f$, the ratio of the total difference in mass
 * of between the heavy and light mass eigenstates of the neutral
 * @f$B_d@f$ system, over the SM difference.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$R(\Delta M_{B_d})@f$. The short distance contributions to this
 * parameter is calculable both in the SM and in NP models. In the SM
 * the short distance contribution is dominated by the top loop in the
 * box diagram.
 */
class RmBd : public ThObservable {
public:

    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    RmBd(const StandardModel& SM_i) : ThObservable(SM_i){};

    /**
     *
     * @return theoretical value of @f$R (\Delta M_{B_s})@f$
     */

    double computeThValue();
};

/**
 * @class CBd
 * @ingroup Flavour
 * @brief Model parameter CBd from the NP model NPDF2
 * @detail The ratio of the absolute value of the $B_d$ mixing amplitude over the Standard Model value
 */
class CBd : public ThObservable {
public:

    CBd(const StandardModel& SM_i);

    double computeThValue();
};

/**
 * @class PhiBd
 * @ingroup Flavour
 * @brief Model parameter PhiBd from the NP model NPDF2
 * @detail Half the relative phase of the $B_d$ mixing amplitude w.r.t. the Standard Model one
 */
class PhiBd : public ThObservable {
public:

    PhiBd(const StandardModel& SM_i);

    double computeThValue();
};

#endif	/* DMBD_H */
