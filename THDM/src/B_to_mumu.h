/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef B_TO_MUMU_H
#define	B_TO_MUMU_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class B_to_mumu
 * @ingroup THDM
 * @brief A class for @f$B^0_{s/d} \to \mu^+\mu^-@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details We take the expressions from @cite Cheng:2015yfu, linking to our implementation of the Standard Model value.
 */
class B_to_mumu : public ThObservable {
public:

    /**
     * @brief Constructor of the B_to_mumu class.
     */
    B_to_mumu(const StandardModel& SM_i);

    /**
     * @brief Destructor of the BDtaunu class.
     */
    virtual ~B_to_mumu();

    /**
     * @brief Empty function.
     */
    double computeThValue();

    /**
     * @brief Determination of the Wilson coefficient @f$C_{10}@f$
     */
    double computeC10();

    /**
     * @brief Determination of the Wilson coefficient @f$C_P@f$
     */
    double computeCP();

    /**
     * @brief Determination of the Wilson coefficient @f$C_S@f$
     */
    double computeCS();

    const THDM * myTHDM;

protected:

private:
};

/**
 * @class BR_BsmumuTHDM
 * @ingroup THDM
 * @brief A class for @f$B^0_s \to \mu^+\mu^-@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Calculates the branching ratio of @f$B^0_s \to \mu^+\mu^-@f$ according to Eq. (2.3) of @cite Cheng:2015yfu.
 */
class BR_BsmumuTHDM: public B_to_mumu {
public:

    /**
     * @brief Constructor of the BR_BsmumuTHDM class.
     */
    BR_BsmumuTHDM(const StandardModel& SM_i);

    /**
     * @brief Branching ratio of @f$B^0_s \to \mu^+\mu^-@f$ decays.
     */
    double computeThValue();
    
private:
};

/**
 * @class BR_BdmumuTHDM
 * @ingroup THDM
 * @brief A class for @f$B^0_d \to \mu^+\mu^-@f$ decays in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Calculates the branching ratio of @f$B^0_d \to \mu^+\mu^-@f$ according to Eq. (2.3) of @cite Cheng:2015yfu.
 */
class BR_BdmumuTHDM: public B_to_mumu {
public:

    /**
     * @brief Constructor of the BR_BdmumuTHDM class.
     */
    BR_BdmumuTHDM(const StandardModel& SM_i);

    /**
     * @brief Branching ratio of @f$B^0_d \to \mu^+\mu^-@f$ decays.
     */
    double computeThValue();
    
private:
};

#endif	/* B_TO_MUMU_H */
