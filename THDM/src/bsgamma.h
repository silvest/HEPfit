/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMATHDM_H
#define	BSGAMMATHDM_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class bsgammaTHDM
 * @ingroup THDM
 * @brief The @f$b\to s\gamma@f$ branching fraction from tabled values.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The branching ratio @f$BR(B\to X_s \gamma)@f$ for the type II is read from a table and interpolated.
 * The tabled values are state-of-the-art calculations at NNLO from @cite Misiak:2015xwa
 * in a grid of 11 @f$\log_{10}(\tan \beta)@f$ values between -4/3 and 2 times 101 @f$\log_{10}(m_{H^+})@f$ values between 1 and 4.
 * The theoretical error is evaluated at every point and found to be almost constant at 7%, 
 * which is accounted for by the nuisance parameter bsgamma_theoryerror.
 */
class bsgammaTHDM : public ThObservable {
public:
    bsgammaTHDM(const StandardModel& SM_i);
    virtual ~bsgammaTHDM();

    /**
     * @brief Calculates the @f$BR(B\to X_s \gamma)@f$ value in the THDM by linear interpolation of tabled NNLO values.
     * @return @f$BR(B\to X_s \gamma)@f$.
     */
    double computeThValue();

protected:
    THDMcache * mycache;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;
};

#endif	/* BSGAMMATHDM_H */
