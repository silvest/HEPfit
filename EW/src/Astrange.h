/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ASTRANGE_H
#define	ASTRANGE_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Astrange
 * @ingroup EW 
 * @brief An observable class for the left-right asymmetry in
 * @f$e^+e^-\to Z\to s\bar{s}@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry for the
 * @f$e^+e^-\to Z\to s\bar{s}@f$ process at the @f$Z@f$ pole.
 *
 * @sa EW_NPZff::Astrange() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Astrange : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Astrange(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The left-right asymmetry for @f$e^+e^-\to Z\to s\bar{s}@f$
     * at the @f$Z@f$ pole, @f$\mathcal{A}_s@f$.
     * @return @f$\mathcal{A}_s@f$
     */
    double computeThValue();

    
private:


};

#endif	/* ASTRANGE_H */

