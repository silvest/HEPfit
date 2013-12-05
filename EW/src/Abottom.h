/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ABOTTOM_H
#define	ABOTTOM_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @addtogroup EW
 * @brief A project for electroweak precision observables.
 * @details This project includes the classes used for the computation of
 * the theory predictions for the W mass and width, as well as electroweak 
 * precision pseudo observables at the Z pole.
 * These pseudo observables cannot be directly measured, but can be extracted from real
 * observables upon substraction of initial-state QED radiation and final-state QED/QCD 
 * corrections.
 * 
 * @{
 */

/**
 * @class Abottom 
 * @brief A class for the left-right asymmetry of @f$Z\to b\bar{b}@f$ at the 
 * @f$Z@f$ pole.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry of the 
 * @f$Z\to b\bar{b}@f$ channel at the @f$Z@f$ pole. 
 */
class Abottom : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    Abottom(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @return The left-right asymmetry of @f$Z\to b\bar{b}@f$ at the @f$Z@f$ pole.
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    const EW& myEW;
};

/** 
 * @}
 */

#endif	/* ABOTTOM_H */

