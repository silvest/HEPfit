/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THOBSERVABLE_H
#define	THOBSERVABLE_H


#include <StandardModel.h>
#include "ThObsType.h"

/**
 * @class ThObservable
 * @ingroup Observable
 * @brief A class for a model prediction of an observable. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class provides a base for the computation of the values
 * of different theory observables.
 */
class ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] ObsType_i a reference to an object of ThObsType() class
     */
    ThObservable(const ThObsType& ObsType_i);
    
    /**
     * @brief The copy constructor.
     */
    ThObservable(const ThObservable& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ThObservable();
    
    /**
     * @brief A memeber to be overloaded by the respective theory observable.
     * class that calculates the value of the observable
     */
    virtual double computeThValue() = 0;

    /**
     * @brief The conversion factor from @f$ GeV^{-2} @f$ to @f$nb@f$.
     */
    static const double GeVminus2_to_nb;

protected:
    
    const ThObsType& ObsType; ///< A reference to an object of ThObsType class.
    const StandardModel& SM; ///< A reference to an object of StandardMode class.
};

#endif	/* THOBSERVABLE_H */
