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
 * @details  
 */
class ThObservable {
public:
    ThObservable(const ThObsType& ObsType_i);
    ThObservable(const ThObservable& orig);
    virtual ~ThObservable();
    virtual double getThValue() = 0;

    /**
     * The conversion factor from GeV^{-2} to nb. 
     */
    static const double GeVminus2_to_nb;

protected:
    
    /**
     * A reference to an object of ThObsType class. 
     */
    const ThObsType& ObsType;
    
    /**
     * A reference to an object of StandardModel class. 
     */
    const StandardModel& SM;
};

#endif	/* THOBSERVABLE_H */
