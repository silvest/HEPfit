/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PHID_H
#define	PHID_H

#include "ThObservable.h"

/**
 * @class phid
 * @ingroup Flavour
 * @brief A class for @f$\phi_{d}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of the phase of the @f$B_d@f$ mixing amplitude relative to the @f$b\to c\bar{c}s@f$ decay amplitude. The phase is given in radians. 
 * @f$\phi_d@f$.
 */
class phid : public ThObservable {
public:
    
    /**
     * @brief Constructor. 
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    
    phid(const StandardModel& SM_i);
    
    /**
     *
     * @return theoretical value of @f$\phi_{d}@f$
     */
    virtual double computeThValue();
    
};

#endif	/* PHID_H */

