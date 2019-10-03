/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef C2BETA_H
#define C2BETA_H

#include "ThObservable.h"
#include "AmpDB2.h"

/**
 * @class C2beta
 * @ingroup Flavour
 * @brief A class for @f$cos2\beta@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$cos2\beta@f$.
 */
class C2beta: public ThObservable, AmpDB2 {
public:
    /**
     * @brief Constructor. 
     * @param[in] SM_i
     */
    
    C2beta(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) {};
    
    /**
     *
     * @return theoretical value of @f$cos2\beta@f$
     */
    virtual double computeThValue();
};

#endif /* C2BETA_H */

