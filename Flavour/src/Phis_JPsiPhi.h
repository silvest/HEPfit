/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BETASJPSIPHI_H
#define	BETASJPSIPHI_H

#include "ThObservable.h"

/**
 * @class Betas_JPsiPhi
 * @ingroup Flavour
 * @brief A class for @f$S_{J/\psi\phi}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$S_{J/\psi\phi}@f$.
 */
class Phis_JPsiPhi : public ThObservable{
public:
    
   /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Phis_JPsiPhi(const StandardModel& SM_i) : ThObservable(SM_i){};
    
    /**
     *
     * @return theoretical value of @f$S_{J/\psi\phi}@f$
     */
    virtual double computeThValue();

};

#endif	/* BETASJPSIPHI_H */

