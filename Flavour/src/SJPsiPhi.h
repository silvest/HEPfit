/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SJPSIPHI_H
#define	SJPSIPHI_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDB2.h"

/**
 * @class SJPsiPhi
 * @ingroup Flavour
 * @brief A class for @f$S_{J/\psi\phi}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$S_{J/\psi\phi}@f$.
 */
class SJPsiPhi : public ThObservable, AmpDB2 {
public:
    
    /**
     * constructor
     * @param Flavour
     */
    
    SJPsiPhi(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};
    
    /**
     *
     * @return theoretical value of @f$S_{J/\psi\phi}@f$
     */
    virtual double computeThValue();

};

#endif	/* SJPSIPHI_H */

