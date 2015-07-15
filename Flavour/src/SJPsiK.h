/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SJPSIK_H
#define	SJPSIK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDB2.h"

/**
 * @class SJPsiK
 * @ingroup Flavour
 * @brief A class for @f$S_{J/\psi K}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$S_{J/\psi K}@f$.
 */
class SJPsiK : public ThObservable, AmpDB2 {
public:
    
    /**
     * @brief Constructor. 
     * @param[in] SM_i
     */
    
    SJPsiK(const StandardModel& SM_i) : ThObservable(SM_i), AmpDB2(SM_i) {};
    
    /**
     *
     * @return theoretical value of @f$S_{J/\psi K}@f$
     */
    virtual double computeThValue();
    
};

#endif	/* SJPSIK_H */

