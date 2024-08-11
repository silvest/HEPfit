/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SJPSIK_H
#define	SJPSIK_H

#include "ThObservable.h"

/**
 * @class SJPsiK
 * @ingroup Flavour
 * @brief A class for @f$S_{J/\psi K}@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$S_{J/\psi K}@f$.
 */
class SJPsiK : public ThObservable {
public:
    
    /**
     * @brief Constructor. 
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    
    SJPsiK(const StandardModel& SM_i);
    
    /**
     *
     * @return theoretical value of @f$S_{J/\psi K}@f$
     */
    virtual double computeThValue();
    
};

#endif	/* SJPSIK_H */

