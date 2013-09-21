/*
 * Copyright (C) 2012 SusyFit Collaboration
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
 * @addtogroup Flavour
 * @brief A project for Flavour observables.
 * @{
 */

/**
 * @class SJPsiK
 * @brief A class for @f$S_{J/\psi K}@f$
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$S_{J/\psi K}@f$.
 */

class SJPsiK : public ThObservable, AmpDB2 {
public:
    
    /**
     * constructor
     * @param Flavour
     */
    
    SJPsiK(Flavour& Flavour) : ThObservable(Flavour), AmpDB2(Flavour) {};
    
    /**
     *
     * @return theoretical value of @f$S_{J/\psi K}@f$
     */
    virtual double getThValue();
    
};

/**
 * @}
 */

#endif	/* SJPSIK_H */

