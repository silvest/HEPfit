/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONP_O_EPSILON_H
#define	EPSILONP_O_EPSILON_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDS1.h"

using namespace gslpp;

/**
 * @addtogroup Flavour
 * @brief A project for Flavour observables.
 * @{
 */

/**
 * @class EpsilonK
 * @brief A class for @f$|\epsilon'_K/\epsilon_K|@f$ that parametrizes
 * direct CPV in the Kaon sector
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$|\epsilon'_K/\epsilon_K|@f$. This parameter gets contributions both
 * from the SM and many NP models.
 */

class EpsilonP_O_Epsilon : public ThObservable, AmpDS1 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonP_O_Epsilon(Flavour& Flavour): ThObservable(Flavour), AmpDS1(Flavour) {};
    
    /**
     * 
     * @return theoretical value of @f$|\epsilon'_K/\epsilon_K|@f$
     */
    double computeThValue();
    
private:
    
};

#endif	/* EPSILONP_O_EPSILON_H */
