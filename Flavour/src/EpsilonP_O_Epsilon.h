/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONP_O_EPSILON_H
#define	EPSILONP_O_EPSILON_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDS1.h"

/**
 * @class EpsilonK
 * @ingroup Flavour
 * @brief A class for @f$|\epsilon'_K/\epsilon_K|@f$ that parametrizes
 * direct CPV in the Kaon sector
 * @author HEPfit Collaboration
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
    EpsilonP_O_Epsilon(const StandardModel& SM_i): ThObservable(SM_i), AmpDS1(SM_i) {};
    
    /**
     * 
     * @return theoretical value of @f$|\epsilon'_K/\epsilon_K|@f$
     */
    double computeThValue();
    
private:
    
};

#endif	/* EPSILONP_O_EPSILON_H */
