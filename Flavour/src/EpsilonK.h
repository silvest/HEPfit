/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EPSILONK_H
#define	EPSILONK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

using namespace gslpp;

/**
 * @class EpsilonK
 * @ingroup Flavour
 * @brief A class for @f$\epsilon_K@f$ that parametrizes
 * indirect CPV in the Kaon sector
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\epsilon_K@f$. This parameter gets contributions both
 * from the SM and many NP models. It is an important
 * constraint on FCNC from NP models, specially those contributing
 * to chirality flipping @f$\Delta F = 2@f$ operators due to large
 * enhancement by QCD running of these operators to Koan mass scales.
 */

class EpsilonK : public ThObservable, AmpDK2 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonK(Flavour& Flavour): ThObservable(Flavour), AmpDK2(Flavour) {};
    
    /**
     * 
     * @return theoretical value of @f$\epsilon_K@f$ 
     */
    double computeThValue();
    
};

#endif	/* EPSILONK_H */