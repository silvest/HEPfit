/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_DRHO_H
#define	MW_DRHO_H

#include <ThObservable.h>


#include "SUSY.h"

/**
 * @class Mw_dRho
 * @ingroup SUSY
 * @brief A class for the W-boson mass in the delta rho approximation. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class Mw_dRho : public ThObservable {
public:

    Mw_dRho(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    double computeThValue()
    {
        return (static_cast<const SUSY*> (&SM))->Mw_dRho();
    };

private:

};

#endif	/* MW_DRHO_H */

