/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGAB_H
#define	DELTAGAB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

/**
 * @addtogroup NewPhysics
 * @brief A project for new physics with general parameterizations. 
 * @{
 */

/**
 * @class deltaGAb
 * @brief A class for @f$\delta g_A^b@f$. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaGAb : public ThObservable {
public:

    deltaGAb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.deltaGAb();
    };
    
private:

};

/** 
 * @}
 */

#endif	/* DELTAGAB_H */

