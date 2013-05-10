/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGRB_H
#define	DELTAGRB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

/**
 * @class deltaGRb
 * @brief A class for @f$\delta g_R^b@f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaGRb : public ThObservable {
public:

    deltaGRb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( (SM.deltaGVq(SM.BOTTOM) - SM.deltaGAq(SM.BOTTOM))/2.0 );
    };

private:

};

#endif	/* DELTAGRB_H */

