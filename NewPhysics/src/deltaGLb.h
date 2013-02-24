/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGLB_H
#define	DELTAGLB_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "NPZbbbar.h"

/**
 * @class deltaGLb
 * @brief A class for @f$\delta g_L^b@f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaGLb : public ThObservable {
public:

    deltaGLb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return ( (SM.deltaGVb() + SM.deltaGAb())/2.0 );
    };
    
private:

};

#endif	/* DELTAGLB_H */

