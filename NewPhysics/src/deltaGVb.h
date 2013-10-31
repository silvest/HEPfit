/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DELTAGVB_H
#define	DELTAGVB_H

#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class deltaGVb
 * @brief A class for @f$\delta g_V^b@f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class deltaGVb : public ThObservable {
public:

    deltaGVb(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double computeThValue()
    {
        return SM.deltaGVq(SM.BOTTOM);
    };
    
private:

};

#endif	/* DELTAGVB_H */

