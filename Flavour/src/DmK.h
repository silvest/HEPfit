/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DMK_H
#define	DMK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

class DmK : public ThObservable, AmpDK2 {
public:
    /**
     * constructor
     * @param ObsType
     */
    DmK(Flavour& ObsType) : ThObservable(ObsType), AmpDK2(ObsType) {};
    
    /**
     * 
     * @return theoretical value of Delta M_K 
     */
    virtual double getThValue();
};

#endif	/* DMK_H */


