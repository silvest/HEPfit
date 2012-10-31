/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDK2_H
#define	AMPDK2_H

#include <gslpp_complex.h>
#include "Flavour.h"

using namespace gslpp;

class AmpDK2 {
public:
    /**
     * 
     * @brief comupte the amplitude for kaon oscillations
     * @param Flavour
     */
    AmpDK2(Flavour& Flavour);

protected:
    complex AmpDK(orders order);
    
private:
    Flavour& myFlavour;

};

#endif	/* AMPDK2_H */


