/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ZFAFBbottom.h"


double ZFAFBbottom::getThValue() {
    return ( 3.0/4.0*myZF.Af(1)*myZF.Af(9) );
}

