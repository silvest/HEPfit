/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmK.h"

using namespace std;

double DmK::getThValue() {
        return(2.*AmpDK(NLO).real() + SM.getDmk());
}
