/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Vcb.h"

double Vcb::getThValue() { 
    return(SM.getCKM().getVcb());
}
