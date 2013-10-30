/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Vud.h"

double Vud::computeThValue() { 
    return((SM.getCKM()).getVud());
}
