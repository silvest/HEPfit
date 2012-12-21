/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Vud.h"

double Vud::getThValue() { 
    return((SM.getCKM()).getVud());
}
