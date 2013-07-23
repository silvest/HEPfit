/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBd.h"
#include <iostream>

using namespace std;

 
 double  DmBd::getThValue() 
 
 {
     return(2. * AmpBd(FULLNLO).abs());
 }