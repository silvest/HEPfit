/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DmBd.h"
#include <iostream>

using namespace std;

 
 double  DmBd::getThValue() 
 
 { 

     
     std::cout << "Delta MB_d = " << 2.*AmpBd(NLO).abs() << std::endl;
     
     return(2.*AmpBd(NLO).abs()); 

 
 }
 
 