/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BernoulliNumbers.h"


BernoulliNumbers::BernoulliNumbers() 
{
    B[0] = 1.0; B[1] = -1.0/2.0; B[2] = 1.0/6.0; B[3] = 0.0; 
    B[4] = -1.0/30.0; B[5] = 0.0; B[6] = 1.0/42.0; B[7] = 0.0;
    B[8] = -1.0/30.0; B[9] = 0.0; B[10] = 5.0/66.0; B[11] = 0.0; 
    B[12] = -691.0/2730.0; B[13] = 0.0; B[14] = 7.0/6.0; B[15] = 0.0;
    B[16] = -3617.0/510.0; B[17] = 0.0; B[18] = 43867.0/798.0;
}


