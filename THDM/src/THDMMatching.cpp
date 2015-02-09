/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMMatching.h"
#include "THDM.h"
#include <math.h>
#include <stdexcept>

THDMMatching::THDMMatching(const THDM & THDM_i) :

    StandardModelMatching(THDM_i),
    myCKM(3, 3, 0.),
    myTHDM(THDM_i)
{}

void THDMMatching::updateTHDMParameters()
{
    myCKM = myTHDM.getVCKM();
    tanb = myTHDM.getTanb();
    v = myTHDM.v();
    v1 = myTHDM.v1();
    v2 = myTHDM.v2();
    gW = sqrt(8. * myTHDM.getGF() / sqrt(2.)) * myTHDM.Mw_tree();; 
}

