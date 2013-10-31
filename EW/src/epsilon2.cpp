/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "epsilon2.h"


double epsilon2::computeThValue() 
{  
    double eps2 = SM.epsilon2();
    
    //if ( myEW.checkSTUVWXY() ) {
    //    double c2 = myEW.cW2_SM();
    //    double s2 = myEW.sW2_SM();
    //    eps2 += myEW.Uhat()
    //            - myEW.V() - myEW.W() + 2.0*sqrt(s2)/sqrt(c2)*myEW.X();
    //}
    
    return eps2; 
}
 

