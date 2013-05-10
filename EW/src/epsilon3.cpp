/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "epsilon3.h"


double epsilon3::getThValue() 
{  
    double eps3 = SM.epsilon3();
    
    //if ( myEW.checkSTUVWXY() ) {
    //    double c2 = myEW.cW2_SM();
    //    double s2 = myEW.sW2_SM();
    //    eps3 += myEW.Shat()
    //            - myEW.W() + myEW.X()/sqrt(s2)/sqrt(c2) - myEW.Y();
    //}
    
    return eps3;
}


