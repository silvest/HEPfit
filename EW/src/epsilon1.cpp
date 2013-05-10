/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "epsilon1.h"


double epsilon1::getThValue() 
{  
    double eps1 = SM.epsilon1();
    
    //if ( myEW.checkSTUVWXY() ) {
    //    double c2 = myEW.cW2_SM();
    //    double s2 = myEW.sW2_SM();
    //    eps1 += myEW.That()
    //            - myEW.W() + 2.0*sqrt(s2)/sqrt(c2)*myEW.X()- s2/c2*myEW.Y();
    //}
    
    return eps1;
}


