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
    
    if ( myEW.checkModelForSTU() )
        eps1 += myEW.That() - myEW.W() 
                + 2.0*sqrt(SM.s02())/sqrt(SM.c02())*myEW.X()
                - SM.s02()/SM.c02()*myEW.Y();
    
    return eps1;
}


