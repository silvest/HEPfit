/* 
 * File:   epsilon3.cpp
 * Author: mishima
 */

#include "epsilon3.h"


double epsilon3::getThValue() {  
    double DeltaRhoPrime = 2.0*( sqrt(SM.rhoZ_l(SM.ELECTRON).abs()) - 1.0 );
    double DeltaKappaPrime = myEW.sin2thetaEff(SM.ELECTRON)/SM.s02() - 1.0;
    double eps3 = SM.c02()*DeltaRhoPrime + (SM.c02() - SM.s02())*DeltaKappaPrime;
    
    if ( myEW.checkModelForSTU() )
        eps3 += myEW.Shat() - myEW.W() 
                + myEW.X()/sqrt(SM.s02())/sqrt(SM.c02()) - myEW.Y();
    
    return eps3;
}


