/* 
 * File:   epsilon3.cpp
 * Author: mishima
 */

#include "epsilon3.h"


double epsilon3::getThValue() {  
    double DeltaRhoPrime = 2.0*( sqrt(SM.rhoZ_l(SM.ELECTRON).abs()) - 1.0 );
    double DeltaKappaPrime = myEW.sin2thetaEff(SM.ELECTRON)/SM.s02() - 1.0;

    return ( SM.c02()*DeltaRhoPrime + (SM.c02() - SM.s02())*DeltaKappaPrime );
}


