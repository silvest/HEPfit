/* 
 * File:   sin2thetaEff.cpp
 * Author: mishima
 */

#include "sin2thetaEff.h"


sin2thetaEff::sin2thetaEff(const EW& EW_i) : ThObservable(EW_i) {
    sin2_theta_eff = EW_i.sin2thetaEff(SM.ELECTRON);
}

double sin2thetaEff::getThValue() {   
    return sin2_theta_eff;
}

