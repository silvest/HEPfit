/* 
 * File:   DmBd.cpp
 * Author: marco
 * 
 * Created on April 29, 2011, 12:45 PM
 */

#include "DmBd.h"

double DmBd::getThValue() { 
    return(getDmBd(1));
}

double DmBd::getDmBd(int I) { 
    gslpp::complex amp = SM.getDBD2Amplitude(I);
    return(amp.abs());
}
