/* 
 * File:   gamma.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:46 PM
 */

#include "gamma.h"

double Gamma::getThValue() { 
    return(SM.getCKM().getGamma()/M_PI*180.);
}
