/* 
 * File:   alpha_2a.cpp
 * Author: silvest
 * 
 * Created on April 4, 2011, 10:58 AM
 */

#include "alpha_2a.h"

double Alpha_2a::getThValue() { 
    double twoa = SM.getCKM().getAlpha()/M_PI*180.;
    return(twoa < 0. ? twoa + 180. : twoa);
}
