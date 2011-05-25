/* 
 * File:   alpha.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:45 PM
 */

#include "alpha.h"

double Alpha::getThValue() {
    return(SM.getCKM().getAlpha()/M_PI*180.);
}
