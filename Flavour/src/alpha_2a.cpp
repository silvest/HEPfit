/* 
 * File:   alpha_2a.cpp
 * Author: silvest
 * 
 * Created on April 4, 2011, 10:58 AM
 */

#include "alpha_2a.h"

alpha_2a::alpha_2a(StandardModel * myModel_i) {
    myModel = myModel_i;
}

alpha_2a::alpha_2a(const alpha_2a& orig) {
    myModel = orig.myModel;
}

alpha_2a::~alpha_2a() {
}

double alpha_2a::getThValue() { 
    double twoa = (myModel->getMyCKM()).getAlpha()/M_PI*180.;
    return(twoa < 0. ? twoa + 180. : twoa);
}
