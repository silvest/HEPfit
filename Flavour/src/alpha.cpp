/* 
 * File:   alpha.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:45 PM
 */

#include "alpha.h"

alpha::alpha(StandardModel * myModel_i) {
    myModel = myModel_i;
}

alpha::alpha(const alpha& orig) {
    myModel = orig.myModel;
}

alpha::~alpha() {
}

double alpha::getThValue() { 
    return((myModel->getMyCKM()).getAlpha());
}
