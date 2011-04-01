/* 
 * File:   gamma.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:46 PM
 */

#include "gamma.h"

gammac::gammac(StandardModel * myModel_i) {
    myModel = myModel_i;
}

gammac::gammac(const gammac& orig) {
    myModel = orig.myModel;
}

gammac::~gammac() {
}

double gammac::getThValue() { 
    return((myModel->getMyCKM()).getGamma());
}
