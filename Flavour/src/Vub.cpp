/* 
 * File:   Vub.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:41 PM
 */

#include "Vub.h"

Vub::Vub(StandardModel * myModel_i) {
    myModel = myModel_i;
}

Vub::Vub(const Vub& orig) {
    myModel = orig.myModel;
}

Vub::~Vub() {
}

double Vub::getThValue() { 
    return((myModel->getMyCKM()).getVub());
}
