/* 
 * File:   Vus.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:17 PM
 */

#include "Vus.h"

Vus::Vus(StandardModel * myModel_i) {
    myModel = myModel_i;
}

Vus::Vus(const Vus& orig) {
    myModel = orig.myModel;
}

Vus::~Vus() {
}

double Vus::getThValue() { 
    return((myModel->getMyCKM()).getVus());
}