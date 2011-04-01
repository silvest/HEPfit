/* 
 * File:   Vcb.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:45 PM
 */

#include "Vcb.h"

Vcb::Vcb(StandardModel * myModel_i) {
    myModel = myModel_i;
}

Vcb::Vcb(const Vcb& orig) {
    myModel = orig.myModel;
}

Vcb::~Vcb() {
}

double Vcb::getThValue() { 
    return((myModel->getMyCKM()).getVcb());
}
