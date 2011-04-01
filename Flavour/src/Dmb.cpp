/* 
 * File:   Dmb.cpp
 * Author: silvest
 * 
 * Created on March 29, 2011, 12:51 PM
 */

#include "Dmb.h"

Dmb::Dmb(StandardModel * myModel_i, const int LE_i) {
    myModel = myModel_i;
    LE = LE_i;
}

Dmb::Dmb(const Dmb& orig) {
    myModel = orig.myModel;
    LE = orig.LE;
}

Dmb::~Dmb() {
}

double Dmb::getThValue() { 
    gslpp::complex amp = myModel->getDBD2Amplitude(LE);
    return(amp.abs());
}