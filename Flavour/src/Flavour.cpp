/* 
 * File:   Flavour.cpp
 * Author: silvest
 * 
 * Created on March 29, 2011, 12:49 PM
 */

#include "Flavour.h"

Flavour::Flavour(StandardModel& Model_i) {
    myModel = &Model_i;
}

Flavour::Flavour(const Flavour& orig) {
    myModel = orig.myModel;
}

Flavour::~Flavour() {
}

