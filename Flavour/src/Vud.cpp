/* 
 * File:   Vud.cpp
 * Author: silvest
 * 
 * Created on April 1, 2011, 2:45 PM
 */

#include "Vud.h"

Vud::Vud(StandardModel * myModel_i) {
    myModel = myModel_i;
}

Vud::Vud(const Vud& orig) {
    myModel = orig.myModel;
}

Vud::~Vud() {
}

double Vud::getThValue() { 
    return((myModel->getMyCKM()).getVud());
}
