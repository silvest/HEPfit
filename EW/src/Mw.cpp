/* 
 * File:   Mw.cpp
 * Author: mishima
 */

#include "Mw.h"


Mw::Mw(const EW& EW_i) : ThObservable(EW_i) {
    myMw = EW_i.getMw();
}

double Mw::getThValue() {   
    return myMw;
}

