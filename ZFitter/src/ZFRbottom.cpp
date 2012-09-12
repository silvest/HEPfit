/* 
 * File:   ZFRbottom.cpp
 * Author: mishima
 */

#include "ZFRbottom.h"


double ZFRbottom::getThValue() {
    return myZF.Gamma_f(9)/myZF.Gamma_had();
}

