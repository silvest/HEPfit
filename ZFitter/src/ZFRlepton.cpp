/* 
 * File:   ZFRlepton.cpp
 * Author: mishima
 */

#include "ZFRlepton.h"


double ZFRlepton::getThValue() {
    return myZF.Gamma_had()/myZF.Gamma_f(1);
}

