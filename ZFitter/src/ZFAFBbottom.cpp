/* 
 * File:   ZFAFBbottom.cpp
 * Author: mishima
 */

#include "ZFAFBbottom.h"


double ZFAFBbottom::getThValue() {
    return ( 3.0/4.0*myZF.Af(1)*myZF.Af(9) );
}

