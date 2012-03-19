/* 
 * File:   DmK.cpp
 * Author: stefano
 *
 * Created on 10 gennaio 2012, 15.58
 */

#include "DmK.h"

using namespace std;

double DmK::getThValue() {
        return(2.*AmpDK(NLO).real() + SM.getDmk());
}
