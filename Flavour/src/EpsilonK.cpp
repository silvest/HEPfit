/* 
 * File:   EpsilonK.cpp
 * Author: stefano
 * Created on 1 dicembre 2011, 10.38
 */

#include "EpsilonK.h"

double EpsilonK::getThValue(){
    std::cout << 1./SM.getDeltaMK() * (AmpDK(NLO).imag()) * SM.getKbarEpsK() / sqrt(2.) << std::endl;
    return(1./SM.getDeltaMK() * (AmpDK(NLO).imag()) * SM.getKbarEpsK() / sqrt(2.));
}
