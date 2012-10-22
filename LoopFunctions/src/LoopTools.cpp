/* 
 * File:   LoopTools.cpp
 * Author: mishima
 */

#include <iostream>
#include <clooptools.h>
#include "LoopTools.h"

static bool LoopToolsInit = false;

LoopTools::LoopTools() {
    if (!LoopToolsInit) {
        //std::cout << std::endl;
        ltini();
        std::cout << std::endl;
        LoopToolsInit = true;
    }
}

LoopTools::~LoopTools() {
    // for debug
    //std::cout << std::endl
    //          << "************* LoopTools ****************" << std::endl;
    //ltexi();
    //std::cout << "****************************************" << std::endl;
}

double LoopTools::PV_A0(const double mu, const double m) const {
    setmudim(mu*mu);
    return ( - real(A0(m*m)) );
}
    
complex LoopTools::PV_B0(const double mu, const double p2, 
                         const double m0, const double m1) const {
    setmudim(mu*mu);
    return complex(real(B0(p2, m0*m0, m1*m1)), imag(B0(p2, m0*m0, m1*m1)), false);
}
    
complex LoopTools::PV_C0(const double p2, 
                         const double m0, const double m1, const double m2) const {
    return complex(real(-C0(0.0, 0.0, p2, m0*m0, m1*m1, m2*m2)),
                   imag(-C0(0.0, 0.0, p2, m0*m0, m1*m1, m2*m2)), false);
}
    
complex LoopTools::PV_D0(const double s, const double t, const double m0, 
                         const double m1, const double m2, const double m3) const {
    return complex(real(D0(0.0, 0.0, 0.0, 0.0, s, t, m0*m0, m1*m1, m2*m2, m3*m3)), 
                   imag(D0(0.0, 0.0, 0.0, 0.0, s, t, m0*m0, m1*m1, m2*m2, m3*m3)),
                   false);
}


