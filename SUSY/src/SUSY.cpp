/* 
 * File:   SUSY.cpp
 * Author: marco
 * 
 * Created on December 2, 2010, 3:32 PM
 */

#include "SUSY.h"
#include <math.h>

SUSY::SUSY() :
        StandardModel(), Ru(6,6,0.), Rd(6,6,0.), Rl(6,6,0.), Rn(6,6,0),
        U(2,2,0.), V(2,2,0.), N(4,4,0.), Msu2(6,0.), Msd2(6,0.), Msl2(6,0.),
        Msn2(6,0.), Mch(2,0.), Mneu(4,0.), muH(muH_i){
}

///////////////////////////////////////////////////////////////////////////

double SUSY::v1() {
    return v()*cosb;
}

double SUSY::v2() {
    return v()*sinb;
}


///////////////////////////////////////////////////////////////////////////

void SUSY::setY(double tanb_i) {

    setTanb(tanb_i);

//    Yd.assign(0,0,md/v1()*sqrt(2.));
//    Yd.assign(1,1,ms/v1()*sqrt(2.));
//    Yd.assign(2,2,mb/v1()*sqrt(2.));
//    Yu.assign(0,0,mu/v2()*sqrt(2.));
//    Yu.assign(1,1,mc/v2()*sqrt(2.));
//    Yu.assign(2,2,mt/v2()*sqrt(2.));
//    Yu = Yu*VCKM;
}

void SUSY::setTanb(double tanb) {
    sinb = tanb * sqrt(1. / (1. + tanb * tanb));
    cosb = sqrt(1. / (1. + tanb * tanb));
    this->tanb = tanb;
}

void SUSY::setSinb(double sinb) {
    cosb = sqrt(1. - sinb * sinb);
    tanb = sinb / cosb;
    this->sinb = sinb;
}

void SUSY::setCosb(double cosb) {
    sinb = sqrt(1. - cosb * cosb);
    tanb = sinb / cosb;
    this->cosb = cosb;
}


///////////////////////////////////////////////////////////////////////////

double SUSY::Mw() const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::Mw() " << std::endl;
    return (80.3613);
}

gslpp::complex SUSY::gZf(const int INDF) const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::gZf() " << std::endl;
    gslpp::complex tmp(0.0738065, -0.0120949, false);
    return (tmp);
}

gslpp::complex SUSY::rhoZf(const int INDF) const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::rhoZf() " << std::endl;
    gslpp::complex tmp(1.00516, -0.00473674, false);
    return (tmp);
}

double SUSY::Delta_r() const {

    /* SM + MSSM */

    std::cout << "Write codes for SUSY::Delta_r() " << std::endl;
    return (0.0378211);
}





