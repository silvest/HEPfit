/* 
 * File:   Meson.cpp
 * Author: silvest
 * 
 * Created on April 12, 2011, 2:17 PM
 */

#include "Meson.h"

Meson::Meson(double mass, double width=0., double decayconst=0.) : bpars(5) {
    this->mass = mass;
    this->width = width;
    this->decayconst = decayconst;
}

Meson::Meson(const Meson& orig) : bpars(orig.getBpars()) {
    Meson(orig.getMass(), orig.getWidth(), orig.getDecayconst());
}

Meson::~Meson() {
}

