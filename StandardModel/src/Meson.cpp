/* 
 * File:   Meson.cpp
 * Author: silvest
 * 
 * Created on April 12, 2011, 2:17 PM
 */

#include "Meson.h"

Meson::Meson(double mass, double width=0., double decayconst=0.) {
    this->mass = mass;
    this->width = width;
    this->decayconst = decayconst;
}

Meson::~Meson() {
}

