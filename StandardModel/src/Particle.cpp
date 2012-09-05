/* 
 * File:   Particle.cpp
 * Author: marco
 * 
 * Created on February 23, 2011, 3:35 PM
 */

#include "Particle.h"

Particle::Particle(double mass, double width = 0., double charge = 0., double isospin = 0.) {
    this->mass = mass;
    this->width = width;
    this->charge = charge;
    this->isospin = isospin;
}

