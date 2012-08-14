/* 
 * File:   Particle.cpp
 * Author: marco
 * 
 * Created on February 23, 2011, 3:35 PM
 */

#include "Particle.h"

Particle::Particle(double mass, double mass_scale, double width, double charge) {
    this->mass = mass;
    this->width = width;
    this->charge = charge;
    this->mass_scale = mass_scale;
}

