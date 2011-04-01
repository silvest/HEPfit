/* 
 * File:   Particle.cpp
 * Author: marco
 * 
 * Created on February 23, 2011, 3:35 PM
 */

#include "Particle.h"

Particle::Particle(double mass, double width) : bpars(5) {
    this->mass = mass;
    this->width = width;
}

Particle::Particle(double mass) : bpars(5) {
    Particle(mass, 0.);
}

Particle::Particle(const Particle& orig) : bpars(orig.getBpars()) {
    Particle(orig.getMass(), orig.getWidth());
}


