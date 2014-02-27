/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Particle.h"

Particle::Particle(double mass, double mass_scale, double width, double charge, double isospin) 
{
    this->mass = mass;
    this->width = width;
    this->charge = charge;
    this->isospin = isospin;
    this->mass_scale = mass_scale;
}

