/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <string>

#include "Particle.h"

Particle::Particle(std::string name, double mass, double mass_scale, double width, double charge, double isospin)
{
    this->name = name;
    this->mass = mass;
    this->width = width;
    this->charge = charge;
    this->isospin = isospin;
    this->mass_scale = mass_scale;
    setIndex();
}

bool Particle::is(std::string name_i) const
{
    if (name_i.compare("LEPTON") == 0) {
        if (index >= 0 && index <= 5)
            return true;
        else
            return false;
    }
    if (name_i.compare("QUARK") == 0) {
        if (index >= 6 && index <= 11)
            return true;
        else
            return false;
    }
    return (name.compare(name_i) == 0);
}

void Particle::setIndex()
{
    if (name.compare("NEUTRINO_1") == 0)
        index = 0;
    else if (name.compare("NEUTRINO_2") == 0)
        index = 2;
    else if (name.compare("NEUTRINO_3") == 0)
        index = 4;
    else if (name.compare("ELECTRON") == 0)
        index = 1;
    else if (name.compare("MU") == 0)
        index = 3;
    else if (name.compare("TAU") == 0)
        index = 5;
    else if (name.compare("UP") == 0)
        index = 6;
    else if (name.compare("CHARM") == 0)
        index = 8;
    else if (name.compare("TOP") == 0)
        index = 10;
    else if (name.compare("DOWN") == 0)
        index = 7;
    else if (name.compare("STRANGE") == 0)
        index = 9;
    else if (name.compare("BOTTOM") == 0)
        index = 11;
    else
        index = -1;
}
