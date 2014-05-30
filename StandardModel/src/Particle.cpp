/* 
 * Copyright (C) 2012 SusyFit Collaboration
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
}

bool Particle::is(std::string name_i) const
{
    if(name_i.compare("QUARK")==0)
        return (this->is("UP") || this->is("CHARM") || this->is("TOP") || this->is("DOWN") || this->is("STRANGE") || this->is("BOTTOM"));
    if(name_i.compare("LEPTON")==0)
        return (this->is("NEUTRINO_1") || this->is("NEUTRINO_2") || this->is("NEUTRINO_3") || this->is("ELECTRON") || this->is("MU") || this->is("TAU"));
    return (name.compare(name_i)==0);
}

int Particle::index() const
{
    if(name.compare("NEUTRINO_1")==0)
        return 0;
    if(name.compare("NEUTRINO_2")==0)
        return 2;
    if(name.compare("NEUTRINO_3")==0)
        return 4;
    if(name.compare("ELECTRON")==0)
        return 1;
    if(name.compare("MU")==0)
        return 3;
    if(name.compare("TAU")==0)
        return 5;
    if(name.compare("UP")==0)
        return 6;
    if(name.compare("CHARM")==0)
        return 8;
    if(name.compare("TOP")==0)
        return 10;
    if(name.compare("DOWN")==0)
        return 7;
    if(name.compare("STRANGE")==0)
        return 9;
    if(name.compare("BOTTOM")==0)
        return 11;
    return -1;
}