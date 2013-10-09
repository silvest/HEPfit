/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Meson.h"

Meson::Meson(double mass, double lifetime=5.e29, double decayconst=0.)
{
    this->mass = mass;
    this->lifetime = lifetime;
    this->decayconst = decayconst;
}

Meson::~Meson() 
{
}

double Meson::getLifetime() const
{
    return lifetime;
}

void Meson::setLifetime(double lifetime)
{
    this->lifetime = lifetime;
}

double Meson::computeWidth() const
{
    return (HCUT / lifetime);
}

double Meson::getDecayconst() const
{
    return decayconst;
}

void Meson::setDecayconst(double decayconst)
{
    this->decayconst = decayconst;
}