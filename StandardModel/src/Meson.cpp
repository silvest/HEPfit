/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Meson.h"

Meson::Meson(double mass, double lifetime = 5.e29, double decayconst = 0., 
        double lambdaM = 0., double gegenalpha1 = 0., double gegenalpha2 = 0.)
{
    this->mass = mass;
    this->lifetime = lifetime;
    this->decayconst = decayconst;
    this->lambdaM = lambdaM;
    gegenalpha[0] = gegenalpha1;
    gegenalpha[1] = gegenalpha2;
}

Meson::~Meson()
{
}

double Meson::computeWidth() const
{
    return (HCUT / lifetime);
}