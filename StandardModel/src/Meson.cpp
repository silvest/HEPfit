/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Meson.h"

Meson::Meson(double mass, double width=0., double decayconst=0.) 
{
    this->mass = mass;
    this->width = width;
    this->decayconst = decayconst;
}

Meson::~Meson() 
{
}

