/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MESON_H
#define	MESON_H

#include "Particle.h"
#include "BParameter.h"

using namespace gslpp;

class Meson : public Particle {
public:

    Meson() {};
    Meson(double mass, double width, double decayconst);
    virtual ~Meson();

    /**
     *
     * @return the particle lifetime in ps
     */
    double Lifetime() const {
        return (HCUT / width);
    }

    double getDecayconst() const {
        return decayconst;
    }

    void setDecayconst(double decayconst) {
        this->decayconst = decayconst;
    }

private:
    double decayconst;

};

#endif	/* MESON_H */

