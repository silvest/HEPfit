/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MESON_H
#define	MESON_H

#include "Particle.h"
#include "BParameter.h"

using namespace gslpp;

/**
 * @class Meson
 * @ingroup StandardModel
 * @brief A class for mesons. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Meson : public Particle {
public:

    Meson() 
    {};
    Meson(double mass, double width, double decayconst);
    virtual ~Meson();

    /**
     *
     * @return The particle lifetime in ps. 
     */
    double Lifetime() const 
    {
        return (HCUT / width);
    }

    double getDecayconst() const 
    {
        return decayconst;
    }

    void setDecayconst(double decayconst)
    {
        this->decayconst = decayconst;
    }

private:
    double decayconst;

};

#endif	/* MESON_H */

