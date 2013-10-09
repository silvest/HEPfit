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
    Meson(double mass, double lifetime, double decayconst);
    virtual ~Meson();
    
    double getLifetime() const;

    /**
     *
     * @return The particle lifetime in ps. 
     */
    void setLifetime(double lifetime);
    
    double computeWidth() const;

    double getDecayconst() const;

    void setDecayconst(double decayconst);

private:
    double decayconst;
    double lifetime;

};

#endif	/* MESON_H */

