/* 
 * Copyright (C) 2012 SusyFit Collaboration
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
 * @details This class is used to define a meson and three of its
 * characteristics: mass, lifetime and decay constant. 
 * This class inherits the public access members of the Particle class.
 */
class Meson : public Particle {
public:

    /**
     * @brief The default constructor.
     */
    Meson()
    {
    };

    /**
     * @brief Constructor.
     * @param[in] mass the mass of the meson in GeV
     * @param[in] lifetime the lifetime of the meson in \f$ \mathrm{ps}^{-1} \f$
     * @param[in] decayconst the decay constant of the meson in GeV
     */
    Meson(double mass, double lifetime, double decayconst);

    /**
     * @brief The default destructor.
     */
    virtual ~Meson();

    /**
     * @brief A get method for the lifetime of the meson.
     * @return the lifetime of the meson in \f$ \mathrm{ps}^{-1} \f$
     */
    double getLifetime() const
    {
        return lifetime;
    }

    /**
     * @brief A set method for the decay constant of the meson.
     * @param[in] lifetime the lifetime of the meson in \f$ \mathrm{ps}^{-1} \f$
     */
    void setLifetime(double lifetime)
    {
        this->lifetime = lifetime;
    }

    /**
     * @brief A get method for the decay constant of the meson.
     * @return the decay constant of the meson in GeV
     */
    const double& getDecayconst() const
    {
        return decayconst;
    }

    /**
     * @brief A set method for the decay constant of the meson.
     * @param[in] decayconst the decay constant of the meson in GeV
     */
    void setDecayconst(double decayconst)
    {
        this->decayconst = decayconst;
    }

    /**
     * @brief A method to compute the width of the meson from its lifetime.
     * @return the width of the meson in GeV
     */
    double computeWidth() const;

private:
    double decayconst; ///< The decay constant of the meson.
    double lifetime; ///< The lifetime of the meson.

};

#endif	/* MESON_H */

