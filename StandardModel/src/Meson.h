/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
 * @details The Meson class is used to define a meson and three of its 
 * characteristics: mass, lifetime and decay constant. All three of these
 * are read by the QCD class from the SomeModel.conf file. The suggested 
 * name for the  Model Parameters in the SomeModel.conf file for a meson 
 * \f$ M_a \f$ is mMa (mass), tMa (lifetime) and FMa (decay constant).
 * Please note that these names have to match the ones defind in the QCD class.
 * The Meson class inherits the public access members of the Particle class.
 */
class Meson : public Particle {
public:
    /**
     * @brief Meson constructor
     */
    Meson() 
    {};
    
    Meson(double mass, double lifetime, double decayconst);
    virtual ~Meson();
    
    /**
     * @brief The get method for the lifetime of a meson as specified in the SomeModel.conf file.
     * @return the lifetime of a meson in \f$ ps^{-1} \f$.
     */
    double getLifetime() const
    {
        return lifetime;
    }
    /**
     * @brief The set method for the decay constant of a meson as specified in the SomeModel.conf file.
     * @param[in] lifetime the lifetime of a meson in \f$ ps^{-1} \f$.
     */
    void setLifetime(double lifetime)
    {
        this->lifetime = lifetime;
    }
    /**
     * @brief The get method for the decay constant of a meson as specified in the SomeModel.conf file.
     * @return the decay constant of a meson in \f$ GeV \f$.
     */
    double getDecayconst() const
    {
        return decayconst;
    }
    /**
     * @brief The set method for the decay constant of a meson as specified in the SomeModel.conf file.
     * @param[in] decayconst the decay constant of a meson in \f$ GeV \f$.
     */
    void setDecayconst(double decayconst)
    {
        this->decayconst = decayconst;
    }
    /**
     * @brief Computes the width of the meson from its lifetime \f$ (ps^{-1}) \f$.
     * @return the width of the meson in \f$ GeV \f$.
     */
    double computeWidth() const;

private:
    double decayconst;
    double lifetime;

};

#endif	/* MESON_H */

