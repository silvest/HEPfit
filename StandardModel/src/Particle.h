/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>

#define HCUT 6.58211899E-13 // GeV psec

/**
 * @class Particle
 * @ingroup StandardModel
 * @brief A class for particles. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining particles and their properties
 * like mass \f$ (GeV) \f$, width \f$ (GeV) \f$, charge, isospin. One can also
 * define the scale, in \f$ GeV \f$, at which the mass of a particle is specified. This class can
 * be used for both fundamental particles and for composite ones like mesons.
 */
class Particle {
public:
    
    /**
     * @brief The default constructor. It sets the scale at which the particle.
     * mass is defined, mass_scale, to 0.
     */
    Particle()
    {
        mass_scale = 0.;
    };
    
    /**
     * @brief The overloaded constructor. It sets the properties of the particle.
     * @param[in] mass the mass of the particle in \f$ GeV \f$
     * @param[in] mass_scale the scale in \f$ GeV \f$ at which the mass is defined, set to 0. by default
     * @param[in] width the decay width of the particle in \f$ GeV \f$, set to 0. by default
     * @param[in] charge the charge of the particle, set to 0. by default
     * @param[in] isospin the isospin of the particle, set to 0. by default
     */
    Particle(double mass, double mass_scale = 0., double width = 0., double charge = 0.,double isospin = 0.);

    /**
     * @return the particle mass in \f$ GeV \f$
     */
    double getMass() const 
    {
        return mass;
    }

    /**
     * set the particle mass
     * @param[in] mass the particle mass in \f$ GeV \f$
     */
    void setMass(double mass) 
    {
        this->mass = mass;
    }

    /**
     * @return the particle width in  \f$ GeV \f$
     */
    double getWidth() const 
    {
        return width;
    }

    /**
     * set the particle width
     * @param[in] width the particle width in  \f$ GeV \f$
     */
    void setWidth(double width) 
    {
        this->width = width;
    }

    /**
     * @return the particle charge
     */
    double getCharge() const 
    {
        return charge;
    }

    /**
     * set the particle charge
     * @param[in] width the particle charge
     */
    void setCharge(double charge) 
    {
        this->charge = charge;
    }        

    /**
     * @return the particle isospin
     */
    double getIsospin() const 
    {
        return isospin;
    }

    /**
     * set the particle isospin
     * @param[in] width the particle isospin
     */
    void setIsospin(double isospin) 
    {
        this->isospin = isospin;
    }
    
    /**
     * @return the scale at which the particle mass is defined
     */
    double getMass_scale() const
    {
        return mass_scale;
    }

    /**
     * set the scale at which the particle mass is defined
     * @param[in] mass_scale the scale at which the particle mass is defined
     */
    void setMass_scale(double mass_scale)
    {
        this->mass_scale = mass_scale;
    }

    
protected:
    double mass;
    double width;
    double charge;
    double mass_scale;
    double isospin;
};

#endif	/* PARTICLE_H */
