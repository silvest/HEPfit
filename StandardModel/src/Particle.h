/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>
#include <string>

#define HCUT 6.58211899E-13 // GeV psec

/**
 * @class Particle
 * @ingroup StandardModel
 * @brief A class for particles. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining particles and their properties
 * like mass (GeV), width (GeV), charge, isospin. One can also specify the scale,
 * in GeV, at which the mass of a particle is defined. 
 */
class Particle {
public:

    /**
     * @brief The default constructor. 
     * @details 
     */
    Particle()
    {
        name = "";
        mass = 0.;
        mass_scale = 0.;
        width = 0.;
        charge = 0.;
        isospin = 0.;
        setIndex();
    };

    /**
     * @brief Constructor. 
     * @details The properties of the particle can be initialized by passing arguments.
     * @param[in] name the name of the particle
     * @param[in] mass the mass of the particle in GeV
     * @param[in] mass_scale the scale in GeV at which the mass is defined, set to 0 by default
     * @param[in] width the decay width of the particle in GeV, set to 0 by default
     * @param[in] charge the charge of the particle, set to 0 by default
     * @param[in] isospin the isospin of the particle, set to 0 by default
     */
    Particle(std::string name, double mass, double mass_scale = 0., double width = 0., double charge = 0., double isospin = 0.);

    virtual ~Particle(){};
    /**
     * @brief A get method to access the particle mass. 
     * @return the particle mass in GeV
     */
    const double& getMass() const
    {
        return mass;
    }

    /**
     * @brief A set method to fix the particle mass. 
     * @param[in] mass the particle mass in GeV
     */
    void setMass(double mass)
    {
        this->mass = mass;
    }

    /**
     * @brief A get method to access the particle width
     * @return the particle width in GeV
     */
    const double& getWidth() const
    {
        return width;
    }

    /**
     * @brief A set method to fix the particle width.
     * @param[in] width the particle width in GeV
     */
    void setWidth(double width)
    {
        this->width = width;
    }

    /**
     * @brief A get method to access the particle charge. 
     * @return the particle charge
     */
    double getCharge() const
    {
        return charge;
    }

    /**
     * @brief A set method to fix the particle charge.
     * @param[in] charge the particle charge
     */
    void setCharge(double charge)
    {
        this->charge = charge;
    }

    /**
     * @brief A get method to access the particle isospin.
     * @return the particle isospin
     */
    double getIsospin() const
    {
        return isospin;
    }

    /**
     * @brief A set method to fix the particle isospin.
     * @param[in] isospin the particle isospin
     */
    void setIsospin(double isospin)
    {
        this->isospin = isospin;
    }

    /**
     * @brief A get method to access the scale at which the particle mass is defined. 
     * @return the scale in GeV at which the particle mass is defined
     */
    double getMass_scale() const
    {
        return mass_scale;
    }

    /**
     * @brief A set method to fix the scale at which the particle mass is defined.
     * @param[in] mass_scale the scale in GeV at which the particle mass is defined
     */
    void setMass_scale(double mass_scale)
    {
        this->mass_scale = mass_scale;
    }

    std::string getName() const
    {
        return name;
    }

    void setName(std::string name)
    {
        this->name = name;
        setIndex();
    }

    bool is(std::string name_i) const;

    int getIndex() const
    {
        return index;
    }

protected:
    double mass; ///< The particle mass in GeV.
    double width; ///< The particle width in GeV.
    double charge; ///< The particle charge.
    double mass_scale; ///< The scale in GeV at which the particle mass is defined.
    double isospin; ///< The particle isospin.
    std::string name; ///< The particle name.
    int index; ///< The index of the particle.

    void setIndex();
};

#endif	/* PARTICLE_H */
