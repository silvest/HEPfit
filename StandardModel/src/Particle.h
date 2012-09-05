/* 
 * File:   Particle.h
 * Author: marco
 *
 * Created on February 23, 2011, 3:35 PM
 */

#ifndef PARTICLE_H
#define	PARTICLE_H

#include <vector>

#define HCUT 6.58211899E-13 // GeV psec

class Particle {
public:
    Particle() {
    mass_scale = 0.;};
    Particle(double mass, double mass_scale = 0., double width = 0., double charge = 0.);
    /**
     * 
     * @return the particle mass in GeV
     */
    double getMass() const {
        return mass;
    }

    /**
     * set the particle mass
     * @param mass the particle mass in GeV
     */
    void setMass(double mass) {
        this->mass = mass;
    }

    /**
     *
     * @return the particle width in GeV
     */
    double getWidth() const {
        return width;
    }

    /**
     * set the particle width
     * @param width the particle width in GeV
     */
    void setWidth(double width) {
        this->width = width;
    }

    /**
     *
     * @return the particle charge
     */
    double getCharge() const {
        return charge;
    }

    /**
     * set the particle charge
     * @param width the particle charge
     */
    void setCharge(double charge) {
        this->charge = charge;
    }    
    

    /**
     *
     * @return the particle isospin
     */
    double getIsospin() const {
        return isospin;
    }

    /**
     * set the particle isospin
     * @param width the particle isospin
     */
    void setIsospin(double isospin) {
        this->isospin = isospin;
    }
    
    double getMass_scale() const {
        return mass_scale;
    }

    void setMass_scale(double mass_scale) {
        this->mass_scale = mass_scale;
    }

    
protected:
    double mass, width, charge, mass_scale, isospin;
};

#endif	/* PARTICLE_H */
