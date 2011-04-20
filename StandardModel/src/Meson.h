/* 
 * File:   Meson.h
 * Author: silvest
 *
 * Created on April 12, 2011, 2:17 PM
 */

#ifndef MESON_H
#define	MESON_H

#include "Particle.h"

class Meson : public Particle {
public:

    Meson() : bpars(5) {};
    Meson(double mass, double width, double decayconst);
    Meson(const Meson& orig);
    virtual ~Meson();

    /**
     *
     * @return the particle lifetime in sec
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

    std::vector<double> getBpars() const {
        return bpars;
    }

    void setBpars(std::vector<double> v) {
        this->bpars = v;
    }

    void setBpars(const int i, const double value) {
        this->bpars.at(i - 1) = value;
    }


private:
    double decayconst;
    std::vector<double> bpars;

};

#endif	/* MESON_H */

