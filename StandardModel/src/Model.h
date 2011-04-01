/* 
 * File:   Model.h
 * Author: marco
 *
 * Created on February 23, 2011, 3:48 PM
 */

#ifndef MODEL_H
#define	MODEL_H

#include "Particle.h"
#include <gslpp_complex.h>
#include <map>

class Model {
public:
    Model() {};
    virtual void update(const std::map<std::string, double>&) = 0;
    virtual bool init(const std::map<std::string, double>&) = 0;
    virtual ~Model() {};
    /**
     * get the @f$\Delta B=\Delta D=2@f$ amplitude
     * @return @f$\langle \bar B_d \vert \mathcal{H}_\mathrm{eff}\vert B_d\rangle@f$ //CHECK!!
     */
    virtual gslpp::complex getDBD2Amplitude(const int LE) const = 0;
private:

};

#endif	/* MODEL_H */

