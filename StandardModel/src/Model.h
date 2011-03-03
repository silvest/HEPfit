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

class Model {
public:
    Model();
    Model(const Model& orig);
    virtual ~Model();
    /**
     * get the @f$\Delta B=\Delta D=2@f$ amplitude
     * @return @f$\langle \bar B_d \vert \mathcal{H}_\mathrm{eff}\vert B_d\rangle@f$ //CHECK!!
     */
    virtual gslpp::complex getDBD2Amplitude(const int LE) const;
private:

};

#endif	/* MODEL_H */

