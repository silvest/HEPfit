/* 
 * File:   AmpDD2.h
 * Author: Mauro_87
 *
 * Created on 24 novembre 2011, 16.30
 */

#ifndef AMPDD2_H
#define	AMPDD2_H

#include <gslpp_complex.h>
#include "Flavour.h"

using namespace gslpp;

class AmpDD2 {
/**
 * @brief \f$ | \Delta C = 2 | \f$ Amplitude Class
 * @details This class is related to the calculation of the \f$ D^{0}-\bar{D}^{0}\f$
 * mixing amplitude in relation to the absorptive part of the mixing;  
 * it is written to consider also the contribution coming from supersymmetric theories
 * in the squark-gluino-quark sector in the mass insertion approximation.
 * @param Flavour an object of Flavour class
 */
public:
/**
 * @brief AmpDD2 constructor  
 * @param Flavour an object of Flavour class
 */
    AmpDD2(Flavour& Flavour);
protected:
/**
 * @brief a method returning the complex amplitude for the \f$ D^{0}-\bar{D}^{0}\f$ mixing
 * @param an enum "orders" for the order of perturbation theory of the evolutor 
 * @return the complex \f$ | \Delta C = 2 | \f$ amplitude 
 */
    complex AmpDD(orders order);
    
private:
    Flavour& myFlavour;
};

#endif	/* AMPDD2_H */

