/* 
 * File:   AmpDK2.h
 * Author: stefano
 *
 * Created on 29 novembre 2011, 15.46
 */

#ifndef AMPDK2_H
#define	AMPDK2_H

#include <gslpp_complex.h>
#include "Flavour.h"

using namespace gslpp;

class AmpDK2 {
public:
    /**
     * 
     * @brief comupte the amplitude for kaon oscillations
     * @param Flavour
     */
    AmpDK2(Flavour& Flavour);

protected:
    complex AmpDK(orders order);
    
private:
    Flavour& myFlavour;

};

#endif	/* AMPDK2_H */


