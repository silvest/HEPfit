/* 
 * File:   EpsilonK.h
 * Author: stefano
 *
 * Created on 1 dicembre 2011, 10.38
 */

#ifndef EPSILONK_H
#define	EPSILONK_H

#include <ThObservable.h>
#include "Flavour.h"
#include "AmpDK2.h"

using namespace gslpp;

class EpsilonK : public ThObservable, AmpDK2 {
public:   
    /**
     * constructor
     * @param Flavour
     */
    EpsilonK(Flavour& Flavour): ThObservable(Flavour), AmpDK2(Flavour) {};
    
    /**
     * 
     * @return theoretical value of Epsilon_K 
     */
    double getThValue();
    
};

#endif	/* EPSILONK_H */



