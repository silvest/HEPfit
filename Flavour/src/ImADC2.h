/* 
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef IMADC2_H
#define IMADC2_H

#include "ThObservable.h"
#include "AmpDD2.h"

class ImADC2 : public ThObservable, AmpDD2 {
/**
 * @brief a class for the imaginary part of the complex amplitude of D^{0} oscillations, given in \f$ps^-1\f$
 */
public:
    /**
     * ImADC2 constructor
     * @param SM_i a reference to a model
     */
    ImADC2(const StandardModel& SM_i) : ThObservable(SM_i), AmpDD2(SM_i) {};
    /**
     * @brief a method returning the imaginary part of the complex amplitude for 
     * the dispersive part of the\f$ | \Delta C = 2 | \f$ mixing
     * @return the imaginary part of the complex amplitude of \f$ D^{0} \f$ oscillations, given in \f$ps^-1\f$
     */
    double computeThValue() {
        return(AmpDD(FULLNLO).imag());
    };

private:

};

#endif /* IMADC2_H */
