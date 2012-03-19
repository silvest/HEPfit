/* 
 * File:   DmD.h
 * Author: Mauro_87
 *
 * Created on 24 novembre 2011, 22.56
 */

#ifndef M12D_H
#define	M12D_H

#include <ThObservable.h>
#include <Flavour.h>
#include <AmpDD2.h>

class M12D : public ThObservable, AmpDD2 {
/**
 * @brief a class for the absolute value of the complex amplitude of D^{0} oscillations
 * @param Flavour an object of Flavour class
 */
public:
    /**
     * @brief M12D constructor
     * @param Flavour an object of Flavour class
     */
    M12D(Flavour& Flavour) : ThObservable(Flavour), AmpDD2(Flavour) {};
    /**
     * @brief a method returning the absolute value of the complex amplitude for 
     * the absorptive part of the\f$ | \Delta C = 2 | \f$ mixing
     * @return the absolute value of the complex amplitude of \f$ D^{0} \f$ oscillations
     */
    double getThValue() {
        return(AmpDD(NLO).abs());
         
    };
};

#endif	/* MD12_H */

