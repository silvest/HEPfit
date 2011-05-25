/* 
 * File:   gamma.h
 * Author: silvest
 *
 * Created on April 1, 2011, 2:46 PM
 */

#ifndef GAMMA_H
#define	GAMMA_H

#include <ThObservable.h>
#include <ThObsType.h>

class Gamma : public ThObservable {
public:
    Gamma(const ThObsType& ObsType) : ThObservable(ObsType) {};
    Gamma(const Gamma& orig) : ThObservable(orig.ObsType) {};
    virtual ~Gamma() {};

    double getThValue();
};

#endif	/* GAMMA_H */
