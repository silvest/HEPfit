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

class CKMGamma : public ThObservable {
public:
    CKMGamma(const ThObsType& ObsType) : ThObservable(ObsType) {};

    double getThValue() { 
        return(SM.getCKM().getGamma()/M_PI*180.);
    };
};

#endif	/* GAMMA_H */
